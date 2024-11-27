/* 
	This do-file complements the paper by Borusyak, Hull, and Jaravel "A Practical Guide to Shift-Share Instruments" (Appx A.9, A.10 and B.1).
	
	In this example we use obtain exposure-robust standard errors for an IV regression with 2 endogenous variables and 2  shift-share IVs. Specifically, it's the original shift-share IV of Autor et al. (2013) and its interaction with any unit-level (i.e., region-level) variable
	
	Inputs: files location_level.dta, industries.dta, Lshares.dta, and shocks.dta from the replication package of Borusyak, Hull, Jaravel (ReStud 2022, henthforth BHJ). Place them in the same folder. See the documentation for that replication package for the meaning of the files.
	
	Author: Kirill Borusyak; November 10, 2024
*/

// Pick any variable to interact the endogenous variable and the instrument with, from the location_level file
global intvar reg_wncen

// List of controls from BHJ Table 4 col.3, plus interactions of the incomplete share controls with the variable of interest 
global controls t2##c.Lsh_manuf##$intvar reg* l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource 

// Benchmark: Uninteracted version (close to BHJ Table 4 col.3 except slightly different controls)
use location_level, clear
ivreg2 y (x=z) $controls [aw=wei], cluster(clus) partial($controls)

* And its shift-level IV equivalent
ssaggregate y x [aw=wei], n(sic87dd) s(ind_share) t(year) sfilename(Lshares) l(czone) controls("$controls") 
merge 1:1 sic87dd year using shocks, assert(3) nogen keepusing(g)
merge m:1 sic87dd using industries, assert(1 3) nogen keepusing(sic3)
ivreg2 y (x=g) [aw=s_n], cluster(sic3)

/* Now interacted with the chosen variable */
// Generate interacted shares
use Lshares, clear
merge m:1 czone year using location_level, assert(3) nogen keepusing($intvar)
gen ind_share_int = ind_share * $intvar
drop ind_share $intvar
save Lshares_int, replace

// Run the interacted regression
use location_level, clear
gen x_int = x * $intvar 
gen z_int = z * $intvar
ivreg2 y (x x_int=z z_int) $controls [aw=wei], cluster(clus) 
predict eps, residuals

// Aggregate it to the industry (=shift) level with each set of shares: in this case, original and interacted
* First with original shares
preserve
	ssaggregate y x x_int eps [aw=wei], n(sic87dd) s(ind_share) t(year) sfilename(Lshares) l(czone) controls("$controls") 
	rename (y x x_int eps s_n) (y1 x1 x_int1 eps1 s_n1)
	tempfile f1
	save `f1'
restore

* Now with interacted shares
ssaggregate y x x_int eps [aw=wei], n(sic87dd) s(ind_share_int) t(year) sfilename(Lshares_int) l(czone) controls("$controls") 
rename (y x x_int eps s_n) (y2 x2 x_int2 eps2 s_n2)
merge 1:1 sic87dd year using `f1', nogen // note: if all interacted shares are zero for some industry, the industry will be dropped from *_2, and the merge will not be perfect
order *2, after(x_int1)

qui {
	// Bring the industry shifts and also industry clustering variable
	merge 1:1 sic87dd year using shocks, assert(3) nogen keepusing(g)
	merge m:1 sic87dd using industries, assert(1 3) nogen keepusing(sic3)

	// Generate the relevant shift-level variables for GMM estimation as in Appendix B.1 of the practical guide (and with similar notation)
	* g_tilde indexed by r (set of weights)
	qui reg g year [aw=s_n1]
	predict g_res1, resid
	qui reg g year [aw=s_n2]
	predict g_res2, resid

	* Elements of the Omega and M matrices: sums of (g_gilde*x_tilde) and (g_tilde*y_tilde)
	* Also psi_k = epsilon_tilde * g_tilde
	* We use the output of ssaggregate: x_tilde1=s_n1*x1, and similar for all other variables
	gen sxg11 = s_n1 * x1 * g_res1 // first index = r (set of weights), second index = j (endog.var)
	gen sxg12 = s_n1 * x_int1 * g_res1 
	gen sxg21 = s_n2 * x2 * g_res2
	gen sxg22 = s_n2 * x_int2 * g_res2 
	gen sy1 = s_n1 * y1 * g_res1
	gen sy2 = s_n2 * y2 * g_res2
	gen psi1 =s_n1 * eps1 * g_res1
	gen psi2 =s_n2 * eps2 * g_res2
	foreach v of varlist sxg* sy* psi* {
		replace `v' = 0 if mi(`v') // missing values may happen if not all industries merged
	}

	* To compute sic3-clustered exposure-robust SE, add up psi by sic3 code
	egen tag_sic3 = tag(sic3)
	egen psi1_sic3 = total(psi1), by(sic3)
	egen psi2_sic3 = total(psi2), by(sic3)

	* Fill the Omega and M matrices
	matrix Omega = J(2,2,.)
	qui sum sxg11
	matrix Omega[1,1] = r(sum)
	qui sum sxg12
	matrix Omega[1,2] = r(sum)
	qui sum sxg21
	matrix Omega[2,1] = r(sum)
	qui sum sxg22
	matrix Omega[2,2] = r(sum)

	matrix M = J(2,1,.)
	qui sum sy1
	matrix M[1,1] = r(sum)
	qui sum sy2
	matrix M[2,1] = r(sum)

	* Compute the mean of the sandwich variance formula, equation (9) in the practical guide
	matrix accum Psi = psi1_sic3 psi2_sic3 if tag_sic3, nocon

	// Compute the coefficient estimates (to double check they match those from the regional regression) and the variance matrix
	matrix b = inv(Omega) * M // coefficient estimates (matches regional level)
	matrix V = inv(Omega) * Psi * inv(Omega)' // variances, without the degree-of-freedom adjustment -- matches ivreg2 but not ivreg
}

// Report the coefficients and standard errors (without the degree-of-freedom adjustment that can be added)
matrix list b
di V[1,1]^0.5 " " V[2,2]^0.5
