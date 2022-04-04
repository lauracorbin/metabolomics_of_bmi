** Code to extract data from source files.

****************************************************************************************************************
** set up a log of activity
capture log close
log using "/--/working/results/logs/01_bmi_rbg_alspac_data_extract.log", replace
****************************************************************************************************************
* Mother BASED files

clear
set maxvar 32767	
use "/--/--/mz_5a.dta", clear
sort aln
gen in_mz=1
merge 1:1 aln using "/--/--/a_3e.dta", nogen
merge 1:1 aln using "/--/--/b_4f.dta", nogen
merge 1:1 aln using "/--/--/c_8a.dta", nogen
merge 1:1 aln using "/--/--/bestgest.dta", nogen

keep aln mz001 mz010a ///
mz013 mz014 mz028b a006 a525 b032 b650 b663 - b667 c645a b720 c755 c765 c800 - c804 ///
bestgest

* Dealing with withdrawal of consent
order aln mz010a, first
order bestgest, last

do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/mother_quest_WoC.do"

save "/--/working/data/analysis01/intermediate/motherQ.dta", replace

*****************************************************************************************************************************************************************************************************************************.
* Child BASED files 

use "/--/--/kz_5c.dta", clear
sort aln qlet
gen in_kz=1
merge 1:1 aln qlet using "/--/--/cp_2b.dta", nogen

keep aln qlet kz011b kz021 kz030 in_core in_alsp in_phase2 in_phase3 in_phase4 tripquad


* Dealing with withdrawal of consent
order aln qlet kz021, first
order in_alsp tripquad, last

do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_based_WoC.do"

drop kz021 tripquad

save "/--/working/data/analysis01/intermediate/childB.dta", replace

*****************************************************************************************************************************************************************************************************************************.
* Child COMPLETED files 

* Phenotype data extraction will be completed in several steps in order to avoid high memory usage.
*** The first step extracts F24 data.
use "/--/--/kz_5c.dta", clear
sort aln qlet
merge 1:1 aln qlet using "/--/--/cp_2b.dta", nogen
merge 1:1 aln qlet using "/--/--/F24_5a.dta", nogen
merge 1:1 aln qlet using "/--/--/YPE_4a.dta", nogen

keep aln qlet kz021 ///
FKMS1040 FKAR0040 FKAR0041 FKAR0002 FKAR0010 FKAR0011 /// bmi, age, data of visit
FKDX1001 FKDX1011 FKDX1021 FKDX1031 FKDX1041 FKDX1051 FKDX1060 FKDX1090 /// fat mass
FKDX1002 FKDX1070 FKDX1080 FKDX1091 /// lean mass
FKMS1052 FKMS1062 /// waist and hip circumference
FKBP1030 FKBP1031 FKBP2020 FKBP2021 /// blood pressure
FKSM1010 FKAL1020 FKAC1110 /// smoke, alcohol, mvpa
FKMS1030 /// weight
FKSA1022 FKSA1023 /// time since last food, Last food/drink was consumed <8 hours
YPE9000-YPE9096 /// food preferences questions
tripquad


* Dealing with withdrawal of consent
order aln qlet kz021, first
order tripquad, last

do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_completed_WoC.do"

* drop kz021 tripquad
save "/--/working/data/analysis01/intermediate/childC_F24.dta", replace

*** The second step extracts TF1 to TF4/F17 data.

use "/--/--/kz_5c.dta", clear
sort aln qlet
merge 1:1 aln qlet using "/--/--/cp_2b.dta", nogen
merge 1:1 aln qlet using "/--/--/tf4_5a.dta", nogen
merge 1:1 aln qlet using "/--/--/tf3_4c.dta", nogen
merge 1:1 aln qlet using "/--/--/tf2_4a.dta", nogen
merge 1:1 aln qlet using "/--/--/tf1_3b.dta", nogen

keep aln qlet kz021 ///
FJMR022a fh3019 fg3139 ff2039 /// bmi
ff2020 fg3120 fh4020 /// waist circumference
FJDX108 FJDX117 FJDX126 FJDX138 FJDX141 FJDX135 /// F17 fat mass
fh2227 fh2236 fh2245 fh2257 fh2260 fh2254 /// TF3 fat mass
fg3227 fg3236 fg3245 fg3257 fg3260 fg3254 /// TF2 fat mass
FJDX127 FJDX118 FJDX109 FJDX139 FJDX142 FJDX136 /// TF4/F17 lean mass
fh2228 fh2237 fh2246 fh2258 fh2261 fh2255 /// TF3 lean mass
fg3228 fg3237 fg3246 fg3258 fg3261 fg3255 /// TF2 lean mass
ff2620 ff2621 ff2622 ff2625 ff2626 ff2627 /// TF1 BP
fg6120 fg6121 fg6122 fg6125 fg6126 fg6127 /// TF2 BP
fh2030 fh2031 fh2032 fh2035 fh2036 fh2037 /// TF3 BP
FJAR019a FJAR019b FJAR019c FJAR020a FJAR020b FJAR020c /// TF4/F17 BP
ff0011a fg0011a fh0011a FJ003a /// age in mth
ff0011b fg0011b fh0011b /// age in wks (no F17)
ff2030 fg3130 fh3010 FJMR022 /// weight
fg1400-fg1601 /// diet diary and food Q's (TF2)
fg1150-fg1167 /// actigraphy (TF2)
fg1200-fg1337 /// actigraphy (TF2)
tripquad


* Dealing with withdrawal of consent
order aln qlet kz021, first
order tripquad, last

do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_completed_WoC.do"

drop kz021 tripquad
save "/--/working/data/analysis01/intermediate/childC_TF1_to_TF4.dta", replace


*** The third step extracts F07 to F11 data.

use "/--/--/kz_5c.dta", clear
sort aln qlet
merge 1:1 aln qlet using "/--/--/cp_2b.dta", nogen
merge 1:1 aln qlet using "/--/--/F11_5d.dta", nogen
merge 1:1 aln qlet using "/--/--/f10_6b.dta", nogen
merge 1:1 aln qlet using "/--/--/f09_4c.dta", nogen
merge 1:1 aln qlet using "/--/--/f08_4d.dta", nogen
merge 1:1 aln qlet using "/--/--/f07_5a.dta", nogen

keep aln qlet kz021 ///
fems026a fdms026a f9ms026a f8lf021a f7ms026a /// bmi
fems018 fdms018 f9ms018 f7ms018 fems020 f9ms020 f7ms020 /// waist and hip circumference
f9dx108 f9dx117 f9dx126 f9dx135 fedx108 fedx117 fedx126 fedx135 /// f9 and f11 fat mass
f9dx109 f9dx118 f9dx127 f9dx136 fedx109 fedx118 fedx127 fedx136 /// f9 and f11 lean mass
f7sa021 f7sa022 f9sa021 f9sa022 fesa021 fesa022 /// f7 f9 and f11 BP
f7003a f8003a f9003a fd003a fe003a /// age in days
f7003b f8003b f9003b fd003b fe003b /// age in wks
f7003c f8003c f9003c fd003c fe003c /// age in mth
f7ms026 f8lf021 f9ms026 fdms026 fems026 /// weight
tripquad


* Dealing with withdrawal of consent
order aln qlet kz021, first
order tripquad, last

do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_completed_WoC.do"

drop kz021 tripquad
save "/--/working/data/analysis01/intermediate/childC_F07_to_F11.dta", replace


*** The final step extracts CIF data.

use "/--/--/kz_5c.dta", clear
sort aln qlet
merge 1:1 aln qlet using "/--/--/cp_2b.dta", nogen
merge 1:1 aln qlet using "/--/--/cif_8b.dta", nogen

keep aln qlet kz021 ///
cf060 - cf069 /// bmi
cf075 - cf079 /// waist circumference
cf123 cf124 cf133 cf134 cf143 cf144 /// 37 49 61 mth 
cf010 cf011 - cf019 /// age in wks
cf040 cf041 cf042 cf043 cf044 cf045 cf046 cf047 cf048 cf049 /// weight
tripquad


* Dealing with withdrawal of consent
order aln qlet kz021, first
order tripquad, last

do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_completed_WoC.do"

drop kz021 tripquad
save "/--/working/data/analysis01/intermediate/childC_CIF.dta", replace


* Blood sample/metabolites data extraction.
use "/--/--/kz_5c.dta", clear
sort aln qlet
merge 1:1 aln qlet using "/--/--/cp_2b.dta", nogen
merge 1:1 aln qlet using "/--/--/Child_bloods_5a.dta", nogen
merge 1:1 aln qlet using "/--/--/child_metabolomics_3a.dta", nogen

keep aln qlet kz021 ///
Glucose_* glucose_* CHOL_* Chol_* chol_* Insulin_* ///
trig_* TRIG_* Trig_* HDL_* hdl_* ///
LDL_* ldl_* VLDL_* vldl_* ///
CRP_* crp_* GGT_* ALT_* AST_* ///
Glc_* ///
Ile_* Leu_* Val_* ///
tripquad


* Dealing with withdrawal of consent
order aln qlet kz021, first
order tripquad, last

do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_completed_WoC.do"

drop kz021 tripquad
save "/--/working/data/analysis01/intermediate/childC_metabolites.dta", replace

*****************************************************************************************************************************************************************************************************************************.
** Matching all data together and saving out the final file*.

use "/--/working/data/analysis01/intermediate/childB.dta", clear
merge 1:1 aln qlet using "/--/working/data/analysis01/intermediate/childC_F24.dta", nogen
merge 1:1 aln qlet using "/--/working/data/analysis01/intermediate/childC_TF1_to_TF4.dta", nogen
merge 1:1 aln qlet using "/--/working/data/analysis01/intermediate/childC_F07_to_F11.dta", nogen
merge 1:1 aln qlet using "/--/working/data/analysis01/intermediate/childC_CIF.dta", nogen
merge 1:1 aln qlet using "/--/working/data/analysis01/intermediate/childC_metabolites.dta", nogen
merge m:1 aln using "/--/working/data/analysis01/intermediate/motherQ.dta", nogen

* Remove non-alspac children.
drop if in_alsp!=1.

* Remove trips and quads.
drop if tripquad==1

drop in_alsp tripquad
save "/--/working/data/analysis01/intermediate/bmi_rbg_meta.dta", replace
*****************************************************************************************************************************************************************************************************************************.
* QC checks*
use "/--/working/data/analysis01/intermediate/bmi_rbg_meta.dta", clear

* Check that there are 15645 records.
count

* Restrict to those who attended Y24 clinic
keep if FKAR0002 == 1
count

* re-save
save "/--/working/data/analysis01/intermediate/bmi_rbg_meta.dta", replace

*****************************************************************************************************************************************************************************************************************************.
* remove intermediate files not needed
rm "/--/working/data/analysis01/intermediate/childB.dta"
rm "/--/working/data/analysis01/intermediate/motherQ.dta"
rm "/--/working/data/analysis01/intermediate/childC_F24.dta"
rm "/--/working/data/analysis01/intermediate/childC_TF1_to_TF4.dta"
rm "/--/working/data/analysis01/intermediate/childC_F07_to_F11.dta"
rm "/--/working/data/analysis01/intermediate/childC_CIF.dta"
rm "/--/working/data/analysis01/intermediate/childC_metabolites.dta"
*
*****************************************************************************************************************************************************************************************************************************.
clear all

log close
