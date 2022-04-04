*** PROJECT: BMI RbG in ALSPAC 
*** Processing data extract created with standard syntax script so it can be read into R

clear all
set more off
version 14.2

cd "/--/working/data/analysis01/"

****************************************************************************************************************
** set up a log of activity
capture log close
log using "../../results/logs/02_bmi_rbg_alspac_data_process.log", replace
****************************************************************************************************************

** read in data file made by extract script
use "intermediate/bmi_rbg_meta.dta", clear

** check count
count


****************************************************************************************************************

** rename variables we need
rename FKAR0010 age_mth_F24
rename FKAR0011 age_years_F24
rename FKAR0040 clinic_month_F24
rename FKAR0041 clinic_year_F24
rename FKAL1020 alc_freq_F24
rename FKSM1010 ever_smoke_F24
rename FKBP1030 seated_SBP_F24
rename FKBP1031 seated_DBP_F24
rename FKBP2020 stand_SBP_F24
rename FKBP2021 stand_DBP_F24
rename FKMS1040 bmi_F24
rename FKMS1052 waist_circ_F24
rename FKMS1062 hip_circ_F24
rename FKDX1001 fat_mass_total_F24
rename FKDX1002 lean_mass_total_F24
rename FKDX1011 fat_mass_arms_F24
rename FKDX1021 fat_mass_legs_F24
rename FKDX1031 fat_mass_trunk_F24
rename FKDX1041 fat_mass_android_F24
rename FKDX1051 fat_mass_gynoid_F24
rename FKDX1060 FMI_F24
rename FKDX1070 LMI_F24
rename FKDX1080 LMI_append_F24
rename FKDX1090 FMI_trunk_legs_F24
rename FKDX1091 LMI_trunk_legs_F24
rename FKAC1110 MVPA_minutes_F24
rename FKSA1022 time_since_last_food_F24 
rename FKSA1023	fasting_status_f24
rename FKMS1030 weight_F24

rename FJMR022a bmi_F17
rename FJ003a age_mth_F17
rename (FJDX108 FJDX109) (fat_mass_arms_F17 lean_mass_arms_F17)
rename (FJDX117 FJDX118) (fat_mass_legs_F17 lean_mass_legs_F17)
rename (FJDX126 FJDX127) (fat_mass_trunk_F17 lean_mass_trunk_F17)
rename (FJDX135 FJDX136) (fat_mass_total_F17 lean_mass_total_F17)
rename (FJDX138 FJDX139) (fat_mass_android_F17 lean_mass_android_F17)
rename (FJDX141 FJDX142) (fat_mass_gynoid_F17 lean_mass_gynoid_f17)
rename (FJAR019a FJAR019b FJAR019c) (right_SBP_F17 right_DBP_F17 right_PBP_F17)
rename (FJAR020a FJAR020b FJAR020c) (left_SBP_F17 left_DBP_F17 left_PBP_F17)
rename FJMR022 weight_F17

rename fh3019 bmi_TF3
rename (fh0011a fh0011b) (age_mth_TF3 age_wks_TF3)
rename (fh2030 - fh2032) (SBP1_TF3 DBP1_TF3 PBP1_TF3)
rename (fh2035 - fh2037) (SBP2_TF3 DBP2_TF3 PBP2_TF3)
rename (fh2227 fh2228) (fat_mass_arms_TF3 lean_mass_arms_TF3)
rename (fh2236 fh2237) (fat_mass_legs_TF3 lean_mass_legs_TF3)
rename (fh2245 fh2246) (fat_mass_trunk_TF3 lean_mass_trunk_TF3)
rename (fh2254 fh2255) (fat_mass_total_TF3 lean_mass_total_TF3)
rename (fh2257 fh2258) (fat_mass_android_TF3 lean_mass_android_TF3)
rename (fh2260 fh2261) (fat_mass_gynoid_TF3 lean_mass_gynoid_TF3)
rename fh4020 waist_circ_TF3
rename fh3010 weight_TF3

rename fg3139 bmi_TF2
rename (fg0011a fg0011b) (age_mth_TF2 age_wks_TF2)
rename fg3120 waist_circ_TF2
rename (fg3227 fg3228) (fat_mass_arms_TF2 lean_mass_arms_TF2)
rename (fg3236 fg3237) (fat_mass_legs_TF2 lean_mass_legs_TF2)
rename (fg3245 fg3246) (fat_mass_trunk_TF2 lean_mass_trunk_TF2)
rename (fg3254 fg3255) (fat_mass_total_TF2 lean_mass_total_TF2)
rename (fg3257 fg3258) (fat_mass_android_TF2 lean_mass_android_TF2)
rename (fg3260 fg3261) (fat_mass_gynoid_TF2 lean_mass_gynoid_TF2)
rename (fg6120 - fg6122) (SBP1_TF2 DBP1_TF2 PBP1_TF2)
rename (fg6125 - fg6127) (SBP2_TF2 DBP2_TF2 PBP2_TF2)
rename fg3130 weight_TF2

rename ff2039 bmi_TF1
rename (ff0011a ff0011b) (age_mth_TF1 age_wks_TF1)
rename ff2020 waist_circ_TF1
rename (ff2620 - ff2622) (SBP1_TF1 DBP1_TF1 PBP1_TF1)
rename (ff2625 - ff2627) (SBP2_TF1 DBP2_TF1 PBP2_TF1)
rename ff2030 weight_TF1

rename fems026a bmi_F11
rename (fe003a fe003b fe003c) (age_days_F11 age_wks_F11 age_mth_F11)
rename fems018 waist_circ_F11
rename fems020 hip_circ_F11
rename (fedx108 fedx109) (fat_mass_arms_F11 lean_mass_arms_F11)
rename (fedx117 fedx118) (fat_mass_legs_F11 lean_mass_legs_F11)
rename (fedx126 fedx127) (fat_mass_trunk_F11 lean_mass_trunk_F11)
rename (fedx135 fedx136) (fat_mass_total_F11 lean_mass_total_F11)
rename (fesa021 fesa022) (SBP_F11 DBP_F11)
rename fems026 weight_F11

rename fdms026a bmi_F10
rename (fd003a fd003b fd003c) (age_days_F10 age_wks_F10 age_mth_F10)
rename fdms018 waist_circ_F10
rename fdms026 weight_F10

rename f9ms026a bmi_F9
rename (f9003a f9003b f9003c) (age_days_F9 age_wks_F9 age_mth_F9)
rename f9ms018 waist_circ_F9
rename f9ms020 hip_circ_F9
rename (f9dx108 f9dx109) (fat_mass_arms_F9 lean_mass_arms_F9)
rename (f9dx117 f9dx118) (fat_mass_legs_F9 lean_mass_legs_F9)
rename (f9dx126 f9dx127) (fat_mass_trunk_F9 lean_mass_trunk_F9)
rename (f9dx135 f9dx136) (fat_mass_total_F9 lean_mass_total_F9)
rename (f9sa021 f9sa022) (SBP_F9 DBP_F9)
rename f9ms026 weight_F9

rename f8lf021a bmi_F8
rename (f8003a f8003b f8003c) (age_days_F8 age_wks_F8 age_mth_F8)
rename f8lf021 weight_F8

rename f7ms026a bmi_F7
rename (f7003a f7003b f7003c) (age_days_F7 age_wks_F7 age_mth_F7)
rename f7ms018 waist_circ_F7
rename f7ms020 hip_circ_F7
rename (f7sa021 f7sa022) (SBP_F7 DBP_F7)
rename f7ms026 weight_F7

rename (cf060 - cf069) (bmi_4mth bmi_8mth bmi_12mth bmi_18mth bmi_25mth bmi_31mth bmi_37mth bmi_43mth bmi_49mth bmi_61mth)
rename (cf010 cf011 - cf019) (age_wks_4mth age_wks_8mth age_wks_12mth age_wks_18mth age_wks_25mth age_wks_31mth age_wks_37mth age_wks_43mth age_wks_49mth age_wks_61mth)
rename (cf075 - cf079) (waist_circ_31mth waist_circ_37mth waist_circ_43mth waist_circ_49mth waist_circ_61mth)
rename (cf123 cf124) (SBP_37mth DBP_37mth)
rename (cf133 cf134) (SBP_49mth DBP_49mth)
rename (cf143 cf144) (SBP_61mth DBP_61mth)
rename (cf040 cf041 cf042 cf043 cf044 cf045 cf046 cf047 cf048 cf049) (weight_4mth weight_8mth weight_12mth weight_18mth weight_25mth weight_31mth weight_37mth weight_43mth weight_49mth weight_61mth)

rename kz021 sex
rename kz030 birthweight
rename b650 mum_ever_smoke
rename b720 mum_alc_freq
rename c645a mum_highest_edu
rename c755 mat_social_class
rename c765 pat_social_class
rename bestgest gestation_at_delivery
rename b032 parity
rename mz028b mum_age_at_delivery

rename (glucose_TF3 glucose_TF4) (Glucose_TF3 Glucose_F17)
rename Glc_TF4 Glc_F17
rename (CHOL_F7 CHOL_F9 chol_TF3 CHOL_TF4) (Chol_F7 Chol_F9 Chol_TF3 Chol_F17)
rename (TRIG_F7 trig_f9 trig_TF3 TRIG_TF4) (Trig_F7 Trig_F9 Trig_TF3 Trig_F17)
rename (HDL_f9 hdl_TF3 HDL_TF4) (HDL_F9 HDL_TF3 HDL_F17)
rename (LDL_f9 ldl_TF3 LDL_TF4) (LDL_F9 LDL_TF3 LDL_F17)
rename (VLDL_f9 vldl_TF3 VLDL_TF4) (VLDL_F9 VLDL_TF3 VLDL_F17)
rename (CRP_f9 crp_TF3 CRP_TF4) (CRP_F9 CRP_TF3 CRP_F17)
rename GGT_TF4 GGT_F17
rename ALT_TF4 ALT_F17
rename AST_TF4 AST_F17

****************************************************************************************************************

** restrict to var we need
keep aln qlet sex birthweight mum_ever_smoke mum_alc_freq mum_highest_edu mat_social_class pat_social_class ///
gestation_at_delivery parity mum_age_at_delivery ///
bmi_4mth bmi_8mth bmi_12mth bmi_18mth bmi_25mth bmi_31mth bmi_37mth bmi_43mth bmi_49mth bmi_61mth ///
SBP_37mth DBP_37mth SBP_49mth DBP_49mth SBP_61mth DBP_61mth ///
bmi_F7 SBP_F7 DBP_F7 ///
bmi_F8 ///
bmi_F9 ///
fat_mass_total_F9 lean_mass_total_F9 SBP_F9 DBP_F9 ///
bmi_F10 ///
bmi_F11 ///
fat_mass_total_F11 lean_mass_total_F11 SBP_F11 DBP_F11 ///
bmi_TF1 SBP1_TF1 DBP1_TF1 PBP1_TF1 SBP2_TF1 DBP2_TF1 PBP2_TF1 ///
bmi_TF2 fat_mass_arms_TF2 lean_mass_arms_TF2 ///
fat_mass_total_TF2 lean_mass_total_TF2 ///
SBP1_TF2 DBP1_TF2 PBP1_TF2 SBP2_TF2 DBP2_TF2 ///
bmi_TF3 SBP1_TF3 DBP1_TF3 PBP1_TF3 SBP2_TF3 DBP2_TF3 ///
fat_mass_total_TF3 lean_mass_total_TF3 ///
waist_circ_TF3 ///
bmi_F17 ///
fat_mass_total_F17 lean_mass_total_F17 ///
right_SBP_F17 right_DBP_F17 left_SBP_F17 left_DBP_F17 ///
clinic_month_F24 clinic_year_F24 alc_freq_F24 ever_smoke_F24 ///
seated_SBP_F24 seated_DBP_F24 stand_SBP_F24 stand_DBP_F24 bmi_F24 ///
fat_mass_total_F24 lean_mass_total_F24 waist_circ_F24 hip_circ_F24 ///
MVPA_minutes_F24 ///
time_since_last_food_F24 fasting_status_f24 ///
age_wks_4mth age_wks_8mth age_wks_12mth age_wks_18mth age_wks_25mth age_wks_31mth age_wks_37mth age_wks_43mth age_wks_49mth age_wks_61mth /// 
age_days_F7 age_wks_F7 age_mth_F7 age_days_F8 age_wks_F8 age_mth_F8 age_days_F9 age_wks_F9 age_mth_F9 ///
age_days_F10 age_wks_F10 age_mth_F10 age_days_F11 age_wks_F11 age_mth_F11 ///
age_mth_TF1 age_wks_TF1 age_mth_TF2 age_wks_TF2 age_mth_TF3 age_wks_TF3 age_mth_F17 age_mth_F24 age_years_F24 ///
weight_4mth weight_8mth weight_12mth weight_18mth weight_25mth weight_31mth weight_37mth weight_43mth weight_49mth weight_61mth ///
weight_F7 weight_F8 weight_F9 weight_F10 weight_F11 weight_TF1 weight_TF2 weight_TF3 weight_F17 weight_F24 ///
Val_F* Val_TF* Leu_F* Leu_TF* Ile_F* Ile_TF* Glucose_F* Glucose_TF* Insulin_* Glc_F* Glc_TF* Chol_CIF* Chol_F* Chol_TF* ///
Trig_CIF* Trig_F* Trig_TF* HDL_CIF* HDL_F* HDL_TF* LDL_CIF* LDL_F* LDL_TF* VLDL_F* VLDL_TF* ///
CRP_F* CRP_TF* ///
YPE9* ///
fg1400 fg1413 fg1414 fg1415 fg1451-fg1460 fg1472-fg1475 fg1510-fg1523 fg1580 fg1581 fg1582 ///



****************************************************************************************************************

** process var where needed


** save out 
save "intermediate/bmi_rbg_meta_processed.dta", replace

*****************************************************************************************************************************************************************************************************************************.

*** add imputed weight data of children
*** script from Katlin Wade

*** phenotype data from growth trajectories ***
use "input/growth data_ long format.dta", clear
codebook aln qlet /*13876*/
*see how many occasions there are
gen num=.
replace num =1 if qlet=="A"
replace num =2 if qlet=="B" 
replace num =3 if qlet=="C" 
replace num =4 if qlet=="D" 
tab qlet
tab num
egen id = concat(aln num)
codebook id /*14068 individuals*/
drop num
destring id, gen(id_num)
bysort id_num: egen numage = count(age_mt)
*generate a number for each clinic for each person
bysort id_num: egen clinic = seq() /*but not everyone were sampled at the same time (i.e., one person's '1' could be another person's '10'*/
generate round_age = round(age_mt)
bysort round_age: egen num = count(round_age)
summ num
bysort id_num: egen mode = mode(round_age)
hist mode
sort mode /*some people have multiple measures taken around the same time, so randomly dropping any duplicates in rounded age*/
duplicates drop aln qlet round_age, force
drop mode
bysort id_num: egen mode = mode(round_age) /*everyone now has unique round_age*/
drop mode clinic num
*reshape
drop source source2
reshape wide weight height bmi age_mt, i(id_num) j(round_age)
summ aln
tab qlet
summ id_num
summ weight* height* bmi* age* 
drop sex

*** apply consent withdrawn
* add dummy var to make wWoC scripts run
gen kz021 = "dummyvar"
gen tripquad = "dummyvar"
gen kz011b = "dummyvar"
gen in_phase4 = "dummyvar"
gen in_alsp = "dummyvar"
* run cw code
do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_completed_WoC.do"
do "/--/working/scripts/bmi_rbg/alspac_woc_20220208/child_based_WoC.do"
drop tripquad kz021 kz011b in_phase4 in_alsp
save "/--/working/data/analysis01/intermediate/growth_clinic_data.dta", replace


*** merge phenotype data together ***
use "/--/working/data/analysis01/intermediate/bmi_rbg_meta_processed.dta", clear
merge 1:1 aln qlet using "/--/working/data/analysis01/intermediate/growth_clinic_data.dta"
drop if _merge==2
drop _merge
order aln qlet
*replace any missing age, weight and BMI values between 4m and 10 years with that from the growth trajectories input (sources: medical records, clinics, quest) data as follows
summ age_wks_* age_mth_*
foreach time in 4m 8m 12m 18m 25m 31m 37m 43m 49m 61m {
	generate age_mth_`time' = age_wks_`time'/4.34524
}
summ age*
summ age_mth_4m
replace age_mth_4m = age_mt4 if age_mth_4m==. | age_mth_4m < 0 
summ age_mth_4m age_mth_8m

egen avage_8_9m = rowmean(age_mt8 age_mt9)
replace age_mth_8m = avage_8_9m if age_mth_8m==. | age_mth_8m < 0 
summ age_mth_8m age_mth_12m

egen avage_12_14m = rowmean(age_mt12 age_mt13 age_mt14)
replace age_mth_12m = avage_12_14m if age_mth_12m==. | age_mth_12m < 0 
summ age_mth_12m age_mth_18m

egen avage_18_20m = rowmean(age_mt18 age_mt19 age_mt20)
replace age_mth_18m = avage_18_20m if age_mth_18m==. | age_mth_18m < 0 
summ age_mth_18m age_mth_25m 

egen avage_25_26m = rowmean(age_mt25 age_mt26)
replace age_mth_25m = avage_25_26m if age_mth_25m==. | age_mth_25m < 0 
summ age_mth_25m age_mth_31m 

egen avage_31_32m = rowmean(age_mt31 age_mt32)
replace age_mth_31m = avage_31_32m if age_mth_31m==. | age_mth_31m < 0 
summ age_mth_31m age_mth_37m 

egen avage_36_38m = rowmean(age_mt36 age_mt37 age_mt38)
replace age_mth_37m = avage_36_38m if age_mth_37m==. | age_mth_37m < 0 
summ age_mth_37m age_mth_43m 

egen avage_43_45m = rowmean(age_mt43 age_mt44 age_mt45)
replace age_mth_43m = avage_43_45m if age_mth_43m==. | age_mth_43m < 0 
summ age_mth_43m age_mth_49m 

egen avage_48_50m = rowmean(age_mt48 age_mt49 age_mt50)
replace age_mth_49m = avage_48_50m if age_mth_49m==. | age_mth_49m < 0 
summ age_mth_49m age_mth_61m 

egen avage_60_67m = rowmean(age_mt60 age_mt61 age_mt62 age_mt63 age_mt64 age_mt65 age_mt66 age_mt67)
replace age_mth_61m = avage_60_67m if age_mth_61m==. |  age_mth_61m < 0 
summ age_mth_61m age_mth_F7 

egen avage_83_105m = rowmean(age_mt83 age_mt84 age_mt85 age_mt86 age_mt87 age_mt88 age_mt89 age_mt90 age_mt91 age_mt92 age_mt93 age_mt94 age_mt95 age_mt96 age_mt97 age_mt98 age_mt99 age_mt100 age_mt101 age_mt102 age_mt103 age_mt104 age_mt105)
replace age_mth_F7 = avage_83_105m if age_mth_F7==. |  age_mth_F7 < 0 
summ age_mth_F7 age_mth_F9 

egen avage_106_120m = rowmean(age_mt106 age_mt107 age_mt108 age_mt109 age_mt110 age_mt111 age_mt112 age_mt113 age_mt114 age_mt115 age_mt116 age_mt117 age_mt118 age_mt119 age_mt120)
replace age_mth_F9 = avage_106_120m if age_mth_F9==. | age_mth_F9 < 0 
summ age_mth_F9  

drop avage*

summ weight_*
summ weight_4mth
replace weight_4mth = weight4 if weight_4mth==. | weight_4mth < 0 
summ weight_4mth weight_8mth

egen avweight_8_9m = rowmean(weight8 weight9)
replace weight_8mth = avweight_8_9m if weight_8mth==. | weight_8mth < 0
summ weight_8mth weight_12mth

egen avweight_12_14m = rowmean(weight12 weight13 weight14)
replace weight_12mth = avweight_12_14m if weight_12mth==.  | weight_12mth < 0
summ weight_12mth weight_18mth

egen avweight_18_20m = rowmean(weight18 weight19 weight20)
replace weight_18mth = avweight_18_20m if weight_18mth==.  | weight_18mth < 0
summ weight_18mth weight_25mth

egen avweight_25_26m = rowmean(weight25 weight26)
replace weight_25mth = avweight_25_26m if weight_25mth==.  | weight_25mth < 0
summ weight_25mth weight_31mth 

egen avweight_31_32m = rowmean(weight31 weight32)
replace weight_31mth = avweight_31_32m if weight_31mth==.  | weight_31mth < 0
summ weight_31mth weight_37mth 

egen avweight_36_38m = rowmean(weight36 weight37 weight38)
replace weight_37mth = avweight_36_38m if weight_37mth==.  | weight_37mth < 0
summ weight_37mth weight_43mth 

egen avweight_43_45m = rowmean(weight43 weight44 weight45)
replace weight_43mth = avweight_43_45m if weight_43mth==.  | weight_43mth < 0
summ weight_43mth weight_49mth 

egen avweight_48_50m = rowmean(weight48 weight49 weight50)
replace weight_49mth = avweight_48_50m if weight_49mth==.  | weight_49mth < 0
summ weight_49mth weight_61mth 

egen avweight_60_67m = rowmean(weight60 weight61 weight62 weight63 weight64 weight65 weight66 weight67)
replace weight_61mth = avweight_60_67m if weight_61mth==.  | weight_61mth < 0
summ weight_61mth weight_F7 

egen avweight_83_105m = rowmean(weight83 weight84 weight85 weight86 weight87 weight88 weight89 weight90 weight91 weight92 weight93 weight94 weight95 weight96 weight97 weight98 weight99 weight100 weight101 weight102 weight103 weight104 weight105)
replace weight_F7 = avweight_83_105m if weight_F7==.  | weight_F7 < 0
summ weight_F7 weight_F9 

egen avweight_106_120m = rowmean(weight106 weight107 weight108 weight109 weight110 weight111 weight112 weight113 weight114 weight115 weight116 weight117 weight118 weight119 weight120)
replace weight_F9 = avweight_106_120m if weight_F9==.  | weight_F9 < 0
summ weight_F9  

drop avweight*

summ bmi_*
summ bmi_4mth
replace bmi_4mth = bmi4 if bmi_4mth==. | bmi_4mth < 0
summ bmi_4mth bmi_8mth
egen avbmi_8_9m = rowmean(bmi8 bmi9)
replace bmi_8mth = avbmi_8_9m if bmi_8mth==. | bmi_8mth < 0
summ bmi_8mth bmi_12mth
egen avbmi_12_14m = rowmean(bmi12 bmi13 bmi14)
replace bmi_12mth = avbmi_12_14m if bmi_12mth==.  | bmi_12mth <0
summ bmi_12mth bmi_18mth
egen avbmi_18_20m = rowmean(bmi18 bmi19 bmi20)
replace bmi_18mth = avbmi_18_20m if bmi_18mth==.  | bmi_18mth < 0
summ bmi_18mth bmi_25mth 
egen avbmi_25_26m = rowmean(bmi25 bmi26)
replace bmi_25mth = avbmi_25_26m if bmi_25mth==.  | bmi_25mth < 0
summ bmi_25mth bmi_31mth 
egen avbmi_31_32m = rowmean(bmi31 bmi32)
replace bmi_31mth = avbmi_31_32m if bmi_31mth==.  | bmi_31mth < 0
summ bmi_31mth bmi_37mth 
egen avbmi_36_38m = rowmean(bmi36 bmi37 bmi38)
replace bmi_37mth = avbmi_36_38m if bmi_37mth==.  | bmi_37mth < 0
summ bmi_37mth bmi_43mth 
egen avbmi_43_45m = rowmean(bmi43 bmi44 bmi45)
replace bmi_43mth = avbmi_43_45m if bmi_43mth==.  | bmi_43mth < 0
summ bmi_43mth bmi_49mth 
egen avbmi_48_50m = rowmean(bmi48 bmi49 bmi50)
replace bmi_49mth = avbmi_48_50m if bmi_49mth==.  | bmi_49mth < 0
summ bmi_49mth bmi_61mth 
egen avbmi_60_67m = rowmean(bmi60 bmi61 bmi62 bmi63 bmi64 bmi65 bmi66 bmi67)
replace bmi_61mth = avbmi_60_67m if bmi_61mth==.  | bmi_61mth < 0
summ bmi_61mth bmi_F7 
egen avbmi_83_105m = rowmean(bmi83 bmi84 bmi85 bmi86 bmi87 bmi88 bmi89 bmi90 bmi91 bmi92 bmi93 bmi94 bmi95 bmi96 bmi97 bmi98 bmi99 bmi100 bmi101 bmi102 bmi103 bmi104 bmi105)
replace bmi_F7 = avbmi_83_105m if bmi_F7==.  | bmi_F7 < 0
summ bmi_F7 bmi_F9 
egen avbmi_106_120m = rowmean(bmi106 bmi107 bmi108 bmi109 bmi110 bmi111 bmi112 bmi113 bmi114 bmi115 bmi116 bmi117 bmi118 bmi119 bmi120)
replace bmi_F9 = avbmi_106_120m if bmi_F9==.  | bmi_F9 < 0
summ bmi_F9  
drop avbmi*

*drop unused data
drop id_num-numage

*make sure all agewk  data is populated
foreach time in 4m 8m 12m 18m 25m 31m 37m 43m 49m 61m F7 F9 F10 F11 TF1 TF2 TF3 {
	replace age_wks_`time' = age_mth_`time'*4.34524 if age_wks_`time' ==.
}
foreach time in F17 F24 {
	generate age_wks_`time' = age_mth_`time'*4.34524
}
summ

* replace '.' in non-character var with -9999
order aln qlet, first
order  sex, last
foreach var of varlist aln- sex {
	capture confirm string var `var'
		if _rc==0 {
		}
		else {
			qui replace `var' = -9999 if  `var' == .
		}
}
		
export delimited using "/--/working/data/analysis01/intermediate/bmi_rbg_meta_processed_with_imputed_weight.csv",  replace nolabel
erase "/--/working/data/analysis01/intermediate/growth_clinic_data.dta"
erase "/--/working/data/analysis01/intermediate/bmi_rbg_meta_processed.dta"

****************************************************************************************************************

clear all

log close
