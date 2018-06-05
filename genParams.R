######################
#March 9, 2015
#Generate parametrizations
#######################

if(!"niter"%in%ls()) {niter = 2.5e3;}

if(!"tox_eff_curves"%in%ls()){
  tox_eff_curves = list(
    #Same as 1 above
    #list(tox=c(0.08,0.15,0.25,0.40,0.58),tox_target = 0.25, eff=c(0.05,0.15,0.25,0.40,0.55), eff_min = 0.05, eff_target = 0.20,avg_pref_dose=3),
    list(tox=c(0.05,0.15,0.25,0.35,0.45),tox_target = 0.25, eff=c(0.05,0.15,0.25,0.40,0.55), eff_min = 0.05, eff_target = 0.20,avg_pref_dose=3,jnci_lab="DTC 10 from JNCI"),
    #Same as 2 above
    #list(tox=c(0.04,0.08,0.13,0.19,0.26),tox_target = 0.25, eff=c(0.05,0.06,0.07,0.08,0.33), eff_min = 0.05, eff_target = 0.20,avg_pref_dose=5),
    list(tox=c(0.07,0.10,0.14,0.19,0.25),tox_target = 0.25, eff=c(0.05,0.055,0.06,0.065,0.33), eff_min = 0.05, eff_target = 0.20,avg_pref_dose=5,jnci_lab="DTC 12 from JNCI"),
    #Replace with doses 1--5 from DTC 11 in JNCI paper?
    #list(tox=c(0.10,0.13,0.16,0.27,0.60),tox_target = 0.25, eff=c(0.20,0.23,0.30,0.70,0.90), eff_min = 0.20, eff_target = 0.40,avg_pref_dose=4),
    list(tox=c(0.07,0.11,0.17,0.25,0.35),tox_target = 0.25, eff=c(0.20,0.35,0.355,0.36,0.365), eff_min = 0.20, eff_target = 0.35,avg_pref_dose=c(2,3,4),jnci_lab="DTC 11 from JNCI"),
    #Same as 2 above
    #list(tox=c(0.04,0.08,0.13,0.19,0.26),tox_target = 0.25, eff=c(0.20,0.21,0.22,0.24,0.48), eff_min = 0.20, eff_target = 0.40,avg_pref_dose=5),
    list(tox=c(0.07,0.10,0.14,0.19,0.25),tox_target = 0.25, eff=c(0.20,0.205,0.21,0.215,0.48), eff_min = 0.20, eff_target = 0.35,avg_pref_dose=5,jnci_lab="DTC 12 from JNCI"),
    #New?
    list(tox=c(0.12,0.17,0.23,0.30,0.38),tox_target = 0.25, eff=c(0.05,0.09,0.23,0.33,0.335), eff_min = 0.05, eff_target = 0.20,avg_pref_dose=c(3,4),jnci_lab="New DTC"),
    #Replace with doses 1--5 from DTC 13 in JNCI paper?
    #list(tox=c(0.19,0.24,0.30,0.37,0.60),tox_target = 0.25, eff=c(0.05,0.305,0.31,0.315,0.32),eff_min = 0.05,  eff_target = 0.20,avg_pref_dose=2))
    list(tox=c(0.17,0.25,0.33,0.41,0.49),tox_target = 0.25, eff=c(0.05,0.30,0.305,0.31,0.315),eff_min = 0.05,  eff_target = 0.20,avg_pref_dose=2,jnci_lab="DTC 13 from JNCI"),
    #Scenario 7
    list(tox=c(0.20,0.20,0.20,0.20,0.20),tox_target = 0.25, eff=c(0.05,0.09,0.23,0.33,0.335),eff_min = 0.05,  eff_target = 0.20,avg_pref_dose=c(3,4,5),jnci_lab="New"),
    #Scenario 8
    list(tox=c(0.05,0.05,0.05,0.05,0.05),tox_target = 0.25, eff=c(0.05,0.09,0.23,0.33,0.335),eff_min = 0.05,  eff_target = 0.20,avg_pref_dose=c(3,4,5),jnci_lab="New"),
    #Scenario 9
    list(tox=c(0.07,0.11,0.17,0.25,0.35),tox_target = 0.25, eff=c(0.35,0.35,0.35,0.35,0.35),eff_min = 0.20,  eff_target = 0.35,avg_pref_dose=c(1,2,3,4),jnci_lab="New"),
    #Scenario 10
    list(tox=c(0.05,0.05,0.05,0.05,0.05),tox_target = 0.25, eff=c(0.35,0.35,0.35,0.35,0.35),eff_min = 0.20,  eff_target = 0.35,avg_pref_dose=1:5,jnci_lab="New"))
}

if(!"skeleton_list"%in%ls()) {
  skeleton_list = list(c(0.05,0.15,0.25,0.35,0.45),
                       c(0.07,0.11,0.17,0.25,0.35));
}

if(!"eff_perturb_sd"%in%ls()){
  eff_perturb_sd = 0;
}
if(!"eff_perturb_sd_null"%in%ls()){
  eff_perturb_sd_null = 0;
}

if(!"dec_description"%in%ls()) {
  dec_description = list(
    #15 patients per DEC
    list(n_dec = 5,two_stage_dec = F,
         skip_first_stage = F,
         pocock_alpha = 0.00185,#For global toxicity stopping rule
         ci_level_stage1 = 0.80,
         ci_level_stage2 = NA),#For efficacy analyses
    #30 patients per DEC with interim analysis after 15 patients in each DEC
    list(n_dec = 5,two_stage_dec = T,
         skip_first_stage = F,
         pocock_alpha=0.0014,
         ci_level_stage1 = 0.80,#For efficacy analyses
         ci_level_stage2 = 0.80), 
    #30 patients per DEC with no interim analysis
    list(n_dec = 5,two_stage_dec = T,
         skip_first_stage = T,
         pocock_alpha=0.0014,
         ci_level_stage1 = NA,#For efficacy analyses
         ci_level_stage2 = 0.80))#For efficacy analyses
}

if(!"arglist"%in%ls()){arglist = list();}
sim_number = length(arglist)+1;

#Initial six scenarios
for(i in 1:length(skeleton_list)) {
  for(null_scenario in T) {
    for(j in 1:length(dec_description)) {
      for(k in 1:6) {
        random_seed = sample(.Machine$integer.max,1);
        assign(paste("sim",sim_number,".params",sep=""),list(niter=niter,
                                                             setting = k, 
                                                             null_scenario = null_scenario,
                                                             tox_curve = tox_eff_curves[[k]]$tox,
                                                             tox_target = tox_eff_curves[[k]]$tox_target,
                                                             eff_curve = tox_eff_curves[[k]]$eff,
                                                             eff_curve_null = rep(tox_eff_curves[[k]]$eff_min,5),
                                                             eff_min = tox_eff_curves[[k]]$eff_min,
                                                             eff_target = tox_eff_curves[[k]]$eff_target,
                                                             eff_perturb_sd = eff_perturb_sd,
                                                             eff_perturb_sd_null = eff_perturb_sd_null,
                                                             n_dec = dec_description[[j]]$n_dec,
                                                             two_stage_dec = dec_description[[j]]$two_stage_dec,
                                                             skip_first_stage = dec_description[[j]]$skip_first_stage,
                                                             pocock_alpha = dec_description[[j]]$pocock_alpha,
                                                             ci_level_stage1 = dec_description[[j]]$ci_level_stage1,
                                                             ci_level_stage2 = dec_description[[j]]$ci_level_stage2,
                                                             skeleton = skeleton_list[[i]],
                                                             random_seed=random_seed,sim_number=sim_number))
        arglist = c(arglist,list(get(paste("sim",sim_number,".params",sep=""))))
        sim_number = sim_number + 1;
      }
    }
  }
}

#Extra four scenarios
for(i in 1:length(skeleton_list)) {
  for(null_scenario in T) {
    for(j in 1:length(dec_description)) {
      for(k in 7:10) {
        random_seed = sample(.Machine$integer.max,1);
        assign(paste("sim",sim_number,".params",sep=""),list(niter=niter,
                                                             setting = k, 
                                                             null_scenario = null_scenario,
                                                             tox_curve = tox_eff_curves[[k]]$tox,
                                                             tox_target = tox_eff_curves[[k]]$tox_target,
                                                             eff_curve = tox_eff_curves[[k]]$eff,
                                                             eff_curve_null = rep(tox_eff_curves[[k]]$eff_min,5),
                                                             eff_min = tox_eff_curves[[k]]$eff_min,
                                                             eff_target = tox_eff_curves[[k]]$eff_target,
                                                             eff_perturb_sd = eff_perturb_sd,
                                                             eff_perturb_sd_null = eff_perturb_sd_null,
                                                             n_dec = dec_description[[j]]$n_dec,
                                                             two_stage_dec = dec_description[[j]]$two_stage_dec,
                                                             skip_first_stage = dec_description[[j]]$skip_first_stage,
                                                             pocock_alpha = dec_description[[j]]$pocock_alpha,
                                                             ci_level_stage1 = dec_description[[j]]$ci_level_stage1,
                                                             ci_level_stage2 = dec_description[[j]]$ci_level_stage2,
                                                             skeleton = skeleton_list[[i]],
                                                             random_seed=random_seed,sim_number=sim_number))
        arglist = c(arglist,list(get(paste("sim",sim_number,".params",sep=""))))
        sim_number = sim_number + 1;
      }
    }
  }
}