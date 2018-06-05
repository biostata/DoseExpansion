#################
#June 10, 2015
#Phil Boonstra
#Code to run a single simulation setting or a suite of simulation settings. 
#If my.work.computer = F, then it is assumed that 
#this is called from a Linux cluster running Slurm
################

rm(list=ls(all=TRUE));
require(dfcrm);require(rstan);require(binom);
my.work.computer = T;
if(my.work.computer){
  setwd("/Users/philb/Desktop/Work/ExpCohortStudy/Efficacy/github");
  titesim_folder = "/users/philb/Desktop/Work/ExpCohortStudy/Efficacy/github/";
  write_to_folder = "/Users/philb/Desktop/Work/ExpCohortStudy/Efficacy/github/";
  array_id = 31;
  niter = 10;
} else {
  setwd("/home/philb/DecEfficacy");
  array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));  
  titesim_folder = "";
  write_to_folder = "";
}

file_name = "DEC_eff";

set.seed(100);
source("genParams.R");

rm(list=setdiff(ls(),c("arglist","my.work.computer","file_name","array_id","offset","titesim_folder")))
source(paste(titesim_folder,"titesimfunctions_phil_V2.R",sep=""));
source("Functions.R");

attach(arglist[[array_id]]);

set.seed(random_seed);
#############################################
###Stage 1a: Dose-escalation + Initial DEC####
#############################################

##3pl3 toxicity data;
dat_3pl3_esc = sim_3pl3_esc(tox_curve,2,niter);
MTD_3pl3 = dat_3pl3_esc$mtd_est;
#Store initialization data for DEC phase
sim_id_forDEC = which(dat_3pl3_esc$mtd_est>0);
dose_num_forDEC = dat_3pl3_esc$mtd_est[dat_3pl3_esc$mtd_est>0];
subj_id_start_forDEC = 1 + dat_3pl3_esc$enrollment[dat_3pl3_esc$mtd_est>0];
avg_n_enroll_3pl3 = mean(dat_3pl3_esc$enrollment[dat_3pl3_esc$mtd_est>0]);
##local stoppage toxicity data;
dat_loc_emp = sim_3pl3_exp(tox_curve = tox_curve,sim_id = sim_id_forDEC,
                           subj_id_start = subj_id_start_forDEC,
                           dose_num = dose_num_forDEC,
                           n_dec = n_dec,
                           dec_size = ifelse(two_stage_dec,30,15),
                           cut_data_at = ifelse(two_stage_dec,15,Inf),
                           deescalation_rule="local",
                           local_thresh = 1/3);
#These take into account stopping for toxicity. They are adjusted
#later in the code to account for stopping for futility. 
n_enroll_loc_emp = 
  n_enroll_loc_mod = dat_loc_emp$enrollment;
#The complete data do not take into account stopping for futility, 
#which is analysis-specific. 
dat_loc_emp = data.frame(dat_loc_emp$all_results);

##global stoppage toxicity data;
dat_glob_emp = sim_3pl3_exp(tox_curve = tox_curve,sim_id = sim_id_forDEC,
                            subj_id_start = subj_id_start_forDEC,
                            dose_num = dose_num_forDEC,
                            n_dec = n_dec,
                            dec_size = ifelse(two_stage_dec,30,15),
                            cut_data_at = ifelse(two_stage_dec,15,Inf),
                            deescalation_rule="global",
                            pocock_alpha = pocock_alpha, pocock_thresh = tox_target);
#These take into account stopping for toxicity. They are adjusted
#later in the code to account for stopping for futility. 
n_enroll_glob_emp = 
  n_enroll_glob_mod = dat_glob_emp$enrollment;
#The complete data do not take into account stopping for futility, 
#which is analysis-specific. 
dat_glob_emp = data.frame(dat_glob_emp$all_results);

##CRM toxicity data;
foo = sim_dec_crm(tox_curve = tox_curve,start = 2,
                  niter = niter,skeleton = skeleton,
                  n_in_esc = avg_n_enroll_3pl3, 
                  n_dec = n_dec, 
                  dec_size = ifelse(two_stage_dec,30,15), 
                  cut_data_at = ifelse(two_stage_dec,15,Inf),
                  tox_target = tox_target);

#The complete data do not take into account stopping for futility, 
#which is analysis-specific. 
dat_crm_emp = data.frame(foo$all_results);
#These take into account stopping for toxicity. They are adjusted
#later in the code to account for stopping for futility. 
n_enroll_crm_emp =
  n_enroll_crm_mod = tapply(dat_crm_emp$subj_id,list(dat_crm_emp$sim_id,dat_crm_emp$dec_id),length);

#Corresponds to the model-based analyses; the beta hats for
#empiric analyses may differ
crm_beta_ests = foo$beta_ests;
crm_beta_scales = foo$beta_scales;
rm(foo);
ptox = sapply(skeleton,"^",exp(crm_beta_ests));
#Calculate MTDs, equal to zero if stopped for toxicity
MTD_crm = apply(cbind(1,abs(ptox-tox_target)/(ptox<=(tox_target+0.05))),1,which.min)-1;
MTD_crm[crm_beta_ests==-Inf] = 0;

#############################################
###Stage 1b: Generate efficacy data###########
#############################################
#Code is written such that each dec in each sim 
#may have a different efficacy curve, which is 
#based on a backbone efficacy curve with 
#random normal perturbations on the logit scale. 
#In the simulations in the manuscript, this feature
#is not used. 
dat_loc_emp[,"eff"] = dat_loc_emp[,"eff_prob"] = 0;
dat_glob_emp[,"eff"] = dat_glob_emp[,"eff_prob"] = 0
dat_crm_emp[,"eff"] = dat_crm_emp[,"eff_prob"] = 0;
all_perturbations = rnorm(niter*n_dec,sd=eff_perturb_sd);
#Matrix (niter by ndec) to keep track of which decs are null
#All are non-null by default. If the scenario is flagged as 
#null, then 3 decs out of 5 are assumed to be 'null', meaning
#that their efficacy curve will be shallow. 
null_decs = matrix(F,n_dec,niter);
if(null_scenario) {
  null_index = rep((0:(niter-1))*n_dec,each=3)+rep(1:3,times=niter);
  all_perturbations[null_index] = rnorm(length(null_index),sd=eff_perturb_sd_null);
  null_decs[null_index] = T;
} 
counter = 0;
store_eff_curves = data.frame(matrix(0,niter*n_dec,3+2*length(eff_curve)));
colnames(store_eff_curves) = c("sim_id","dec_id",paste0("dose",1:length(eff_curve)),paste0("dose",1:length(eff_curve),"_admiss"),"pref_dose");
for(i in 1:niter) {
  for(d in 1:n_dec) {
    counter = counter+1;
    if(!null_decs[d,i]) {
      curr_efficacy_curve = expit(logit(eff_curve)+all_perturbations[counter]);
    } else {
      curr_efficacy_curve = expit(logit(eff_curve_null)+all_perturbations[counter]);
    }
    curr_admissibility = (curr_efficacy_curve>=eff_target)&(tox_curve<=tox_target+0.05);
    if(any(curr_admissibility)) {
      curr_pref = max(which(curr_admissibility));
      other_cands = setdiff(which(curr_admissibility),curr_pref);
      while(length(other_cands)>0){
        if(curr_efficacy_curve[curr_pref]-max(curr_efficacy_curve[other_cands])<0.05) {
          curr_pref = other_cands[length(other_cands)];
          other_cands = other_cands[-length(other_cands)];
        } else {
          break;
        }
      }
    } else {
      curr_pref = 0;
    }
    store_eff_curves[counter,] = c(i,d,curr_efficacy_curve,curr_admissibility,curr_pref);
    curr_loc_index = which(dat_loc_emp[,"sim_id"]==i&dat_loc_emp[,"dec_id"]==d);
    dat_loc_emp[curr_loc_index,"eff_prob"] = curr_efficacy_curve[dat_loc_emp[curr_loc_index,"dose_num"]];
    curr_glob_index = which(dat_glob_emp[,"sim_id"]==i&dat_glob_emp[,"dec_id"]==d);
    dat_glob_emp[curr_glob_index,"eff_prob"] = curr_efficacy_curve[dat_glob_emp[curr_glob_index,"dose_num"]];
    curr_crm_index = which(dat_crm_emp[,"sim_id"]==i&dat_crm_emp[,"dec_id"]==d);
    dat_crm_emp[curr_crm_index,"eff_prob"] = curr_efficacy_curve[dat_crm_emp[curr_crm_index,"dose_num"]];  
  }
}
store_eff_curves = data.frame(store_eff_curves);

dat_loc_emp[,"eff"] = rbinom(nrow(dat_loc_emp),1,dat_loc_emp[,"eff_prob"]);
dat_loc_mod = dat_loc_emp;
dat_glob_emp[,"eff"] = rbinom(nrow(dat_glob_emp),1,dat_glob_emp[,"eff_prob"]);
dat_glob_mod = dat_glob_emp;
dat_crm_emp[,"eff"] = rbinom(nrow(dat_crm_emp),1,dat_crm_emp[,"eff_prob"]);
dat_crm_mod = dat_crm_emp;

######################################################
###Stage 1c: First efficacy analysis#######
###Two approaches: model efficacy curves,
###or estimate proportions to look for 
###evidence of inactivity
###Every analyzed DEC gets a code based on its recommendation:
#101 = no recommendation (toxic during dose escalation)
#201 = no recommendation (toxic during first stage)
#301 = no recommendation (failed first efficacy analysis)
#401 = no recommendation (toxic during second stage; only possible with 30 patient DECs)
#501 = no recommendation (failed second efficacy analysis; only possible with 30 patient DECs)
#601 = recommend dose that is beneath the acceptability range
#701 = recommend dose that is acceptable
#801 = recommend dose that is above the acceptability range
######################################################
#How many bytes large should a stan file be
#in order to trigger a garbage cleanup?
gc_trigger = .1 * (2^30);
oneGrp_stan_compiled = F;

for(curr_design in c("crm","glob","loc")) {
  curr_dat_emp = get(paste0("dat_",curr_design,"_emp"));
  curr_dat_mod = get(paste0("dat_",curr_design,"_mod"));
  curr_n_enroll_emp = get(paste0("n_enroll_",curr_design,"_emp"));
  curr_n_enroll_mod = get(paste0("n_enroll_",curr_design,"_mod"));
  #Restrict only to data from the first stage 
  #(when there is only the first stage, this will just be everything)
  curr_dat_emp_first_stage = curr_dat_emp[curr_dat_emp$first_stage==1,];
  curr_eff_est = cbind(store_eff_curves[,c("sim_id","dec_id")],0,0,101,NA,0,0,101,NA);
  colnames(curr_eff_est) = c("sim_id","dec_id","finalDose_mod","rp2d_mod","rp2dCode_mod","effProbAtMTD_mod","finalDose_emp","rp2d_emp","rp2dCode_emp","effProbAtMTD_emp");
  
  #Which trials did not open DECs are necessarily identical
  #for the empiric and model-based analyses, so we can
  #make both empiric- and model-based decisions based 
  #on the current empiric data 
  #and then adjust the second-stage data separately
  #for each analysis based on the results
  for(curr_sim in unique(curr_dat_emp_first_stage$sim_id)) {
    ####################################
    curr_dat_first_stage_index = which(curr_dat_emp_first_stage$sim_id==curr_sim);
    #Which decs enrolled all 15 patients without triggering toxicity rule
    decs_completed = which(tapply(curr_dat_emp_first_stage[curr_dat_first_stage_index,"dec_id"],curr_dat_emp_first_stage[curr_dat_first_stage_index,"dec_id"],length) == 15);
    ####################################
    
    #RP2D will be the MTD if it meets the efficacy target,
    #otherwise it will be 0. 
    
    ####################################
    #MTD is the same up until this point (difference
    #between analyses is in whether the MTD will be 
    #recommended or not)
    MTD = rep(0,n_dec);
    effProbAtMTD_emp = effProbAtMTD_mod = rep(NA,n_dec);
    rp2d_emp = rp2d_mod = numeric(n_dec);
    rp2dCode_emp = rp2dCode_mod = 201 + numeric(n_dec);
    #Which decs enrolled 15 patients are necessarily 
    #the same for both types of analyses
    for(curr_dec in decs_completed) {
      curr_dat_first_stage_index = which(curr_dat_emp_first_stage$sim_id==curr_sim&curr_dat_emp_first_stage$dec_id==curr_dec);
      curr_admiss_dose = which(store_eff_curves[store_eff_curves[,"sim_id"]==curr_sim&store_eff_curves[,"dec_id"]==curr_dec,grep("_admiss",colnames(store_eff_curves))]==1);
      if(length(curr_admiss_dose)==0) {curr_admiss_dose = -Inf;}
      n = length(curr_dat_first_stage_index);
      y = as.integer(curr_dat_emp_first_stage[curr_dat_first_stage_index,"eff"]);
      x = curr_dat_emp_first_stage[curr_dat_first_stage_index,"dose_num"];
      MTD[curr_dec] = x[n];  
      if(!skip_first_stage) {
        #Model-based Analysis
        if(!oneGrp_stan_compiled) {
          assign("sepMod_stanFit",stan(file = "Bayes2PLogit_OneGrp.stan", data = c("n","y","x"), 
                                       warmup = 250, iter = 3250, chains = 1, thin = 1));
          oneGrp_stan_compiled = T;
        } else { 
          assign("sepMod_stanFit",stan(fit = sepMod_stanFit, data = c("n","y","x"), 
                                       warmup = 250, iter = 3250, chains = 1, thin = 1));
        }  
        sepMod_parEst = extract(sepMod_stanFit);
        fitted_probs = expit(sepMod_parEst$mu + x[n] * sepMod_parEst$beta);
        effProbAtMTD_mod[curr_dec] = mean(fitted_probs);
        rp2d_mod[curr_dec] = ifelse(quantile(fitted_probs,ci_level_stage1)>=eff_min,x[n],0);
        rp2dCode_mod[curr_dec] = 601*(rp2d_mod[curr_dec]<min(curr_admiss_dose)&rp2d_mod[curr_dec]>0)+
          701*(rp2d_mod[curr_dec]%in%curr_admiss_dose)+
          801*(rp2d_mod[curr_dec]>max(curr_admiss_dose)&rp2d_mod[curr_dec]>0)+
          301*(rp2d_mod[curr_dec]==0&MTD[curr_dec]!=0);
        rm(sepMod_parEst);
        
        #Empiric Analysis
        effProbAtMTD_emp[curr_dec] = mean(y[x==x[n]]);
        #Minimum number of patients to continue
        eff_thresh = max(min(which(binom.confint(0:sum(x==x[n]),sum(x==x[n]),2*ci_level_stage1-1,"wilson")[,"upper"]>=eff_min))-1,1);
        rp2d_emp[curr_dec] = ifelse(sum(y[x==x[n]])>=eff_thresh,x[n],0);
        rp2dCode_emp[curr_dec] = 601*(rp2d_emp[curr_dec]<min(curr_admiss_dose)&rp2d_emp[curr_dec]>0)+
          701*(rp2d_emp[curr_dec]%in%curr_admiss_dose)+
          801*(rp2d_emp[curr_dec]>max(curr_admiss_dose)&rp2d_emp[curr_dec]>0)+
          301*(rp2d_emp[curr_dec]==0&MTD[curr_dec]!=0);
      } else {
        effProbAtMTD_mod[curr_dec] = Inf;
        rp2d_mod[curr_dec] = Inf;
        rp2dCode_mod[curr_dec] = Inf;
        effProbAtMTD_emp[curr_dec] = Inf;
        rp2d_emp[curr_dec] = Inf;
        rp2dCode_emp[curr_dec] = Inf;
      }
    }
    #If stopped due to futility at first interim analysis, 
    #the final 15 patients would not have been enrolled. 
    curr_n_enroll_mod[as.character(curr_sim),which(rp2dCode_mod==301)] = 15
    #If stopped due to futility at first interim analysis, 
    #the final 15 patients would not have been enrolled. 
    curr_n_enroll_emp[as.character(curr_sim),which(rp2dCode_emp==301)] = 15;
    ####################################
    
    #Adjust DEC enrollment at second stage:
    #Any DECs that stopped at for futility after 15 patients 
    #would not have enrolled patients in 2nd stage. Based on this, 
    #we shift all patients in the remaining DECs up, and determine their
    #new dose assignment as if the other patients had never been seen
    #This requires determining a revised efficacy probablility and observed outcome
    remove_emp = sum(curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0&(curr_dat_emp$dec_id%in%which(rp2dCode_emp==301)));
    if(remove_emp) {
      new_dec_assignments = which(curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0&(!curr_dat_emp$dec_id%in%which(rp2dCode_emp==301)));
      
      if(length(new_dec_assignments)) {
        curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,][1:length(new_dec_assignments),"dec_id"] = curr_dat_emp[new_dec_assignments,"dec_id"];
        curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,][1:length(new_dec_assignments),"eff_prob"] = 
          store_eff_curves[store_eff_curves$sim_id==curr_sim,][cbind(curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,][1:length(new_dec_assignments),"dec_id"],2+curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,][1:length(new_dec_assignments),"dose_num"])];
        curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,][1:length(new_dec_assignments),"eff"] = 
          rbinom(length(new_dec_assignments),1,curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,][1:length(new_dec_assignments),"eff_prob"]);
        
        curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,"dec_id"][-(1:length(new_dec_assignments))] = Inf;
      } else {
        curr_dat_emp[curr_dat_emp$sim_id==curr_sim&curr_dat_emp$first_stage==0,"dec_id"] = Inf;
      }
      rm(new_dec_assignments);
    }
    rm(remove_emp);
    remove_mod = sum(curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0&(curr_dat_mod$dec_id%in%which(rp2dCode_mod==301)));
    if(remove_mod) {
      new_dec_assignments = which(curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0&(!curr_dat_mod$dec_id%in%which(rp2dCode_mod==301)));
      
      if(length(new_dec_assignments)) {
        curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,][1:length(new_dec_assignments),"dec_id"] = curr_dat_mod[new_dec_assignments,"dec_id"];
        curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,][1:length(new_dec_assignments),"eff_prob"] = 
          store_eff_curves[store_eff_curves$sim_id==curr_sim,][cbind(curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,][1:length(new_dec_assignments),"dec_id"],2+curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,][1:length(new_dec_assignments),"dose_num"])];
        curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,][1:length(new_dec_assignments),"eff"] = 
          rbinom(length(new_dec_assignments),1,curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,][1:length(new_dec_assignments),"eff_prob"]);
        
        curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,"dec_id"][-(1:length(new_dec_assignments))] = Inf;
      } else {
        curr_dat_mod[curr_dat_mod$sim_id==curr_sim&curr_dat_mod$first_stage==0,"dec_id"] = Inf;
      }
      rm(new_dec_assignments);
    }
    rm(remove_mod);
    
    if(exists("sepMod_stanFit")&&(object.size(sepMod_stanFit)>gc_trigger|curr_sim%%249==0)) {
      oneGrp_stan_compiled =  F;
      rm(sepMod_stanFit);
      cat(paste0("gc tiggered: design ", curr_design,"; array_id ",array_id,"; sim ",curr_sim,"\n"),file="gc_list.txt",append=T);
    }
    curr_eff_est[curr_eff_est[,"sim_id"]==curr_sim,c("finalDose_mod","rp2d_mod","rp2dCode_mod","effProbAtMTD_mod","finalDose_emp","rp2d_emp","rp2dCode_emp","effProbAtMTD_emp")] = 
      cbind(MTD,rp2d_mod,rp2dCode_mod,effProbAtMTD_mod,MTD,rp2d_emp,rp2dCode_emp,effProbAtMTD_emp); 
    rm(MTD,rp2d_mod,rp2dCode_mod,effProbAtMTD_mod,rp2d_emp,rp2dCode_emp,effProbAtMTD_emp);
  }
  curr_dat_emp = curr_dat_emp[curr_dat_emp$dec_id<Inf,];
  curr_dat_mod = curr_dat_mod[curr_dat_mod$dec_id<Inf,];
  
  assign(paste0("dat_",curr_design,"_emp"),curr_dat_emp);
  assign(paste0("dat_",curr_design,"_mod"),curr_dat_mod);
  assign(paste0("eff_est_",curr_design),data.frame(curr_eff_est,row.names=NULL));
  assign(paste0("n_enroll_",curr_design,"_emp"),curr_n_enroll_emp);
  assign(paste0("n_enroll_",curr_design,"_mod"),curr_n_enroll_mod);
  rm(curr_dat_emp,curr_dat_mod,curr_dat_first_stage_index,
     decs_completed,curr_eff_est,
     curr_n_enroll_emp,curr_n_enroll_mod,
     n,x,y);
}

####################################
#Garbage cleanup####################
rm(curr_sim,curr_design,curr_dec,sepMod_stanFit);
####################################

if(two_stage_dec) {
  ######################################################
  ###Stage 2: Second efficacy analysis#######
  ###Two approaches: model efficacy curves,
  ###or estimate proportions to look for 
  ###evidence of inactivity
  ######################################################
  #How many bytes large should a stan file be
  #in order to trigger a garbage cleanup?
  gc_trigger = .1 * (2^30);
  oneGrp_stan_compiled = F;
  
  for(curr_design in c("crm","glob","loc")) {
    curr_dat_emp = get(paste0("dat_",curr_design,"_emp"));
    curr_dat_mod = get(paste0("dat_",curr_design,"_mod"));
    #Now go through and revise results for decs that
    #progressed through to the second stage
    #Here the two analysis approaches potentially branch, 
    #as one may have indicated to stop for futility after
    #the first stage whereas the other indicated to proceed
    curr_eff_est = get(paste0("eff_est_",curr_design))
    
    #Model-based
    #Look at sims in which at least one DEC proceeded to the 2nd stage
    for(curr_sim in which(tapply(curr_eff_est$rp2dCode_mod>301,curr_eff_est$sim_id,any))) {
      ####################################
      curr_dat_mod_index = which(curr_dat_mod$sim_id==curr_sim);
      curr_eff_est_index = which(curr_eff_est$sim_id==curr_sim);
      #Which decs enrolled all 30 patients without triggering toxicity rule or interim efficacy rule
      decs_completed = which(tapply(curr_dat_mod[curr_dat_mod_index,"dec_id"],curr_dat_mod[curr_dat_mod_index,"dec_id"],length) == 30 &
                               curr_eff_est[curr_eff_est_index,"rp2dCode_mod"]>301);
      decs_stopped = which(tapply(curr_dat_mod[curr_dat_mod_index,"dec_id"],curr_dat_mod[curr_dat_mod_index,"dec_id"],length) < 30 &
                             curr_eff_est[curr_eff_est_index,"rp2dCode_mod"]>301);
      ####################################
      
      ####################################
      #Separate efficacy models (no sharing between DECs)
      MTD_mod = curr_eff_est[curr_eff_est_index,"finalDose_mod"];
      effProbAtMTD_mod = curr_eff_est[curr_eff_est_index,"effProbAtMTD_mod"];
      rp2d_mod = curr_eff_est[curr_eff_est_index,"rp2d_mod"];
      rp2dCode_mod = curr_eff_est[curr_eff_est_index,"rp2dCode_mod"];
      
      if(length(decs_stopped)) {
        MTD_mod[decs_stopped] = 0;
        effProbAtMTD_mod[decs_stopped] = NA;
        rp2d_mod[decs_stopped] = 0; 
        rp2dCode_mod[decs_stopped] = 401;
      }
      
      for(curr_dec in decs_completed) {
        curr_dat_mod_index = which(curr_dat_mod$sim_id==curr_sim&curr_dat_mod$dec_id==curr_dec);
        curr_admiss_dose = which(store_eff_curves[store_eff_curves[,"sim_id"]==curr_sim&store_eff_curves[,"dec_id"]==curr_dec,grep("_admiss",colnames(store_eff_curves))]==1);
        if(length(curr_admiss_dose)==0) {curr_admiss_dose = -Inf;}
        n = length(curr_dat_mod_index);
        y = as.integer(curr_dat_mod[curr_dat_mod_index,"eff"]);
        x = curr_dat_mod[curr_dat_mod_index,"dose_num"];
        MTD_mod[curr_dec] = x[n];  
        if(!oneGrp_stan_compiled) {
          assign("sepMod_stanFit",stan(file = "Bayes2PLogit_OneGrp.stan", data = c("n","y","x"), 
                                       warmup = 250, iter = 3250, chains = 1, thin = 1));
          oneGrp_stan_compiled = T;
        } else { 
          assign("sepMod_stanFit",stan(fit = sepMod_stanFit, data = c("n","y","x"), 
                                       warmup = 250, iter = 3250, chains = 1, thin = 1));
        }  
        sepMod_parEst = extract(sepMod_stanFit);
        fitted_probs = expit(sepMod_parEst$mu + x[n] * sepMod_parEst$beta);
        effProbAtMTD_mod[curr_dec] = mean(fitted_probs);
        rp2d_mod[curr_dec] = ifelse(quantile(fitted_probs,1-ci_level_stage2)>=eff_min,x[n],0);
        rp2dCode_mod[curr_dec] = 601*(rp2d_mod[curr_dec]<min(curr_admiss_dose)&rp2d_mod[curr_dec]>0)+
          701*(rp2d_mod[curr_dec]%in%curr_admiss_dose)+
          801*(rp2d_mod[curr_dec]>max(curr_admiss_dose)&rp2d_mod[curr_dec]>0)+
          501*(rp2d_mod[curr_dec]==0&MTD_mod[curr_dec]!=0);
      }
      curr_eff_est[curr_eff_est[,"sim_id"]==curr_sim,c("finalDose_mod","rp2d_mod","rp2dCode_mod","effProbAtMTD_mod")] = cbind(MTD_mod,rp2d_mod,rp2dCode_mod,effProbAtMTD_mod); 
      
      if(exists("sepMod_stanFit")&&(object.size(sepMod_stanFit)>gc_trigger|curr_sim%%249==0)) {
        oneGrp_stan_compiled =  F;
        rm(sepMod_stanFit);
        cat(paste0("gc tiggered: design ", curr_design,"; array_id ",array_id,"; sim ",curr_sim,"\n"),file="gc_list.txt",append=T);
      }
    }
    
    rm(rp2d_mod,rp2dCode_mod,effProbAtMTD_mod,sepMod_parEst,
       n,x,y,decs_stopped,decs_completed,
       curr_dat_mod,curr_dat_mod_index,curr_eff_est_index);
    
    #Empirical proportion-based
    #Look at sims in which at least one DEC proceeded to the 2nd stage
    for(curr_sim in which(tapply(curr_eff_est$rp2dCode_emp>301,curr_eff_est$sim_id,any))) {
      ####################################
      curr_dat_emp_index = which(curr_dat_emp$sim_id==curr_sim);
      curr_eff_est_index = which(curr_eff_est$sim_id==curr_sim);
      #Which decs enrolled all 30 patients without triggering toxicity rule or interim efficacy rule
      decs_completed = which(tapply(curr_dat_emp[curr_dat_emp_index,"dec_id"],curr_dat_emp[curr_dat_emp_index,"dec_id"],length) == 30 &
                               curr_eff_est[curr_eff_est_index,"rp2dCode_emp"]>301);
      decs_stopped = which(tapply(curr_dat_emp[curr_dat_emp_index,"dec_id"],curr_dat_emp[curr_dat_emp_index,"dec_id"],length) < 30 &
                             curr_eff_est[curr_eff_est_index,"rp2dCode_emp"]>301);
      ####################################
      
      ####################################
      #Separate efficacy models (no sharing between DECs)
      MTD_emp = curr_eff_est[curr_eff_est_index,"finalDose_emp"];
      effProbAtMTD_emp = curr_eff_est[curr_eff_est_index,"effProbAtMTD_emp"];
      rp2d_emp = curr_eff_est[curr_eff_est_index,"rp2d_emp"];
      rp2dCode_emp = curr_eff_est[curr_eff_est_index,"rp2dCode_emp"];
      
      if(length(decs_stopped)) {
        MTD_emp[decs_stopped] = 0;
        effProbAtMTD_emp[decs_stopped] = NA;
        rp2d_emp[decs_stopped] = 0; 
        rp2dCode_emp[decs_stopped] = 401;
      }
      
      for(curr_dec in decs_completed) {
        curr_dat_emp_index = which(curr_dat_emp$sim_id==curr_sim&curr_dat_emp$dec_id==curr_dec);
        curr_admiss_dose = which(store_eff_curves[store_eff_curves[,"sim_id"]==curr_sim&store_eff_curves[,"dec_id"]==curr_dec,grep("_admiss",colnames(store_eff_curves))]==1);
        if(length(curr_admiss_dose)==0) {curr_admiss_dose = -Inf;}
        n = length(curr_dat_emp_index);
        y = as.integer(curr_dat_emp[curr_dat_emp_index,"eff"]);
        x = curr_dat_emp[curr_dat_emp_index,"dose_num"];
        MTD_emp[curr_dec] = x[n];  
        
        effProbAtMTD_emp[curr_dec] = mean(y[x==x[n]]);
        #Minimum number of patients to continue
        eff_thresh = min(which(binom.confint(0:sum(x==x[n]),sum(x==x[n]),2*ci_level_stage2-1,"wilson")[,"lower"]>=eff_min))-1;
        rp2d_emp[curr_dec] = ifelse(sum(y[x==x[n]])>=eff_thresh,x[n],0);
        rp2dCode_emp[curr_dec] = 601*(rp2d_emp[curr_dec]<min(curr_admiss_dose)&rp2d_emp[curr_dec]>0)+
          701*(rp2d_emp[curr_dec]%in%curr_admiss_dose)+
          801*(rp2d_emp[curr_dec]>max(curr_admiss_dose)&rp2d_emp[curr_dec]>0)+
          501*(rp2d_emp[curr_dec]==0&MTD_emp[curr_dec]!=0);
      }
      
      ####################################
      curr_eff_est[curr_eff_est[,"sim_id"]==curr_sim,c("finalDose_emp","rp2d_emp","rp2dCode_emp","effProbAtMTD_emp")] = cbind(MTD_emp,rp2d_emp,rp2dCode_emp,effProbAtMTD_emp); 
      
    }
    
    rm(rp2d_emp,rp2dCode_emp,effProbAtMTD_emp,
       n,x,y,decs_stopped,decs_completed,
       curr_dat_emp,curr_dat_emp_index,curr_eff_est_index);
    
    assign(paste0("eff_est_",curr_design),data.frame(curr_eff_est,row.names=NULL));
  }
  
}

assign(paste0("sim",array_id),list(params = arglist[[array_id]],
                                   avg_n_enroll_3pl3 = avg_n_enroll_3pl3, 
                                   n_enroll_loc_emp = n_enroll_loc_emp,
                                   n_enroll_loc_mod = n_enroll_loc_mod,
                                   n_enroll_glob_emp = n_enroll_glob_emp,
                                   n_enroll_glob_mod = n_enroll_glob_mod,
                                   n_enroll_crm_emp = n_enroll_crm_emp,
                                   n_enroll_crm_mod = n_enroll_crm_mod,
                                   eff_curves = store_eff_curves,
                                   null_decs = null_decs,
                                   MTD_3pl3 = MTD_3pl3, 
                                   MTD_crm = MTD_crm,
                                   dat_loc_mod = dat_loc_mod,
                                   dat_glob_mod = dat_glob_mod,
                                   dat_crm_mod = dat_crm_mod,
                                   dat_loc_emp = dat_loc_emp,
                                   dat_glob_emp = dat_glob_emp,
                                   dat_crm_emp = dat_crm_emp,
                                   eff_est_loc = eff_est_loc,
                                   eff_est_glob = eff_est_glob,
                                   eff_est_crm = eff_est_crm));


detach();

save(list=paste0("sim",array_id),file=paste(write_to_folder,file_name,array_id,".RData",sep=""));


