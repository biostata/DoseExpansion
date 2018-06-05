
expit = function(x) {1/(1+exp(-x));}
logit = function(x) {log(1/(1/x-1));}

#Simulate the escalation phase of a 3plus3
#The strategy is to efficiently create all possible
#dose assignments first, then go through and
#prune and reorder as necessary. 

sim_3pl3_esc = function(tox_curve,start=2,niter,seed=sample(.Machine$integer.max,1)) {
  
  set.seed(seed);
  #Create all possible dose assignments a priori
  all_results = matrix(0,nrow=6*length(tox_curve)*niter,ncol=5);
  colnames(all_results) = c("sim_id","subj_id","dose_num","tox_prob","tox");
  
  reorder_all_rows = rep(NA,nrow(all_results));
  mtd_est = rep(NA,niter);
  
  all_results[,"sim_id"] = rep(1:niter,each=6*length(tox_curve));
  all_results[,"dose_num"] = rep(rep(c(start,(1:length(tox_curve))[-start]),each=6),times=niter);
  all_results[,"tox_prob"] = rep(rep(tox_curve[c(start,(1:length(tox_curve))[-start])],each=6),times=niter);
  all_results[,"tox"] = rbinom(nrow(all_results),1,all_results[,"tox_prob"]);
  
  for(i in 1:niter) {
    sim_index = which(all_results[,"sim_id"]==i);
    curr_results = all_results[sim_index,];
    curr_dose = start;
    curr_dose_index = row_order = 1:3;
    max_safe_dose = length(tox_curve);
    #Keep track of whether six have already been at this dose or not
    six_at_dose = rep(F,length(tox_curve));
    finished = F;
    while(!finished) {
      #If no toxicity at all
      if(sum(curr_results[curr_dose_index,"tox"])==0) {
        #escalate if the next dose is safe
        if(curr_dose<max_safe_dose) {
          curr_dose = curr_dose + 1;
          curr_dose_index = which(curr_results[,"dose_num"]==curr_dose)[1:3];
          row_order = c(row_order,curr_dose_index);
          next;
        } else {
          #if six are at the current dose, this is the MTD
          if(six_at_dose[curr_dose]) {
            finished = T;
            next;
            #otherwise enroll three more
          } else {
            curr_dose_index = which(curr_results[,"dose_num"]==curr_dose);
            row_order = c(row_order,curr_dose_index[4:6]);
            six_at_dose[curr_dose]=T;
            next;    
          }
        }
        #If one toxicity
      } else if(sum(curr_results[curr_dose_index,"tox"])==1) {
        #If already six at the current dose
        if(six_at_dose[curr_dose]) {
          #escalate if safe to
          if(curr_dose<max_safe_dose) {
            curr_dose = curr_dose + 1;
            curr_dose_index = which(curr_results[,"dose_num"]==curr_dose)[1:3];
            row_order = c(row_order,curr_dose_index);
            next;  
            #otherwise this is the MTD
          } else {
            finished = T;
            next;
          }
          #Otherwise enroll three more  
        } else {
          curr_dose_index = which(curr_results[,"dose_num"]==curr_dose);
          row_order = c(row_order,curr_dose_index[4:6]);
          six_at_dose[curr_dose]=T;
          next;    
        }
        #Otherwise there are 2 or more toxicities, so de-escalate  
      } else {
        curr_dose = max_safe_dose = curr_dose - 1;
        #stop if everything is too toxic
        if(curr_dose==0) {
          finished = T;
          next;
        } else {
          #if six are already at the lower dose, this is the MTD
          if(six_at_dose[curr_dose]) {
            finished = T;
            next;
          } else {
            #if you have already visited this dose with three patients
            if(curr_dose%in%curr_results[row_order,"dose_num"]) {
              curr_dose_index = which(curr_results[,"dose_num"]==curr_dose);
              row_order = c(row_order,curr_dose_index[4:6]);
              six_at_dose[curr_dose]=T;
              next;
            } else {
              #otherwise you havent ever visited this dose
              curr_dose_index = which(curr_results[,"dose_num"]==curr_dose)[1:3];              
              row_order = c(row_order,curr_dose_index);
              next;
            }
          }
        }
      }
    }        
    mtd_est[i] = curr_dose;
    reorder_all_rows[min(sim_index)-1+(1:length(row_order))] = min(sim_index)-1 + row_order;
  }
  all_results = all_results[reorder_all_rows[which(!is.na(reorder_all_rows))],];
  all_results[,"subj_id"] = unlist(sapply(diff(c(0,which(!duplicated(all_results[,"sim_id"],fromLast=T)))),seq,from=1,simplify=F))
  
  enrollment = as.numeric(tapply(all_results[,"subj_id"],all_results[,"sim_id"],length))
  
  return(list(all_results = all_results, mtd_est = mtd_est, enrollment = enrollment));
}


#Simulate the expansion cohort of a 3plus3
sim_3pl3_exp = function(tox_curve,sim_id,subj_id_start,dose_num,n_dec,
                        dec_size,deescalation_rule="local",
                        local_thresh = 1/3,
                        pocock_alpha = 0.0012, pocock_thresh = 0.30,
                        seed=sample(.Machine$integer.max,1),
                        cut_data_at = Inf) {
  
  set.seed(seed);
  #Check that arguments are equal length
  if(!isTRUE(all.equal(rep(length(sim_id),2),c(length(subj_id_start),length(dose_num))))) {
    stop("lengths of 'sim_id', 'subj_id_start', and 'dose_num' must be equal");
  }
  
  #First create all outcomes assuming that no de-escalation rules are in place
  all_results = matrix(0,nrow=length(sim_id)*n_dec*dec_size,ncol=7);
  colnames(all_results) = c("sim_id","subj_id","dose_num","tox_prob","tox","dec_id","first_stage");
  
  all_results[,"sim_id"] = rep(sim_id, each = n_dec * dec_size);
  all_results[,"subj_id"] = unlist(sapply(subj_id_start-1,rep,each=n_dec*dec_size,simplify=F))+unlist(sapply(diff(c(0,which(!duplicated(all_results[,"sim_id"],fromLast=T)))),seq,from=1,simplify=F))
  if(cut_data_at<dec_size) {
    all_results[,"dec_id"] = as.numeric(rbind(replicate(length(sim_id),sample(rep(1:n_dec,cut_data_at))),
                                              replicate(length(sim_id),sample(rep(1:n_dec,dec_size-cut_data_at)))));
  } else {
    all_results[,"dec_id"] = as.numeric(replicate(length(sim_id),sample(rep(1:n_dec,dec_size))));
  }
  
  all_results[,"dose_num"] = rep(dose_num, each = n_dec * dec_size);
  all_results[,"tox_prob"] = rep(tox_curve[dose_num], each = n_dec * dec_size);
  all_results[,"tox"] = rbinom(nrow(all_results),1,all_results[,"tox_prob"]);
  all_results[,"first_stage"] = 1;
  
  #Now implement de-escalation rules, which 
  #are assumed to act locally on each DEC or globally over
  #all DECs. 
  remove_rows = NULL;
  if(deescalation_rule=="local"&dec_size>3) {
    #Loop over DECs
    for(k in 1:n_dec) {
      #SimID-specific list of numbers of subjects to look over
      #Initially equal to dec_size
      list_length = rep(dec_size,length(sim_id));
      #The code uses `index' to look down each DEC 
      #and determine when dose would have been de-escalated. 
      index = which(all_results[,"dec_id"]==k);
      #Calculate cumulative toxicity rates after each patient
      tox_rates = lapply(tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),cumsum),"/",1:dec_size);
      #Only look at rates after the third patient
      tox_rates = lapply(tox_rates,"[",-(1:2))
      #If rate exceeds 'local_thresh', will de-escalate at the next patient. 
      #This equals Inf if there is never a need to de-escalate
      decrease_at = 3 + unlist(lapply(lapply(lapply(tox_rates,">",local_thresh),which),min,Inf));
      #Keep looping through until no de-escalation occurs
      while(any(decrease_at<=list_length)) {
        #Sims to look at
        sim_set = as.numeric(names(decrease_at)[which(decrease_at<=list_length)])
        #Determine whether this is the first run through the loop
        if(all(list_length==dec_size)) {
          #Start will all sim ids that need de-escalating
          index = which(all_results[,"dec_id"]==k&
                          all_results[,"sim_id"]%in%sim_set);              
        } else {
          #Otherwise prune away the sim ids that have been handled by an earlier run
          index = setdiff(index,which(!all_results[,"sim_id"]%in%sim_set));
        }
        for(s in sim_set) {
          #For a given sim id, decrease index gives the set of rows, 
          #the dose assignments of which would have been de-escalated
          #before the subjects are enrolled. It will be as if their
          #toxicity outcomes at the un-de-escalated dose were never observed
          #keep index gives the complement, the rows that are left alone
          decrease_index = intersect(index,
                                     which(all_results[,"dec_id"]==k&all_results[,"sim_id"]==s));
          keep_index = decrease_index[which(all_results[decrease_index,"subj_id"]<all_results[decrease_index,"subj_id"][decrease_at[as.character(s)]])];
          decrease_index = decrease_index[which(all_results[decrease_index,"subj_id"]>=all_results[decrease_index,"subj_id"][decrease_at[as.character(s)]])];
          #if(sum(duplicated(all_results[decrease_index,"dose_num"]))!=length(decrease_index)-1) {stop();}
          #De-escalate
          all_results[decrease_index,"dose_num"] = all_results[decrease_index,"dose_num"] - 1;
          if(all_results[decrease_index,"dose_num"][1]>0) {
            #assign new toxicity prob
            all_results[decrease_index,"tox_prob"] = tox_curve[all_results[decrease_index,"dose_num"]];
            #Determine toxicities under the de-escalated dose. 
            all_results[decrease_index,"tox"] = rbinom(length(decrease_index),1,all_results[decrease_index,"tox_prob"]);
            #Prune away the keep index rows from the global index. 
            index = setdiff(index,keep_index);
          } else {
            #Keep track of rows that would never have occurred because
            #de-escalation went down to zero
            remove_rows = c(remove_rows,decrease_index);
            #Remove this entire sim id (for this given dec id) from 
            #consideration, as it has been addressed
            index = setdiff(index,which(all_results[,"sim_id"]==s));
          }
        }
        #Recalculate the list length to reflect the patients/sims that have been addresed
        #This should be monotonically decreasing for each run of the while loop
        list_length = tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),length)
        #Recalculate the toxicity rates 
        tox_rates = mapply("/",tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),cumsum),sapply(list_length,seq,from=1,simplify=F),SIMPLIFY=F);
        decrease_at = rep(Inf,length(list_length));
        names(decrease_at) = names(list_length);
        if(any(list_length>3)) {
          eligible_to_decrease = which(list_length>3);
          tox_rates = lapply(tox_rates,"[",-(1:2));
          #recalculate when to de-escalate the next time, if applicable
          decrease_at[eligible_to_decrease] = 3 + unlist(lapply(lapply(lapply(tox_rates[eligible_to_decrease],">",local_thresh),which),min,Inf));
        }
      }
    }
  } else if(deescalation_rule=="global"&((n_dec*dec_size)>3)) {
    #If the de-escalation rule is not local to each dec,
    #assume it is global over all decs
    #SimID-specific list of numbers of subjects to look over
    #Initially equal to n_dec*dec_size
    list_length = rep(n_dec*dec_size,length(sim_id));
    #The code uses `index' to
    #determine when dose would have been de-escalated. 
    #Calculate cumulative toxicity rates after each patient
    tox_nums = tapply(all_results[,"tox"],list(all_results[,"sim_id"]),cumsum);
    #Only look at rates after the third patient
    tox_nums = lapply(tox_nums,"[",-(1:2));
    #If pocock-type boundary is crossed, de-escalate at the next patient. 
    #This equals Inf if there is never a need to de-escalate
    decrease_at = 3 + unlist(lapply(lapply(lapply(mapply("pbinom",tox_nums,lapply(list_length,"seq",from=3),MoreArgs=list(p=pocock_thresh),SIMPLIFY=F),">",1-pocock_alpha),"which"),min,Inf));
    #Keep looping through until no de-escalation occurs
    while(any(decrease_at<=list_length)) {
      #Sims to look at
      sim_set = as.numeric(names(decrease_at)[which(decrease_at<=list_length)])
      #Determine whether this is the first run through the loop
      if(all(list_length==n_dec*dec_size)) {
        #Start will all sim ids that need de-escalating
        index = which(all_results[,"sim_id"]%in%sim_set);              
      } else {
        #Otherwise prune away the sim ids that have been handled by an earlier run
        index = setdiff(index,which(!all_results[,"sim_id"]%in%sim_set));
      }
      for(s in sim_set) {
        #For a given sim id, decrease index gives the set of rows, 
        #the dose assignments of which would have been de-escalated
        #before the subjects are enrolled. It will be as if their
        #toxicity outcomes at the un-de-escalated dose were never observed
        #keep index gives the complement, the rows that are left alone
        decrease_index = intersect(index,which(all_results[,"sim_id"]==s));
        keep_index = decrease_index[which(all_results[decrease_index,"subj_id"]<all_results[decrease_index,"subj_id"][decrease_at[as.character(s)]])];
        decrease_index = decrease_index[which(all_results[decrease_index,"subj_id"]>=all_results[decrease_index,"subj_id"][decrease_at[as.character(s)]])];
        #De-escalate
        all_results[decrease_index,"dose_num"] = all_results[decrease_index,"dose_num"] - 1;
        if(all_results[decrease_index,"dose_num"][1]>0) {
          #assign new toxicity prob
          all_results[decrease_index,"tox_prob"] = tox_curve[all_results[decrease_index,"dose_num"]];
          #Determine toxicities under the de-escalated dose. 
          all_results[decrease_index,"tox"] = rbinom(length(decrease_index),1,all_results[decrease_index,"tox_prob"]);
          #Prune away the keep index rows from the global index. 
          index = setdiff(index,keep_index);
        } else {
          #Keep track of rows that would never have occurred because
          #de-escalation went down to zero
          remove_rows = c(remove_rows,decrease_index);
          #Remove this entire sim id (for this given dec id) from 
          #consideration, as it has been addressed
          index = setdiff(index,which(all_results[,"sim_id"]==s));
        }
      }
      #Recalculate the list length to reflect the patients/sims that have been addresed
      #This should be monotonically decreasing for each run of the while loop
      list_length = tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),length)
      #Recalculate the toxicity rates 
      tox_nums = tapply(all_results[index,"tox"],list(all_results[index,"sim_id"]),cumsum);
      decrease_at = rep(Inf,length(list_length));
      names(decrease_at) = names(list_length);
      if(any(list_length>3)) {
        eligible_to_decrease = which(list_length>3);
        tox_nums = lapply(tox_nums,"[",-(1:2));
        #recalculate when to de-escalate the next time, if applicable
        decrease_at[eligible_to_decrease] = 3 + unlist(lapply(lapply(lapply(mapply("pbinom",tox_nums[eligible_to_decrease],lapply(list_length[eligible_to_decrease],"seq",from=3),MoreArgs=list(p=pocock_thresh),SIMPLIFY=F),">",1-pocock_alpha),"which"),min,Inf));
      }
    }
  }
  #Remove those rows that would never have been enrolled because the
  #DEC(s) stopped for toxicity
  if(length(remove_rows)) {
    all_results = all_results[-remove_rows,,drop=F];
  }
  #Determine whether this is 'first stage' data,
  #meaning data collected before the first interim 
  #efficacy analysis
  if(cut_data_at<dec_size) {
    all_results[which(all_results[,"subj_id"]>unlist(mapply("rep",x=subj_id_start+n_dec*cut_data_at-1,each=tapply(all_results[,"sim_id"],all_results[,"sim_id"],length),SIMPLIFY=F))),"first_stage"] = 0;
  } 
  #Renumber the subjects in 'all_results' to reflect that the patients
  #in the `remove rows' set would never have been enrolled in the first place
  #Subtle note: it is crucial that the 'first_stage' indicator
  #be set before the subject ids are relabeled: 
  skip_id = 1 + which(diff(all_results[,"subj_id"])>1);
  while(length(skip_id>0)) {
    all_results[skip_id,"subj_id"] = all_results[skip_id,"subj_id"] - 1;
    skip_id = 1 + which(diff(all_results[,"subj_id"])>1);
  }
  
  list(#enrollment =  1+tapply(all_results[,"subj_id"],list(all_results[,"sim_id"]),max)-tapply(all_results[,"subj_id"],list(all_results[,"sim_id"]),min),
       enrollment =  tapply(all_results[,"subj_id"],list(all_results[,"sim_id"],all_results[,"dec_id"]),length),
       all_results = all_results[order(all_results[,"sim_id"],all_results[,"subj_id"]),]);
}


sim_dec_crm = function(tox_curve, start = 2, niter, 
                   skeleton, n_in_esc, n_dec, dec_size,
                   cut_data_at = Inf,
                   seed = sample(.Machine$integer.max,1), 
                   restrict = T, 
                   tox_target, 
                   scale = 0.6, 
                   method = "bayes",
                   no.exceed = 0.05, 
                   cohort.size = 3, 
                   first.cohort.only = T, 
                   n.at.MTD = Inf) {
  set.seed(seed);
  
  n_in_exp = n_dec*dec_size;
  all_results = matrix(0,nrow=(ceiling(n_in_esc)+n_in_exp)*niter,ncol=7);
  colnames(all_results) = c("sim_id","subj_id","dose_num","tox_prob","tox","dec_id","first_stage");
  beta_ests = beta_scales = numeric(niter);
  index = 0;
  all_n = n_in_exp + floor(n_in_esc) + rbinom(niter,1,n_in_esc-floor(n_in_esc));
  for(k in 1:niter) {
    curr_n = all_n[k];
    curr_seed = sample(.Machine$integer.max,1);
    big_crm = titesim_phil(PI = tox_curve,
                           prior = skeleton,  
                           target = tox_target,
                           n = curr_n,
                           x0=start,nsim=1,seed = curr_seed,count=F,
                           conf.level = 0.5,
                           method=method,model="empiric",scale=scale,
                           no.exceed=no.exceed,cohort.size=cohort.size,
                           first.cohort.only=first.cohort.only,n.at.MTD=n.at.MTD,
                           obswin=1,rate=1,accrual="fixed")
    
    if(big_crm$last_sim$stop.for.tox>0&big_crm$last_sim$stop.for.tox<=curr_n-n_in_exp) {
      beta_ests[k] = -Inf;
      beta_scales[k] = Inf;
    } else {
      new_dat = cbind(big_crm$last_sim$level,(big_crm$last_sim$PI)[big_crm$last_sim$level],big_crm$last_sim$tox)[-(1:(curr_n-n_in_exp)),,drop=F];
      colnames(new_dat) = c("level","dose","tox");
      stop_at = nrow(new_dat);
      index = max(index) + (1:stop_at);
      if(cut_data_at<dec_size) {
        dec_assignments = c(sample(rep(1:n_dec,cut_data_at)),sample(rep(1:n_dec,dec_size-cut_data_at)));
        first_stage = c(rep(1,cut_data_at*n_dec),rep(0,(dec_size-cut_data_at)*n_dec))[1:stop_at];
      } else {
        dec_assignments = sample(rep(1:n_dec,dec_size));
        first_stage = 1;
      }
      dec_assignments = dec_assignments[1:stop_at];
      all_results[index,] = as.matrix(cbind(k,cbind(curr_n-n_in_exp + (1:stop_at),new_dat),dec_assignments,first_stage)); 
      if(any(c(diff(all_results[index,"dose_num"]),0)==-1&all_results[index,"tox"]==0)) {
        cat(paste0("De-escalation following no DLT in CRM with seed ", curr_seed,"\n")); 
      }
      beta_ests[k] = big_crm$last_sim$final.est;
      beta_scales[k] = sqrt(big_crm$last_sim$post.var);
    }
    rm(big_crm);
    if(k%%100==0) {cat(k,"\n");}
  }    
  all_results = all_results[1:max(index),];
  
  list(all_results = all_results,
       beta_ests=beta_ests,beta_scales=beta_scales,
       seed=seed);
  
}

