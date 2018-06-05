#################
#June 10, 2015
#Phil Boonstra
#Code to analyze the results of running 'DECEfficacy.R'
################

rm(list=ls(all=TRUE));
require(dfcrm);require(rstan);
require(lattice);require(RColorBrewer);
my.work.computer=T;
if(my.work.computer){
  setwd("/Users/philb/Desktop/Work/ExpCohortStudy/Efficacy");
  titesim_folder = "/users/philb/Desktop/Work/ExpCohortStudy/Efficacy/github";
  write_to_folder = "/Users/philb/Desktop/Work/ExpCohortStudy/Efficacy";
} 

#sim_set should be any subset of 1:60 (assuming the default settings were used in genParams.R) and corresponds to the simulation labels constructed by genParams.R
#scenario_set should be any subset of 1:10 and corresponds to the scenario labels used in the paper

#sim_set is constructed according to the 'arglist' variable constructed by genParams
#1:18 = scenarios 1-6 from the main AOC paper (3 DEC types for each scenario)
#19:36 = scenarios 1-6 from the main AOC paper using a different skeleton for the crm
#37:48 = scenarios 1-6 from the supplement to the AOC paper
#49:60 = scenarios 1-6 from the supplement to the AOC paper using a different skeleton for the crm
sim_set = c(1:18,37:48);#load the first skeleton sets of skims
file.name = "DEC_eff";

#################################################################
#################################################################
scenario_set = 1:6;#choose this (and comment out the below line) for the figures in the main text 
#scenario_set = 7:10;#choose this (and comment out the above line) for the figures in the supplement 


nonnull_category_codes = ((1:8)*100+1);
nonnull_category_codes_short = c(1,3,6,7,8)*100+1;
onestage_nonnull_category_labels = c("StopToxEsc","StopTox(15)","StopEff(15)","        ","       ","RecBelow","RecOK","RecAbove")
twostage_nonnull_category_labels = c("StopToxEsc","StopTox(15)","StopEff(15)","StopTox(30)","StopEff(30)","RecBelow","RecOK","RecAbove")
twostage_nointer_nonnull_category_labels = c("StopToxEsc","StopTox(15)","        ","StopTox(30)","StopEff(30)","RecBelow","RecOK","RecAbove")
########
onestage_nonnull_category_labels_short = 
  twostage_nonnull_category_labels_short = 
  twostage_nointer_nonnull_category_labels_short = c("No Recom.\n(Toxicity)","No Recom.\n(Futility)","Recom. Sub-\ntherapeutic","Recom.\nAcceptable","Recom.\nToxic")

null_category_labels = c("sens","spec");
probability_distance = probability_bias = 
  TPR = FPR = 
  nPerNullDEC_enrolled = 
  nPerNonNullDEC_enrolled = data.frame();
onestage_nonnull_props = 
  twostage_nonnull_props = 
  twostage_nointer_nonnull_props =
  onestage_nonnull_props_short = 
  twostage_nonnull_props_short = 
  twostage_nointer_nonnull_props_short = data.frame(numeric(0),character(0),character(0),character(0),numeric(0),stringsAsFactors = F);
colnames(onestage_nonnull_props) = 
  colnames(twostage_nonnull_props) = 
  colnames(twostage_nointer_nonnull_props) =
  colnames(onestage_nonnull_props_short) = 
  colnames(twostage_nonnull_props_short) = 
  colnames(twostage_nointer_nonnull_props_short) = c("setting","design","analysis","code","prop");
#null_props = data.frame(numeric(0),character(0),character(0),character(0),numeric(0),stringsAsFactors = F);
#colnames(null_props) = c("setting","design","analysis","type","prop");

for(curr_sim in sim_set) {
  load(file=paste(file.name,curr_sim,".RData",sep=""));
}

null_counter = 1;
onestage_nonnull_counter = twostage_nonnull_counter = 
  twostage_nointer_nonnull_counter = 1;
onestage_nonnull_counter_short = twostage_nonnull_counter_short = 
  twostage_nointer_nonnull_counter_short = 1;

for(curr_sim in sim_set) {
  attach(eval(as.name(paste("sim",curr_sim,sep=""))));
  ndose = length(params$tox_curve);
  
  if(params$null_scenario) {
    for(curr_design in c("loc","glob","crm")) {
      curr_eff_est = get(paste0("eff_est_",curr_design));
      for(curr_analysis in c("mod","emp")) {
        if(!params$two_stage_dec) {
          curr_props = numeric(length(nonnull_category_codes));
          names(curr_props) = nonnull_category_codes;
          #Assumes first three decs are non-null
          foo = summary(factor(curr_eff_est[rep(4:5,times=params$niter) + 5*rep(0:(params$niter-1),each=2),paste0("rp2dCode_",curr_analysis)]))/(2*params$niter);
          curr_props[names(foo)] = foo;
          onestage_nonnull_props[onestage_nonnull_counter+(0:7),] = data.frame(params$setting,curr_design,curr_analysis,onestage_nonnull_category_labels,curr_props,stringsAsFactors = F,row.names=NULL);
          onestage_nonnull_counter = onestage_nonnull_counter + 8;  
          #Combined codes
          curr_props_short = numeric(length(nonnull_category_codes_short));
          names(curr_props_short) = nonnull_category_codes_short;
          foo["101"] = sum(foo[names(foo)%in%c("101","201","401")]);
          foo["301"] = sum(foo[names(foo)%in%c("301","501")]);
          foo = foo[names(foo)%in%nonnull_category_codes_short];
          curr_props_short[names(foo)] = foo;
          onestage_nonnull_props_short[onestage_nonnull_counter_short+(0:4),] = data.frame(params$setting,curr_design,curr_analysis,onestage_nonnull_category_labels_short,curr_props_short,stringsAsFactors = F,row.names=NULL);
          onestage_nonnull_counter_short = onestage_nonnull_counter_short + 5;
          ###
        } else {
          if(params$skip_first_stage) {
            curr_props = numeric(length(nonnull_category_codes));
            names(curr_props) = nonnull_category_codes;
            foo = summary(factor(curr_eff_est[rep(4:5,times=params$niter) + 5*rep(0:(params$niter-1),each=2),paste0("rp2dCode_",curr_analysis)]))/(2*params$niter);
            curr_props[names(foo)] = foo;
            twostage_nointer_nonnull_props[twostage_nointer_nonnull_counter+(0:7),] = data.frame(params$setting,curr_design,curr_analysis,twostage_nointer_nonnull_category_labels,curr_props,stringsAsFactors = F,row.names=NULL);
            twostage_nointer_nonnull_counter = twostage_nointer_nonnull_counter + 8;  
            #Combined codes
            curr_props_short = numeric(length(nonnull_category_codes_short));
            names(curr_props_short) = nonnull_category_codes_short;
            foo["101"] = sum(foo[names(foo)%in%c("101","201","401")]);
            names(foo)[names(foo)=="501"] = "301";
            foo = foo[names(foo)%in%nonnull_category_codes_short];
            curr_props_short[names(foo)] = foo;
            twostage_nointer_nonnull_props_short[twostage_nointer_nonnull_counter_short+(0:4),] = data.frame(params$setting,curr_design,curr_analysis,twostage_nointer_nonnull_category_labels_short,curr_props_short,stringsAsFactors = F,row.names=NULL);
            twostage_nointer_nonnull_counter_short = twostage_nointer_nonnull_counter_short + 5;  
            ###
          } else {
            curr_props = numeric(length(nonnull_category_codes));
            names(curr_props) = nonnull_category_codes;
            foo = summary(factor(curr_eff_est[rep(4:5,times=params$niter) + 5*rep(0:(params$niter-1),each=2),paste0("rp2dCode_",curr_analysis)]))/(2*params$niter);
            curr_props[names(foo)] = foo;
            twostage_nonnull_props[twostage_nonnull_counter+(0:7),] = data.frame(params$setting,curr_design,curr_analysis,twostage_nonnull_category_labels,curr_props,stringsAsFactors = F,row.names=NULL);
            twostage_nonnull_counter = twostage_nonnull_counter + 8;  
            #Combined codes
            curr_props_short = numeric(length(nonnull_category_codes_short));
            names(curr_props_short) = nonnull_category_codes_short;
            foo["101"] = sum(foo[names(foo)%in%c("101","201","401")]);
            foo["301"] = sum(foo[names(foo)%in%c("301","501")]);
            foo = foo[names(foo)%in%nonnull_category_codes_short];
            curr_props_short[names(foo)] = foo;
            twostage_nonnull_props_short[twostage_nonnull_counter_short+(0:4),] = data.frame(params$setting,curr_design,curr_analysis,twostage_nonnull_category_labels_short,curr_props_short,stringsAsFactors = F,row.names=NULL);
            twostage_nonnull_counter_short = twostage_nonnull_counter_short + 5;  
            ###
          }
        }
      }
    }    
    curr_sens = c(sum(eff_est_crm$rp2dCode_emp>=601&eff_curves$pref_dose>0)/sum(eff_curves$pref_dose>0),
                  sum(eff_est_crm$rp2dCode_mod>=601&eff_curves$pref_dose>0)/sum(eff_curves$pref_dose>0),
                  sum(eff_est_glob$rp2dCode_emp>=601&eff_curves$pref_dose>0)/sum(eff_curves$pref_dose>0),
                  sum(eff_est_glob$rp2dCode_mod>=601&eff_curves$pref_dose>0)/sum(eff_curves$pref_dose>0),
                  sum(eff_est_loc$rp2dCode_emp>=601&eff_curves$pref_dose>0)/sum(eff_curves$pref_dose>0),
                  sum(eff_est_loc$rp2dCode_mod>=601&eff_curves$pref_dose>0)/sum(eff_curves$pref_dose>0))
    TPR = rbind(TPR,c(params$setting,15+15*params$two_stage_dec,params$skip_first_stage,curr_sens));
    curr_spec = c(sum(eff_est_crm$rp2dCode_emp<601&eff_curves$pref_dose==0)/sum(eff_curves$pref_dose==0),
                  sum(eff_est_crm$rp2dCode_mod<601&eff_curves$pref_dose==0)/sum(eff_curves$pref_dose==0),
                  sum(eff_est_glob$rp2dCode_emp<601&eff_curves$pref_dose==0)/sum(eff_curves$pref_dose==0),
                  sum(eff_est_glob$rp2dCode_mod<601&eff_curves$pref_dose==0)/sum(eff_curves$pref_dose==0),
                  sum(eff_est_loc$rp2dCode_emp<601&eff_curves$pref_dose==0)/sum(eff_curves$pref_dose==0),
                  sum(eff_est_loc$rp2dCode_mod<601&eff_curves$pref_dose==0)/sum(eff_curves$pref_dose==0))
    FPR = rbind(FPR,c(params$setting,15+15*params$two_stage_dec,params$skip_first_stage,1-curr_spec));
    curr_nPerNullDEC_enrolled = c(sum(n_enroll_crm_emp*t(null_decs[,as.numeric(rownames(n_enroll_crm_emp))]),na.rm=T)/sum(null_decs[,as.numeric(rownames(n_enroll_crm_emp))]),
                                  sum(n_enroll_crm_mod*t(null_decs[,as.numeric(rownames(n_enroll_crm_mod))]),na.rm=T)/sum(null_decs[,as.numeric(rownames(n_enroll_crm_mod))]),
                                  sum(n_enroll_glob_emp*t(null_decs[,as.numeric(rownames(n_enroll_glob_emp))]),na.rm=T)/sum(null_decs[,as.numeric(rownames(n_enroll_glob_emp))]),
                                  sum(n_enroll_glob_mod*t(null_decs[,as.numeric(rownames(n_enroll_glob_mod))]),na.rm=T)/sum(null_decs[,as.numeric(rownames(n_enroll_glob_mod))]),
                                  sum(n_enroll_loc_emp*t(null_decs[,as.numeric(rownames(n_enroll_loc_emp))]),na.rm=T)/sum(null_decs[,as.numeric(rownames(n_enroll_loc_emp))]),
                                  sum(n_enroll_loc_mod*t(null_decs[,as.numeric(rownames(n_enroll_loc_mod))]),na.rm=T)/sum(null_decs[,as.numeric(rownames(n_enroll_loc_mod))]));
    nPerNullDEC_enrolled = rbind(nPerNullDEC_enrolled,c(params$setting,15+15*params$two_stage_dec,params$skip_first_stage,curr_nPerNullDEC_enrolled));
    curr_nPerNonNullDEC_enrolled = c(sum(n_enroll_crm_emp*t(1-null_decs[,as.numeric(rownames(n_enroll_crm_emp))]),na.rm=T)/sum(1-null_decs[,as.numeric(rownames(n_enroll_crm_emp))]),
                                     sum(n_enroll_crm_mod*t(1-null_decs[,as.numeric(rownames(n_enroll_crm_mod))]),na.rm=T)/sum(1-null_decs[,as.numeric(rownames(n_enroll_crm_mod))]),
                                     sum(n_enroll_glob_emp*t(1-null_decs[,as.numeric(rownames(n_enroll_glob_emp))]),na.rm=T)/sum(1-null_decs[,as.numeric(rownames(n_enroll_glob_emp))]),
                                     sum(n_enroll_glob_mod*t(1-null_decs[,as.numeric(rownames(n_enroll_glob_mod))]),na.rm=T)/sum(1-null_decs[,as.numeric(rownames(n_enroll_glob_mod))]),
                                     sum(n_enroll_loc_emp*t(1-null_decs[,as.numeric(rownames(n_enroll_loc_emp))]),na.rm=T)/sum(1-null_decs[,as.numeric(rownames(n_enroll_loc_emp))]),
                                     sum(n_enroll_loc_mod*t(1-null_decs[,as.numeric(rownames(n_enroll_loc_mod))]),na.rm=T)/sum(1-null_decs[,as.numeric(rownames(n_enroll_loc_mod))]));
    nPerNonNullDEC_enrolled = rbind(nPerNonNullDEC_enrolled,c(params$setting,15+15*params$two_stage_dec,params$skip_first_stage,curr_nPerNonNullDEC_enrolled));
    
    true_mtd = which.min(abs(params$tox_curve-params$tox_target));
    curr_prob_dist = c(mean(abs(eff_est_crm$effProbAtMTD_emp-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean(abs(eff_est_crm$effProbAtMTD_mod-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean(abs(eff_est_glob$effProbAtMTD_emp-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean(abs(eff_est_glob$effProbAtMTD_mod-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean(abs(eff_est_loc$effProbAtMTD_emp-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean(abs(eff_est_loc$effProbAtMTD_mod-eff_curves[,paste0("dose",true_mtd)]),na.rm=T))
    probability_distance = rbind(probability_distance,c(params$setting,15+15*params$two_stage_dec,params$skip_first_stage,params$null_scenario,curr_prob_dist));
    curr_prob_bias = c(mean((eff_est_crm$effProbAtMTD_emp-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean((eff_est_crm$effProbAtMTD_mod-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean((eff_est_glob$effProbAtMTD_emp-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean((eff_est_glob$effProbAtMTD_mod-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean((eff_est_loc$effProbAtMTD_emp-eff_curves[,paste0("dose",true_mtd)]),na.rm=T),
                       mean((eff_est_loc$effProbAtMTD_mod-eff_curves[,paste0("dose",true_mtd)]),na.rm=T))
    probability_bias = rbind(probability_bias,c(params$setting,15+15*params$two_stage_dec,params$skip_first_stage,params$null_scenario,curr_prob_bias));
  }
  detach();
}

colnames(probability_distance) = colnames(probability_bias) = 
  c("setting","dec_size","skip_first_stage","null_scenario","crm_emp","crm_mod","glob_emp","glob_mod","loc_emp","loc_mod");

colnames(TPR) = 
  colnames(FPR) = 
  colnames(nPerNullDEC_enrolled) = 
  colnames(nPerNonNullDEC_enrolled) = 
  c("setting","dec_size","skip_first_stage","crm_emp","crm_mod","glob_emp","glob_mod","loc_emp","loc_mod");

onestage_colors = twostage_colors = twostage_nointer_colors = 
  c(brewer.pal(5,"Reds"),brewer.pal(3,"Paired"));
onestage_colors_short = twostage_colors_short = twostage_nointer_colors_short = 
  brewer.pal(5,"Blues");
#  rev(grey.colors(5,start=0,end=0.95));
#c("#F1F1F1",rev(grey.colors(4,start=0,end=0.7)));
onestage_colors[(4:5)] = "#FFFFFF";
twostage_nointer_colors[3] = "#FFFFFF";


design_names = c("LOCAL","GLOBAL","CRM");
names(design_names) = c("loc","glob","crm");
analysis_names = c("MODEL","EMPIRIC");
names(analysis_names) = c("mod","emp");
onestage_nonnull_props$x = paste(design_names[onestage_nonnull_props$design],analysis_names[onestage_nonnull_props$analysis],sep="/\n");
onestage_nonnull_props$x_hutch = design_names[onestage_nonnull_props$design];
onestage_nonnull_props_short$x = paste(design_names[onestage_nonnull_props_short$design],analysis_names[onestage_nonnull_props_short$analysis],sep="/\n");
onestage_nonnull_props_short$x_hutch = design_names[onestage_nonnull_props_short$design];
twostage_nonnull_props$x = paste(design_names[twostage_nonnull_props$design],analysis_names[twostage_nonnull_props$analysis],sep="/\n");
twostage_nonnull_props$x_hutch = design_names[twostage_nonnull_props$design];
twostage_nonnull_props_short$x = paste(design_names[twostage_nonnull_props_short$design],analysis_names[twostage_nonnull_props_short$analysis],sep="/\n");
twostage_nonnull_props_short$x_hutch = design_names[twostage_nonnull_props_short$design];
twostage_nointer_nonnull_props$x = paste(design_names[twostage_nointer_nonnull_props$design],analysis_names[twostage_nointer_nonnull_props$analysis],sep="/\n");
twostage_nointer_nonnull_props$x_hutch = design_names[twostage_nointer_nonnull_props$design];
twostage_nointer_nonnull_props_short$x = paste(design_names[twostage_nointer_nonnull_props_short$design],analysis_names[twostage_nointer_nonnull_props_short$analysis],sep="/\n");
twostage_nointer_nonnull_props_short$x_hutch = design_names[twostage_nointer_nonnull_props_short$design];
###

onestage_nonnull_props$code = factor(onestage_nonnull_props$code,levels=unique(onestage_nonnull_props$code),ordered=T);
onestage_nonnull_props_short$code = factor(onestage_nonnull_props_short$code,levels=unique(onestage_nonnull_props_short$code),ordered=T);
twostage_nonnull_props$code = factor(twostage_nonnull_props$code,levels=unique(twostage_nonnull_props$code),ordered=T);
twostage_nonnull_props_short$code = factor(twostage_nonnull_props_short$code,levels=unique(twostage_nonnull_props_short$code),ordered=T);
twostage_nointer_nonnull_props$code = factor(twostage_nointer_nonnull_props$code,levels=unique(twostage_nointer_nonnull_props$code),ordered=T);
twostage_nointer_nonnull_props_short$code = factor(twostage_nointer_nonnull_props_short$code,levels=unique(twostage_nointer_nonnull_props_short$code),ordered=T);

#Normalized specificify
norm_FPR = FPR;
norm_FPR[,-(1:3)] = matrix(0.25,nrow=18,ncol=6);

norm_TPR = TPR;
norm_TPR[,-(1:3)] = ((1-((1-norm_FPR)/(1-FPR)*(1-TPR)))*(FPR<=norm_FPR))[,-(1:3)]+
  (TPR*norm_FPR/FPR*(FPR>norm_FPR))[,-(1:3)];



my.settings <- list(
  superpose.polygon=list(col=onestage_colors, border="transparent"),
  strip.background=list(col="#AAAAAA"),
  strip.border=list(col="black"),
  panel.background = list(col="#E5E5E5"),
  layout.heights = list(top.padding = 0,
                        main.key.padding = 0,
                        key.axis.padding = 0,
                        axis.xlab.padding = 0,
                        xlab.key.padding = 0,
                        key.sub.padding = 0,
                        bottom.padding = 0),
  layout.widths =
    list(left.padding = 0,
         key.ylab.padding = 0,
         ylab.axis.padding = 0,
         axis.key.padding = 0,
         right.padding = 0),
  add.line = list(col.line="#909090")
)
#9D9D9D

#################################################################
#################################################################
#Figure 3 in paper
#Top part of Figure 3 in paper
my.settings$superpose.polygon$col = onestage_colors_short;
#trellis.device("pdf",file="recommendations_dec15_mod.pdf",height=5,width=7,family="serif");
trellis.device("jpeg",file="recommendations_dec15_mod.jpg",height=5,width=7,family="serif",units="in",res=500);
barchart(prop~x|factor(setting),groups=code,
         data=onestage_nonnull_props_short[onestage_nonnull_props_short$analysis=="mod" & onestage_nonnull_props_short$setting %in% scenario_set,],
         main = "# in DEC = 15",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=5, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 1.75),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();

#Bottom part of Figure 3 in paper
my.settings$superpose.polygon$col = twostage_colors_short;
#trellis.device("pdf",file="recommendations_dec30_mod.pdf",height=5,width=7,family="serif");
trellis.device("jpeg",file="recommendations_dec30_mod.jpg",height=5,width=7,family="serif",units="in",res=500);
barchart(prop~x|factor(setting),groups=code,
         data=twostage_nonnull_props_short[twostage_nonnull_props_short$analysis=="mod"& twostage_nonnull_props_short$setting %in% scenario_set,],
         main = "# in DEC = 30 (With Futility Analysis at 15)",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=5, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 1.75),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();

#Figure 3 in paper (Black and white version)
#Top part of Figure 3 in paper
my.settings$superpose.polygon$col = rev(grey.colors(5,start=0,end=0.8));
#trellis.device("pdf",file="recommendations_dec15_mod.pdf",height=5,width=7,family="serif");
trellis.device("jpeg",file="recommendations_dec15_mod_bw.jpg",height=5,width=7,family="serif",units="in",res=600);
barchart(prop~x|factor(setting),groups=code,
         data=onestage_nonnull_props_short[onestage_nonnull_props_short$analysis=="mod" & onestage_nonnull_props_short$setting %in% scenario_set,],
         main = "# in DEC = 15",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=5, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 1.75),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();

#Bottom part of Figure 3 in paper (Black and white version)
my.settings$superpose.polygon$col = rev(grey.colors(5,start=0,end=0.8));
#trellis.device("pdf",file="recommendations_dec30_mod.pdf",height=5,width=7,family="serif");
trellis.device("jpeg",file="recommendations_dec30_mod_bw.jpg",height=5,width=7,family="serif",units="in",res=600);
barchart(prop~x|factor(setting),groups=code,
         data=twostage_nonnull_props_short[twostage_nonnull_props_short$analysis=="mod"& twostage_nonnull_props_short$setting %in% scenario_set,],
         main = "# in DEC = 30 (With Futility Analysis at 15)",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=5, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 1.75),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();

#################################################################
#################################################################


#################################################################
#################################################################
#Figure SX in paper supplement
#Top part of Figure SX in paper
my.settings$superpose.polygon$col = onestage_colors_short;
trellis.device("pdf",file="recommendations_dec15_emp.pdf",height=5,width=7,family="serif");
barchart(prop~x|factor(setting),groups=code,
         data=onestage_nonnull_props_short[onestage_nonnull_props_short$analysis=="emp" & onestage_nonnull_props_short$setting %in% scenario_set,],
         main = "# in DEC = 15",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=5, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 1.75),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();

#Bottom part of Figure SX in paper
my.settings$superpose.polygon$col = twostage_colors_short;
trellis.device("pdf",file="recommendations_dec30_emp.pdf",height=5,width=7,family="serif");
barchart(prop~x|factor(setting),groups=code,
         data=twostage_nonnull_props_short[twostage_nonnull_props_short$analysis=="emp"& twostage_nonnull_props_short$setting %in% scenario_set,],
         main = "# in DEC = 30 (With Futility Analysis at 15)",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=5, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 1.75),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();
#################################################################
#################################################################


#################################################################
#################################################################
#Auxiliary figures not used
my.settings$superpose.polygon$col = onestage_colors;
trellis.device("pdf",file="recommendations_dec15.pdf",height=5,width=7,family="serif");
barchart(prop~x|factor(setting),groups=code,
         data=onestage_nonnull_props[onestage_nonnull_props$setting %in% scenario_set,],
         main = "# in DEC = 15",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75)),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=4, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 0.1),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0); 
           panel.barchart(x,y,subscripts=subscripts,...)
           panel.text(x=x[subscripts%%8==7],y=tapply(y,rep(1:6,each=8),sum)-y[subscripts%%8==0]-0.5*y[subscripts%%8==7],labels=formatC(100*y[subscripts%%8==7],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();

my.settings$superpose.polygon$col = twostage_colors;
trellis.device("pdf",file="recommendations_dec30.pdf",height=5,width=7,family="serif");
barchart(prop~x|factor(setting),groups=code,
         data=twostage_nonnull_props[twostage_nonnull_props$setting %in% scenario_set,],
         main = "# in DEC = 30",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75)),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=4, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 0.1),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%8==7],y=tapply(y,rep(1:6,each=8),sum)-y[subscripts%%8==0]-0.5*y[subscripts%%8==7],labels=formatC(100*y[subscripts%%8==7],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();

my.settings$superpose.polygon$col = twostage_nointer_colors;
trellis.device("pdf",file="recommendations_nointer_dec30.pdf",height=5,width=7,family="serif");
barchart(prop~x|factor(setting),groups=code,
         data=twostage_nointer_nonnull_props[twostage_nointer_nonnull_props$setting %in% scenario_set,],
         main = "# in DEC = 30 (No Efficacy Analysis at 15)",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75)),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=4, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 0.1),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%8==7],y=tapply(y,rep(1:6,each=8),sum)-y[subscripts%%8==0]-0.5*y[subscripts%%8==7],labels=formatC(100*y[subscripts%%8==7],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();


my.settings$superpose.polygon$col = twostage_nointer_colors_short;
trellis.device("pdf",file="recommendations_nointer_dec30_mod.pdf",height=5,width=7,family="serif");
barchart(prop~x|factor(setting),groups=code,
         data=twostage_nointer_nonnull_props_short[twostage_nonnull_props_short$analysis=="mod"&twostage_nonnull_props_short$setting %in% scenario_set,],
         main = "# in DEC = 30 (No Efficacy Analysis at 15)",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,ceiling(length(scenario_set)/2)),as.table=T,ylab="",
         auto.key=list(space="top", columns=5, 
                       points=FALSE, rectangles=TRUE,
                       title="", cex.title=1,
                       between = 0, between.columns = 1.75),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"));
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="Scenario",strip.names=T,sep=" ",...);
         });
dev.off();
#################################################################
#################################################################


#################################################################
######################
##########Plots for 10/26/15 Talk at Hutch
hutch1_data = rbind(cbind(twostage_nonnull_props_short[twostage_nonnull_props_short$analysis=="mod"&twostage_nonnull_props_short$setting==4,],futility="With Interim Futility Analysis"),
                    cbind(twostage_nointer_nonnull_props_short[twostage_nointer_nonnull_props_short$analysis=="mod"&twostage_nointer_nonnull_props_short$setting==4,],futility="No Interim Futility Analysis"))
hutch2_data = rbind(cbind(twostage_nonnull_props_short[twostage_nonnull_props_short$analysis=="mod"&twostage_nonnull_props_short$setting==5,],futility="With Interim Futility Analysis"),
                    cbind(twostage_nointer_nonnull_props_short[twostage_nointer_nonnull_props_short$analysis=="mod"&twostage_nointer_nonnull_props_short$setting==5,],futility="No Interim Futility Analysis"))
hutch1_data$code = sub("Recom.\nToxic","\n",hutch1_data$code);
hutch1_data$code = factor(hutch1_data$code,levels=unique(hutch1_data$code),ordered=T)

my.settings$superpose.polygon$col = twostage_colors_short;
my.settings$superpose.polygon$col[1:2] = c("#FC9272","#DE2D26");
my.settings$superpose.polygon$col[5] = "#FFFFFF";
my.settings$add.text = trellis.par.get()$add.text;
my.settings$add.text$cex = 0.8;
trellis.device("pdf",file="hutch1a.pdf",height=10/4,width=21/4,family="serif");
barchart(prop~x_hutch|factor(futility),groups=code,
         data=hutch1_data[grep("With",hutch1_data$futility),],
         main = "",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,1),as.table=T,ylab="Proportion",
         auto.key=F,
         #auto.key=list(space="top", columns=5, reverse=F,
         #             points=FALSE, rectangles=TRUE,
         #             title="", cex.title=1,cex=0.77,
         #             between = 0, between.columns = 0.5),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"),cex=1.2);
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="",strip.names=T,sep=" ",...);
         }
);
dev.off();

trellis.device("pdf",file="hutch1b.pdf",height=10/4,width=21/4,family="serif");
barchart(prop~x_hutch|factor(futility),groups=code,
         data=hutch1_data,
         main = "",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,1),as.table=T,ylab="Proportion",
         auto.key=F,
         #auto.key=list(space="top", columns=5, reverse=F,
         #               points=FALSE, rectangles=TRUE,
         #              title="", cex.title=1,cex=0.77,
         #             between = 0, between.columns = 0.5),
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"),cex=1.2);
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="",strip.names=T,sep=" ",...);
         }
);
dev.off();

my.settings$superpose.polygon$col = twostage_colors_short;
my.settings$superpose.polygon$col[1:2] = c("#FC9272","#DE2D26");
trellis.device("pdf",file="hutch2a.pdf",height=10/4,width=21/4,family="serif");
barchart(prop~x_hutch|factor(futility),groups=code,
         data=hutch2_data[grep("With",hutch2_data$futility),],
         main = "",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,1),as.table=T,ylab="Proportion",
         auto.key=F,
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"),cex=1.2);
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="",strip.names=T,sep=" ",...);
         }
);
dev.off();

trellis.device("pdf",file="hutch2b.pdf",height=10/4,width=21/4,family="serif");
barchart(prop~x_hutch|factor(futility),groups=code,
         data=hutch2_data,
         main = "",
         stack = T,
         par.settings=my.settings,ylim=c(-0.001,1),
         scales = list(x=list(cex=0.75),y=list(at=c(0.2,0.4,0.6,0.8),labels=c(20,40,60,80))),
         layout=c(2,1),as.table=T,ylab="Proportion",
         auto.key=F,
         panel=function(x,y,subscripts,...){
           panel.grid(h=-1, v=0, col=trellis.par.get()$add.line$col.line); 
           panel.barchart(x,y,subscripts=subscripts,...);
           panel.text(x=x[subscripts%%5==4],y=tapply(y,rep(1:3,each=5),sum)-y[subscripts%%5==0]-0.5*y[subscripts%%5==4],labels=formatC(100*y[subscripts%%5==4],digits=0,format="f"),cex=1.2);
         }, 
         strip = function(var.name,...) {
           strip.default(var.name="",strip.names=T,sep=" ",...);
         }
);
dev.off();
#################################################################
#################################################################


#################################################################
#################################################################
###TPR/FPR
row1_name = "Model";
row2_name = "Empiric";
file_counter = 0;

for(size_skip in list(c(15,0),c(30,1),c(30,0))) {
  for(setting_start in c(1,3,5)) {
    file_counter = file_counter + 1;
    curr_index = which(TPR$setting%in%c(setting_start,setting_start+1)&TPR$dec_size==size_skip[1]&TPR$skip_first_stage==size_skip[2]);
    
    #TPR/FPR
    file_name = paste0(write_to_folder,"/sensspec",LETTERS[file_counter],".txt");
    write.table("",file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\n",append=F)
    row1a = formatC(100*as.matrix(cbind(TPR[curr_index[1],c("crm_mod","glob_mod","loc_mod")],
                                        TPR[curr_index[2],c("crm_mod","glob_mod","loc_mod")])),format="f",digits=0);
    row2a = formatC(100*as.matrix(cbind(TPR[curr_index[1],c("crm_emp","glob_emp","loc_emp")],
                                        TPR[curr_index[2],c("crm_emp","glob_emp","loc_emp")])),format="f",digits=0);
    row1b = formatC(100*as.matrix(cbind(FPR[curr_index[1],c("crm_mod","glob_mod","loc_mod")],
                                        FPR[curr_index[2],c("crm_mod","glob_mod","loc_mod")])),format="f",digits=0);
    row2b = formatC(100*as.matrix(cbind(FPR[curr_index[1],c("crm_emp","glob_emp","loc_emp")],
                                        FPR[curr_index[2],c("crm_emp","glob_emp","loc_emp")])),format="f",digits=0);
    row1 = cbind(row1a,"&",row1b);#paste0("(",row1a,",",row1b,")",row1c);
    row2 = cbind(row2a,"&",row2b);#paste0("(",row2a,",",row2b,")",row2c);
    write.table(matrix(c(row1_name,row1),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
    write.table(matrix(c(row2_name,row2),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
    #Normalized TPR/FPR
    file_name = paste0(write_to_folder,"/norm_sensspec",LETTERS[file_counter],".txt");
    write.table("",file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\n",append=F)
    row1a = formatC(100*as.matrix(cbind(norm_TPR[curr_index[1],c("crm_mod","glob_mod","loc_mod")],
                                        norm_TPR[curr_index[2],c("crm_mod","glob_mod","loc_mod")])),format="f",digits=0);
    row2a = formatC(100*as.matrix(cbind(norm_TPR[curr_index[1],c("crm_emp","glob_emp","loc_emp")],
                                        norm_TPR[curr_index[2],c("crm_emp","glob_emp","loc_emp")])),format="f",digits=0);
    row1b = formatC(100*as.matrix(cbind(norm_FPR[curr_index[1],c("crm_mod","glob_mod","loc_mod")],
                                        norm_FPR[curr_index[2],c("crm_mod","glob_mod","loc_mod")])),format="f",digits=0);
    row2b = formatC(100*as.matrix(cbind(norm_FPR[curr_index[1],c("crm_emp","glob_emp","loc_emp")],
                                        norm_FPR[curr_index[2],c("crm_emp","glob_emp","loc_emp")])),format="f",digits=0);
    row1 = cbind(row1a,"&",row1b);#paste0("(",row1a,",",row1b,")",row1c);
    row2 = cbind(row2a,"&",row2b);#paste0("(",row2a,",",row2b,")",row2c);
    write.table(matrix(c(row1_name,row1),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
    write.table(matrix(c(row2_name,row2),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
    
  }
}
#################################################################
#################################################################


#################################################################
#################################################################
#Sample size
row1_name = "Model";
row2_name = "Empiric";
file_counter = 0;

for(setting_start in c(1,3,5)) {
  file_counter = file_counter + 1;
  no_interim_index = which(nPerNonNullDEC_enrolled$setting%in%c(setting_start,setting_start+1)&nPerNonNullDEC_enrolled$dec_size==30&nPerNonNullDEC_enrolled$skip_first_stage==1);
  with_interim_index = which(nPerNonNullDEC_enrolled$setting%in%c(setting_start,setting_start+1)&nPerNonNullDEC_enrolled$dec_size==30&nPerNonNullDEC_enrolled$skip_first_stage==0);
  
  file_name = paste0(write_to_folder,"/sampsize",LETTERS[file_counter],".txt");
  write.table("",file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\n",append=F)
  row1a = formatC(as.matrix(cbind(nPerNonNullDEC_enrolled[no_interim_index[1],c("crm_mod","glob_mod","loc_mod")]-nPerNonNullDEC_enrolled[with_interim_index[1],c("crm_mod","glob_mod","loc_mod")],
                                  nPerNonNullDEC_enrolled[no_interim_index[2],c("crm_mod","glob_mod","loc_mod")]-nPerNonNullDEC_enrolled[with_interim_index[2],c("crm_mod","glob_mod","loc_mod")])),format="f",digits=1);
  row2a = formatC(as.matrix(cbind(nPerNonNullDEC_enrolled[no_interim_index[1],c("crm_emp","glob_emp","loc_emp")]-nPerNonNullDEC_enrolled[with_interim_index[1],c("crm_emp","glob_emp","loc_emp")],
                                  nPerNonNullDEC_enrolled[no_interim_index[2],c("crm_emp","glob_emp","loc_emp")]-nPerNonNullDEC_enrolled[with_interim_index[2],c("crm_emp","glob_emp","loc_emp")])),format="f",digits=1);
  row1b = formatC(as.matrix(cbind(nPerNullDEC_enrolled[no_interim_index[1],c("crm_mod","glob_mod","loc_mod")]-nPerNullDEC_enrolled[with_interim_index[1],c("crm_mod","glob_mod","loc_mod")],
                                  nPerNullDEC_enrolled[no_interim_index[2],c("crm_mod","glob_mod","loc_mod")]-nPerNullDEC_enrolled[with_interim_index[2],c("crm_mod","glob_mod","loc_mod")])),format="f",digits=1);
  row2b = formatC(as.matrix(cbind(nPerNullDEC_enrolled[no_interim_index[1],c("crm_emp","glob_emp","loc_emp")]-nPerNullDEC_enrolled[with_interim_index[1],c("crm_emp","glob_emp","loc_emp")],
                                  nPerNullDEC_enrolled[no_interim_index[2],c("crm_emp","glob_emp","loc_emp")]-nPerNullDEC_enrolled[with_interim_index[2],c("crm_emp","glob_emp","loc_emp")])),format="f",digits=1);
  row1 = cbind(row1a,"&",row1b);#paste0("(",row1a,",",row1b,")",row1c);
  row2 = cbind(row2a,"&",row2b);#paste0("(",row2a,",",row2b,")",row2c);
  write.table(matrix(c(row1_name,row1),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
  write.table(matrix(c(row2_name,row2),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
  
}
#################################################################
#################################################################

#################################################################
#################################################################
#Distance
row1_name = "Model";
row2_name = "Empiric";
file_counter = 0;
for(size_skip in list(c(15,0),c(30,1),c(30,0))) {
  for(setting_start in c(1,3,5)) {
    file_counter = file_counter + 1;
    curr_index = which(probability_distance$setting%in%c(setting_start,setting_start+1)&probability_distance$dec_size==size_skip[1]&probability_distance$skip_first_stage==size_skip[2]);

    file_name = paste0(write_to_folder,"/distance",LETTERS[file_counter],".txt");
    write.table("",file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\n",append=F)
    row1 = formatC(as.matrix(cbind(probability_distance[curr_index[1],c("crm_mod","glob_mod","loc_mod")],
                                   probability_distance[curr_index[2],c("crm_mod","glob_mod","loc_mod")])),format="f",digits=3);
    row2 = formatC(as.matrix(cbind(probability_distance[curr_index[1],c("crm_emp","glob_emp","loc_emp")],
                                   probability_distance[curr_index[2],c("crm_emp","glob_emp","loc_emp")])),format="f",digits=3);
    write.table(matrix(c(row1_name,row1),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
    write.table(matrix(c(row2_name,row2),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\n",append=T)
  } 
}
#################################################################
#################################################################
