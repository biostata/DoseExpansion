####
#Oct 3, 2014
#Modification of titesim function that adds various desirable properties
##########

titesim_phil = function (PI, prior, target, n, x0, nsim = 1, restrict = TRUE, 
                         obswin = 1, tgrp = obswin, rate = 1, accrual = "fixed", surv = "uniform", surv_rate = obswin/10,
                         scheme = "polynomial", scheme_args = list(scheme_power = 1), count = TRUE, method = "bayes", model = "empiric", 
                         intcpt = 3, scale = sqrt(1.34), seed = 1009, 
                         conf.level = 0.5, no.exceed = Inf, cohort.size = 1, first.cohort.only = T, n.at.MTD = Inf, followup_b4_esc = obswin) 
{
  
  set.seed(seed)
  
  nexpt <- ntox <- matrix(0,nsim,length(prior));
  sel <- matrix(0, nsim, length(prior)+1);
  BETAHAT <- matrix(rep(NA, nsim * n), nrow = nsim)
  final.est <- DURATION <- n.enrolled <- rep(NA, nsim)
  ###Begin Phil's modification
  stop.for.tox = numeric(nsim);
  first.cohort.size = cohort.size;
  ###End Phil's modification	
  for (r in 1:nsim) {
    cohort.size = first.cohort.size;###Phil's modification
    if (count) {
      cat("simulation number:", r, "\n")
    }
    if (accrual == "fixed") {
      next.arrival <- obswin/rate
    } else if (accrual == "poisson") {
      next.arrival <- rexp(1, rate/obswin)
    }
    if (length(x0) > 1) {
      stop("Phil's modifications only made for standard 1-stage TITE-CRM")          
    } else {
      if (method == "mle") {
        stop(" Require an initial design for mle-CRM!")
      }
      bethat <- 0
      u <- y <- level <- arrival <- NULL
      cur <- x0
      i=1;
      while(i<=n) {
        arrival <- c(arrival, next.arrival)
        level <- c(level, cur)
          ynew <- rbinom(1, 1, PI[cur])
          if (surv == "uniform") {
            unew <- runif(1, 0, obswin)/ynew;
          } else if (surv == "exponential") {
            unew <- min(obswin,rexp(1,1/surv_rate))/ynew;
          }
          y <- c(y, ynew)
          u <- c(u, unew)
          utox <- u + arrival
        
        if(i<n) {
          if (accrual == "fixed") {
            next.arrival <- next.arrival + obswin/rate
          } else if (accrual == "poisson") {
            next.arrival <- next.arrival + rexp(1, rate/obswin)
          }
          B <- rep(0, length(y))
          B[utox <= next.arrival] <- 1
          censor <- pmin(next.arrival, utox) - arrival
          followup <- pmin(censor, obswin)
          obj <- titecrm_phil(prior, target, B, level, followup = followup, 
                         obswin = obswin, scheme = scheme, scheme_args = scheme_args, 
                         method = method, 
                         model = model, intcpt = intcpt, scale = scale, 
                         var.est = T, conf.level = conf.level)
        } else {
          followup <- pmin(utox-arrival,obswin)
          obj <- titecrm_phil(prior, target, y, level, weights = rep(1,n), method = method, model = model, intcpt = intcpt, 
                         scale = scale, var.est = T, conf.level = conf.level);
        }
        
        ###Begin Phil's modification
        if(i%%cohort.size==0) {#Only modify the dose at the end of each cohort
          if(restrict) {
            max.possible = max(c(0,which((obj$ptox-target)<=no.exceed)));#To indicate whether to stop early due to toxicity
            if(max.possible==0){
              if(i>=6) {              
                stop.for.tox[r] = i;#Variable to indicate at what patient the trial stopped
                bethat <- c(bethat, rep(-Inf,n-i));
                cur = 0;
                break;
              } else {
                max.possible = 1;#Continue the trial at the lowest dose if fewer than 6 patients have been enrolled
              }    
            }
            if((!any(followup[which(level==level[i])]>=followup_b4_esc))) {
              #Do not escalate until you've followed at least one patient for a certain period
              max.possible = min(max.possible,cur)
            }
            cur <- min(obj$mtd, (cur + 1), max.possible)#Do not skip a dose, do not exceed all of the previous restrictions
          } else {
            cur <- obj$mtd
          }
        }
        
        if(i%%cohort.size==0&first.cohort.only == T) {cohort.size=1;}#If cohorts only apply to the first batch of patients, set subsequent cohort sizes to 1. 
        ###End Phil's modification
        if(sum(level==cur)>=n.at.MTD|i==n) {
          bethat = c(bethat, rep(Inf,n-i));
          break;
        }
        i=i+1;
        bethat <- c(bethat, obj$est);
      }
      BETAHAT[r, ] <- bethat;
      est <- obj$est;
      msg <- "Okay"
    }
    
    sel[r,cur+1] <- 1;#first index of 'sel' corresponds to stopping for toxicity
    final.est[r] <- est;
    DURATION[r] <- max(arrival) + obswin;
    n.enrolled[r] = i;
    for (k in 1:length(prior)) {
      nexpt[r,k] <- length(which(level == k))			
      ntox[r,k] <- length(which(y == 1 & level == k))
    }
  }
  if (length(x0) == 1) {
    design <- paste("TITE-CRM starting at dose", x0);
  }
  else {
    design <- "Two-stage TITE-CRM";
  }
  foo <- list(all_sim = list(PI = PI, prior = prior, target = target, n = n, n.enrolled = n.enrolled,
                             x0 = x0, nsim = nsim, MTD = colMeans(sel), level = colMeans(nexpt), tox = colMeans(ntox), 
                             beta.hat = BETAHAT, final.est = final.est, Duration = DURATION, 
                             design = design, method = method, prior.var = scale^2, 
                             model = model, intcpt = intcpt, restriction = restrict, 
                             seed = seed, tite = TRUE, dosescaled = obj$dosescaled, 
                             msg = msg, obswin = obswin, tgrp = tgrp, rate = rate, 
                             accrual = accrual, scheme = scheme, scheme_args = scheme_args,
                             no.exceed = no.exceed, cohort.size = first.cohort.size, first.cohort.only = first.cohort.only, 
                             stop.for.tox = stop.for.tox, n.at.MTD = n.at.MTD, followup_b4_esc = followup_b4_esc, 
                             sel = sel, nexpt = nexpt, ntox = ntox),
              last_sim = list(PI = PI, prior = prior, target = target, n = n, n.enrolled = n.enrolled,
                              x0 = x0, nsim = nsim, MTD = cur, level = level, tox = y, 
                              beta.hat = bethat, final.est = est, arrival = arrival, 
                              toxicity.time = u, toxicity.study.time = utox, design = design, 
                              method = method, prior.var = scale^2, model = model, 
                              intcpt = intcpt, restriction = restrict, seed = seed, 
                              tite = TRUE, dosescaled = obj$dosescaled, msg = msg, 
                              obswin = obswin, tgrp = tgrp, rate = rate, accrual = accrual, 
                              scheme = scheme, scheme_args = scheme_args,
                              post.var = obj$post.var, ptox = obj$ptox, 
                              ptoxL = obj$ptoxL, ptoxU = obj$ptoxU, conf.level = obj$conf.level, 
                              no.exceed = no.exceed, cohort.size = first.cohort.size, first.cohort.only = first.cohort.only, 
                              stop.for.tox = stop.for.tox, n.at.MTD = n.at.MTD, followup_b4_esc = followup_b4_esc))
  
  class(foo) <- "philsim"
  foo
}

print.philsim = function (x, dgt = 3, patient.detail = TRUE, ...) 
{
  
  if (x$all_sim$nsim == 1) {
    x = x$last_sim
    n.enrolled = n = x$n.enrolled;
    PI <- x$PI
    prior <- round(x$prior, digits = dgt)
    target <- x$target
    K <- length(prior)
    y <- x$tox
    level <- x$level
    bethat <- signif(x$beta.hat, digits = 1)
    est <- round(x$final.est, digits = dgt)
    ptox <- round(x$ptox, digits = dgt)
    ptoxL <- round(x$ptoxL, digits = dgt)
    ptoxU <- round(x$ptoxU, digits = dgt)
    if (patient.detail) {
      if (x$tite) {
        arrival <- signif(x$arrival, digits = dgt)
        utox <- round(x$toxicity.study.time, digits = dgt)
        u <- round(x$toxicity.time, digits = dgt)
        tevent <- round(c(arrival, utox), digits = dgt)
        pid <- rep(1:n, 2)
        event <- c(rep("enrol", n), rep("TOX", n))
        level2 <- c(level, level)
        est2 <- round(c(bethat, rep(NA, n)), digits = dgt)
        o <- order(tevent)
        tevent <- tevent[o]
        pid <- pid[o]
        event <- event[o]
        level2 <- level2[o]
        est2 <- est2[o]
        ind <- which(tevent < Inf)
        tevent <- tevent[ind]
        pid <- pid[ind]
        event <- event[ind]
        level2 <- level2[ind]
        est2 <- est2[ind]
        m <- length(ind)
        cat("Trial summary on study time\n")
        cat("Time \t PID \t Event \t Level \t Beta \n")
        for (j in 1:m) {
          cat(tevent[j], "\t", pid[j], "\t", event[j], 
              "\t", level2[j], "\t", est2[j], "\n")
        }
        cat("\nPatient summary (TITE-CRM) \n")
        cat("PID \t Arrive \t Beta \t Level \t Tox \t Tox.time \n")
        for (i in 1:n) {
          cat(i, "\t", arrival[i], "\t\t", bethat[i], 
              "\t", level[i], "\t", y[i], "\t", u[i], "\n")
        }
      }
      else {
        cat("Patient summary (CRM) \n")
        cat("PID", "\t", "Beta", "\t", "Level", "\t", 
            "Toxicity", "\n")
        for (i in 1:n) {
          cat(i, "\t", bethat[i], "\t", level[i], "\t", 
              y[i], "\n")
        }
      }
    }
    cat("\nToxicity probability summary (with", x$conf.level * 
          100, "percent probability interval):", "\n")
    cat("Level", "\t", "Ptrue", "Prior", "\t", "n", "\t", 
        "ntox", "\t", "Posterior", "\t", "LoLmt", "\t", "UpLmt", 
        "\n")
    ntox <- nexpt <- rep(0, K)
    for (k in 1:K) {
      nexpt[k] <- length(which(level == k))
      ntox[k] <- length(which(level == k & y == 1))
      cat(k, "\t", PI[k], "\t", prior[k], "\t", nexpt[k], 
          "\t", ntox[k], "\t", ptox[k], "\t\t", ptoxL[k], 
          "\t", ptoxU[k], "\n")
    }
    PIjit <- PI
    PIjit[1] <- PI[1] - (1e-05)
    PIjit[K] <- PI[K] + (1e-05)
    cat("True MTD:", order(abs(PIjit - target))[1], "\tEstimated MTD:", 
        x$MTD, "\tTarget DLT rate:", target, "\n")
    cat("\nThis trial is generated by a", x$design, "\n")
    if (length(x$x0) > 1) {
      cat("Dose escalation proceeds as follows before any toxicity is seen:")
      xtab <- cbind(1:K, rep(NA, K))
      for (k in 1:K) {
        xtab[k, 2] = length(which(x$x0 == k))
      }
      colnames(xtab) <- c("dose.level", "cohort.size")
      rownames(xtab) <- rep("", K)
      print(t(xtab))
    }
    if (x$restrict) {
      ###Begin Phil's modification  
      cat("\nSafety constraints implemented:\n")
      cat("\t (1) No skipping doses in escalation;\n")
      cat("\t (2) No escalation before followup of",x$followup_b4_esc, "at current dose.\n")
      cat("\t (3) No assignment to dose with estimated DLT rate beyond ", x$no.exceed + x$target, ".\n",sep="")
      cat("\t (4) Stopping trial altogether if at any point after patient 6 the estimated DLT rate at all dose levels exceeds ", x$no.exceed + x$target, ".\n",sep="")
    }
    if(x$first.cohort.only) {
      cat("\nThe first", x$cohort.size, "patients were enrolled at the starting dose level, with patients assigned individually thereafter")
    } else {
      cat("\nPatients were enrolled in cohorts of", x$cohort.size)
    }
    if (x$stop.for.tox>0) {
      cat("\n-->IMPORTANT: At patient ", x$stop.for.tox, ", the lowest dose had an estimated toxicity rate \n\texceeding the pre-specified safety bound of ",x$no.exceed + x$target,", and the trial stopped. Final estimates of toxicity correspond to the patient immediately prior\n",sep="");
    }
    ###End Phil's modification
    cat("\nThe working model is", x$model, "\n")
    if (x$model == "empiric") {
      cat("\tptox = dose^{exp(beta)} with doses =", round(x$dosescaled, 
                                                          digits = dgt), "\n")
    }
    else {
      cat("\tlogit(ptox) = a + exp(beta)*dose, with a =", 
          x$intcpt, "\n\tand doses =", signif(x$dosescaled, 
                                              digits = dgt), "\n")
    }
    if (x$method == "bayes") {
      cat("\tand beta is estimated by its posterior mean \n\tassuming a normal prior with mean 0 and variance", 
          x$prior.var, "\n")
    }
    else if (x$method == "mle") {
      cat("\tand beta is estimated by its mle\n")
    }
    cat("\nThe final estimate of beta", round(x$final.est, 
                                              digits = dgt))
    if (x$method == "bayes") {
      cat(" with posterior variance", round(x$post.var, 
                                            digits = dgt), "\n")
    }
    else if (x$method == "mle") {
      cat(" with variance", round(x$post.var, digits = dgt), 
          "\n")
      if (x$msg != "Okay") {
        print(x$msg)
      }
    }
    if (x$tite) {
      cat("\nThe", x$scheme, "function is used to assign weights to patients, using parameters\n",unlist(x$scheme_args),".\n\n")
      cat("Patient arrival is modeled as a", x$accrual, 
          "process\n")
      cat("\twith rate", x$rate, "patients per", x$obswin, 
          "time units (= observation window).\n")
      if (length(x$x0) > 1) {
        cat("\tA minimum waiting time of", x$tgrp, "time units is imposed\n")
        cat("\tbetween two dose cohorts in the initial stage.\n")
      }
    }
  } else {
    x = x$all_sim
    n.success <- mean(x$n.enrolled[which(x$stop.for.tox==0)]);
    n.fail <- mean(x$n.enrolled[which(x$stop.for.tox>0)]);
    PI <- x$PI
    prior <- round(x$prior, digits = dgt)
    target <- x$target
    K <- length(prior)
    oc <- t(cbind(PI, prior, x$MTD[-1], x$level, x$tox))
    colnames(oc) <- as.character(1:K)
    rownames(oc) <- c("Truth", "Prior", "Selected", "Nexpt", 
                      "Ntox")
    cat("\nNumber of simulations:\t", x$nsim, "\n")
    cat("Patients accrued in successfully completed trials:\t", n.success, "\n")
    cat("Patients accrued in trials stopped for toxicity:\t", n.fail, "\n")
    cat("Target DLT rate:\t", target, "\n")
    print(round(oc, digits = dgt))
    if (x$tite) {
      cat("\nThe distribution of trial duration:\n")
      print(summary(x$Dur))
    }
    cat("\nThe trials are generated by a", x$design, "\n")
    if (length(x$x0) > 1) {
      cat("Dose escalation proceeds as follows before any toxicity is seen:")
      xtab <- cbind(1:K, rep(NA, K))
      for (k in 1:K) {
        xtab[k, 2] = length(which(x$x0 == k))
      }
      colnames(xtab) <- c("dose.level", "cohort.size")
      rownames(xtab) <- rep("", K)
      print(t(xtab))
    }
    if (x$restrict) {
      cat("\nSafety constraints implemented:\n")
      ###Begin Phil's modification
      cat("\t (1) No skipping doses in escalation;\n")
      cat("\t (2) No escalation before followup of",x$followup_b4_esc, "at current dose.\n")
      cat("\t (3) No assignment to dose with estimated DLT rate beyond ", x$no.exceed + x$target, ".\n",sep="")
      cat("\t (4) Stopping trial altogether if at any point after patient 6 the estimated DLT rate at all dose levels exceeds ", x$no.exceed + x$target, ".\n",sep="")
    }
    if(x$first.cohort.only) {
      cat("\nThe first", x$cohort.size, "patients were enrolled at the starting dose level, with patients assigned individually thereafter")
    } else {
      cat("\nPatients were enrolled in cohorts of", x$cohort.size)
    }
    cat("\nThe proportion of trials for which at some point during the trial the lowest dose had an estimated rate of toxicity exceeding ", x$no.exceed + x$target, ", which is the pre-specified safety bound, was ", mean(x$stop.for.tox>0),".\n",sep="");
    ###End Phil's modification
    cat("\nThe working model is", x$model, "\n")
    if (x$model == "empiric") {
      cat("\tptox = dose^{exp(beta)} with doses =", round(x$dosescaled, 
                                                          digits = dgt), "\n")
    }
    else {
      cat("\tlogit(ptox) = a + exp(beta)*dose, with a =", 
          x$intcpt, "\n\tand doses =", signif(x$dosescaled, 
                                              digits = dgt), "\n")
    }
    x0 <- x$x0
    if (x$method == "bayes") {
      cat("\tand beta is estimated by its posterior mean \n\tassuming a normal prior with mean 0 and variance", 
          x$prior.var, "\n")
    }
    else if (x$method == "mle") {
      cat("\tand beta is estimated by its mle\n")
    }
    if (x$tite) {
      cat("\nThe", x$scheme, "function is used to assign weights to patients, using parameters\n",unlist(x$scheme_args),".\n\n")
      cat("Patient arrival is modeled as a", x$accrual, 
          "process\n")
      cat("\twith rate", x$rate, "patients per", x$obswin, 
          "time units (= observation window).\n\n")
      if (length(x$x0) > 1) {
        cat("\tA minimum waiting time of", x$tgrp, "time units is imposed\n")
        cat("\tbetween two dose cohorts in the initial stage.\n")
      }
    }
  }
}

crm_phil = function (prior, target, tox, level, n = length(level), dosename = NULL, 
                     include = 1:n, pid = 1:n, conf.level = 0.9, method = "bayes", 
                     model = "empiric", intcpt = 3, scale = sqrt(1.34), model.detail = TRUE, 
                     patient.detail = TRUE, var.est = TRUE) 
{
  y1p <- tox[include]
  w1p <- rep(1, length(include))
  if(model!="empiric"|method!="bayes") {
    stop("only 'bayes' method and 'empiric' model implemented");
  }
  
  dosescaled <- prior
  x1p <- prior[level[include]]
  
  den <- integrate(crmh_phil, -Inf, Inf, x1p, y1p, w1p, 
                   scale, abs.tol = 0)[[1]]
  est <- integrate(crmht_phil, -10, 10, x1p, y1p, w1p, scale, 
                   abs.tol = 0)[[1]]/den
  if (var.est) {
    e2 <- integrate(crmht2_phil, -10, 10, x1p, y1p, w1p, 
                    scale, abs.tol = 0)[[1]]/den
  }
  
  ptox <- 1-(1-prior)^exp(est)
  if (var.est) {
    post.var <- e2 - est^2
    crit <- qnorm(0.5 + conf.level/2)
    lb <- est - crit * sqrt(post.var)
    ub <- est + crit * sqrt(post.var)
    ptoxL <- 1-(1-prior)^exp(lb)
    ptoxU <- 1-(1-prior)^exp(ub)
  }
  
  if (all(ptox <= target)) {
    rec <- length(prior)
  }
  else if (all(ptox >= target)) {
    rec <- 1
  }
  else {
    rec <- order(abs(ptox - target))[1]
  }
  if (!var.est) {
    post.var <- ptoxU <- ptoxL <- NA
  }
  foo <- list(prior = prior, target = target, tox = tox, level = level, 
              dosename = dosename, subset = pid[include], estimate = est, 
              model = model, prior.var = scale^2, post.var = post.var, 
              method = method, mtd = rec, include = include, pid = pid, 
              model.detail = model.detail, intcpt = intcpt, ptox = ptox, 
              ptoxL = ptoxL, ptoxU = ptoxU, conf.level = conf.level, 
              patient.detail = patient.detail, tite = FALSE, dosescaled = dosescaled, 
              var.est = var.est)
  class(foo) <- "mtd"
  foo
}

crmh_phil = function (a, x, y, w, s) 
{
  v = exp(-a^2/2/s^2);
  for (i in 1:length(x)) {
    v = v * (1-(1-x[i])^exp(a))^y[i] * (1 - w[i] * (1-(1-x[i])^exp(a)))^(1 - y[i])
  }
  return(v)
}
crmht_phil = function (a, x, y, w, s) 
{
  v = a * exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    v = v * (1-(1-x[i])^exp(a))^y[i] * (1 - w[i] * (1-(1-x[i])^exp(a)))^(1 - y[i])
  }
  return(v)
}
crmht2_phil = function (a, x, y, w, s) 
{
  v = a^2 * exp(-a^2/2/s^2)
  for (i in 1:length(x)) {
    v = v * (1-(1-x[i])^exp(a))^y[i] * (1 - w[i] * (1-(1-x[i])^exp(a)))^(1 - y[i])
  }
  return(v)
}

titecrm_phil = function (prior, target, tox, level, n = length(level), weights = NULL, 
          followup = NULL, entry = NULL, exit = NULL, obswin = NULL, 
          scheme = "polynomial", scheme_args = list(scheme_power=1), conf.level = 0.9, dosename = NULL, include = 1:n, 
          pid = 1:n, method = "bayes", model = "empiric", var.est = TRUE, 
          scale = sqrt(1.34), intcpt = 3, model.detail = TRUE, patient.detail = TRUE, 
          tite = TRUE) 
{
  if (is.null(weights)) {
    if (is.null(followup)) {
      followup <- exit - entry
    }
    if (scheme == "polynomial") {
      weights <- (followup^scheme_args$scheme_power)/(obswin^scheme_args$scheme_power);
    } else if (scheme == "logistic") {
      weights <- ifelse(followup>=obswin,1,1/(1+exp(-scheme_args$scheme_int - followup * scheme_args$scheme_slope)));
    } else if (scheme == "adaptive") {
      support <- sort(followup[tox == 1])
      z <- length(support)
      if (z) {
        for (i in 1:n) {
          m <- length(support[support <= followup[i]])
          if (!m) 
            weights[i] <- followup[i]/support[1]/(z + 
                                                    1)
          else if (m == z) 
            weights[i] <- (z + (followup[i] - support[z])/(obswin - 
                                                             support[z]))/(z + 1)
          else weights[i] <- (m + (followup[i] - support[m])/(support[m + 
                                                                        1] - support[m]))/(z + 1)
        }
      }
      else {
        weights <- followup/obswin
      }
    } else {
      stop(" Weighting scheme undefined!")
    }
    weights <- pmin(weights, 1)
  }
  if (any(weights > 1) | any(weights < 0)) 
    stop(" Weights have to be between 0 and 1!")
  if (is.null(pid)) {
    if (!(length(tox) == length(level) & length(tox) == length(weights))) 
      stop(" tox, level, and weights are of different lengths!")
  }
  else {
    if (!(length(tox) == length(level) & length(tox) == length(weights) & 
          length(tox) == length(pid))) 
      stop(" pid, tox, level, and weights are of different lengths!")
  }
  weights[tox == 1] <- 1;
  y1p <- tox[include]
  w1p <- weights[include]
  if (model == "empiric") {
    dosescaled <- prior
    x1p <- prior[level[include]]
    if (method == "mle") {
      if (sum(y1p) == 0 | sum(y1p) == length(y1p)) 
        stop(" mle does not exist!")
      est <- optimize(lcrm, c(-10, 10), x1p, y1p, w1p, 
                      tol = 1e-04, maximum = TRUE)$max
      if (var.est) {
        e2 <- integrate(crmht2, -10, 10, x1p, y1p, w1p, 
                        500, abs.tol = 0)[[1]]/integrate(crmh, -10, 
                                                         10, x1p, y1p, w1p, 500, abs.tol = 0)[[1]]
      }
    }
    else if (method == "bayes") {
      den <- integrate(crmh, -Inf, Inf, x1p, y1p, w1p, 
                       scale, abs.tol = 0)[[1]]
      est <- integrate(crmht, -10, 10, x1p, y1p, w1p, scale, 
                       abs.tol = 0)[[1]]/den
      if (var.est) {
        e2 <- integrate(crmht2, -10, 10, x1p, y1p, w1p, 
                        scale, abs.tol = 0)[[1]]/den
      }
    }
    else {
      stop(" unknown estimation method")
    }
    ptox <- prior^exp(est)
    if (var.est) {
      post.var <- e2 - est^2
      crit <- qnorm(0.5 + conf.level/2)
      lb <- est - crit * sqrt(post.var)
      ub <- est + crit * sqrt(post.var)
      ptoxL <- prior^exp(ub)
      ptoxU <- prior^exp(lb)
    }
  }
  else if (model == "logistic") {
    dosescaled <- log(prior/(1 - prior)) - intcpt
    if (!all(dosescaled < 0)) {
      stop("intercept parameter in logit model is too small: scaled doses > 0!")
    }
    x1p <- dosescaled[level[include]]
    if (method == "mle") {
      if (sum(y1p) == 0 | sum(y1p) == length(y1p)) 
        stop(" mle does not exist!")
      est <- optimize(lcrmlgt, c(-10, 10), x1p, y1p, w1p, 
                      intcpt, tol = 1e-04, maximum = TRUE)$max
      if (var.est) {
        e2 <- integrate(crmht2lgt, -10, 10, x1p, y1p, 
                        w1p, 500, intcpt, abs.tol = 0)[[1]]/integrate(crmhlgt, 
                                                                      -10, 10, x1p, y1p, w1p, 500, abs.tol = 0)[[1]]
      }
    }
    else if (method == "bayes") {
      den <- integrate(crmhlgt, -Inf, Inf, x1p, y1p, w1p, 
                       scale, intcpt, abs.tol = 0)[[1]]
      est <- integrate(crmhtlgt, -10, 10, x1p, y1p, w1p, 
                       scale, intcpt, abs.tol = 0)[[1]]/den
      if (var.est) {
        e2 <- integrate(crmht2lgt, -10, 10, x1p, y1p, 
                        w1p, scale, intcpt, abs.tol = 0)[[1]]/den
      }
    }
    else {
      stop(" unknown estimation method")
    }
    ptox <- (1 + exp(-intcpt - exp(est) * dosescaled))^{
      -1
    }
    if (var.est) {
      post.var <- e2 - est^2
      crit <- qnorm(0.5 + conf.level/2)
      lb <- est - crit * sqrt(post.var)
      ub <- est + crit * sqrt(post.var)
      ptoxL <- (1 + exp(-intcpt - exp(ub) * dosescaled))^{
        -1
      }
      ptoxU <- (1 + exp(-intcpt - exp(lb) * dosescaled))^{
        -1
      }
    }
  }
  else {
    stop(" model specified not available.")
  }
  if (all(ptox <= target)) {
    rec <- length(prior)
  }
  else if (all(ptox >= target)) {
    rec <- 1
  }
  else {
    rec <- order(abs(ptox - target))[1]
  }
  if (!var.est) {
    post.var <- ptoxL <- ptoxU <- NA
  }
  foo <- list(prior = prior, target = target, tox = tox, level = level, 
              dosename = dosename, subset = pid[include], estimate = est, 
              weights = weights, followup = followup, entry = entry, 
              exit = exit, obswin = obswin, scheme = scheme, scheme_args = scheme_args, model = model, 
              prior.var = scale^2, post.var = post.var, method = method, 
              mtd = rec, include = include, pid = pid, model.detail = model.detail, 
              intcpt = intcpt, ptox = ptox, ptoxL = ptoxL, ptoxU = ptoxU, 
              conf.level = conf.level, patient.detail = patient.detail, 
              tite = tite, dosescaled = dosescaled)
  class(foo) <- "mtd"
  foo
}

