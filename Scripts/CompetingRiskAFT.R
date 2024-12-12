#One stage competing risks AFT model
library(lme4)
library(tidyverse)
library(ggpubr)
library(geeM)

expit = function(x) exp(x)/(1 + exp(x))

#aft function with linear mixed model approach
aft_lmm = function(data, treat.mod, cens.mod, outcome.mod, cause){
  
  #Treatment model and censoring model 
  
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")

  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data$a - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  weights = weights[index]
  
  #Fitting linear mixed model
  
  model = lmer(formula = outcome.mod, data = data, weights = weights)
  fixef(model)
  
}

#aft function with GEE's
aft= function(data, treat.mod, cens.mod, outcome.mod, cause, corstr = "exchangeable", treat = "DonorHCV"){
  
  #Treatment model and censoring model 
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")
  
  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data[[treat]] - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  w = weights[index]
  
  data = data %>% mutate(w = w)
  
  #Fitting GEE with exchangeable correlation structure, geeM package requires data.frame
  
  data = data %>% arrange(group) %>% as.data.frame()
  model = geem(formula = outcome.mod, data = data, weights = w, id = group, corstr = corstr, sandwich = F, maxit=1)
  
  #Get main effect and interactions with treatment (all named elements of vector containing name of treatment)
  
  coefs = coef(model)
  
  blip = coefs[grepl(treat, names(coefs), fixed = T)]
  
  #list(mod = model, blip = blip)
  list(blip = blip, mod = model)
}


#One step GEE

aft1step= function(data, treat.mod, cens.mod, outcome.mod, cause, treat = "DonorHCV"){
  
  #Treatment model and censoring model 
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")
  
  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data[[treat]] - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  w = weights[index]
  
  data = data %>% mutate(w = w)
  
  #One step GEE
  
  data = data %>% arrange(group) %>% as.data.frame()
  model = gee1step(formula = outcome.mod, data = data, cluster = "group", family = "gaussian", w = w)
  
  #Get main effect and interactions with treatment (all named elements of vector containing name of treatment)
  
  coefs = coef(model)
  
  blip = coefs[grepl(treat, names(coefs), fixed = T)]
  
  #list(mod = model, blip = blip)
  list(blip = blip, mod = model)
}

bias_SE = function(psi_rep, psi){
  tab = matrix(0, nrow =2, ncol = length(psi))
  tab[1,] = apply(psi_rep, 2, mean) - psi
  tab[2,] = apply(psi_rep, 2, sd)
  
  tab
}

eval_DTR = function(data, psi1_hat, psi2_hat, cause_prob, regime){
  
  #Get the optimal regime depending on true underlying cause
  #Note : have to make sure this works the way i think it does
  blip_vars = c("RecHCV", "DON_TY", "age_n")
  blip_matrix = cbind(1, data %>% select(blip_vars) %>% as.matrix())
  
  
  data = data %>% mutate(blip1 = blip_matrix %*% psi1_hat, blip2 = blip_matrix %*% psi2_hat, blip = ifelse(epsilon ==1, blip1, blip2), opt= as.numeric(blip > 0))
  
  
  if(regime == "weighted"){
    data  = data %>% mutate(blip_hat = blip1*cause_prob + blip2*(1-cause_prob),opt_hat = as.numeric(blip_hat > 0))
  }
  else if(regime == "greedy"){
    data = data %>% mutate(blip_hat = ifelse(cause_prob > 0.5, blip1, blip2), opt_hat = as.numeric(blip_hat > 0))
  }
  else if(regime == "other"){
    data = data %>% mutate(blip_hat = blip1, opt_hat = as.numeric(blip_hat > 0))
  }
  else if(regime == "oracle"){
    data = data %>% mutate(blip_hat = blip, opt_hat  = opt)
  }
  
  #Identify proportion of optimal treatment (POT)
  n_row = nrow(data)
  pot = sum(data$opt == data$opt_hat)/n_row

  
  list(pot = pot, blip = data$blip_hat)
  
}


#Helper function to plot graph for each cause

get_plot_data = function(test, blips, blips.est){
  
  #Get 2.5%, 50% and 97.5% percentiles for blips of each individual in the test set
  blip_quantiles = t(apply(blips, 2, function(x) quantile(x, probs = c(0.025, 0.975)))) %>% as_tibble() 
  colnames(blip_quantiles) <- c("lower", "upper")
  
  dat = bind_cols(test, blip_quantiles) %>% mutate(blip = blips.est)%>% arrange(blip) #oracle

  dat
 
}


#Should turn this into two different functions to decouple estimation and plotting

#The fact that this returns cause_prob is not good -> too large

AFT.boot = function(n_boot, data, models,save = F, file = "default", regime = "joint", ...){
  
  
  #Estimate AFT models and cause model
  if(regime == "joint"){
    m1 = aft(data, models$treat.mod, models$cens.mod, models$out.mod,1)
    gc()
    m2 = aft(data, models$treat.mod, models$cens.mod, models$out.mod,2)
    gc()
    
    cause_mod.est = glm(formula = models$cause.mod, data = data %>% filter(delta == 1) , family = binomial)
    cause_prob.est = predict(cause_mod.est,newdata = data, type = "response")
    
    est = list(psi1_hat = m1$blip, psi2_hat = m2$blip, cause_prob = cause_prob.est, m1 = m1$mod, m2 = m2$mod)
  }
  else if(regime == "censor"){
    m = aft(data %>% mutate(delta = ifelse(delta == 1 & epsilon == 1, 1, 0)), 
            models$treat.mod, models$cens.mod, models$out.mod,1)
    gc()
    
    est = list(psi_hat = m$blip, m = m)
  }
  else if(regime == "composite"){
    m = aft(data %>% mutate(epsilon = 1), 
            models$treat.mod, models$cens.mod, models$out.mod,1)
    gc()
    
    est = list(psi_hat = m$blip, m = m)
  }
    
    
    
  
  
  #### Bootstrap ####
  
  #Parallel processing
  cl = makeCluster(3)
  clusterEvalQ(cl,{
    library(tidyverse)
    source("Scripts/CompetingRiskAFT.R")
  })
  
  #All unique cluster IDs
  groups = unique(data$group)
  
  #Main loop
      
  boot.est = function(i, data, models){
    
      boot = NULL
      
      #Cluster bootstrap
      sample_groups = sample(groups, length(groups), replace = T)
      for(j in 1:length(groups)){
        boot = rbind(boot, data %>% filter(group == sample_groups[j])) #Note this is probably very slow
      }
      
      #Basic bootstrap
      #n = nrow(data)
      #boot = data[sample(1:n, n, replace = T),]
      
      #Stratified bootstrap
      #case = data %>% filter(DonorHCV == 1); ncase = nrow(case)
      #control = data %>% filter(DonorHCV == 0); ncontrol = nrow(control)
      #boot = rbind(case[sample(1:ncase, ncase, replace = T),], control[sample(1:ncontrol, ncontrol, replace = T),])
      
      #Try catch in case estimation fails
      res = tryCatch({
          
         if(regime == "joint"){
           psi1 = aft(boot, models$treat.mod, models$cens.mod, models$out.mod,1)$blip
           gc()
           psi2 = aft(boot, models$treat.mod, models$cens.mod, models$out.mod,2)$blip
           gc()
           
           cause_mod = glm(formula = models$cause.mod, data = boot %>% filter(delta == 1) , family = binomial)
           cause_prob = predict(cause_mod,newdata = data, type = "response")
           
           list(psi1_hat = psi1, psi2_hat = psi2, cause_prob = cause_prob)
         }
         else if(regime == "censor"){
           psi_censor = aft(boot %>% mutate(delta = ifelse(delta == 1 & epsilon == 1, 1, 0)), models$treat.mod, models$cens.mod, models$out.mod,1)$blip
           gc()
           
           list(psi_hat = psi_censor)
         }
         else if(regime == "composite"){
           psi_composite = aft(boot %>% mutate(epsilon = 1), models$treat.mod, models$cens.mod, models$out.mod,1)$blip
           gc()
           
           list(psi_hat = psi_composite)
         }
       }, error = function(c){NULL})
      
      res
  }
  
  
  boot.est = parLapply(cl,1:n_boot, boot.est, data = data, models = models)
  stopCluster(cl)
  
  #### Results ####
  results = list(est = est, boot.est = boot.est, regime = regime)
  
  #Save raw results
  if(save){
    
    if(regime == "joint"){
      saveRDS(results, file = paste("Results/",file,".rds", sep=""))
    }
    else if(regime == "censor"){
      saveRDS(results, file = paste("Results/Censor/",file,".rds", sep=""))
    }
    else if(regime == "composite"){
      saveRDS(results, file = paste("Results/Composite/",file,".rds", sep=""))
    }
    
  }
  
  results

}

summaryAFT= function(data, results){
  
  #Initialize final results
  #Final results
  res = list()
  
  #Attach estimated blips from main analysis to data
  est = results$est
  
  #Collect raw bootstrap results
  raw = results$boot.est

  n_boot = length(raw)
  n = nrow(data)
  
  if(results$regime == "joint"){
    psi1_res = matrix(0, nrow = n_boot, ncol = 4)
    psi2_res = matrix(0, nrow = n_boot, ncol = 4)
    
    weighted_blip= matrix(0, nrow = n_boot, ncol = n)
    greedy_blip = matrix(0, nrow = n_boot, ncol = n)
    
    pot_w = rep(0, n_boot)
    pot_g = rep(0, n_boot)
    
    for(i in 1:n_boot){
      psi1_res[i,] = raw[[i]]$psi1_hat
      psi2_res[i,] = raw[[i]]$psi2_hat
      
      eval_w = eval_DTR(data, raw[[i]]$psi1_hat, raw[[i]]$psi2_hat, raw[[i]]$cause_prob, regime = "weighted")
      eval_g = eval_DTR(data, raw[[i]]$psi1_hat, raw[[i]]$psi2_hat, raw[[i]]$cause_prob, regime = "greedy")
      
      weighted_blip[i,]= eval_w$blip
      greedy_blip[i,] = eval_g$blip
      
      pot_w[i] = eval_w$pot
      pot_g[i] = eval_g$pot
      
    }
    
    eval.est_w = eval_DTR(data, est$psi1_hat, est$psi2_hat, est$cause_prob, regime = "weighted")
    eval.est_g = eval_DTR(data, est$psi1_hat, est$psi2_hat, est$cause_prob, regime = "greedy")
    
    dat.plot = data %>% mutate(oracle = eval_DTR(data, est$psi1_hat, est$psi2_hat, est$cause_prob, regime = "oracle")$blip)
    
    #Blip plots
    data_w = get_plot_data(dat.plot, weighted_blip, eval.est_w$blip)
    data_g = get_plot_data(dat.plot, greedy_blip, eval.est_g$blip)
    
    #POT and value
    res$measures =  c(weightedPOT = c(eval.est_w$pot,quantile(pot_w, probs = c(0.025, 0.975))), 
                      greedyPOT = c(eval.est_g$pot,quantile(pot_g, probs = c(0.025, 0.975))))
    
    #bootstrap CIs for parameters
    res$psi1_hat = psi1_res
    res$psi2_hat = psi2_res
    
    res$psi1_CI = t(apply(psi1_res, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
    res$psi2_CI = t(apply(psi2_res, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
    
    #Plotting data
    res$data_w = data_w 
    res$data_g = data_g 
    
  }else{
    psi_res = matrix(0, nrow = n_boot, ncol = 4)
    
    blip= matrix(0, nrow = n_boot, ncol = n)

    pot= rep(0, n_boot)
    
    
    for(i in 1:n_boot){
      psi_res[i,] = raw[[i]]$psi_hat

      
      eval = eval_DTR(data, raw[[i]]$psi_hat, raw[[i]]$psi_hat, NULL, regime = "other")
      
      blip[i,]= eval$blip
      
      pot[i] = eval$pot
      
    }
    
    eval.est = eval_DTR(data, est$psi_hat, est$psi_hat, NULL, regime = "other")
    
    #Blip plots
    data_comp = get_plot_data(data, blip, eval.est$blip)
    
    #POT and value
    res$measures =  c(POT = c(eval.est$pot,quantile(pot, probs = c(0.025, 0.975))))
    
    #bootstrap CIs for parameters
    res$psi_hat = psi_res
    
    res$psi_CI = t(apply(psi_res, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
    
    #Plotting data
    res$data_comp = data_comp 
    
  }
  
  
  #Return result list
  res
}

plotAFT = function(dat, name){
  
  dat1 = dat %>% filter(epsilon == 1)
  dat2 = dat %>% filter(epsilon == 2)
  
  colors <- c("Chosen Regime" = "black", "Oracle Regime" = "red") # geom_line(aes(y = oracle, color = "Oracle Regime"), linewidth = 1.2)  
  
  p1 =ggplot(dat1, aes(x = 1:nrow(dat1))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = blip, color = "Chosen Regime"), linewidth= 1.2) + labs(x="Observation number (Ordered)", y = "Estimated benefit",title ="Cause 1", color = "Legend",caption = name) + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  p2 =ggplot(dat2, aes(x = 1:nrow(dat2))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = blip,color = "Chosen Regime"), linewidth= 1.2) + labs(x="Observation number (Ordered)", y = "Estimated benefit",title ="Cause 2",color = "Legend", caption = name) + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  p3 =ggplot(dat, aes(x = 1:nrow(dat))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = blip,color = "Chosen Regime"), linewidth= 1.2) + labs(x="Observation number (Ordered)", y = "Estimated benefit",title ="Overall",color = "Legend", caption = name) + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  

  list(p1 = p1, p2 = p2, p3 = p3)
}