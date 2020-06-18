require(DEoptim) # DEoptimR ?
require(tidyverse)

orig = read.table("x6-tc_reorder/orig_order4x4.txt", header=F)
x6_1 = read.table("x6-tc_reorder/reord_orig_1.txt", header=F) # avg 1 Successive Repetition
x6_2 = read.table("x6-tc_reorder/reord_orig_2.txt", header=F) # avg 1.41 SRs
x6_3 = read.table("x6-tc_reorder/reord_orig_3.txt", header=F)
x6_4 = read.table("x6-tc_reorder/reord_orig_4.txt", header=F)

tc_ords = list(tc0=orig, x6_1=x6_1, x6_2=x6_2, x6_3=x6_3, x6_4=x6_4)
#save(tc_ords, file="tc_trial_orders.RData")

x13_1olap2tr = read.table("x13-temp-cont/1olap2tr_contr.txt", header=F)
x13_1olap3tr = read.table("x13-temp-cont/1olap3tr_contr.txt", header=F)
x13_2olap2tr = read.table("x13-temp-cont/2olap2tr_contr.txt", header=F)

#x8_1 = read.table("x8-maxTC-wwo-spatTC/1_max_temp_spat_cont_orig.txt", header=F) 
x8_4 = read.table("x8-maxTC-wwo-spatTC/4_no_spat_orig_max_tc.txt", header=F) 

# Exp 1
#ords = list(tc0=orig, x6_1=x6_1, x6_2=x6_2, x6_3=x6_3, x6_4=x6_4) 
# Exp 2
#ords = list(x8_4=x8_4, x13_1=x13_1olap2tr, x13_2=x13_1olap3tr, x13_3=x13_2olap2tr)

ords = list(tc0=orig, x6_1=x6_1, x6_2=x6_2, x6_3=x6_3, x6_4=x6_4, 
            x8_4=x8_4, x13_1=x13_1olap2tr, x13_2=x13_1olap3tr, x13_3=x13_2olap2tr)

load("merged_temporal_contiguity_data.RData")

hum = bdi %>% group_by(cond, item, itemSR, condSR) %>% 
  summarize(accuracy = mean(Correct)) 

mnames = c("model","fazly") # ..
# for(mname in mnames) {}
mname = "model"
mname = "model_WM"
source(paste("models/",mname,".R", sep=''))
sse = c()
for(o in names(ords)) {
  m = model(c(.1, 2.0, .95, .8), ord=ords[[o]])$perf
  sse = c(sse, sum((m - subset(hum, cond==o)$accuracy)^2))
  print(cor(subset(hum, cond==o)$accuracy, m))
}
sum(sse)
# model
# 0.4306615 -0.5034796 0.3476981 0.1449971 0.06054025
# sum(sse) 1.375584

#evalSSE <- function(par, ord, hum) {
#  return( sum((model(par, ord)$perf - hum$accuracy)^2) )
#}

simulate_subjects <- function(ord, param, n=500) {
  voc_sz = max(unlist(ord))
  #lex = matrix(0, nrow=n, ncol=voc_sz)
  totlex = matrix(0, nrow=voc_sz, ncol=voc_sz)
  for(s in 1:n) {
    mod = model(param, ord)
    #lex[s,] = mod$perf
    totlex = totlex + mod$matrix
  }
  pcor = diag(totlex) / n
  #print(totlex)
  #print(pcor)
  #print(cor(pcor, colMeans(lex)))
  return(pcor)
}

evalSSE <- function(par, ord, hum, stochastic=F) {
  if(stochastic) {
    mod = simulate_subjects(orig, par) 
  } else {
    mod = model(par, ord)$perf 
  }
  SSE = sum((mod - hum$accuracy)^2)
  return( SSE )
}


model_perf_df <- function(par, ord, condname, mname) {
  source(paste("models/",mname,".R", sep=''))
  if(grepl("stochastic", mname)==1) {
    mod = simulate_subjects(orig, par)  # need to run and average 100s of Ss
  } else {
    mod = model(par, ord)$perf 
  }
  df = data.frame(model=mname, cond=condname, item=1:length(mod), accuracy=mod)
  return(df)
}

# should also make a fit_by_exp
fit_per_cond <- function(ords, hum, mname, lowerp, upperp) {
  if(grepl("stochastic", mname)==1) {
    stochastic = T # need to run and average 100s of Ss
  } else {
    stochastic = F
  }
  source(paste("models/",mname,".R", sep=''))
  df = data.frame(model=mname, order=names(ords), fit=NA) 
  parcols = paste("par",1:length(lowerp),sep='')
  df[,parcols] = NA
  for(oname in names(ords)) {
  #oname = "x6_1"
    hum_c = subset(hum, cond==oname)
    best <- DEoptim(fn=evalSSE, ord=ords[[oname]], hum=hum_c, stochastic=stochastic, lower=lowerp, upper=upperp, 
                  DEoptim.control(reltol=.001, steptol=50, itermax=200, trace=10)) # was itermax=5000, trace=250))
    df[which(df$order==oname),"fit"] = best$optim$bestval
    df[which(df$order==oname),parcols] = best$optim$bestmem 
  }
  return(df)
  #return(list(par=goodpar, fit=goodfit))
}


# find best-fitting parameters for each condition separately
set.seed(123)

fit_mod = fit_per_cond(ords, hum, "model", c(.001, .05, .7), c(2, 7, 1))
mean(fit_mod$fit) # .17 / .197

fitWM_mod = fit_per_cond(ords, hum, "model_WM", c(.001, .1, .7, 0), c(2, 7, 1, 1))
mean(fitWM_mod$fit) # .17 / .19

fit_bay = fit_per_cond(ords, hum, "Bayesian_decay", c(.01, .01, .01), c(1, 1, 1))
mean(fit_bay$fit) # .21 / .239

fit_tru = fit_per_cond(ords, hum, "stochastic/trueswell2012_model_detailed_uncertain", c(.001, .001), c(1,1))
mean(fit_tru$fit) # .17 / .185

fit_hyp = fit_per_cond(ords, hum, "stochastic/hypoth_model_detailed", c(.001, .001), c(1,1))
mean(fit_hyp$fit) # .16 / .178

fitPur_mod = fit_per_cond(ords, hum, "stochastic/pursuit_detailed", c(.001, .001, .00001), c(1, 1, .2))
mean(fitPur_mod$fit) # .61 / .65 # was run with 500 subjects

fit_faz = fit_per_cond(ords, hum, "fazly", c(1e-7, 10), c(1,10000))
mean(fit_faz$fit) # 0.23 / .252

fit_str = fit_per_cond(ords, hum, "strength", c(1e-5, .7), c(2,1))
mean(fit_str$fit) # .19 / .226
# can also fit strength, novelty, and uncertainty models

fit_mod_samp = fit_per_cond(ords, hum, "stochastic/model_detailed_sampling", c(.001, .05, .7), c(2, 7, 1))
mean(fit_mod_samp$fit) # .177

two_par = rbind(fit_hyp, fit_tru, fit_faz, fit_str)
two_par$par3 = NA
model_pars = rbind(two_par, fitPur_mod, fit_mod, fit_bay, fit_mod_samp)
model_pars$par4 = NA
model_pars = rbind(model_pars, fitWM_mod)
save(model_pars, file="temporal_contiguity_model_ALL_fits.RData")


# need to get model predicted item accuracies per condition in a similar dataframe
mod = data.frame()
for(i in 1:nrow(model_pars)) { # 
  condname = as.character(model_pars[i,]$order)
  par = model_pars[i , paste("par", 1:4, sep='') ]
  par = par[ , colSums(is.na(par)) == 0]
  mod = rbind(mod, model_perf_df(as.numeric(par), ords[[condname]], condname, as.character(model_pars[i,]$model)))
}

save(mod, file="exp1n2_fitted_models_item_perf.RData")
#model_perf_df(c(.1, 1, .95), tc_ords[["tc0"]], "tc0", "model")
