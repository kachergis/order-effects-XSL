require(DEoptim) # DEoptimR ?

orig = read.table("x6-tc_reorder/orig_order4x4.txt", header=F)
x6_1 = read.table("x6-tc_reorder/reord_orig_1.txt", header=F) # avg 1 Successive Repetition
x6_2 = read.table("x6-tc_reorder/reord_orig_2.txt", header=F) # avg 1.41 SRs
x6_3 = read.table("x6-tc_reorder/reord_orig_3.txt", header=F)
x6_4 = read.table("x6-tc_reorder/reord_orig_4.txt", header=F)

x13_1olap2tr = read.table("x13-temp-cont/1olap2tr_contr.txt", header=F)
x13_1olap3tr = read.table("x13-temp-cont/1olap3tr_contr.txt", header=F)
x13_2olap2tr = read.table("x13-temp-cont/2olap2tr_contr.txt", header=F)

#x8_1 = read.table("x8-maxTC-wwo-spatTC/1_max_temp_spat_cont_orig.txt", header=F) 
x8_4 = read.table("x8-maxTC-wwo-spatTC/4_no_spat_orig_max_tc.txt", header=F) 

ords = list(tc0=orig, x6_1=x6_1, x6_2=x6_2, x6_3=x6_3, x6_4=x6_4) 
# x8_4=x8_4, x13_1=x13_1olap2tr, x13_2=x13_1olap3tr, x13_3=x13_2olap2tr)

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

evalSSE <- function(par, ord, hum) {
  return( sum((model(par, ord)$perf - hum$accuracy)^2) )
}

fit_per_cond <- function(ords, hum, mname, lowerp, upperp) {
  source(paste("models/",mname,".R", sep=''))
  df = data.frame(model=mname, order=names(ords), fit=NA) 
  parcols = paste("par",1:length(lowerp),sep='')
  df[,parcols] = NA
  for(oname in names(ords)) {
  #oname = "x6_1"
    hum_c = subset(hum, cond==oname)
    best <- DEoptim(fn=evalSSE, ord=ords[[oname]], hum=hum_c, lower=lowerp, upper=upperp, 
                  DEoptim.control(reltol=.001, steptol=50, itermax=200, trace=10)) # was itermax=5000, trace=250))
    df[which(df$order==oname),"fit"] = best$optim$bestval
    df[which(df$order==oname),parcols] = best$optim$bestmem 
  }
  return(df)
  #return(list(par=goodpar, fit=goodfit))
}

fit_mod = fit_per_cond(ords, hum, "model", c(.001, .05, .7), c(2, 7, 1))
mean(fit_mod$fit) # .17

fitWM_mod = fit_per_cond(ords, hum, "model_WM", c(.001, .1, .7, 0), c(2, 7, 1, 1))
mean(fitWM_mod$fit) # .17

fit_bay = fit_per_cond(ords, hum, "Bayesian_decay", c(.01, .01, .01), c(1, 1, 1))
mean(fit_bay$fit) # .21