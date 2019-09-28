# got all the temporal contiguity experiment data and 4x4 TC=0 (0 SR) data
require(tidyverse)
require(ggplot2)

orig = read.table("x6-tc_reorder/orig_order4x4.txt", header=F)
x6_1 = read.table("x6-tc_reorder/reord_orig_1.txt", header=F) # avg 1 Successive Repetition
x6_2 = read.table("x6-tc_reorder/reord_orig_2.txt", header=F) # avg 1.41 SRs
x6_3 = read.table("x6-tc_reorder/reord_orig_3.txt", header=F)
x6_4 = read.table("x6-tc_reorder/reord_orig_4.txt", header=F)
x6_data = read.table("x6-tc_reorder/x6-tc_data9-9-2009.txt", header=T)

load("x6-tc_reorder/freqCD_item_data.Rdata")
tc0 = subset(raw, Cond=="4x4") # 76 subjects, mean acc = .32

x8_1 = read.table("x8-maxTC-wwo-spatTC/1_max_temp_spat_cont_orig.txt", header=F) 
# 2_x8_39.txt 3_x8_369.txt - from freqCD paper
x8_4 = read.table("x8-maxTC-wwo-spatTC/4_no_spat_orig_max_tc.txt", header=F) 
x8_data = read.table("x8-maxTC-wwo-spatTC/x8_9-22-08.txt", header=T)


# not useful for evaluating effects of order on particular items (because these are not just shuffled from orig)
x13_data = read.table("x13-temp-cont/x13_data_1-21-2009.txt", header=T)
x13_1olap2tr = read.table("x13-temp-cont/1olap2tr_contr.txt", header=F)
x13_1olap3tr = read.table("x13-temp-cont/1olap3tr_contr.txt", header=F)
x13_2olap2tr = read.table("x13-temp-cont/2olap2tr_contr.txt", header=F)
# block4_369_39mx.txt - from freqCD paper

get_trial_item_successive_repetitions <- function(ord, condName='') {
  voc_sz = length(unique(unlist(ord)))
  sr_item = rep(0, voc_sz) # track # of times each item is repeated in consecutive trials
  sr_tr = rep(0, nrow(ord)) # track # of repeated items in each trial
  for(t in 2:nrow(ord)) {
    repeated = intersect(unlist(ord[t,]), unlist(ord[t-1,]))
    if(length(repeated)>0) {
      sr_tr[t] = length(repeated)
      sr_item[repeated] = sr_item[repeated] + 1
    } 
  }
  df = data.frame(cond=condName, item=1:voc_sz, itemSR=sr_item, condSR=mean(sr_tr))
  #return(list(trial=sr_tr, item=sr_item))
  return(df)
}

mean_SR <- function (ord, verbose=T) {
  sr = get_trial_item_successive_repetitions(ord)
  if(verbose) cat("trial: ", mean(sr$condSR), "item: ", mean(sr$itemSR), "\n")
  return(list(cond=mean(sr$condSR), item=mean(sr$itemSR)))
}

# item-level data frames for each experiment (built from trial orderings)
x6sr = rbind(get_trial_item_successive_repetitions(x6_1, condName="x6_1"),
      get_trial_item_successive_repetitions(x6_2, condName="x6_2"),
      get_trial_item_successive_repetitions(x6_3, condName="x6_3"),
      get_trial_item_successive_repetitions(x6_4, condName="x6_4"))

x8sr = rbind(get_trial_item_successive_repetitions(x8_1, condName="x8_1"),
           get_trial_item_successive_repetitions(x8_4, condName="x8_4"))

x13sr = rbind(get_trial_item_successive_repetitions(x13_1olap2tr, condName="x13_1"),
           get_trial_item_successive_repetitions(x13_1olap3tr, condName="x13_2"),
           get_trial_item_successive_repetitions(x13_2olap2tr, condName="x13_3"))
# 1_ is 1olap2tr
#	2_ is 1olap3tr
#	3_ is 3_2olap2tr

tc0sr = get_trial_item_successive_repetitions(orig, cond="tc0")


# x6 has shuffles of orig (SR=0)
mean_SR(x6_1) # tr: 1 SR, item: 1.5
mean_SR(x6_2) # tr: 1.41 item: 2.11
mean_SR(x6_3) # tr: .67  item: 1
mean_SR(x6_4) # tr: .33  item: .5

# can start by looking at differences in average performance across these 18 items 
# (and predictions from the models for the reordered items)

# now get binary accuracy data for logistic regression
x13_ia = x13_data %>% filter(Condition!="4_") %>%
  group_by(Condition, item=CorrectAns, Subject=Subject)

length(unique(x8_data$Subject)) # 33 subjects
x8_ia = x8_data %>% filter(Condition=="1_" | Condition=="4_") %>%
  group_by(Condition, item=CorrectAns, Subject=Subject) 

length(unique(x6_data$Subject)) # 36 subjects (2 Ss missed condition 4)
x6_ia = x6_data %>% 
  group_by(Condition, item=CorrectAns, Subject=Subject) 

x6_ia$cond = with(x6_ia, ifelse(Condition=="1_", "x6_1", 
                                  ifelse(Condition=="2_", "x6_2", 
                                         ifelse(Condition=="3_", "x6_3", "x6_4"))))

x8_ia$cond = with(x8_ia, ifelse(Condition=="1_", "x8_1", "x8_4"))
x8_ia$Subject = x8_ia$Subject + 100 # unique subject IDs

x13_ia$cond = with(x13_ia, ifelse(Condition=="1_", "x13_1", ifelse(Condition=="2_", "x13_2", "x13_3")))
x13_ia$Subject = x13_ia$Subject + 200 # unique subject IDs

tc0$cond = "tc0"
tc0$item = tc0$CorrectAns

# could look at accuracy as a function of TestIndex to look for memory decay effects during test
# (lack thereof may be evidence for somewhat durable knowledge! at least not quickly-fading) -- looks good in x6 & x8
x6_ia = left_join(x6_ia, x6sr, by=c("item","cond"))
x8_ia = left_join(x8_ia, x8sr, by=c("item","cond"))
x13_ia = left_join(x13_ia, x13sr, by=c("item","cond"))
tc0_ia = left_join(tc0, tc0sr, by=c("item","cond"))
wanted_cols = c("Subject","Correct","Response","RT","CorrectAns","item","cond","itemSR","condSR")
bdi = rbind(x6_ia[,wanted_cols], x8_ia[,wanted_cols], x13_ia[,wanted_cols], tc0_ia[,wanted_cols])

save(bdi, file="merged_temporal_contiguity_data.RData")


# look for each trial of ord1 in ord2
are_orders_shuffled <- function(ord1, ord2) {
  if(nrow(ord1)!=nrow(ord2)) return(FALSE)
  found = rep(FALSE, nrow(ord1))
  for(r in 1:nrow(ord1)) {
    for(s in 1:nrow(ord2)) {
      if(setequal(unlist(ord1[r,]), unlist(ord2[s,]))) found[r] = TRUE
    }
  }
  return(found)
}

are_orders_shuffled(x6_1, x6_2) # yes
are_orders_shuffled(x6_1, x6_3) # yes 
are_orders_shuffled(x6_1, x6_4) # yes
are_orders_shuffled(x6_1, x6_2) # yes
are_orders_shuffled(orig, x6_1) # yes
# x6 orders are all shuffles of each other, and of original 4x4

are_orders_shuffled(x8_1, x8_4) # yes
are_orders_shuffled(orig, x8_4) # yes
are_orders_shuffled(x8_1, x6_1) # yes
# x8_1 and x8_4 are shuffles of original 

# x13 is all distinct trials (not shuffled of any others)
are_orders_shuffled(x13_1olap2tr, x13_1olap3tr) # no
are_orders_shuffled(x13_1olap2tr, orig) # no
are_orders_shuffled(x13_2olap2tr, orig) # no
are_orders_shuffled(x13_2olap2tr, x13_1olap2tr) # no (but largely)
are_orders_shuffled(x13_2olap2tr, x13_1olap3tr) # no