# source("preprocessing.R") 
require(lme4)
require(tidyverse)
require(ggplot2)
require(tidyboot)
load("merged_temporal_contiguity_data.RData")

bd = bdi %>% group_by(cond, item, itemSR, condSR) %>% 
  summarize(accuracy = mean(Correct)) 

bd %>% group_by(cond) %>% 
  summarize(mean_acc = mean(accuracy)) %>% 
  arrange(mean_acc)

ggplot(data=bd, aes(x=jitter(itemSR), y=accuracy, group=cond, colour=cond)) + 
  geom_point(alpha=.7) + theme_bw() + ylim(0,1) + xlab("Item SR")

cor.test(bd$accuracy, bd$itemSR) # .35 p<.001
cor.test(bd$accuracy, bd$condSR) # .35 p<.001

bdi$itemSRsc = scale(bdi$itemSR)
bdi$condSRsc = scale(bdi$condSR)
summary(glmer(Correct ~ itemSRsc*condSRsc + (itemSR|Subject/cond), data=bdi, family=binomial)) 
summary(glmer(Correct ~ itemSRsc*condSRsc + (1|Subject), data=bdi, family=binomial)) 
# main effects of itemSR (.19) and itemSR:condSR (-.12), marginal of condSRsc .10

# exclude spatial+temporal contiguity condition? subset(bdi, cond!="x8_1")
summary(glmer(Correct ~ itemSR*condSR + (itemSR|Subject/cond), data=bdi, family=binomial)) # /cond or not?
summary(glmer(Correct ~ itemSR*condSR + (1|Subject), data=bdi, family=binomial)) # /cond or not?
#               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)    -1.01296    0.11150  -9.085  < 2e-16 ***
# itemSR         0.25965    0.07403   3.507 0.000452 ***
# condSR         0.33487    0.12803   2.616 0.008908 ** 
# itemSR:condSR -0.14077    0.05732  -2.456 0.014055 *  
# Significant positive effects of item SR (beta=0.26, z=3.51, p<.001) and
# of average trial-to-trial SR in the condition (beta=0.33, z=2.62, p<.01), showing that more
# successive repetitions increase accuracy of both individual items and all items in the condition.
# However, there was also a significant negative interaction of item SR and condition SR 
# (beta=-.14, z=2.46, p<.05), showing that in conditions with more repeated items there was 
# less of a benefit for being an item with many repetitions.


# graph the comparable items (from reshuffled orders)
ggplot(data=subset(bd, cond=="tc0" | cond=="x6_1" | cond=="x6_2" | cond=="x6_3" | cond=="x6_4" | cond=="x8_1" | cond=="x8_4"), aes(x=item, y=accuracy, group=cond, colour=cond)) + 
  geom_point() + theme_bw() + ylim(0,1)

# sum of squared accuracy differences between equivalent items in different shufflings;
var(subset(bd, cond=='tc0')$accuracy) # .004
var(subset(bd, cond=='x6_1')$accuracy) # .02
var(subset(bd, cond=='x6_2')$accuracy) # .014
sum((subset(bd, cond=='tc0')$accuracy - subset(bd, cond=='x6_1')$accuracy)^2) # .57
sum((subset(bd, cond=='tc0')$accuracy - subset(bd, cond=='x6_2')$accuracy)^2) # .46
sum((subset(bd, cond=='tc0')$accuracy - subset(bd, cond=='x6_3')$accuracy)^2) # .24
sum((subset(bd, cond=='tc0')$accuracy - subset(bd, cond=='x6_4')$accuracy)^2) # .16
sum((subset(bd, cond=='x6_1')$accuracy - subset(bd, cond=='x6_2')$accuracy)^2) # .58
sum((subset(bd, cond=='x6_1')$accuracy - subset(bd, cond=='x6_3')$accuracy)^2) # .64
sum((subset(bd, cond=='x6_2')$accuracy - subset(bd, cond=='x6_3')$accuracy)^2) # .51
sum((subset(bd, cond=='x6_2')$accuracy - subset(bd, cond=='x6_4')$accuracy)^2) # .53

# could use these shufflings together to test strength of recency, primacy, context familiarity

# can compare to bootstrapped shufflings of the items, e.g.:
sum((subset(bd, cond=='tc0')$accuracy - sample(subset(bd, cond=='x6_2')$accuracy))^2) # .32...76

# x8 has 2 maxTC conditions -- 1_ spatial (+temporal) contiguity, 4_ temporal (no spatial) contiguity
ggplot(data=subset(bd, cond=="x8_1" | cond=="x8_4"), aes(x=item, y=accuracy, group=cond, colour=cond)) + 
  geom_point() + theme_bw() + ylim(0,1)
cor.test(subset(bd, cond=="x8_1")$accuracy, subset(bd, cond=="x8_4")$accuracy) # cor=.03 p=.89
# item-level accuracy is uncorrelated (wow)
# Q: does spatial contiguity matter beyond temporal contiguity?
t.test(subset(bd, cond=="x8_1")$accuracy - subset(bd, cond=="x8_4")$accuracy)
# one-sample t(17) = 2.77, diff=.084, p=.01 
# Temp+Spat Contiguity = .49, Temp Contiguity = .4
# (could look at difference between TC and non-TC items)