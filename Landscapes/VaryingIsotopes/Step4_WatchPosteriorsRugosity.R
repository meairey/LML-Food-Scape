library(tidyverse)
cat = read.csv("trying.csv")
rugosity.summary = cat[-1,] %>%
  pivot_longer(sdrtrough:sds_dat, names_to = "metric", values_to = "value")

rugosity.summary$pos %>% unique() %>% length()
bayes_cri <- function(x) {
  data.frame(
    y = mean(x),  # Mean of posterior samples
    ymin = quantile(x, 0.025),  # 2.5% quantile (Lower bound)
    ymax = quantile(x, 0.975)   # 97.5% quantile (Upper bound)
  )
}

rugosity.summary 

rugosity.summary %>% 
  ggplot(aes(x = as.factor(comm), y = value))+
  stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  )  + 
  facet_wrap(~metric, scales = "free_y")




cat = read.csv("trying.csv")
rugosity.summary = cat[-1,]

sds1 = rugosity.summary %>%
  filter(pos %in% c(1:100)) %>%
  filter(comm == 1)
sds2 = rugosity.summary %>%
  filter(pos %in% c(1:100)) %>%
  filter(comm == 2)
sds3 = rugosity.summary %>%
  filter(pos %in% c(1:100)) %>%
  filter(comm == 3)

mean(sds1$sdrtrough > sds3$sdrtrough)
mean(sds1$sdrtrough > sds2$sdrtrough)
mean(sds2$sdrtrough > sds3$sdrtrough)


mean(sds1$sds_dat > sds3$sds_dat)
mean(sds1$sds_dat > sds2$sds_dat)


mean(sds1$s10z_dat > sds2$s10z_dat) ## significant

