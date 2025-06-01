## Libraries and functions -----------------------------
library(tidyverse)
library(wesanderson)
bayes_cri <- function(x) {
  data.frame(
    y = mean(x),  # Mean of posterior samples
    ymin = quantile(x, 0.025),  # 2.5% quantile (Lower bound)
    ymax = quantile(x, 0.975)   # 97.5% quantile (Upper bound)
  )
}

## Data formatting


cat = read.csv("trying.csv") #%>% 
#  rugosity.summary %>%
 # rename("sdr" = "sdrtrough",
  #       "s10z" = "sds_dat",
   #      "sds" = "sds_dat")
rugosity.summary = cat[-1,] %>%
  pivot_longer(3:6, names_to = "metric", values_to = "value")


## How 

rugosity.summary$pos %>% unique() %>% length()

## Plot Landscape Metrics

rugosity.summary %>% 
  ggplot(aes(x = as.factor(comm), y = value, col = as.factor(comm)))+
  stat_summary(
    fun.data = bayes_cri,   # Use Bayesian credible interval function
    geom = "pointrange"
  )  + 
  facet_wrap(~metric, scales = "free_y", labeller = labeller(metric = c("s10z" = "Ten-point height",
                                                                        "sdr" = "Rugosity",
                                                                        "sds" = "Distinct \npeaks"))) + 
  theme_minimal(base_size = 15) +
  ylab("Landscape Heterogeneity Metric") + xlab("Period") +
  scale_x_discrete(labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O.")) +
  scale_color_manual( values = wes_palette("Darjeeling1", type = "discrete", n = 4),labels = c("1" = "P.R.", "2" = "P.I.", "3" = "M.O."))+
  theme(legend.position = "none")




#cat = read.csv("trying.csv")
rugosity.summary = cat[-1,]

sds1 = rugosity.summary %>%
  filter(pos %in% c(1:500)) %>%
  filter(comm == 1)
sds2 = rugosity.summary %>%
  filter(pos %in% c(1:500)) %>%
  filter(comm == 2)
sds3 = rugosity.summary %>%
  filter(pos %in% c(1:500)) %>%
  filter(comm == 3)

## SDR
mean(sds1$sdr_dat > sds3$sdr_dat)
mean(sds1$sdr_dat > sds2$sdr_dat)
mean(sds2$sdr_dat > sds3$sdr_dat)


## SDS

mean(sds1$sds_dat > sds3$sds_dat)
mean(sds1$sds_dat > sds2$sds_dat)
mean(sds3$sds_dat > sds2$sds_dat)

## s10z
mean(sds1$s10z > sds2$s10z) 
mean(sds2$s10z > sds3$s10z)
mean(sds1$s10z > sds3$s10z)

## scl
mean(sds1$scl > sds2$scl) 
mean(sds2$scl > sds3$scl)
mean(sds1$scl > sds3$scl)








 
