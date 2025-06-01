library(SIBER)
h1 <- 1
load( file = "Data/VaryingIsotopesData/data.siber.RData")
load( file = "Data/VaryingIsotopesData/data.siber.RData")
legend.ellips = legend
dat1 <- data_setup(data_siber %>% filter(group != 12), h1)
spp1 <- length(names(dat1[[2]]))

dat.corr1 <- dat1[[3]] %>%
  left_join(baseline_simmr %>%
              mutate(community = as.numeric(as.factor(community))) %>%
              filter(GROUP == "B")) %>%
  mutate(D13C_c = iso1 - mean_d13C,
         D15N_c = (iso2 - mean_d15N)/3.4 + 1) %>%
  left_join(legend)

ellip1 <- array(0, dim = c(n.points, 2, spp1))
for(i in 1:spp1){
  ellip1[,,i] <- ellipse::ellipse(
    x = matrix(c(mean(dat1[[2]][[i]][,1]),
                 mean(dat1[[2]][[i]][,2]),
                 mean(dat1[[2]][[i]][,3]),
                 mean(dat1[[2]][[i]][,4])), 2, 2),
    centre = c(mean(dat1[[2]][[i]][,5]), mean(dat1[[2]][[i]][,6])),
    level = .4,
    npoints = n.points
  )
}

mean_d13C1 <- baseline_simmr %>%
  mutate(community = as.numeric(as.factor(community))) %>%
  filter(GROUP == "B", community == h1) %>%
  pull(mean_d13C)

mean_d15N1 <- baseline_simmr %>%
  mutate(community = as.numeric(as.factor(community))) %>%
  filter(GROUP == "B", community == h1) %>%
  pull(mean_d15N)

x.vec1 <- as.vector(ellip1[,1,]) - mean_d13C1
y.vec1 <- (as.vector(ellip1[,2,]) - mean_d15N1)/3.4 + 1

p1 <- ggplot() +
  geom_point(data = dat.corr1, aes(x = D13C_c, y = D15N_c, color = as.factor(group)), alpha = 0.25) +
  geom_path(aes(x = x.vec1, y = y.vec1, color = rep(as.factor(unique(dat1[[3]]$group)), each = n.points)),
            linewidth = 1) +
  scale_color_manual(
    values = (legend.ellips %>% filter(CODE %in% dat.corr1$CODE, community == h1))$color,
    labels = (legend.ellips %>% filter(CODE %in% dat.corr1$CODE, community == h1))$scientific,
    name = "Species"
  ) +
  theme_minimal(base_size = 14) +
  ylab("Trophic Position") + xlab("δ13C (corrected)") +
  ylim(1.75, 4.25) + xlim(-11, 5)



## Second graph

h2 <- 2
dat2 <- data_setup(data_siber %>% filter(group != 12), h2)
spp2 <- length(names(dat2[[2]]))

dat.corr2 <- dat2[[3]] %>%
  left_join(baseline_simmr %>%
              mutate(community = as.numeric(as.factor(community))) %>%
              filter(GROUP == "B")) %>%
  mutate(D13C_c = iso1 - mean_d13C,
         D15N_c = (iso2 - mean_d15N)/3.4 + 1) %>%
  left_join(legend)

ellip2 <- array(0, dim = c(n.points, 2, spp2))
for(i in 1:spp2){
  ellip2[,,i] <- ellipse::ellipse(
    x = matrix(c(mean(dat2[[2]][[i]][,1]),
                 mean(dat2[[2]][[i]][,2]),
                 mean(dat2[[2]][[i]][,3]),
                 mean(dat2[[2]][[i]][,4])), 2, 2),
    centre = c(mean(dat2[[2]][[i]][,5]), mean(dat2[[2]][[i]][,6])),
    level = .4,
    npoints = n.points
  )
}

mean_d13C2 <- baseline_simmr %>%
  mutate(community = as.numeric(as.factor(community))) %>%
  filter(GROUP == "B", community == h2) %>%
  pull(mean_d13C)

mean_d15N2 <- baseline_simmr %>%
  mutate(community = as.numeric(as.factor(community))) %>%
  filter(GROUP == "B", community == h2) %>%
  pull(mean_d15N)

x.vec2 <- as.vector(ellip2[,1,]) - mean_d13C2
y.vec2 <- (as.vector(ellip2[,2,]) - mean_d15N2)/3.4 + 1

p2 <- ggplot() +
  geom_point(data = dat.corr2, aes(x = D13C_c, y = D15N_c, color = as.factor(group)), alpha = 0.25) +
  geom_path(aes(x = x.vec2, y = y.vec2, color = rep(as.factor(unique(dat2[[3]]$group)), each = n.points)),
            linewidth = 1) +
  scale_color_manual(
    values = (legend.ellips %>% filter(CODE %in% dat.corr2$CODE, community == h2))$color,
    labels = (legend.ellips %>% filter(CODE %in% dat.corr2$CODE, community == h2))$scientific,
    name = "Species"
  ) +
  theme_minimal(base_size = 14) +
  ylab("Trophic Position") + xlab("δ13C (corrected)") +
  ylim(1.75, 4.25) + xlim(-11, 5)

p2


grid.arrange(p1, p2, ncol = 2)





























































