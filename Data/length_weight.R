'%nin%' = Negate(`%in%`)
## Read in the data file with lengths -----------

fish = read.csv("../AFRP/Data/FISH_MEASUREMENT_LML.csv") %>%
  filter(!(SPECIES == "RS" & LENGTH == 8 & WEIGHT == 2))

## Manually enter the slimy sculpin data -----------

SS = data.frame(SPECIES = rep("SS", 6), 
                LENGTH = c(67,55,48,67,68,71),
                WEIGHT= c(3,2,1,3,4,5))


## Filter out just Little Moose fish ----------------
### And select for only those that have recorded weights
LML.fish =  fish %>%
  filter(grepl("LML",FISH_N), 
         is.na(WEIGHT)==F ) %>%
  select(SPECIES, LENGTH, WEIGHT) %>%
  rbind(SS) %>%
  mutate(LENGTH = LENGTH / 10) ## Convert mm to cm because fishbase coefficients are calculated with cm not mm

## Set a vector that contains all LML fish in it
lml.species = unique(LML.fish$SPECIES) 
## Set vectors to fill with the coefficients from the model
a = vector()
b = vector()

fish %>% 
  filter(SPECIES == "RS") %>%
  select(LENGTH, WEIGHT) %>%
  na.omit() %>%
  arrange(WEIGHT)


## For loop ------------
## Calculates the coefficients for each of the species in LML using the log log relationship

for(i in 1:length(lml.species)){
  
  taxa.data = LML.fish %>% 
    filter(SPECIES == lml.species[i]) %>%
    mutate(LENGTH = LENGTH / 1) 
  
  model = lm(log(taxa.data$WEIGHT) ~ log(taxa.data$LENGTH ))
  
  a[i] =   exp(coef(model)[1]) 
  b[i] <- coef(model)[2] 
  
}

## Create a data frame with the coefficients and the species names for weiht estimation

coef = data.frame(SPECIES= lml.species, coef_a = a, coef_b = b) %>% 
  filter(SPECIES != "BB")
coef = rbind(coef, data.frame(SPECIES = "BB", coef_a =0.02820 , coef_b = 2.760))

## Weight estimation formula ------------
weight.estimated = function(x, a, b){a * x^b}

weight.frame = fish %>% # go back to using the frame that has all fish w & w/o weights
  filter(SPECIES %in% lml.species) %>%
  filter(SPECIES %nin% c( "SPL","NRD", "RT")) %>% # filter out low abundance fish
  filter(grepl("LML",YSAMP_N)) %>%
  separate(YSAMP_N, into = c("gear", "water","YEAR", "ysamp") ) %>%
  select(SPECIES, LENGTH, WEIGHT, YEAR) %>%
  #rbind(SS) %>% 
  left_join(coef) %>%
  mutate(LENGTH_cm = LENGTH / 10) %>%
  mutate(weight_e = weight.estimated(x = LENGTH_cm, a = coef_a, b = coef_b)) 

weight.frame %>% filter(SPECIES == "CC")

#write.csv(weight.frame, "Data/weight_frame.csv")

weight.frame %>%
  ggplot(aes(y = LENGTH, x = as.numeric(weight_e), col = SPECIES)) + 
  geom_point(aes(y = LENGTH, x = WEIGHT), col = "black", alpha = .5, key_glyph = "rect") +
  geom_smooth() +
  facet_wrap(~SPECIES, scales = "free") + 
  theme_minimal() +
  theme(legend.position = "none")
