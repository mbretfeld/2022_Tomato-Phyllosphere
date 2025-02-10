#### Non-microbial Tomato Data Models ####
#### Models that previously used random effects now include a second "nonrandom
#### model that treats Plant ID as a fixed effect 

require("tidyverse"); require("ggsignif"); require("lmerTest")
require("lubridate"); require("rstatix"); require("GLMMadaptive"); require("RColorBrewer")

ggtheme <- theme_bw() +
  theme(panel.grid = element_blank())

####### LOAD DATA #####################

# Trial 2 Fluorescence Data 
raw_fluorescence_t2 <- read.csv(file = "round2multispecdata.csv", stringsAsFactors = T, fileEncoding="UTF-8-BOM") %>%  filter(Series == "Whole Data Set") %>%
  mutate(DateTime = parse_date_time(time, orders = c("mdy H:M")),
         Sample = as.factor(paste0(Row, Pot.ID)),
         Salinity = as.factor(case_when(Sample == "A1" | Sample == "A3" | Sample == "A5" | Sample == "A7" |
                                          Sample == "B2" | Sample == "B4" | Sample == "B6" | Sample == "B8" |
                                          Sample == "C1" | Sample == "C3" | Sample == "C5" | Sample == "C7" |
                                          Sample == "D2" | Sample == "D4" | Sample == "D6" | Sample == "D8"
                                        ~ "No Salt",
                                        TRUE ~ "Salt")),
         Inoculation = as.factor(case_when(Row == "C" | Row == "D" ~ "Inoculated",
                                           TRUE ~ "Not Inoculated")), 
         Treatment = as.factor(paste(Inoculation, Salinity)),
         Date = as.Date(DateTime))

data_fluor_t2 <- raw_fluorescence_t2 %>% 
  filter(FvP_over_FmP > 0.45) %>%
  filter(DateTime >= "2022-09-13") %>%
  mutate(DOY = yday(DateTime),
         DateFac = as.factor(DateTime),
         LogitPhi2 = qlogis(Phi2),  # you can backtransform using plogis()
         LogitFvP_over_FmP = qlogis(FvP_over_FmP),
         scaledlight = (Light.Intensity..PAR./1000),
         FluorType = case_when(Leaflet.Age == "New (Top 3 Leaves)" 
                               ~ "Dark",
                               TRUE ~ "Light"))


#### Trial 2 Tomato Harvest Data - includes fruit weight and sugar data 
raw_yield_t2 <- read.csv(file = "Round2TomatoFruitData.csv", header = TRUE, stringsAsFactors = TRUE) %>% 
  mutate(Pot = as.factor(paste0(Row, Pot)), 
         Salinity = as.factor(case_when(Pot == "A1" | Pot == "A3" | Pot == "A5" | Pot == "A7" |
                                          Pot == "B2" | Pot == "B4" | Pot == "B6" | Pot == "B8" |
                                          Pot == "C1" | Pot == "C3" | Pot == "C5" | Pot == "C7" |
                                          Pot == "D2" | Pot == "D4" | Pot == "D6" | Pot == "D8"
                                        ~ "No Salt",
                                        TRUE ~ "Salt")),
         Inoculation = as.factor(case_when(Row == "C" | Row == "D" ~ "Inoculated",
                                           TRUE ~ "Not Inoculated")), 
         Treatment = as.factor(paste(Inoculation, Salinity))) %>%
  rename(BER = BER.)


Sugar<-raw_yield_t2 %>% filter(BER == "N" & SugarAverage > 0)


### I created a dataframes for analyzing all tomatoes regardless of BER status 
#### and the proportion of BER tomatoes per plant ###

data_yield_t2_total <- raw_yield_t2 %>% 
  group_by(Salinity, Inoculation, Treatment, Pot) %>% 
  tally()

data_yield_t2 <- raw_yield_t2 %>% 
  group_by(Salinity, Inoculation, BER, Treatment, Pot) %>% 
  tally()

data_yield_t2_wide <- data_yield_t2 %>% 
  pivot_wider(names_from = BER, values_from = n, values_fill = 0) %>%
  mutate(Prop_BER = Y/(Y+N))

### Trial 2 Takedown of Tomato Plants Data - Includes plant shoot biomass 
#### and remaining fruits on plant
raw_yield_t2_takedown <- read.csv(file = "Round2Biomass&HarvestData.csv", header = TRUE, stringsAsFactors = T) %>% 
  mutate(Salinity = as.factor(case_when(Plant == "A1" | Plant == "A3" | Plant == "A5" | Plant == "A7" |
                                          Plant == "B2" | Plant == "B4" | Plant == "B6" | Plant == "B8" |
                                          Plant == "C1" | Plant == "C3" | Plant == "C5" | Plant == "C7" |
                                          Plant == "D2" | Plant == "D4" | Plant == "D6" | Plant == "D8"
                                        ~ "No Salt",
                                        TRUE ~ "Salt")),
         Inoculation = as.factor(case_when(Inoculation == "Y" ~ "Inoculated",
                                           TRUE ~ "Not Inoculated")), 
         Treatment = as.factor(paste(Inoculation, Salinity))) %>%
  rename(n = Fruit.Count)




### Phi 2 Models ##### 
summary(mod_Phi2_t2_random <- (lmer(
  Phi2 ~ Inoculation * Salinity + Light.Intensity..PAR. + (1|Sample) ,
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

summary(mod_Phi2_t2_norandom <- (lm(
  Phi2 ~ Inoculation * Salinity + Light.Intensity..PAR. + Sample ,
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))


plot(density(residuals(mod_Phi2_t2_random)))
plot(density(residuals(mod_Phi2_t2_norandom)))

anova(mod_Phi2_t2.2, type = "III")


### Fv'/Fm' Models #### - Removing random effect makes inoculation non-significant 
summary(mod_FvP_FmP_t2_random <- (lmer(
  FvP_over_FmP ~ Inoculation * Salinity + Light.Intensity..PAR. + (1|Sample),
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

summary(mod_FvP_FmP_t2_norandom <- (lm(
  FvP_over_FmP ~ Inoculation * Salinity + Light.Intensity..PAR. + Sample,
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

plot(density(residuals(mod_FvP_FmP_t2_random)))
plot(density(residuals(mod_FvP_FmP_t2_norandom)))

#### NPQt ###### - Removing random effect makes inoculation non-significant

summary(mod_NPQt_t2_random <- (glmer(
  NPQt ~ Inoculation * Salinity + scaledlight + (1|Sample), family = gaussian(link = "log"),
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

summary(mod_NPQt_t2_norandom <- (glm(
  NPQt ~ Inoculation * Salinity + scaledlight + Sample, family = gaussian(link = "log"),
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

plot(density(residuals(mod_NPQt_t2_random)))
plot(density(residuals(mod_NPQt_t2_norandom)))

                                                                                                                                              
##SPAD - nothing changes 

summary(mod_SPAD_t2_random <- (lmer(
  SPAD ~ Inoculation * Salinity + (1|Sample),
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

summary(mod_SPAD_t2_norandom <- (lm(
  SPAD ~ Inoculation * Salinity + Sample,
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

plot(density(residuals(mod_SPAD_t2_random)))
plot(density(residuals(mod_SPAD_t2_norandom)))


## LEF Test - nothing changes 

summary(mod_LEF_t2_random <- (lmer(
  LEF ~ Inoculation * Salinity + Light.Intensity..PAR. + (1|Sample),
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

summary(mod_LEF_t2_norandom <- (lm(
  LEF ~ Inoculation * Salinity + Light.Intensity..PAR. + Sample,
  data = data_fluor_t2 %>% filter(FluorType == "Light"))))

plot(density(residuals(mod_LEF_t2_random)))
plot(density(residuals(mod_LEF_t2_norandom)))

### qL ###### - Removing random effect makes inoculation non-significant

summary(mod_qL_t2_random <- (lmer(
  (qL) ~ Inoculation * Salinity + Light.Intensity..PAR.+ (1|Sample),
  data = data_fluor_t2 %>% filter(FluorType == "Light") )))

summary(mod_qL_t1_norandom <- (lm(
  (qL) ~ Inoculation * Salinity + Light.Intensity..PAR.+ Sample,
  data = data_fluor_t2 %>% filter(FluorType == "Light") )))

#####Total Number of Fruits (Included BER and no BER)########
#### This wasn't a model that used random effects #####

data_yield_t2_total

summary(aov_totalfruit_t2.1 <- aov(n ~ Salinity * Inoculation, 
                                   data = data_yield_t2_total))

anova(aov_totalfruit_t2.1)


####Fruit Sugar - I am unsure how to go about this since I used two random effects 
#### the original model. Both pot and "logged" (the date we weighed our ripe tomatoes)
#### were fairly small in number of entries. However, I have created a model that keeps 
#### "logged" as a random effect

#### I get an error running the second model - so maybe were over-parameterizing? ###

summary(sugarmod_t2_random <-lmer(SugarAverage ~ Salinity * Inoculation + (1|Pot) + (1|Logged), 
                        data = Sugar))

summary(sugarmod_t2_random <-lmer(SugarAverage ~ Salinity * Inoculation + Pot + (1|Logged), 
                                  data = Sugar))



## Proportion of BER per plant - this was not a model I used random effects for ###

summary(propmod1<-glm(Prop_BER ~ Salinity * Inoculation, data = data_yield_t2_wide))


#  (plant weights) - significant
### Plant weight take a lognormal distribution hence the transformation 
#### This was not a model I used random effects for ##### 


summary(shootmass1<- glm(PlantWeightgrams ~ Salinity * Inoculation, family = gaussian(link = "log"),
                         data = raw_yield_t2_takedown))

anova(shootmass1)


