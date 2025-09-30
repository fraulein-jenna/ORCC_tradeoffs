# Analytical code to accompany the manuscript 
#"Morphological plasticity in a reef-building coral is context-dependent and trades off with resistance to thermal stress"
# by Dilworth et al. (2025)

# Author: Jenna Dilworth

#Set up ####
# set working directory and load necessary packages 
setwd("~/Desktop/pathtoyourdirectory")

library(raster)
library(sf)
library(geodata)
library(mapproj)
library(ggrepel)
library(ggspatial)
library(tidyverse)
library(readxl)
library(janitor)
library(gtools)
library(survminer)
library(RColorBrewer)
library(survival)
library(coxme)
library(forestmodel)
library(car)
library(lubridate)
library(lme4)
library(lmerTest)
library(emmeans)
library(cowplot)
library(vegan)
library(ggthemes)

# Background information ####

##site map #####
sites<-read_xlsx("site_coords.xlsx")
shapefile_path <- "~/Desktop/USC/Research_Projects/ORCC_Outplant/Analysis/gadm41_USA_shp/gadm41_USA_1.shp"

# Read the shapefile using sf
gadm_states <- st_read(shapefile_path)
# Filter for a specific state (e.g., California)
florida <- gadm_states[gadm_states$NAME_1 == "Florida", ]
Fl_sf <- st_as_sf(florida)

#create map for Figure 1a
site_map <- ggplot(Fl_sf) +
  geom_sf(fill = "#F7F7F7") +
  geom_point(data = sites, aes(x = Lon, y = Lat),
             color = ifelse(sites$Name == "Looe Key Nursery", "#808080", "#246BAE"),
             shape = ifelse(sites$Name == "Looe Key Nursery", 8, 16),
             size =5)+
  geom_text_repel(data = sites, aes(x = Lon, y = Lat, label = Name), size = 5)+
  xlim(-81.7, -81.22) +
  ylim(24.38, 24.75) +
  scale_color_identity() +
  scale_size_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#D1EDF2"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
  ) +
  annotation_scale(location = "bl", width_hint = 0.3,
                   height = unit(0.15, "in"),
                   text_cex = 1.3
  ) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(0.3, "in"),
                         height = unit(0.8, "in"),
                         width = unit(0.8, "in"),
                         style = north_arrow_fancy_orienteering())
site_map

##temperature data ####
looe.temps <-read.csv("Carly_Looe_2023.09.12.csv")%>%
  clean_names()%>%
  dplyr::rename(datetime=date_time_gmt_04_00)%>%
  dplyr::rename(temperature=temp_c_lgr_s_n_20265724_sen_s_n_20265724)%>%
  dplyr::select(datetime, temperature)%>%
  separate(datetime, into=c("date", "time"), sep= " ")%>%
  dplyr::select(-time)%>%
  mutate(date = gsub("/", "-", date))%>%
  mutate(date=mdy(date))%>%
  filter(date < as.Date("2023-07-30"))%>%
  group_by(date)%>%
  dplyr::summarise(daily_average=mean(temperature), 
                   sd = sd(temperature),
                   .groups="keep")%>%
  na.omit()%>%
  mutate(site="Looe Key")

daves.temps <-read.csv("CK_Dave's_Ledge_2023.09.12.csv")%>%
  clean_names()%>%
  dplyr::rename(datetime=date_time_gmt_04_00)%>%
  dplyr::rename(temperature=temp_c_lgr_s_n_20361634_sen_s_n_20361634)%>%
  dplyr::select(datetime, temperature)%>%
  separate(datetime, into=c("date", "time"), sep= " ")%>%
  dplyr::select(-time)%>%
  mutate(date = gsub("/", "-", date))%>%
  mutate(date=mdy(date))%>%
  filter(date < as.Date("2023-07-30"))%>%
  group_by(date)%>%
  dplyr::summarise(daily_average=mean(temperature), 
                   sd = sd(temperature),
                   .groups="keep")%>%
  na.omit()%>%
  mutate(site="Dave's Ledge")

field_temps <- full_join(looe.temps, daves.temps)

# ambient temperatures for figure 3a
overall_temperatures <- ggplot ()+
  scale_color_manual(values=c("tomato2","steelblue2"))+
  scale_fill_manual(values=c("tomato2","steelblue2"))+
  geom_line(data = field_temps, aes(x=date, y = daily_average, color = site))+
  geom_ribbon(data= field_temps, aes(x=date,y = daily_average, ymin =daily_average - sd, ymax = daily_average + sd, fill = site), alpha = .2) +
  theme_light()+
  scale_x_date(lim = c(as.Date("2022-10-14"),as.Date("2023-06-29")), breaks = "1 month", date_labels = "%b", expand = c(0,0))+
  scale_y_continuous(breaks = c(20, 22, 24, 26, 28,30,32), expand = c(0,0))+
  ylab("Temperature (C)")+
  xlab("Month")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))
overall_temperatures

# bleaching event temperatures for figure 4a
bleaching_temps <- field_temps%>%
  filter(date < as.Date("2023-07-31") & date > as.Date("2023-06-01"))

bleaching_temperatures <- ggplot ()+
  scale_color_manual(values=c("tomato2","steelblue2"))+
  scale_fill_manual(values=c("tomato2","steelblue2"))+
  geom_hline(aes(yintercept = 29.63),color = "grey")+
  geom_hline(aes(yintercept = 30.63), linetype = 2,color = "grey")+
  geom_vline(aes(xintercept = as.Date("2023-06-29")), color = "black")+
  geom_vline(aes(xintercept = as.Date("2023-07-26")), color = "black")+
  geom_line(data = bleaching_temps, aes(x=date, y = daily_average, color = site))+
  geom_ribbon(data= bleaching_temps, aes(x=date,y = daily_average, ymin =daily_average - sd, ymax = daily_average + sd, fill = site), alpha = .2) +
  theme_bw()+
  scale_x_date(breaks = "2 weeks", date_labels = "%b %d", expand = c(0,0))+
  scale_y_continuous(lim = c(27,33.5),breaks = c(20, 22, 24, 26, 28,30,32), expand = c(0,0))+
  ylab("Temperature (C)")+
  xlab("Date")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))
bleaching_temperatures

#calculate experimental degree heating weeks accumulated at each site
#following Leggat et al. 2022
# first, calculate daily hotspots
daily_HS <- bleaching_temps%>%
  mutate(eHS = daily_average - 29.63)%>%
  filter(eHS >0)

looe <- daily_HS%>%
  filter(site == "Looe Key")%>%
  filter(eHS >=1)

looe_DHWs = sum(looe$eHS)/7
# 10.32354

daves <- daily_HS%>%
  filter(site == "Dave's Ledge")%>%
  filter(eHS >=1)

daves_DHWs = sum(daves$eHS)/7
# 10.9567

##survivorship ####

#read in mortality data
surv<-read_xlsx("Mortality/ORCC_bleaching_mortality_Jul23_QCed.xlsx", sheet = "Mortality_binary")%>%
  clean_names()%>%
  dplyr::select(site, genotype, tag_number, replicate, jun_23_surv, jul_23_surv, tod_jun, tod_jul, nov_23, tod_nov)%>%
  dplyr::rename(jun_23=jun_23_surv)%>%
  dplyr::rename(jul_23=jul_23_surv)%>%
  mutate_at("site", as.factor)%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("tag_number", as.factor)%>%
  mutate_at("replicate", as.factor)%>%
  dplyr::mutate(across(ends_with("23"), as.numeric))%>%
  dplyr::mutate(across(starts_with("tod"), as.numeric))%>%
  na.omit()

#AMBIENT SURVIVAL
ambient.died=Surv(surv$tod_jun, surv$jun_23)

#Kaplan-Meier survival curves to look at overall trends
ambient.scurveG <- survfit(ambient.died~genotype, data=surv) #function cannot handle mixed effects or interactions

colorG <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)
ggsurvplot(ambient.scurveG, data = surv,
           palette=colorG,fun = "pct", conf.int = FALSE,legend="right",legend.title="Genotype",
           ggtheme = theme_bw(),
           ylim = c(40,100), 
           break.x.by = 3, xlab = "Months")

ambient.res<-data.frame(summary(ambient.scurveG)$table)
genet<-str_split(rownames(ambient.res),"=",simplify=TRUE)[,2]
df<-data.frame(ambient.res,genet,stringsAsFactors=TRUE)

df$n.surv<-df$n.start-df$events
df$PercSurv<-df$n.surv/df$n.start
rownames(df)<-seq_len(nrow(df))
ambient.surv.sort<-df[order(df$genet),]

#reorder survival dataframe according to Kaplan-Meier results:
surv$genotype<-factor(surv$genotype,levels=c("36","41","3","50","13","1","31","44","7","62"))

## Cox proportional hazards - is there an effect of genotype under ambient temperature conditions?
ambient.genotype.cox<-coxph(ambient.died~genotype, data=surv)

ambient.cox.rank <- data.frame(ambient.genotype.cox$coefficients)%>%
  rownames_to_column()%>%
  dplyr::rename(genotype=rowname)%>%
  rename(ambient_cox = ambient.genotype.cox.coefficients)

# plot results for figure 3b
print(forest_model(coxph(ambient.died~genotype, data=surv)))

#creating survival curve for Nov - illustrating that most corals ended up dying - supplemental figure
# need to reorder genotypes in ascending order for color palette
surv$genotype<-factor(surv$genotype,levels=c("1","3","7","13","31","36","41","44","50","62"))

overall.died=Surv(surv$tod_nov, surv$nov_23)

overall.scurveG <- survfit(overall.died~genotype, data=surv) #function cannot handle mixed effects or interactions

colorG <-  c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF")
ggsurvplot(overall.scurveG, data = surv,
           palette=colorG,fun = "pct", conf.int = FALSE,legend="right",legend.title="Genotype",
           ggtheme = theme_bw(),
           #ylim = c(50,100), 
           break.x.by = 3, xlab = "Months")

#BLEACHING SURVIVAL - only until July time point - after that we have too much mortality for these models
bleach.died=Surv(surv$tod_jul, surv$jul_23)

###### Plotting survival curve
bleach.scurveG <- survfit(bleach.died~genotype, data=surv) #function cannot handle mixed effects or interactions

colorG <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)
ggsurvplot(bleach.scurveG, data = surv,
           palette=colorG,fun = "pct", conf.int = FALSE,legend="right",legend.title="Genotype",
           ggtheme = theme_bw(),
           #ylim = c(50,100), 
           break.x.by = 3, xlab = "Months")

bleach.res<-data.frame(summary(bleach.scurveG)$table)
genet<-str_split(rownames(bleach.res),"=",simplify=TRUE)[,2]
df<-data.frame(bleach.res,genet,stringsAsFactors=TRUE)

df$n.surv<-df$n.start-df$events
df$PercSurv<-df$n.surv/df$n.start
rownames(df)<-seq_len(nrow(df))
bleach.surv.sort<-df[order(df$genet),]

# growth data ####
# read in mortality datasheet for outplant date info and other frag metadata
frag_info <- read_xlsx("Mortality/ORCC_bleaching_mortality_Jul23_QCed.xlsx", sheet = "Mortality_binary")%>%
  clean_names()%>%
  dplyr::select(site, genotype, replicate, tag_number, outplant_date)%>%
  mutate_at("genotype", as.factor)

##October 2022 outplants ####
oct_dl_data <- read_xlsx("Morphology/Oct2022_AcerMorphology.xlsx", sheet = "Daves Ledge")
oct_lk_data <- read_xlsx("Morphology/Oct2022_AcerMorphology.xlsx", sheet = "Looe")

oct_data <- full_join(oct_dl_data, oct_lk_data)%>%
  clean_names()%>%
  dplyr::select(timepoint,array, genotype, number_of_branches,true_sa, true_volume, true_tle)%>%
  mutate_at("genotype", as.factor)%>%
  rename(number_branches = number_of_branches)%>%
  rename(sa = true_sa)%>%
  rename(volume = true_volume)%>%
  rename(tle = true_tle)%>%
  left_join(., frag_info, by = c("array"="replicate","genotype"="genotype"))%>%
  mutate_at("tag_number", as.factor)%>%
  filter(outplant_date == "Oct")

##January 2023 outplants####
jan_dl_data <- read_xlsx("Morphology/Jan2023_AcerMorphology_RACKS.xlsx", sheet = "Daves Ledge")
jan_lk_data <- read_xlsx("Morphology/Jan2023_AcerMorphology_RACKS.xlsx", sheet = "Looe")

# filter out super tiny polyps that were counted as branches and calculate new TLE accounting for rack holder
jan_data <- full_join(jan_dl_data, jan_lk_data)%>%
  clean_names()%>%
  filter(!is.na(site))%>%
  dplyr::mutate(across(starts_with("tle"), ~replace(., . <0.5, "NA")))%>%
  dplyr::mutate(across(starts_with("tle"), as.numeric))%>%
  dplyr::mutate(tle_corrected = (rowSums(across(starts_with("tle")), na.rm=TRUE)+1.33))%>%
  dplyr::select(array_replicate, genotype, number_of_branches, true_sa, true_volume, tle_corrected)%>%
  rename(number_branches = number_of_branches)%>%
  rename(sa = true_sa)%>%
  rename(volume = true_volume)%>%
  rename(tle = tle_corrected)%>%
  rename(array = array_replicate)%>%
  mutate_at("genotype", as.factor)%>%
  left_join(., frag_info, by = c("array"="replicate","genotype"="genotype"))%>%
  mutate_at("tag_number", as.factor)%>%
  mutate(timepoint = "T1")

##January 2023 monitoring ####
jan2_dl_data <-read_xlsx("Morphology/Jan2023_AcerMorphology_v2.xlsx", sheet = "Daves Ledge")
jan2_lk_data <-read_xlsx("Morphology/Jan2023_AcerMorphology_v2.xlsx", sheet = "Looe Reef")

jan2_data <- full_join(jan2_dl_data, jan2_lk_data)%>%
  clean_names()%>%
  dplyr::mutate(across(starts_with("tle"), ~replace(., . <0.5, "NA")))%>%
  dplyr::mutate(across(starts_with("tle"), as.numeric))%>%
  dplyr::mutate(tle_new = rowSums(across(starts_with("tle")), na.rm=TRUE))%>%
  dplyr::select(site, timepoint, genotype, array, tag_number, breakage, status, number_branches, sa, volume, tle_new)%>%
  rename(tle = tle_new)%>%
  mutate_at("timepoint", as.factor)%>%
  mutate_at("genotype", as.factor)%>%
  left_join(., frag_info)%>%
  dplyr::select(-replicate)

##April 2023 monitoring ####
apr_dl_data <-read_xlsx("Morphology/April2023_AcerMorphology_v2.xlsx", sheet = "Daves Ledge")
apr_lk_data <-read_xlsx("Morphology/April2023_AcerMorphology_v2.xlsx", sheet = "Looe Reef")

apr_data <- full_join(apr_dl_data, apr_lk_data)%>%
  clean_names()%>%
  dplyr::mutate(across(starts_with("tle"), ~replace(., . <0.5, "NA")))%>%
  dplyr::mutate(across(starts_with("tle"), as.numeric))%>%
  dplyr::mutate(tle_new = rowSums(across(starts_with("tle")), na.rm=TRUE))%>%
  dplyr::select(site, timepoint, genotype, array, tag_number, breakage, status, number_branches, sa, volume, tle_new)%>%
  rename(tle = tle_new)%>%
  mutate_at("timepoint", as.factor)%>%
  mutate_at("genotype", as.factor)%>%
  left_join(., frag_info)%>%
  dplyr::select(-replicate)

## June 2023 monitoring ####
june_dl_data <-read_xlsx("Morphology/Jun23_AcerMorphology_v2.xlsx", sheet = "Daves Ledge")
june_lk_data <-read_xlsx("Morphology/Jun23_AcerMorphology_v2.xlsx", sheet = "Looe Reef")

june_data <- full_join(june_dl_data, june_lk_data)%>%
  clean_names()%>%
  dplyr::mutate(across(starts_with("tle"), ~replace(., . <0.5, "NA")))%>%
  dplyr::mutate(across(starts_with("tle"), as.numeric))%>%
  dplyr::mutate(tle_new = rowSums(across(starts_with("tle")), na.rm=TRUE))%>%
  dplyr::select(site, timepoint, genotype, replicate, tag_number, breakage, status, number_branches, sa, volume, tle_new)%>%
  rename(tle = tle_new)%>%
  rename(array = replicate)%>%
  mutate_at("timepoint", as.factor)%>%
  mutate_at("genotype", as.factor)%>%
  left_join(., frag_info)%>%
  dplyr::select(-replicate)

## all growth ####
# combine datasets and set initial timepoint as October or January depending on reset status
growth_working <- june_data%>%
  full_join(apr_data)%>%
  full_join(jan2_data)%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("tag_number", as.factor)

#split into Jan and Oct outplanting subsets to standardize time outplanted
oct_outplant <- growth_working %>%
  filter(outplant_date == "Oct")%>%
  full_join(oct_data)%>%
  mutate(months =case_when(
    timepoint =="T0"~0,
    timepoint == "T1"~3,
    timepoint =="T2"~6,
    timepoint == "3"~9))

jan_outplant <- growth_working %>%
  filter(outplant_date == "Jan")%>%
  full_join(jan_data)%>%
  mutate(months =case_when(
    timepoint =="T1"~0,
    timepoint =="T2"~3,
    timepoint == "3"~6))

#join back together to plot all of the growth data
growth_months <- full_join(oct_outplant, jan_outplant)%>%
  mutate_at("months", as.factor)%>%
  mutate(status = case_when(
    timepoint == "T0" ~"A",
    TRUE ~ as.character(status)))%>%
  filter(status == "A")%>%
  filter(!is.na(sa))%>%
  filter(outplant_date == "Oct")%>% # filter to remove those that were reset in Jan
  mutate(sav = sa/volume)%>%
  mutate_at("array", as.factor)

# initial sizes
initial.size <- growth_months%>%
  filter(timepoint == "T0")

mean(initial.size$tle)
sd(initial.size$tle)

#ANOVAs to check if there are initial diffs in any of the traits of interest
#TLE
res_aov <-aov(tle ~ genotype*site, data = initial.size)

#checking normality
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
#checking homogeneity of variance
leveneTest(tle ~ genotype*site, data = initial.size) # p = 0.4699

summary(res_aov) # significant effect of genotype

#SA
res_aov <-aov(sa ~ genotype*site, data = initial.size)

#checking normality
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
#checking homogeneity of variance
leveneTest(sa ~ genotype*site, data = initial.size) # p = 0.6725

summary(res_aov) # significant effect of genotype

#Volume
res_aov <-aov(volume ~ genotype*site, data = initial.size)

#checking normality
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
#checking homogeneity of variance
leveneTest(volume ~ genotype*site, data = initial.size) # p = 0.1847

summary(res_aov) # significant effect of genotype

#SA:V
res_aov <-aov(sav ~ genotype*site, data = initial.size)

#checking normality
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
#checking homogeneity of variance
leveneTest(sav ~ genotype*site, data = initial.size) # p =0.2101

summary(res_aov) # significant effect of genotype

# growth - change in size over each time interval
# for each of the 3-month sampling intervals
growth_TLE <- growth_months%>%
  dplyr::select(site, genotype, array, tag_number, tle, timepoint)%>%
  pivot_wider(names_from = timepoint, values_from = tle)%>%
  rename(T3="3")%>%
  mutate(growth_3 = (T1-T0))%>%
  mutate(growth_6 = (T2-T1))%>%
  mutate(growth_9 = (T3-T2))%>%
  pivot_longer(cols = 9:11, names_to = "timepoint", values_to = "growth")%>%
  dplyr::select(site, genotype, timepoint, tag_number, timepoint, growth)

growth_TLE$timepoint <- sub("growth_3", "T1", growth_TLE$timepoint)
growth_TLE$timepoint <- sub("growth_6", "T2", growth_TLE$timepoint)
growth_TLE$timepoint <- sub("growth_9", "3", growth_TLE$timepoint)

# breakage - info on severity/branch ordination
breakage <- growth_months%>%
  select(site, timepoint, genotype, tag_number,breakage)%>%
  separate_rows(breakage)%>%
  mutate(value = 1)%>%
  pivot_wider(names_from = breakage, values_fill = 0)%>%
  mutate(secondary = rowSums(across(starts_with("S", ignore.case=FALSE))))%>%
  select(site, timepoint, genotype, tag_number, "NA", P, secondary, T1, C)%>%
  rename("none"="NA")%>%
  mutate(breakage = case_when(
    none == 1 ~ 0,
    (P + secondary + T1) == 1 ~ 1,
    (P + secondary + T1) == 2 ~ 2,
    (P + secondary + T1) == 3 ~ 3,
    (P + secondary + T1) == 4 ~ 4,
    C == 1 ~ 5
  ))%>%
  select(site, genotype, timepoint, tag_number, breakage)

# create data frame which includes growth rates, breakage, and sizes over all timepoints
all_morphology <- growth_months%>%
  select(-breakage, -status, -months, - outplant_date)%>%
  left_join(breakage)%>%
  left_join(growth_TLE)

## plot growth over time ####
#TLE by genotype
tle.data.geno<- plyr::ddply(growth_months, c("genotype","site","months"), summarise,
                            N    = length(tle[!is.na(tle)]),
                            mean = mean(tle, na.rm=TRUE),
                            sd   = sd(tle, na.rm=TRUE),
                            se   = sd / sqrt(N))

tle.data.geno <- tle.data.geno%>%
  mutate(timepoint = case_when(
    months == 0 ~ "October",
    months == 3 ~ "January",
    months == 6 ~ "April",
    months == 9 ~ "June"
  ))

tle.data.geno$genotype<-factor(tle.data.geno$genotype,levels=c("1","3","7","13","31","36","41","44","50","62"))
tle.data.geno$timepoint<-factor(tle.data.geno$timepoint,levels=c("October","January","April", "June"))

# plot for figure 1
tle.plot.geno<-ggplot(data=tle.data.geno, aes(x=timepoint, y=mean, colour=genotype)) +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  theme_light()+ 
  geom_line(aes(group=genotype), alpha = 0.75) +
  geom_point(aes(color=genotype, shape = site), size =3) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=genotype), alpha = 0.5, width=0.05, size=0.5, linetype=1)+
  facet_wrap(~site)+
  ylab("Total Linear Extension (cm)")+
  xlab("Timepoint")+
  guides(color=guide_legend(ncol=2))+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15))
tle.plot.geno

#volume by genotype
vol.data.geno<- plyr::ddply(growth_months, c("genotype","months", "site"), summarise,
                            N    = length(volume[!is.na(volume)]),
                            mean = mean(volume, na.rm=TRUE),
                            sd   = sd(volume, na.rm=TRUE),
                            se   = sd / sqrt(N))

vol.plot.geno<-ggplot(data=vol.data.geno, aes(x=months, y=mean, colour=genotype, shape = site)) +
  scale_color_simpsons()+
  geom_line(aes(group=genotype)) +
  geom_point(aes(group=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=genotype), width=0.05, size=0.5, linetype=1)+
  theme_light() +
  facet_wrap(~site)+
  ylab("Volume (cm3)")+
  xlab("Months Outplanted")+
  theme(legend.position = "none")+
  scale_x_discrete(expand = c(0.01,0.01))
vol.plot.geno


#surface area by genotype
sa.data.geno<- plyr::ddply(growth_months, c("genotype","months", "site"), summarise,
                           N    = length(sa[!is.na(sa)]),
                           mean = mean(sa, na.rm=TRUE),
                           sd   = sd(sa, na.rm=TRUE),
                           se   = sd / sqrt(N))

sa.plot.geno<-ggplot(data=sa.data.geno, aes(x=months, y=mean, colour=genotype, shape = site)) +
  scale_color_simpsons()+
  geom_line(aes(group=genotype)) +
  geom_point(aes(group=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=genotype), width=0.05, size=0.5, linetype=1)+
  theme_light() +
  facet_wrap(~site)+
  ylab("SA (cm2)")+
  xlab("Months Outplanted")+
  theme(legend.position = "none")+
  scale_x_discrete(expand = c(0.01,0.01))
sa.plot.geno

#surface area:volume by genotype
sav.data.geno<- plyr::ddply(growth_months, c("genotype","months", "site"), summarise,
                            N    = length(sav[!is.na(sav)]),
                            mean = mean(sav, na.rm=TRUE),
                            sd   = sd(sav, na.rm=TRUE),
                            se   = sd / sqrt(N))

sav.plot.geno<-ggplot(data=sav.data.geno, aes(x=months, y=mean, colour=genotype, shape = site)) +
  scale_color_simpsons()+
  geom_line(aes(group=genotype)) +
  geom_point(aes(group=genotype)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=genotype), width=0.05, size=0.5, linetype=1)+
  theme_light() +
  facet_wrap(~site)+
  ylab("SA:V")+
  xlab("Months Outplanted")+
  scale_x_discrete(expand = c(0.01,0.01))
sav.plot.geno

#create composite figure for supplement
other_growth <-plot_grid(vol.plot.geno, sa.plot.geno, sav.plot.geno, labels = c("a", "b","c"),rel_widths = c(0.8,0.8,1.0), nrow=1)
other_growth

## linear models ####

#need to include initial size as a random factor for size metrics
initial_size <-initial.size %>% 
  select(tag_number, tle, sa, volume, sav)%>%
  rename(initial_tle = tle)%>%
  rename(initial_sa = sa)%>%
  rename(initial_vol = volume)%>%
  rename(initial_sav = sav)

size_data <- all_morphology%>%
  filter(timepoint != "T0")%>%
  left_join(., initial_size)

all_morphology <- all_morphology%>%
  left_join(., initial_size)
write.csv(all_morphology, file = "ORCC_morphology.csv")

# linear growth
#in this case, we want only positive growth to separate out breakage
model_data <- size_data%>%
  filter(growth >0)

growth.model = lmer(growth ~ timepoint*genotype*site + (1|array) + (1|initial_tle), data=model_data)

#check assumptions
qqPlot(residuals(growth.model)) #normality
leveneTest(residuals(growth.model)~model_data$timepoint*model_data$genotype*model_data$site)

summary(growth.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
rand(growth.model) # array does play a role - but we are accounting for this in the model
anova(growth.model, type =3) # use this view significant effects and get p-values

boxplot(growth ~ timepoint, dat = model_data)

#TLE
# density plot to see how we might need to transform the data
ggplot(aes(x = tle), data = size_data)+
  geom_density()

tle.model = lmer(1/tle ~ timepoint*genotype*site + (1|array) + (1|initial_tle), data=size_data)

rand(tle.model)
#check assumptions
qqPlot(residuals(tle.model)) #normality
leveneTest(residuals(tle.model)~size_data$timepoint*size_data$genotype*size_data$site)

summary(tle.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(tle.model, type =3) # use this view significant effects and get p-values

# Volume
# density plot to see how we might need to transform the data
ggplot(aes(x = volume), data = size_data)+
  geom_density()

vol.model = lmer(1/volume ~ timepoint*genotype*site + (1|array) + (1|initial_vol), data=size_data)

rand(vol.model)
#check assumptions
qqPlot(residuals(vol.model))
leveneTest(residuals(vol.model)~size_data$timepoint*size_data$genotype*size_data$site)

summary(vol.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(vol.model, type=3) #use this to view significant effects/interactions and get p-values
# sig effect of genotype and timepoint

# SA
# density plot to see how we might need to transform the data
ggplot(aes(x = sa), data = size_data)+
  geom_density()

sa.model = lmer(1/sa ~ timepoint*genotype*site + (1|array) + (1|initial_sa), data=size_data)

#check assumptions
qqPlot(residuals(sa.model))
leveneTest(residuals(sa.model)~size_data$timepoint*size_data$genotype*size_data$site)

summary(sa.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(sa.model, type=3) #use this to view significant effects/interactions and get p-values
# sig effect of genotype and timepoint

# SA:V
# density plot to see how we might need to transform the data
ggplot(aes(x = sav), data = size_data)+
  geom_density()

sav.model = lmer(sav ~ timepoint*genotype*site + (1|array) + (1|initial_sav), data=size_data)

#check assumptions
qqPlot(residuals(sav.model))
leveneTest(residuals(sav.model)~size_data$timepoint*size_data$genotype*size_data$site)

summary(sav.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(sav.model, type=3) #use this to view significant effects/interactions and get p-values
# sig effect of genotype and timepoint


## time seems to be sucking up all of the variance in these models - makes sense as they are getting bigger with each timestep! 
# let's try modeling size/growth rates for just the June timepoint
# can't include T0 size in model bc of same number of observations as grouping factor - so need to subtract initial size
june_size <- size_data%>%
  filter(timepoint =="3")%>%
  mutate(model_tle = tle - initial_tle)%>%
  mutate(model_vol = volume - initial_vol)%>%
  mutate(model_sa = sa - initial_sa)%>%
  mutate(model_sav = sav - initial_sav)

# TLE
tle.jun.model = lmer(1/model_tle ~ genotype*site + (1|array), data=june_size) # cant include T0 size in model bc same number of observations in y

#check assumptions
qqPlot(residuals(tle.jun.model))
leveneTest(residuals(tle.jun.model)~june_size$genotype*june_size$site)

summary(tle.jun.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(tle.jun.model, type=3) #use this to view significant effects/interactions and get p-values
# no sig effects

# volume
vol.jun.model = lmer(model_vol ~ genotype*site + (1|array), data=june_size) # cant include T0 size in model bc same number of observations in y

#check assumptions
qqPlot(residuals(vol.jun.model))
leveneTest(residuals(vol.jun.model)~june_size$genotype*june_size$site)

summary(vol.jun.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(vol.jun.model, type=3) #use this to view significant effects/interactions and get p-values
# sig effect of genotype, site, genotype x site

# SA
sa.jun.model = lmer(model_sa ~ genotype*site + (1|array), data=june_size) # cant include T0 size in model bc same number of observations in y

#check assumptions
qqPlot(residuals(sa.jun.model))
leveneTest(residuals(sa.jun.model)~june_size$genotype*june_size$site)

summary(sa.jun.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(sa.jun.model, type=3) #use this to view significant effects/interactions and get p-values
# sig effect of genotype, genotype x site

# SA:V
sav.jun.model = lmer(model_sav ~ genotype*site + (1|array), data=june_size) # cant include T0 size in model bc same number of observations in y

#check assumptions
qqPlot(residuals(sav.jun.model))
leveneTest(residuals(sav.jun.model)~june_size$genotype*june_size$site)

summary(sav.jun.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
anova(sav.jun.model, type=3) #use this to view significant effects/interactions and get p-values
# sig effect of genotype

# qPCR data ####
options(scipen=999) # To avoid scientific notation

#read in and filter out neg and pos controls from raw qPCR data
rawdata<-read_excel("qPCR/ORCC_bleaching_allAssays.xlsx")%>%
  clean_names()%>%
  filter(well_type!= "NTC")%>%
  filter(well_name!= "+C")

copy.n.A<-9                                              #from Cunning et al. 2017 Supplemental Materials
copy.n.Acer<-1 
copy.n.ratioAAcer<-(copy.n.A/copy.n.Acer)
fluo.A<-0
fluo.Acer<--6.305

data<-rawdata%>%  
  dplyr::select(well_name, target, cq_rn)%>%
  mutate_at("cq_rn", as.numeric)%>%
  na.replace(41)%>%
  group_by(well_name, target)%>%                         ###group technical replicates
  dplyr::summarise(meanCq = mean(cq_rn, na.rm = TRUE),             ###take mean of technical replicates
                   stdev = sd(cq_rn, na.rm = TRUE), 
                   .groups = "keep") 

A<-dplyr::filter(data,target=="FAM")
Acer<-dplyr::filter(data,target=="SYBR")

#final formating of data
final_data<-left_join(A,Acer,by="well_name")%>%
  dplyr::select(-target.x,-target.y)%>%                                              
  dplyr::rename(a_mean=meanCq.x)%>%                                           
  dplyr::rename(a_sd=stdev.x)%>%
  dplyr::rename(acer_mean=meanCq.y)%>%
  dplyr::rename(acer_sd=stdev.y)%>%
  mutate(acer_mean=acer_mean-fluo.Acer)%>%
  mutate(aacer_ratio=(2^(acer_mean-a_mean))/copy.n.ratioAAcer)%>%
  mutate(sd_warning=case_when((a_sd>1&acer_sd>1)~"Ahost*",a_sd>1~"a*",acer_sd>1~"host*"))%>%                                                                                   ###set warnings where standard deviation of tech replicates is >1
  mutate(ct_warning=case_when((a_mean>38&acer_mean>38)~"Ahost*",a_mean>38~"a*",acer_mean>38~"host*"))                                                                             ###set warnings where ct values is later than 38

write.csv(final_data, file = "qPCR/ORCC_sh_ratios.csv")

## compare the two extraction methods - filter just the ones that were extracted via both methods
extraction_data <- final_data%>%
  mutate(extraction = case_when(grepl("JD", well_name) ~ "JD",
                                TRUE ~ "NV"))%>%
  mutate(well_name = gsub("-JD", "", well_name))%>%
  filter(well_name == "15" | well_name == "7" | well_name == "88" |
           well_name == "32" | well_name == "14" |well_name == "97")
# they are veryyyyyyy different - so will have to filter out the DNA/RNA dual extracted samples for actual analysis

## adding in metadata
bleaching_data <- final_data%>%
  mutate(extraction = case_when(grepl("JD", well_name) ~ "JD",
                                TRUE ~ "NV"))%>%
  mutate(well_name = gsub("-JD", "", well_name))%>%
  filter(extraction != "JD")%>%
  dplyr::select(well_name, aacer_ratio)%>%
  mutate_at("well_name", as.numeric)%>%
  left_join(frag_info, by=c("well_name"="tag_number"))%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("replicate", as.factor)%>%
  filter(outplant_date == "Oct")%>%
  filter(well_name != "46")%>%
  mutate_at("well_name", as.factor)
str(bleaching_data)  

#reorder by survivorship 
bleaching_data$genotype<-factor(bleaching_data$genotype,levels=c("36","41","3","50","13","44","1","31","7","62"))

write.csv(bleaching_data, file ="qPCR/sh_ratios_meta.csv")

# exploratory figures - S:H by site/geno - for supplement
nolog_bleaching_plot<-ggplot(bleaching_data) +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  geom_boxplot(aes(x = genotype, y = aacer_ratio))+
  geom_jitter(aes(x = genotype, y = aacer_ratio ,color = genotype), width = 0.2)+
  facet_wrap(~site)+
  theme_bw()
nolog_bleaching_plot

#checking prevalence of C/D in a subset of samples
#read in and filter out neg and pos controls from raw qPCR data
rawdata<-read_excel("qPCR/ORCC_symCD_allAssays.xlsx")%>%
  clean_names()%>%
  filter(well_type!= "NTC")%>%
  filter(well_name!= "+C")

copy.n.A<-9                                              #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioAAcer<-(copy.n.A/copy.n.Acer)
fluo.A<-0
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
copy.n.Acer<-1                                           
copy.n.ratioCD<-(copy.n.C/copy.n.D)
copy.n.ratioCAcer<-(copy.n.C/copy.n.Acer)
copy.n.ratioDAcer<-(copy.n.D/copy.n.Acer)
copy.n.ratioAAcer<-(copy.n.A/copy.n.Acer)
fluo.C<-0
fluo.D<-1.238
fluo.Acer<--5.3

data<-rawdata%>%  
  dplyr::select(well_name, target, cq_rn)%>%
  mutate_at("cq_rn", as.numeric)%>%
  group_by(well_name, target)%>%                         ###group technical replicates
  dplyr::summarise(meanCq = mean(cq_rn, na.rm = TRUE),             ###take mean of technical replicates
                   stdev = sd(cq_rn, na.rm = TRUE), 
                   .groups = "keep") 

C<-dplyr::filter(data,target=="VIC")
D<-dplyr::filter(data,target=="FAM")
Acer<-dplyr::filter(data,target=="SYBR")

#final formating of data
final_data<-left_join(C,D,by="well_name")%>%
  dplyr::select(-target.x,-target.y)%>%                                              
  dplyr::rename(c_mean=meanCq.x)%>%                                           
  dplyr::rename(c_sd=stdev.x)%>%
  dplyr::rename(d_mean=meanCq.y)%>%
  dplyr::rename(d_sd=stdev.y)%>%
  left_join(.,Acer,by="well_name")%>%
  dplyr::select(-target)%>%                                              
  dplyr::rename(Acer_mean=meanCq)%>%                                           
  dplyr::rename(Acer_sd=stdev)%>%
  mutate(Acer_mean=Acer_mean-fluo.Acer)%>%
  left_join(.,A,by="well_name")%>%
  dplyr::select(-target)%>%                                              
  dplyr::rename(a_mean=meanCq)%>%                                           
  dplyr::rename(a_sd=stdev)%>%
  mutate(d_mean=d_mean-fluo.D)%>%
  mutate(dacer_ratio=(2^(Acer_mean-d_mean))/copy.n.ratioDAcer)%>%                                #from Cunning et al. 2012
  mutate(cacer_ratio=(2^(Acer_mean-c_mean))/copy.n.ratioCAcer)%>%     #from Cunning et al. 2012
  mutate(aacer_ratio=(2^(Acer_mean-a_mean))/copy.n.ratioAAcer)%>%
  replace_na(list(dacer_ratio = 0, cacer_ratio = 0, aacer_ratio = 0))%>%
  mutate(prop_c = cacer_ratio/(cacer_ratio+dacer_ratio+aacer_ratio))%>%
  mutate(prop_d = dacer_ratio/(cacer_ratio+dacer_ratio+aacer_ratio))%>%
  mutate(prop_a = aacer_ratio/(cacer_ratio+dacer_ratio+aacer_ratio))%>%
  mutate(ct_warning=case_when((c_mean>38&d_mean>38)~"cd*",d_mean>38~"d*", c_mean>38~"c*")) 

write.csv(final_data, file = "qPCR/ORCC_CD_ratios.csv")

CD_data <- final_data%>%
  mutate(extraction = case_when(grepl("JD", well_name) ~ "JD",
                                TRUE ~ "NV"))%>%
  mutate(well_name = gsub("-JD", "", well_name))%>%
  dplyr::select(well_name, dacer_ratio, cacer_ratio, aacer_ratio, prop_d)%>%
  mutate_at("well_name", as.numeric)%>%
  left_join(frag_info, by=c("well_name"="tag_number"))%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("replicate", as.factor)%>%
  filter(outplant_date == "Oct")%>%
  mutate_at("well_name", as.factor)

CD_plot<-ggplot(CD_data) +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  geom_boxplot(aes(x = genotype, y = sh_ratio))+
  geom_jitter(aes(x = genotype, y = sh_ratio ,color = genotype), width = 0.2)+
  facet_wrap(~site)+
  theme_bw()
CD_plot

CA_plot<- ggplot(CD_data, aes(x = log10(cacer_ratio), y = log10(aacer_ratio), colour = genotype, shape = site))+
  scale_color_manual(values =c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  geom_point(aes(group = site),size = 2)+
  theme_light()+
  geom_abline(linetype = "dashed")+
  ylim(-6.8, 1)+
  xlim(-6.8, 1)+
  coord_fixed()+
  ylab("log(A:H)")+
  xlab("log(C:H)")+
  guides(color=guide_legend(ncol=2))+
  theme(legend.position = "none"
  )
CA_plot

DA_plot<- ggplot(CD_data, aes(x = log10(dacer_ratio), y = log10(aacer_ratio), colour = genotype, shape = site))+
  scale_color_manual(values =c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  geom_point(aes(group = site),size = 2)+
  theme_light()+
  geom_abline(linetype = "dashed")+
  ylim(-6.5, 1)+
  xlim(-6.5, 1)+
  coord_fixed()+
  ylab("log(A:H)")+
  xlab("log(D:H)")+
  guides(color=guide_legend(ncol=2))
DA_plot

# create composite plot for supplement
symplot <- plot_grid(CA_plot, DA_plot, labels = c("a","b"), rel_widths = c(0.73,1))
symplot

# check if proportion D correlates with color score:
D_colscore <- CD_data%>%
  select(well_name, prop_d)%>%
  left_join(., bleached_scores, join_by(well_name == id))%>%
  na.omit()

cor(D_colscore$prop_d, D_colscore$bleach_score)
#-0.05643657

#plot for supplement
ggplot(D_colscore, aes(x = prop_d, y = bleach_score))+
  geom_point(aes(color = genotype, shape = site), size = 3)+
  scale_color_manual(values =c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  theme_light()+
  ylab("Color Score in July")+
  xlab("Proportion Durusdinium")+
  facet_wrap(~site)+
  guides(color=guide_legend(ncol=2))+
  theme(axis.text=element_text(size=14),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15))

# color scores ####
# read in and wrangle data
color_info <- read_xlsx("Mortality/ORCC_bleaching_mortality_Jul23_QCed.xlsx", sheet = "Mortality_binary")%>%
  clean_names()%>%
  dplyr::select(site, genotype, replicate, tag_number, outplant_date, jun_2023_bleach, jul_23)%>%
  mutate_at("jun_2023_bleach", as.numeric)%>%
  mutate_at("jul_23", as.numeric)

bleached_scores <- read.csv("ColorScore/bleached_scores.csv", header = TRUE)%>%
  clean_names()%>%
  dplyr::select(id, scs)%>%
  dplyr::rename(bleach_score = scs)%>%
  full_join(., color_info, join_by(id == tag_number))%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("replicate", as.factor)%>%
  mutate_at("id", as.factor)%>%
  unique()

baseline_scores <- read.csv("ColorScore/baseline_scores.csv", header = TRUE)%>%
  clean_names()%>%
  dplyr::select(id, scs)%>%
  dplyr::rename(base_score = scs)%>%
  full_join(., color_info, join_by(id == tag_number))%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("site", as.factor)%>%
  mutate_at("replicate", as.factor)%>%
  mutate_at("id", as.factor)%>%
  unique()

color_scores <- bleached_scores%>%
  mutate_at("replicate", as.factor)%>%
  full_join(.,baseline_scores)%>%
  mutate_at("genotype", as.factor)%>%
  mutate_at("site", as.factor)%>%
  mutate(score_change = base_score-bleach_score)%>%
  filter(outplant_date =="Oct")
str(color_scores)

# baseline color score boxplot by site and genotype
ggplot(color_scores) +
  geom_boxplot(aes(x = genotype, y = base_score, color=genotype),outlier.shape =NA)+
  geom_jitter(aes(x = genotype, y = base_score, color=genotype), alpha=0.5, width = 0.2)+
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  facet_wrap(~site)+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()

# test for baseline differences in color between genotypes at each site
daves_color <- color_scores%>%
  filter(site == "Daves")

daves_color_anova <- lm(base_score ~ genotype, data = daves_color)

anova(daves_color_anova)
#Df Sum Sq Mean Sq F value  Pr(>F)  
#genotype   9 23.178  2.5754  2.3294 0.02273 *
# Residuals 74 81.813  1.1056                  

TukeyHSD(aov(daves_color_anova))
# sig diff between 7 and 31

looe_color <- color_scores%>%
  filter(site == "Looe")

looe_color_anova <- lm(base_score ~ genotype, data = looe_color)

anova(looe_color_anova)
#Response: base_score
#Df Sum Sq Mean Sq F value Pr(>F)  
#genotype    9 16.554 1.83931  2.3928 0.02267 *
#Residuals 56 43.046 0.76867                  

TukeyHSD(aov(looe_color_anova))
# not seeing significant differences between specific genotypes

# checking on correlation between color score and S:H 
bleaching_measures <- color_scores%>%
  mutate_at("id", as.factor)%>%
  left_join(bleaching_data, by = c("id"="well_name", "genotype"="genotype", "site"="site"))%>%
  na.omit()

shvcol<-ggplot(bleaching_measures, aes(x=bleach_score, y=aacer_ratio)) +
  scale_color_simpsons()+
  geom_smooth(method = "lm")+
  geom_point(size=3, aes(color = genotype))+
  theme_classic()
shvcol

cor(bleaching_measures$bleach_score, bleaching_measures$aacer_ratio)
# 0.3692167

## bleaching stress index ####
# formula from https://www.nature.com/articles/s41467-024-52895-1
## formula: BSI = 1 - ((0xc1 + 1xc2 + 2xc3 + 3xc4 + 4xc5)/N-1)

#binning color socres and mortality into BSI categories
bleach_index <- color_scores%>%
  select(id, site, genotype, bleach_score, jul_23, jun_2023_bleach)%>%
  dplyr::mutate(score = case_when(bleach_score > 4.99 ~ "6",
                                  bleach_score > 3.99 ~ "5",
                                  bleach_score > 2.99 ~ "4",
                                  bleach_score > 1.99 ~ "3",
                                  bleach_score > 0.99 ~ "2",
                                  bleach_score > 0 ~ "1",
                                  is.na(bleach_score)&jul_23 ==1 ~ "0",
                                  is.na(bleach_score)&is.na(jul_23) ~"NA"))%>%
  mutate_at("score", as.numeric)

combined_scores <- bleach_index%>%
  group_by(genotype, site, score, .drop = FALSE)%>%
  dplyr::count(score)%>%
  mutate_at("score", as.factor)%>%
  na.omit()

combined_scores$score <- factor(combined_scores$score ,levels = c("0","1","2","3","4","5","6"))

# plot for figure 4b
col=c("grey30","white","grey80","burlywood3","tan3","tan4","darkorange4")
stacked.plot <-ggplot(data=combined_scores, aes(fill=score, x=genotype, y=n))+
  geom_bar(position="fill", stat="identity", colour="black")+
  facet_wrap(~site, ncol=2)+
  scale_fill_manual(values=col)+
  scale_color_manual(values=col)+
  theme_classic()+
  ylab("Proportion")+
  xlab("Genotype")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))
stacked.plot

# calculating BSI per genotype at each site
## formula: BSI = 1 - ((0xc1 + 1xc2 + 2xc3 + 3xc4 + 4xc5)/N-1)

#daves ledge
daves_genotype_totals <- bleach_index%>%
  filter(site == "Daves")%>%
  select(genotype, id, jun_2023_bleach)%>%
  na.replace(0)%>%
  group_by(genotype, .drop = FALSE)%>%
  summarise_each(funs(sum(.==0)))%>%
  select(genotype, jun_2023_bleach)%>%
  mutate_at("genotype", as.factor)%>%
  rename("n"="jun_2023_bleach")

daves_bsi<- bleach_index%>%
  filter(site == "Daves")%>%
  group_by(genotype, .drop = FALSE)%>%
  dplyr::count(score)%>%
  pivot_wider(names_from = score, values_from = n)%>%
  mutate_at("genotype", as.factor)%>%
  left_join(daves_genotype_totals)%>%
  na.replace(0)%>%
  clean_names()%>%
  mutate(bsi = 1-((6*(x0/n) + 5*(x1/n) + 4*(x2/n) + 3*(x3/n) +2*(x4/n))/6))%>%
  mutate(site = "Daves")%>%
  mutate(x5 = 0)%>%
  mutate(x6 = 0)

#looe
looe_genotype_totals <- bleach_index%>%
  filter(site == "Looe")%>%
  select(genotype, id, jun_2023_bleach)%>%
  na.replace(0)%>%
  group_by(genotype, .drop = FALSE)%>%
  summarise_each(funs(sum(.==0)))%>%
  select(genotype, jun_2023_bleach)%>%
  mutate_at("genotype", as.factor)%>%
  rename("n"="jun_2023_bleach")

looe_bsi<- bleach_index%>%
  filter(site == "Looe")%>%
  group_by(genotype, .drop = FALSE)%>%
  dplyr::count(score)%>%
  pivot_wider(names_from = score, values_from = n)%>%
  mutate_at("genotype", as.factor)%>%
  left_join(looe_genotype_totals)%>%
  na.replace(0)%>%
  clean_names()%>%
  mutate(bsi = 1-((6*(x0/n) + 4*(x2/n) + 3*(x3/n) + 2*(x4/n)  + 1*(x5/n) + 0*(x6/n))/6))%>%
  mutate(site = "Looe")%>%
  mutate(x1 = 0)

# Multivariate analyses ####
#investigating overall morphological plasticty 

# read in and format morphology data - create matrix and explanatory variables for RDA
# new growth variable - positive growth only
# new breakage variable - if there was negative growth, transform to be positive
morpho.data <- read.csv("ORCC_morphology.csv", header=TRUE)%>%
  dplyr::select(-X)%>%
  filter(timepoint != "T0")%>%
  mutate(timepoint = case_when(
    timepoint == "3" ~"T3",
    TRUE ~ as.character(timepoint)
  ))%>%
  mutate(growth_pos = case_when(
    growth >0 ~ as.numeric(growth),
    TRUE ~ 0
  ))%>%
  mutate(breakage_net= case_when(
    growth <0 ~ abs(growth),
    TRUE ~ 0
  ))%>%
  mutate(breakage_type = case_when(
    breakage == 0 & breakage_net < 0 ~ 1,
    TRUE ~ as.numeric(breakage)
  ))%>%
  mutate(model_tle = tle - initial_tle)%>%
  mutate(model_vol = volume - initial_vol)%>%
  mutate(model_sa = sa - initial_sa)%>%
  mutate(model_sav = sav - initial_sav)

# create matrix with relevant traits and explanatory variables
morpho <- morpho.data[c("model_tle", "model_sa", "model_vol", "model_sav", "growth_pos","breakage_net", "breakage_type")] 
time <- (morpho.data$timepoint)
site <- morpho.data$site
genotype <- as.character(morpho.data$genotype)

#RDA
morpho_rda <- rda(morpho ~ genotype + site + genotype:site + Condition(time), na.action=na.exclude, scale = TRUE)
head(summary(morpho_rda))
#timepoint explains 11.15% of the variance
# other variables explain 8.82% of the variance

RsquareAdj(morpho_rda)
# 0.08826938

anova.cca(morpho_rda, by = "terms")
# genotype: p = 0.001
# genotype x site: p = 0.001

anova.cca(morpho_rda, by = "axis")
# RDA1 = 0.001
# RDA2 = 0.503

# visualize
coords <- as.matrix(scores(morpho_rda)$sites)
coef <- data.frame(scores(morpho_rda, display="bp", scaling = 1))

morpho2 <- cbind(morpho.data, coords)
centroids <- morpho2 %>% group_by(genotype, site, timepoint)   %>% 
  summarise(mean_PC1 = mean(RDA1),
            mean_PC2 = mean(RDA2))

morpho2 <- morpho2%>%
  left_join(centroids)

# get scores for arrows
arrows <- as.data.frame(vegan::scores(morpho_rda, display = "species"))
arrows$names <- rownames(arrows)

morpho2$genotype<-factor(morpho2$genotype,levels=c("1","3","7","13","31","36","41","44","50","62"))

morpho2$timepoint<-factor(morpho2$timepoint,levels=c("T3", "T2", "T1"))

# plot all points for figure 2a
morpho_plot <- ggplot(morpho2)+
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  theme_light()+ 
  geom_point(aes(x=RDA1, y = RDA2, color=as.factor(genotype), shape = site), size = 2.5) +
  ylab("RDA2 (20.66%)")+
  xlab("RDA1 (59.08%)")+
  stat_ellipse(aes(x=RDA1, y = RDA2, alpha = timepoint), color = "darkblue", linewidth = 1)+
  geom_segment(data=arrows, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), color = "grey30", linewidth =0.8,arrow = arrow(length = unit(0.25, "cm")))+
  guides(color=guide_legend(ncol=2))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
morpho_plot

# get separate plot with vectors for inset
vector_panel <- ggplot(morpho2)+
  geom_segment(data=arrows, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel(data=arrows, aes(x = RDA1, y = RDA2, label = names))+
  ylab("RDA2 (20.66%)")+
  xlab("RDA1 (59.08%)")+
  scale_x_continuous(lim = c(-4,2))+
  scale_y_continuous(lim = c(-6,3))+
  theme_classic()
vector_panel

# export as multipanel to edit in Affinity Designer
multipanel <- plot_grid(morpho_plot, vector_panel, rel_widths = c(2.5,1))
multipanel

# visualize site/genotype averages over time
coords <- as.matrix(scores(morpho_rda)$sites)
coef <- data.frame(scores(morpho_rda, display="bp", scaling = 1))

morpho3 <- cbind(morpho.data, coords)
centroids <- morpho3 %>% group_by(genotype, site, timepoint)   %>% 
  summarise(mean_RDA1 = mean(RDA1),
            mean_RDA2 = mean(RDA2))

morpho3 <- morpho3%>%
  left_join(centroids)

morpho3$timepoint<-factor(morpho3$timepoint,levels=c("T3", "T2", "T1"))

# create example plasticity calculation plots for genotype 62
example <- morpho3%>%
  mutate(timepoint = case_when(
    timepoint == "T1"~1,
    timepoint =="T2"~2,
    timepoint == "T3"~3
  ))%>%
  arrange(timepoint)%>%
  filter(genotype == 62) 

example_labels <- example%>%
  dplyr::select(genotype, timepoint, site, mean_RDA1, mean_RDA2)%>%
  unique()

# figure 2b - genotype 62 site centroids connected through time
example62 <- ggplot(example)+
  scale_color_manual(values = c("#71D0F5FF"))+
  theme_light()+ 
  geom_path(aes(x = mean_RDA1, y = mean_RDA2, group = interaction(genotype,site)), lwd = 0.75, col = "grey30") +
  geom_point(aes(x = mean_RDA1, y = mean_RDA2, color=as.factor(genotype), shape = site), size = 5.5) +
  geom_point(aes(x = RDA1, y = RDA2, color=as.factor(genotype), shape = site), size = 2, alpha=0.5,show.legend = FALSE) +
  geom_text_repel(data = individual_labels, aes(x = mean_RDA1, y = mean_RDA2, label = timepoint), size=5)+
  ylab("RDA2 (20.66%)")+
  xlab("RDA1 (59.08%)")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")
example62

# for supplement - each timepoint faceted out, to show distance calculated at each interval
supplement62 <- ggplot(example)+
  scale_color_manual(values = c("#71D0F5FF"))+
  theme_light()+ 
  geom_path(aes(x = mean_RDA1, y = mean_RDA2, group = interaction(genotype,site)), lwd = 0.75, col = "grey30") +
  geom_point(aes(x = mean_RDA1, y = mean_RDA2, color=as.factor(genotype), shape = site), size = 5.5) +
  geom_point(aes(x = RDA1, y = RDA2, color=as.factor(genotype), shape = site), size = 2, alpha=0.5,show.legend = FALSE) +
  facet_wrap(~timepoint)+
  ylab("RDA2 (20.66%)")+
  xlab("RDA1 (59.08%)")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
supplement62

# for Figure S2, change out numbers with correct genotype number and replace correct color from scale:
toPlot <- morpho3%>%
  mutate(timepoint = case_when(
    timepoint == "T1"~1,
    timepoint =="T2"~2,
    timepoint == "T3"~3
  ))%>%
  arrange(timepoint)%>%
  filter(genotype == 1)

individual_labels <- toPlot%>%
  dplyr::select(genotype, timepoint, site, mean_RDA1, mean_RDA2)%>%
  unique()

# values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF")
plot1 <- ggplot(toPlot)+
  scale_color_manual(values = c("#FED439FF"))+
  theme_light()+ 
  geom_path(aes(x = mean_RDA1, y = mean_RDA2, group = interaction(genotype,site)), lwd = 0.75, col = "grey30") +
  geom_point(aes(x = mean_RDA1, y = mean_RDA2, color=as.factor(genotype), shape = site), size = 5.5) +
  geom_point(aes(x = RDA1, y = RDA2, color=as.factor(genotype), shape = site), size = 2, alpha=0.5,show.legend = FALSE) +
 geom_text_repel(data = individual_labels, aes(x = mean_RDA1, y = mean_RDA2, label = timepoint), size=5)+
  ylab("RDA2 (20.66%)")+
  xlab("RDA1 (59.08%)")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")
plot1

#combine all into one supplemental figure
genoplots <- plot_grid(plot1, plot3, plot7, plot13, plot31, plot36, plot41, plot44, plot50, nrow = 3)
genoplots

##calculating plasticity ####
#euclidian distances between sites over time = morphological plasticity
coords <- scores(morpho_rda, display="sites", choices=c(1:3))
data.coords <- cbind(morpho.data[,c(1,2,3)], coords)

mean.coords <- data.coords %>% group_by(genotype, site, timepoint)   %>% 
  summarise(RD1 = mean(RDA1),
            RD2 = mean(RDA2)) 

plasticity <- mat.or.vec(0,0)

genets <- levels(as.factor(paste(mean.coords$genotype)))
times <- levels(as.factor(mean.coords$timepoint))
effet.env <- mat.or.vec(0,0)

for (j in 1:length(times)){
  time.sep <- subset(mean.coords, timepoint==times[j])
  
  for (i in 1:length(genets)){
    geno <- genets[i]
    tp <- subset(time.sep, genotype == geno)
    morpho <- tp[c("RD1","RD2")]
    dist.mean <- vegdist(morpho, method="euclidean")
    tp2 <- data.frame(tp, 
                      dist = dist.mean[1])
    effet.env <- rbind(effet.env, tp2)
  }
}

plasticity_morpho <- effet.env%>%
  filter(site =="Daves")%>%
  mutate_at("genotype", as.factor)

plasticity_morpho$genotype<-factor(plasticity_morpho$genotype,levels=c("1","3","7","13","31","36","41","44","50","62"))

#plot changes in plasticity over time - Figure 2c
ggplot(data=plasticity_morpho, aes(x=timepoint, y=dist, col=genotype)) + 
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  geom_point(size =3.5, shape = 15)+
  geom_line(aes(group = genotype), linetype = 3)+
  ylab("Plasticity (Euclidean distance)")+
  xlab("Timepoint")+
  theme_light()+
  scale_x_discrete(expand = c(0.01,0.01))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.position = "none")

# pull out overall genotype plasticity for each timepoint
plasticity_morpho <- plasticity_morpho%>%
  dplyr::select(genotype, timepoint, dist)%>%
  mutate_at("timepoint", as.factor)

bleaching_plasticity <- plasticity_morpho%>%
  filter(timepoint == "T3")

# visualize differences between high and low plasticity genotypes
#create new factor that assigns genotypes to high or low plasticity (top 5 = high, bottom 5 = low)
morpho4 <- morpho3%>%
  mutate(plast_level = case_when(
    genotype == 62 ~ "high",
    genotype == 31 ~ "high",
    genotype == 13 ~ "high",
    genotype == 50 ~ "high",
    genotype == 44 ~ "high",
    genotype == 7 ~ "low",
    genotype == 3 ~ "low",
    genotype == 41 ~ "low",
    genotype == 1 ~ "low",
    genotype == 36 ~ "low",
  ))

# high vs low plasticity plot - figure 2d
contrast_plot <-ggplot(morpho4)+
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  theme_light()+ 
 geom_point(aes(x = RDA1, y = RDA2, color=as.factor(genotype), shape = site), alpha = 0.5,size = 2) +
  stat_ellipse(aes(x=RDA1, y = RDA2, alpha = timepoint), color = "darkblue", linewidth = 1)+
  guides(color = guide_legend(ncol =2))+
  ylab("RDA2 (20.66%)")+
  xlab("RDA1 (59.08%)")+
  facet_wrap(~plast_level)+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
contrast_plot

##calculating noise ####
#distance between centroids and each indvidual datapoint at each site = noise
noise <- morpho3 %>%
  dplyr::select(site, timepoint, genotype, tag_number, RDA1, RDA2, mean_RDA1, mean_RDA2)%>%
  mutate(distance = sqrt((RDA1 - mean_RDA1)^2 + (RDA2 - mean_RDA2)^2))

# for supplement: boxplot of noise for each genotype for each timepoint, ordered by plasticity at that timepoint
geno_colors <- c("1" = "#FED439FF", "3" = "#709AE1FF", "7" = "#8A9197FF", "13" = "#D2AF81FF", "31" = "#FD7446FF", 
                 "36" = "#D5E4A2FF", "41" = "#197EC0FF", "44" = "#C80813FF", "50" = "#46732EFF", "62" = "#71D0F5FF")

T1_noise <- noise%>%
  filter(timepoint == "T1")

T1_noise$genotype<-factor(T1_noise$genotype,levels=c("62", "44", "31", "7", "13", "41", "1", "3", "36", "50"))

T1_noiseplot <- ggplot(T1_noise, aes(x = genotype, y = distance))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color = genotype, shape = site))+
  scale_color_manual(values = geno_colors)+
  theme_light()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.position = "none")+
  ylab("Noise (Distance to Centroid)")+
  xlab("")
T1_noiseplot

T2_noise <- noise%>%
  filter(timepoint == "T2")

T2_noise$genotype<-factor(T2_noise$genotype,levels=c("44", "31", "62", "13", "50", "3", "1", "41", "7", "36"))

T2_noiseplot <- ggplot(T2_noise, aes(x = genotype, y = distance))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color = genotype, shape = site))+
  scale_color_manual(values = geno_colors)+
  theme_light()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.position = "none")+
  ylab("")+
  xlab("Genotype")
T2_noiseplot

T3_noise <- noise%>%
  filter(timepoint == "T3")

T3_noise$genotype<-factor(T3_noise$genotype,levels=c("62", "31", "13", "50", "44", "7", "3", "41", "1", "36"))

T3_noiseplot <- ggplot(T3_noise, aes(x = genotype, y = distance))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color = genotype, shape = site))+
  scale_color_manual(values = geno_colors)+
  theme_light()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=15))+
  ylab("")+
  xlab("")
T3_noiseplot

noise_composite <- plot_grid(T1_noiseplot, T2_noiseplot, T3_noiseplot, nrow = 1, rel_widths = c(1,1,1.3))
noise_composite

## ambient plasticity effects ####
# Time-dependent Cox model for ambient temps 
# take into account genotype, site, and changing plasticity values over time

# first, get data in correct format: code time in discrete blocks, with mortality outcomes and plasticity values for each timestep
# ie T0 - T1: mortality and plasticity for each fragment, then T1-T2 etc.
full_mortality <- read.csv("Mortality/survivorship_time.csv")%>%
  dplyr::select(genotype, site, tag_number, timepoint, mortality)%>%
  filter(timepoint != "october")%>%
  filter(timepoint != "july")

#rename timepoints to match plasticity dataset
full_mortality$timepoint <- sub("january", "T1", full_mortality$timepoint)
full_mortality$timepoint <- sub("april", "T2", full_mortality$timepoint)
full_mortality$timepoint <- sub("june", "T3", full_mortality$timepoint)

interval_mortality <-full_mortality%>%
  group_by(genotype, site, tag_number) %>%
  arrange(tag_number, timepoint) %>%  # Ensure the data is ordered by tag_number and timepoint
  mutate(
    start_interval = lag(timepoint, order_by = timepoint),  # Start of interval (previous timepoint)
    end_interval = timepoint  # End of interval (current timepoint)
  ) %>%
  filter(!is.na(start_interval))%>%  # Remove the first observation for each individual (no previous interval)
  group_by(genotype, site, tag_number, start_interval, end_interval) %>%
  summarise(mortality = max(mortality), .groups = 'drop')%>%  # max ensures we capture mortality during the interval
  mutate_at("genotype", as.factor)%>%
  left_join(plasticity_morpho, join_by(genotype == genotype, start_interval == timepoint))

#rename timepoints to be number of months, turn into numeric 
interval_mortality$start_interval <- sub("T1", "3", interval_mortality$start_interval)
interval_mortality$start_interval <- sub("T2", "6", interval_mortality$start_interval)

interval_mortality$end_interval <- sub("T2", "6", interval_mortality$end_interval)
interval_mortality$end_interval <- sub("T3", "9", interval_mortality$end_interval)

interval_mortality$start_interval <- as.numeric(interval_mortality$start_interval)
interval_mortality$end_interval <- as.numeric(interval_mortality$end_interval)

fit <- coxph(Surv(start_interval, end_interval, mortality) ~ dist, data = interval_mortality)
summary(fit)

#plot along with our genotype model - fogure 3
forest_model(model_list = list("1. Genotype" =ambient.genotype.cox, "2. Dynamic Plasticity" = fit),
             factor_separate_line = FALSE)

## investigate tradeoffs ####
#Bleaching stress index, join with plasticity data
bsi_data <- read.csv("ColorScore/ORCC_bleachindex.csv")%>%
  dplyr::select(genotype, site, bsi)%>%
  mutate_at("genotype", as.factor)%>%
  left_join(., bleaching_plasticity)

bsi_data$genotype<-factor(bsi_data$genotype,levels=c("1","3","7","13","31","36","41","44","50","62"))
#model 
bsi_model <- glm(bsi~dist, data = bsi_data)
summary(bsi_model) 
# not significant

#what about for each site individually?
looe <- bsi_data%>%
  filter(site == "Looe")

looe_model <- glm(bsi~dist, data = looe)
summary(looe_model)
# p = 0.0104

RsquareAdj(looe_model)
# R2 =0.5807892

daves <- bsi_data%>%
  filter(site == "Daves")

daves_model <- glm(bsi~dist, data = daves)
summary(daves_model)
# not signficant, p = 0.38340

RsquareAdj(daves_model)
#R2 = 0.09608979

# plot the relationship - Figure 4b
plast_bsi <-ggplot(bsi_data, aes(x=dist, y=bsi)) +
  scale_color_manual(values = c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF","#197EC0FF", "#C80813FF", "#46732EFF","#71D0F5FF"))+
  geom_smooth(method = "lm", alpha = 0.2, color = "grey")+
  geom_jitter(aes(color=genotype, shape =site), size=3) +
  ylab("Bleaching Stress Index")+
  xlab("Plasticity in June")+
  theme_light()+
  facet_wrap(~site)+
  scale_x_continuous(expand=c(0.01,0.01))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "right")
plast_bsi

