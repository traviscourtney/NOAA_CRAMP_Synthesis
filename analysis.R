#####################################################
# Quantifying census-based net carbonate budgets and chemistry-based net calcification for Pacific Ocean coral reefs following the 2014-2017 global coral bleaching event 
# Travis Courtney (traviscourtney@gmail.com)
#####################################################

###########Load libraries##########
library(tidyverse)
library(ggridges)
library(patchwork)
library(FSA)
library(ggpubr)
library(readxl)
library(nlme)
library(sjPlot)
library(lubridate)
library(lutz)
library(sf)
library(data.table)
library(maps)
library(ggsignif)

###########Import data##########
#load CoralNet calcification rate data
rates=read_csv("Data/CoralNet_Calcification_With_Bioerosion_Rates_1.csv")
#load NOAA CRED labels
CRED_labels=read_csv("Data/NOAA_CRED_CoralNet_labelset2.csv")
#load benthic cover data
cramp3=read_csv("Data/BenthicCover_2010-2020_Tier3_SITE_OCC.csv")
cramp1=read_csv("Data/BenthicCover_2010-2020_Tier1_SITE_OCC.csv")
#import carbonate chemistry data
data=read.csv("Data/V_NCEI_OCC_OCADS_H2O_Pacific_2006-2019.csv")
#load max daily DHW data extracted for each site
DHW_OCC=read_csv("Data/DHW_OCC_2.csv")
#load ReefBudget parrotfish bioerosion rates
RB_parrotfish=read_csv("Data/ReefBudgetParrotfishBioerosionRates.csv")
#load stratified random survey data for each region
WAKE_fishes=read_csv("Data/V0_FISH_REA_PRIAs_2017.csv")
PRIA_fishes=read_csv("Data/V0_FISH_REA_PRIAs_2018.csv")
NWHI_fishes=read_csv("Data/CRCP_Reef_Fish_Surveys_Hawaii_a697_0026_3955.csv")
SAMOA_fishes=read_csv("Data/V0_FISH_REA_SAMOA_2018.csv")
MARIANAS_fishes=read_csv("Data/V0_FISH_REA_MARIAN_2017.csv")
#load stratified random parrotfish survey data
wd_OCC <- read.csv("Data/parrot_raw.csv")
#load stratified random survey area scaling
strata.OCC <- read.csv("Data/strata_area.csv")

###########Calculate gross carbonate production with macro and micro bioerosion##########
#here we will apply calcification and endolithic bioerosion rates to the benthic communities at each site to calculate G and uncertainty for each site
#remove asterisks from CRED labels
CRED_labels$`Short Code`=str_sub(CRED_labels$`Short Code`, 2, str_length(CRED_labels$`Short Code`))
#merge CRED labels with CoralNet calcification rate data and rename columns
calc_rates=merge(CRED_labels,rates,by="Name",all.x=TRUE) %>% 
  rename("name"=`Name`,"label"=`Short Code`,"calc"=`Mean`,"calc_lower"=`Lower bound`,"calc_upper"=`Upper bound`) %>% 
  select(-c(`Functional Group`,"Region"))

#calculate calcifier cover from tier 1 and merge with tier 3 data
cramp1$calcifier_cover=cramp1$CCA+cramp1$CORAL #set calcifier cover as corals + CCA
cramp1sub=bind_cols(SITEVISITID=cramp1$SITEVISITID,CORAL=cramp1$CORAL,CCA=cramp1$CCA,CALCIFIERS=cramp1$calcifier_cover,ENCRUSTING=cramp1$EMA,HALIMEDA=cramp1$HAL,SOFTCORAL=cramp1$SC,INVERT=cramp1$I,MACROALGAE=cramp1$MA,TURF=cramp1$TURF,SEDIMENT=cramp1$SED)
cramp_data=merge(cramp3,cramp1sub,by="SITEVISITID")

#select only climate sites
crampOCC = subset(cramp_data,OCC_SITEID!="NA")

#convert benthic cover data to long format
cramp_long=gather(crampOCC,label,percent,'ACAS':'TAB',factor_key=TRUE)

#merge cover and rates
CRAMP_rates=merge(cramp_long,calc_rates,by="label",all.x=TRUE)
CRAMP_rates$g=(CRAMP_rates$percent/100)*CRAMP_rates$calc
CRAMP_rates$g_lwr=(CRAMP_rates$percent/100)*CRAMP_rates$calc_lower
CRAMP_rates$g_upr=(CRAMP_rates$percent/100)*CRAMP_rates$calc_upper

#calculate carbonate production as the sum of all calcification for each station id in each site name for each year surveyed
CRAMP_siteG = CRAMP_rates %>%
  group_by(REGION,ISLAND,MISSIONID,OCC_SITEID,SITE,OBS_YEAR,CALCIFIERS,CORAL,CCA,HALIMEDA,SOFTCORAL,ENCRUSTING,MACROALGAE,TURF,SEDIMENT,INVERT,LATITUDE,LONGITUDE,REEF_ZONE,DEPTH_BIN) %>%
  summarize(siteG=sum(g,na.rm=TRUE),siteG_lwr=sum(g_lwr,na.rm=TRUE),siteG_upr=sum(g_upr,na.rm=TRUE))

###########Calculate Bioerosion##########
#here we will calculate parrotfish and urchin bioerosion for each island
#fish surveys were stratified random surveys around each island so site-level data
#does not match up with benthic survey data and instead we will calculate island level data

#separate Bolbometopon muricatum and Calotomus spp. before calculating a and b following Lange et al (2020) https://doi.org/10.3390/d12100379 where Bioerosion rates (kg/ind/yr) = a * (TL^b)
RB_parrots_others=RB_parrotfish[1:3,]
RB_parrots_bioeroders=RB_parrotfish[4:nrow(RB_parrotfish),]

#determine a and b for Scaridae and Scarus sp to match survey data
Scaridae=RB_parrots_bioeroders[11,]
Scaridae[1,1]="Scaridae"
Scarus_sp=RB_parrots_bioeroders[11,]
Scarus_sp[1,1]="Scarus sp"

#add Scaridae and Scarus sp means to RB_parrots_bioeroders dataframe
RB_parrots_bioeroders2=rbind(RB_parrots_bioeroders,Scaridae,Scarus_sp)

#convert from wide to long format to account for initial and terminal phase fish
parrotfish_long=gather(RB_parrots_bioeroders2, size_bin, bioerosion, 4:12, factor_key=TRUE)

#extract mid-point of size bins to determine mean length to assosociate with the respective bioerosion rate
parrotfish_long$cm=(as.numeric(substr(parrotfish_long$size_bin,1,2))+
                      as.numeric(substr(parrotfish_long$size_bin,4,5)))/2

#create a linear fit between log(BioerosionRate) and log(FishLength), extract the slope and intercept of the respective fits, and then convert them to a nd b following Bioerosion rates (kg/ind/yr) = a * (TL^b)
bioerosion_scaling=
  parrotfish_long %>% 
  group_by(Species) %>% 
  do(model = lm(log(bioerosion) ~ log(cm), data=.)) %>% 
  mutate(a=exp(coef(model)[1]), b=coef(model)[2]) %>%
  select(-model)

#determine a and b for Bolbometopon muricatum and Calotomus spp.
#Bolbometopon muricatum a is the fixed 3428 rate for all size classes and b is 0 since the rate applies regardless of size
#Rates of zero are applied to grazing Calatomus spp.
other_scaling=RB_parrots_others[1:3,c(1,4,5)]
other_scaling[1,3]=0
colnames(other_scaling)=c("Species","a","b")

#merge all parrotfish a and b data
parrotfish_a_b=rbind(other_scaling,bioerosion_scaling)

#merge parrotfish a and b data with the originating ReefBudget size class data
RB_parrotfish_a_b=merge(rbind(RB_parrots_bioeroders2,RB_parrots_others),parrotfish_a_b,by="Species")

#here we will calculate parrotfish and urchin bioerosion for each island
#fish surveys were stratified random surveys around each island so site-level data
#does not match up with benthic survey data and instead we will calculate island level data

#separate Bolbometopon muricatum and Calotomus spp. before calculating a and b following Lange et al (2020) https://doi.org/10.3390/d12100379 where Bioerosion rates (kg/ind/yr) = a * (TL^b)
RB_parrots_others=RB_parrotfish[1:3,]
RB_parrots_bioeroders=RB_parrotfish[4:nrow(RB_parrotfish),]

#determine a and b for Scaridae and Scarus sp to match survey data
Scaridae=RB_parrots_bioeroders[11,]
Scaridae[1,1]="Scaridae"
Scarus_sp=RB_parrots_bioeroders[11,]
Scarus_sp[1,1]="Scarus sp"

#add Scaridae and Scarus sp means to RB_parrots_bioeroders dataframe
RB_parrots_bioeroders2=rbind(RB_parrots_bioeroders,Scaridae,Scarus_sp)

#convert from wide to long format to account for initial and terminal phase fish
parrotfish_long=gather(RB_parrots_bioeroders2, size_bin, bioerosion, 4:12, factor_key=TRUE)

#extract mid-point of size bins to determine mean length to assosociate with the respective bioerosion rate
parrotfish_long$cm=(as.numeric(substr(parrotfish_long$size_bin,1,2))+
                      as.numeric(substr(parrotfish_long$size_bin,4,5)))/2

#create a linear fit between log(BioerosionRate) and log(FishLength), extract the slope and intercept of the respective fits, and then convert them to a nd b following Bioerosion rates (kg/ind/yr) = a * (TL^b)
bioerosion_scaling=
  parrotfish_long %>% 
  group_by(Species) %>% 
  do(model = lm(log(bioerosion) ~ log(cm), data=.)) %>% 
  mutate(a=exp(coef(model)[1]), b=coef(model)[2]) %>%
  select(-model)

#determine a and b for Bolbometopon muricatum and Calotomus spp.
#Bolbometopon muricatum a is the fixed 3428 rate for all size classes and b is 0 since the rate applies regardless of size
#Rates of zero are applied to grazing Calatomus spp.
other_scaling=RB_parrots_others[1:3,c(1,4,5)]
other_scaling[1,3]=0
colnames(other_scaling)=c("Species","a","b")

#merge all parrotfish a and b data
parrotfish_a_b=rbind(other_scaling,bioerosion_scaling)

#merge parrotfish a and b data with the originating ReefBudget size class data
RB_parrotfish_a_b=merge(rbind(RB_parrots_bioeroders2,RB_parrots_others),parrotfish_a_b,by="Species")

## SITE ESTIMATES
wsd_bioer <- wd_OCC %>% left_join(parrotfish_a_b, by = c("TAXONNAME" = "Species")) %>% # attach a and b factors to dataset
  # calculate bioerosion per record
  mutate(eros_m2 = COUNT*a*(SIZE_^b)/(pi*7.5^2)) %>% 
  mutate(eros_m2 = ifelse(is.na(eros_m2), 0, eros_m2)) %>% # replace NAs with 0s 
  # pool up to site estimates
  mutate(Fish.grp = ifelse(SPECIES == "NONFOCAL.SP", "NONFOCAL.SPP", "FOCAL.SPP")) %>% 
  group_by(REGION, ISLAND, SEC_NAME, REEF_ZONE, DEPTH_BIN, OBS_YEAR, SITEVISITID, METHOD, REP, REPLICATEID, Fish.grp) %>%
  dplyr::summarize(eros.repIDtot = sum(eros_m2)) %>% as.data.frame() %>% # totals per repID
  spread(Fish.grp, eros.repIDtot, fill = 0) %>% gather("Fish.grp", "eros.repIDtot", FOCAL.SPP, NONFOCAL.SPP) %>%   # expand to include instances when no observations of focal/nonfocal spp
  group_by(REGION, ISLAND, SEC_NAME, REEF_ZONE, DEPTH_BIN, OBS_YEAR, SITEVISITID, METHOD, REP, Fish.grp) %>% dplyr::summarize(eros.rep = mean(eros.repIDtot)) %>% as.data.frame() %>% # avg across repIDs
  group_by(REGION, ISLAND, SEC_NAME, REEF_ZONE, DEPTH_BIN, OBS_YEAR, SITEVISITID, METHOD, Fish.grp) %>% dplyr::summarize(eros.site = mean(eros.rep)) %>% as.data.frame() %>% # avg across reps
  # focal species only
  filter(Fish.grp == "FOCAL.SPP") %>% select(-Fish.grp) %>%
  # attach strata area
  left_join(strata.OCC %>% select(SEC_NAME, DEPTH_BIN, REEF_ZONE, AREA_HA), by = c("SEC_NAME", "REEF_ZONE", "DEPTH_BIN"))

## STRATA ESTIMATES
strata_bioer <- wsd_bioer %>% group_by(REGION, ISLAND, SEC_NAME, REEF_ZONE, DEPTH_BIN, OBS_YEAR, AREA_HA) %>%
  dplyr::summarize(strata.mean = mean(eros.site), # average each metric across sites per strata
                   strata.var = var(eros.site), # calculate variance
                   strata.N = n_distinct(SITEVISITID), # calculate # of sites per strata
                   strata.se = sqrt(strata.var)/sqrt(strata.N)) %>% as.data.frame() # calculate SE

# REMOVE strata with N<2 --> just GUA_TUMON
strata_drop <- strata_bioer %>% filter(strata.N>1)

# ISLAND ESTIMATES
pool.totals <- strata_drop %>% group_by(REGION, ISLAND, REEF_ZONE, DEPTH_BIN, OBS_YEAR) %>% dplyr::summarize(TOT_AREA_HA = sum(AREA_HA), tot.N = sum(strata.N)) %>% as.data.frame()

# calculate weighted values per strata
strata_wt <- strata_drop %>% left_join(pool.totals, by=c("REGION", "ISLAND", "REEF_ZONE", "DEPTH_BIN", "OBS_YEAR")) %>% # attach totals
  mutate(grid_size = 50*50, area_units = 100*100) %>% # add grid cell size and units of area = hectares
  mutate(propSampled = (strata.N*grid_size) / (AREA_HA*area_units)) %>%  # calculate proportion of total area sampled per strata
  mutate(wt = AREA_HA/TOT_AREA_HA) %>% # calculate proportional wt per strata: AREA_HA / total AREA_HA per grouping (= wts per strata)
  # calculate weighted values
  mutate(wt.mean = strata.mean*wt) %>% # weighted means
  mutate(wt.var = strata.var*((wt^2)/strata.N)*(1-propSampled))  # weighted var: var*(N^2/N), then adjust for proportion of total area sampled (using Krebs formula 8.18 p276)

# CHECK for propSampled > 1 (*NOTE: if any values >1, then cap at 1)
strata_wt %>% filter(propSampled > 1) # 0 rows --> good to go!

# pool across strata to desired groupings
island_parrot_bioer <- strata_wt %>% group_by(REGION, ISLAND, REEF_ZONE, DEPTH_BIN, OBS_YEAR, TOT_AREA_HA, tot.N) %>% 
  filter(REEF_ZONE=="Forereef") %>% # select only Forereef sites
  dplyr::summarize(pool.mean = sum(wt.mean), # sum across strata: pooled means
                   pool.var = sum(wt.var),  # sum across strata: this is effectively pooled variance of sample means for entire domain
                   pool.se = sqrt(pool.var)) %>% as.data.frame() # SD of sample means for entire domain (=equiv to sample SE)  

# calculate island level parrotfish bioerosion as the mean ± pool.se
parrotfish = with(island_parrot_bioer,cbind(island=ISLAND,islandB=pool.mean,islandB_lwr=pool.mean-pool.se,islandB_upr=pool.mean+pool.se)) %>% as.data.frame() 

#regrettably, there was inadequate sea urchin survey data for each region
#instead, find relative abundance of urchins for each region
(table(MARIANAS_fishes$URCHIN_DACOR)/sum(table(MARIANAS_fishes$URCHIN_DACOR)))*100 #93.6% rare
(table(PRIA_fishes$URCHIN_DACOR)/sum(table(PRIA_fishes$URCHIN_DACOR)))*100 #82.8% rare
(table(SAMOA_fishes$URCHIN_DACOR)/sum(table(SAMOA_fishes$URCHIN_DACOR)))*100 #99.8% rare
(table(WAKE_fishes$URCHIN_DACOR)/sum(table(WAKE_fishes$URCHIN_DACOR)))*100 #87.5% rare
(table(NWHI_fishes$urchin_dacor)/sum(table(NWHI_fishes$urchin_dacor)))*100 #74.8% rare
#conclude that urchins are rare across all surveys and therefore likely had minimal contribution to budget

###########Calculate NET CARBONATE BUDGET##########
#merge CRAMP_siteG data with parrotfish bioerosion data and calculate net budget
CRAMP_G=merge(CRAMP_siteG,parrotfish,by.x="ISLAND",by.y="island")

#we assumed that 50% of mechanical bioerosion by parrotfish is reincoporated into the reef framework following 
#Hubbard et al (1990) https://doi.org/10.1306/212F9197-2B24-11D7-8648000102C1865D
#Perry et al (2018) https://doi.org/10.1038/s41586-018-0194-z
#we therefore multiply the parrotfish bioerosion loss term by 50% to account for this 50% reincorporation in the calculation of our net carbonate budgets

#G is site carbonate production minus island parrotfish bioerosion * 50% reincorporation
CRAMP_G$G=with(CRAMP_G,siteG-as.numeric(islandB)*0.5)
#G_lwr is lower site carbonate production minus upper island parrotfish bioerosion * 50% reincorporation
CRAMP_G$G_lwr=with(CRAMP_G,siteG_lwr-as.numeric(islandB_upr)*0.5)
#G_upr is upper site carbonte production minus lower island parrotfish bioerosion * 50% reincorporation
CRAMP_G$G_upr=with(CRAMP_G,siteG_upr-as.numeric(islandB_lwr)*0.5) 

###########Calculate dnTA##########

#isolate climate_sites data and endmembers data
climate_sites=data[data$ASSOC_OCC_SITEID!="",]
endmembers=data[data$SURVEY_DESIGN=="End_member",]
#calculate mean TA, DIC, and S endmembers for each cruise and location
mean_em=endmembers %>%
  group_by(CRUISE_ID,LOCATIONCODE) %>%
  summarize(meanTAem=mean(TALK_UMOL_KG),meanDICem=mean(DIC_UMOL_KG),meanSem=mean(SALINITY_PSS78_BEST))
#merge mean endmembers with climate site data
climate_sites_em=merge(climate_sites,mean_em,by=c("CRUISE_ID","LOCATIONCODE"))

#majority of reef sites have higher reef salinity than ocean salinity so may be reasonable to assume that evaporation is the dominant process and TA=0@S=0
Srangeplot=
  ggplot()+
  geom_point(data=climate_sites_em,aes(x=LOCATIONCODE,y=SALINITY_PSS78_BEST-meanSem,color=LOCATIONCODE),show.legend=FALSE,size=3,alpha=0.7)+
  geom_hline(yintercept=0)+
  xlab("Location")+ylab(expression(S['reef']~-S['offshore']))+
  theme_minimal()+theme(text = element_text(size=20))

#calculate dnTA assuming mean precip, river, and runoff TA@S=0 from Courtney et al (2021) https://doi.org/10.1371/journal.pone.0261210
#calculate dnTA assuming TA=15@S=0
climate_sites_em$dnta_ocean_15=with(climate_sites_em,meanTAem-((TALK_UMOL_KG-15)*(meanSem/SALINITY_PSS78_BEST)+15))
climate_sites_em$dnta_reef_15=with(climate_sites_em,((meanTAem-15)*(SALINITY_PSS78_BEST/meanSem)+15)-TALK_UMOL_KG)
climate_sites_em$dnta_15=(climate_sites_em$dnta_reef_15+climate_sites_em$dnta_ocean_15)/2
#calculate dnTA assuming TA=1298@S=0
climate_sites_em$dnta_ocean_1298=with(climate_sites_em,meanTAem-((TALK_UMOL_KG-1298)*(meanSem/SALINITY_PSS78_BEST)+1298))
climate_sites_em$dnta_reef_1298=with(climate_sites_em,((meanTAem-1298)*(SALINITY_PSS78_BEST/meanSem)+1298)-TALK_UMOL_KG)
climate_sites_em$dnta_1298=(climate_sites_em$dnta_reef_1298+climate_sites_em$dnta_ocean_1298)/2
#calculate dnTA assuming TA=586@S=0
climate_sites_em$dnta_ocean_586=with(climate_sites_em,meanTAem-((TALK_UMOL_KG-586)*(meanSem/SALINITY_PSS78_BEST)+586))
climate_sites_em$dnta_reef_586=with(climate_sites_em,((meanTAem-586)*(SALINITY_PSS78_BEST/meanSem)+586)-TALK_UMOL_KG)
climate_sites_em$dnta_586=(climate_sites_em$dnta_reef_586+climate_sites_em$dnta_ocean_586)/2
OCC_sitelevel_dnta=climate_sites_em %>%
  group_by(CRUISE_ID,ASSOC_OCC_SITEID) %>%
  summarize(lat=LATITUDE_DEC[1],
            lon=LONGITUDE_DEC[1],
            time=TIMESTAMP_UTC[1],
            dnta=mean(c(dnta_15,dnta_586,dnta_1298)),
            dnta_lwr=min(c(dnta_ocean_15,dnta_ocean_586,dnta_ocean_1298,dnta_reef_15,dnta_reef_586,dnta_reef_1298)),
            dnta_upr=max(c(dnta_ocean_15,dnta_ocean_586,dnta_ocean_1298,dnta_reef_15,dnta_reef_586,dnta_reef_1298)))

#calculate time of data of dnta samples in local time
local_time=OCC_sitelevel_dnta %>% 
  summarise(timezone = tz_lookup_coords(lat = lat, lon = lon, method = "accurate"),
            local_time = map2(.x = time, .y = timezone, 
                              .f = function(x, y) {with_tz(time = x, tzone = y)})) %>% 
  unnest(local_time) %>% 
  select(local_time)
OCC_sitelevel_dnta$local_time=local_time$local_time

###########Calculate CDP,BCC,CVI########
#rename OCC_sitelevel_dnta site level index prior to merging data
colnames(OCC_sitelevel_dnta)[colnames(OCC_sitelevel_dnta) == "ASSOC_OCC_SITEID"] <- "OCC_SITEID"

#add mission columns from CRUISE IDs prior to merging
OCC_sitelevel_dnta$MISSION=substr(OCC_sitelevel_dnta$CRUISE_ID, 1, 6)
CRAMP_G$MISSION=substr(CRAMP_G$MISSIONID, 1, 6)

#merge CRAMP_G and OCC_sitelevel_dnta dataframes for matching G and dnTA data and select only mid depth sites
mergedCRAMP=merge(CRAMP_G,OCC_sitelevel_dnta,by=c('MISSION','OCC_SITEID')) %>% filter(DEPTH_BIN=="Mid")
#remove extraneous columns
mergedCRAMP <- mergedCRAMP[ , ! names(mergedCRAMP) %in% c("MISSIONID", "CRUISE_ID")]

#Calculate CDP, BCC, and CVI indices
CRAMP_OCC=select(mergedCRAMP,mission=MISSION,OCC=OCC_SITEID,region=REGION,island=ISLAND,site=SITE,year=OBS_YEAR,lat=LATITUDE,long=LONGITUDE,
                 zone=REEF_ZONE,depth=DEPTH_BIN,calcifiers=CALCIFIERS,coral=CORAL,cca=CCA,macroalgae=MACROALGAE,encrusting=ENCRUSTING,halimeda=HALIMEDA,
                 softcoral=SOFTCORAL,turf=TURF,sediment=SEDIMENT,invertebrate=INVERT,siteG,siteG_lwr,siteG_upr,islandB,islandB_lwr,islandB_upr,G,G_lwr,G_upr,dnta,dnta_lwr,dnta_upr)

CRAMP_OCC$BCC=if_else(CRAMP_OCC$G_upr < 0, "negative", 
                      if_else(CRAMP_OCC$G_lwr > 0, "positive", "neutral"))

CRAMP_OCC$CDP=if_else(CRAMP_OCC$dnta_upr < 0, "negative", 
                      if_else(CRAMP_OCC$dnta_lwr > 0, "positive", "neutral"))

CRAMP_OCC$CVI=if_else(CRAMP_OCC$CDP == "positive" & CRAMP_OCC$BCC == "positive", "positive",
                      if_else(CRAMP_OCC$CDP == "negative" & CRAMP_OCC$BCC == "negative", "negative","neutral"))

###########MAX DHW 2014-2017########
#The below commented out code chunk will download and compile daily DHW from NOAA servers for all of the CRAMP sites. It has been commented out because it takes a very long time to download all of the data from the servers and compile into the resulting DHW data frame for each site.Instead of running these time intensive steps, this script imports and uses the resulting DHW_OCC_2.csv file created by this code chunk

#load required libraries
#library(tidync)
#library(ncdf4)
#library(RNetCDF)

#setup range of dates to download
#day = format(as.Date(seq(as.Date("2014/1/1"), as.Date("2017/12/31"), "days"), format = "%d/%m/%Y"), "%Y%m%d")

#step through file structure for each year to download the daily DHW netCDF files from NOAA Coral Reef Watch
#for(i in 1:length(day)){
#  url <- paste0("ftp://ftp.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1/nc/v1.0/daily/dhw/", substr(day,1,4)[i], "/")
#  filenamelist = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE, crlf = TRUE) 
#  filenames = paste(url, strsplit(filenamelist, "\r*\n")[[1]], sep = "") 
#  filename = paste0("ct5km_dhw_v3.1_",day[i],".nc")
#  download.file(url = paste0(url,filename), destfile = filename, mode="wb")
#}

#DHW_OCC=NULL

#day = format(as.Date(seq(as.Date("2014/1/1"), as.Date("2017/12/31"), "days"), format = "%d/%m/%Y"), "%Y%m%d")
#for(i in 1:length(day)){
#  for(j in 1:nrow(CRAMP_OCC)){
#    date=day[i]
#    DHW=tidync(paste0("ct5km_dhw_v3.1_",day[i],".nc")) %>%  
#      hyper_filter(lon = index == which.min(abs(lon - CRAMP_OCC$long[j])), lat = index == which.min(abs(lat - CRAMP_OCC$lat[j]))) %>% 
#      hyper_tibble() %>% select(degree_heating_week) %>% first()
#    DHW = if (length(DHW)==0) {tidync(paste0("ct5km_dhw_v3.1_",day[i],".nc")) %>%  
#        hyper_filter(lon = index == which.min(abs(lon - CRAMP_OCC$long[j]+0.05)), lat = index == which.min(abs(lat - CRAMP_OCC$lat[j]+0.05))) %>% 
#        hyper_tibble() %>% select(degree_heating_week) %>% first()} else {DHW}
#    OCC=CRAMP_OCC$OCC[j]
#    DHW_OCC <-rbind(DHW_OCC, data.frame(date,DHW,OCC))
#  }
#}

#write.csv(DHW_OCC,'DHW_OCC_2.csv',row.names=FALSE)

#create year column for annual summaries
DHW_OCC$year=substr(DHW_OCC$date,1,4)
#annual max DHW experienced for each site
maxDHWannual=DHW_OCC %>%
  group_by(OCC,year) %>%
  summarize(maxDHW=max(DHW))
#max DHW experienced for each site
maxDHW=DHW_OCC %>%
  group_by(OCC) %>%
  summarize(maxDHW=max(DHW))

#calculate bleaching stress years and most recent bleaching stress year
years_annual_bleaching = maxDHWannual %>% spread(year, maxDHW) %>% 
  rename("DHW_2014"="2014","DHW_2015"="2015","DHW_2016"="2016","DHW_2017"="2017") %>% 
  mutate(last_year_bleached=ifelse(DHW_2017>4,2017,ifelse(DHW_2016>4,2016,ifelse(DHW_2015>4,2015,ifelse(DHW_2014>4,2014,0)))))

#merge bleaching stress experienced for each site with CRAMP_OCC dataframe
CRAMP_bleaching_years=arrange(merge(CRAMP_OCC,years_annual_bleaching,by="OCC"),region,island)

#merge max DHW experienced for each site with CRAMP_OCC dataframe
CRAMP=arrange(merge(CRAMP_bleaching_years,maxDHW,by="OCC"),region,island)

#add column for reefs that experienced ecologically significant bleaching (DHW≥4) and ecologically severe bleaching (DHW≥8) heat stress
CRAMP$bleaching=ifelse(CRAMP$maxDHW>=8,"severe",ifelse(CRAMP$maxDHW>=4,"significant","none"))

###########Summary of metrics########
mean_se(CRAMP$G,mult=1.96)
min(CRAMP$G)
max(CRAMP$G)

mean_se(CRAMP$dnta,mult=1.96)
min(CRAMP$dnta)
max(CRAMP$dnta)

length(CRAMP$BCC) #56 sites total

length(CRAMP$bleaching[CRAMP$bleaching=="severe"]) #44 sites experienced ecologically severe bleaching heat stress
length(CRAMP$bleaching[CRAMP$bleaching=="severe"])/length(CRAMP$BCC) #79%
length(CRAMP$bleaching[CRAMP$bleaching=="significant"]) #5 sites experienced ecologically significant bleaching heat stress
length(CRAMP$bleaching[CRAMP$bleaching=="significant"])/length(CRAMP$BCC) #9%
length(CRAMP$bleaching[CRAMP$bleaching=="none"]) #7 sites did not experience ecologically significant bleaching heat stress
length(CRAMP$bleaching[CRAMP$bleaching=="none"])/length(CRAMP$BCC) #13%

length(CRAMP$BCC[CRAMP$BCC=="positive"]) #43 sites positive BCC
length(CRAMP$BCC[CRAMP$BCC=="positive"])/length(CRAMP$BCC) #77% sites positive BCC
length(CRAMP$BCC[CRAMP$BCC=="neutral"]) #9 sites neutral BCC
length(CRAMP$BCC[CRAMP$BCC=="neutral"])/length(CRAMP$BCC) #16% sites neutral BCC
length(CRAMP$BCC[CRAMP$BCC=="negative"]) #4 sites negative BCC
length(CRAMP$BCC[CRAMP$BCC=="negative"])/length(CRAMP$BCC) #7% sites negative BCC

length(CRAMP$CDP[CRAMP$CDP=="positive"]) #47 sites positive CDP
length(CRAMP$CDP[CRAMP$CDP=="positive"])/length(CRAMP$BCC) #84% sites positive CDP
length(CRAMP$CDP[CRAMP$CDP=="neutral"]) #6 sites neutral CDP
length(CRAMP$CDP[CRAMP$CDP=="neutral"])/length(CRAMP$BCC) #11% sites neutral CDP
length(CRAMP$CDP[CRAMP$CDP=="negative"]) #3 sites negative CDP
length(CRAMP$CDP[CRAMP$CDP=="negative"])/length(CRAMP$BCC) #5% sites negative CDP

length(CRAMP$CVI[CRAMP$CVI=="positive"]) #38 sites positive CVI
length(CRAMP$CVI[CRAMP$CVI=="positive"])/length(CRAMP$BCC) #68% sites positive CVI
length(CRAMP$CVI[CRAMP$CVI=="neutral"]) #18 sites neutral CVI
length(CRAMP$CVI[CRAMP$CVI=="neutral"])/length(CRAMP$BCC) #32% sites neutral CVI
length(CRAMP$CVI[CRAMP$CVI=="negative"]) #0 sites negative CVI
length(CRAMP$CVI[CRAMP$CVI=="negative"])/length(CRAMP$BCC) #0% sites negative CVI
      
###########Correlations between parameters########
#change lme optimizer and increase interations to improve model convergence
ctrl=lmeControl(maxIter=100000,msMaxIter=100000,opt='optim')

#lme of of dnta~G with random slopes and intercepts by island by maximum likelihood
m1=lme(dnta~G,random=~G|island,data=CRAMP,control=ctrl,method="ML")
#lme of of G~coral with random slopes and intercepts by island by maximum likelihood
m2=lme(G~coral,random=~coral|island,data=CRAMP,control=ctrl,method="ML")
#lme of of coral~G with random slopes and intercepts by island by maximum likelihood
m3=lme(coral~G,random=~G|island,data=CRAMP,control=ctrl,method="ML")
#lme of of dnta~coral with random slopes and intercepts by island by maximum likelihood
m4=lme(dnta~coral,random=~coral|island,data=CRAMP,control=ctrl,method="ML") 
#lme of of coral~maxDHW with random slopes and intercepts by island by maximum likelihood
m5=lme(coral~maxDHW,random=~maxDHW|island,data=CRAMP,control=ctrl,method="ML")
#lme of of G~maxDHW with random slopes and intercepts by island by maximum likelihood
m6=lme(G~maxDHW,random=~maxDHW|island,data=CRAMP,control=ctrl,method="ML") 
#lme of of dnta~maxDHW with random slopes and intercepts by island by maximum likelihood
m7=lme(dnta~maxDHW,random=~maxDHW|island,data=CRAMP,control=ctrl,method="ML")

summary(m1) #G does not predict dnta
summary(m2) #coral cover predicts G
summary(m3) #G predicts coral cover
summary(m4) #coral does not predict dnta
summary(m5) #maxDHW predicts coral cover
summary(m6) #maxDHW does not predict G
summary(m7) #maxDHW predicts dnta

qqnorm(resid(m2)) #visual assessment of normality
plot(resid(m2)) #visual assessment of homoscedasticity

#create predicted data frame with 95% CI to plot
newdat.m2 = data.frame(coral = CRAMP$coral)
newdat.m2$predlme = predict(m2, newdata = newdat.m2, level = 0)
des.m2 = model.matrix(formula(m2)[-2], newdat.m2)
predvar = diag(des.m2 %*% vcov(m2) %*% t(des.m2) )
newdat.m2$lower = with(newdat.m2, predlme - 1.96*sqrt(predvar) )
newdat.m2$upper = with(newdat.m2, predlme + 1.96*sqrt(predvar) )

qqnorm(resid(m3)) #visual assessment of normality
plot(resid(m3)) #visual assessment of homoscedasticity

#create predicted data frame with 95% CI at G = 0 to calculate threshold for maintaining net positive budget
newdat.m3 = data.frame(G = 0)
newdat.m3$predlme = predict(m3, newdata = newdat.m3, level = 0)
des.m3 = model.matrix(formula(m3)[-2], newdat.m3)
predvar = diag(des.m3 %*% vcov(m3) %*% t(des.m3) )
newdat.m3$lower = with(newdat.m3, predlme - 1.96*sqrt(predvar) )
newdat.m3$upper = with(newdat.m3, predlme + 1.96*sqrt(predvar) )
newdat.m3$predlme #8.6 percent coral cover at G=0
newdat.m3$upper-newdat.m3$predlme #95% CI = 6.0 percent coral cover
#0 G threshold @ 8.6±6.0% coral cover

qqnorm(resid(m5)) #visual assessment of normality
plot(resid(m5)) #visual assessment of homoscedasticity

#create predicted data frame with 95% CI to plot
newdat.m5 = data.frame(maxDHW = CRAMP$maxDHW)
newdat.m5$predlme = predict(m5, newdata = newdat.m5, level = 0)
des.m5 = model.matrix(formula(m5)[-2], newdat.m5)
predvar = diag(des.m5 %*% vcov(m5) %*% t(des.m5) )
newdat.m5$lower = with(newdat.m5, predlme - 1.96*sqrt(predvar) )
newdat.m5$upper = with(newdat.m5, predlme + 1.96*sqrt(predvar) )

qqnorm(resid(m7)) #visual assessment of normality
plot(resid(m7)) #visual assessment of homoscedasticity

#create predicted data frame with 95% CI to plot
newdat.m7 = data.frame(maxDHW = CRAMP$maxDHW)
newdat.m7$predlme = predict(m7, newdata = newdat.m7, level = 0)
des.m7 = model.matrix(formula(m7)[-2], newdat.m7)
predvar = diag(des.m7 %*% vcov(m7) %*% t(des.m7) )
newdat.m7$lower = with(newdat.m7, predlme - 1.96*sqrt(predvar) )
newdat.m7$upper = with(newdat.m7, predlme + 1.96*sqrt(predvar) )

###########maxDHW vs. CDP,BCC,CVI Indices########

#calculate mean±SE maxDHW for +,±,- classification of each calcification index
BCC_coral_means_DHW=CRAMP %>% 
  group_by(BCC) %>% 
  summarize(mean_se(maxDHW,mult=1.96))
CDP_coral_means_DHW=CRAMP %>% 
  group_by(CDP) %>% 
  summarize(mean_se(maxDHW,mult=1.96))
CVI_coral_means_DHW=CRAMP %>% 
  group_by(CVI) %>% 
  summarize(mean_se(maxDHW,mult=1.96))

kruskal.test(maxDHW~as.factor(BCC),data=CRAMP) #detectable differences in maxDHW between BCC
dunnTest(maxDHW~as.factor(BCC),data=CRAMP,method="bonferroni") #BCC max DHW negative>positive, negative>neutral, neutral=positive

kruskal.test(maxDHW~as.factor(CDP),data=CRAMP) #no detectable differences in maxDHW between CDP

kruskal.test(maxDHW~as.factor(CVI),data=CRAMP) #no detectable differences in maxDHW between CVI

###########Calcifier Cover vs. CDP,BCC,CVI Indices########
#calculate mean±SE % calcifier cover for +,±,- of each calcification index
BCC_coral_means=CRAMP %>% 
  group_by(BCC) %>% 
  summarize(mean_se(coral,mult=1.96))
CDP_coral_means=CRAMP %>% 
  group_by(CDP) %>% 
  summarize(mean_se(coral,mult=1.96))
CVI_coral_means=CRAMP %>% 
  group_by(CVI) %>% 
  summarize(mean_se(coral,mult=1.96))

kruskal.test(coral~as.factor(BCC),data=CRAMP) #coral cover differs between BCC classifications
dunnTest(coral~as.factor(BCC),data=CRAMP,method="bonferroni") #BCC coral cover positive>negative, positive=neutral, neutral=negative

kruskal.test(coral~as.factor(CDP),data=CRAMP) #no detectable differences in coral cover between CDP classifications
dunnTest(coral~as.factor(CDP),data=CRAMP,method="bonferroni") #no detectable differences in coral cover between CDP classifications

kruskal.test(coral~as.factor(CVI),data=CRAMP) #coral cover differs between CVI classifications
dunnTest(coral~as.factor(CVI),data=CRAMP,method="bonferroni") #CVI coral cover positive>neutral and no negative classifications

###########Figure 2 Display map of study locations########
map <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))
shapes <- c(15,18,16,17,19)
map_plot=ggplot() +
  geom_polygon(data = map, aes(x=long, y = lat, group = group)) +
  geom_point(data = CRAMP, aes(x = ifelse(long < -25, long + 360, long), y = lat, color=region,shape=region),size=3)+
  scale_color_manual(values = c("MARIAN"="#a1dab4","NWHI"="#41b6c4","PRIAs"="#253494","SAMOA"="#2c7fb8"))+
  scale_shape_manual(values=shapes)+
  coord_fixed(ratio=1,xlim = c(135,210), ylim = c(-20,35),expand=TRUE)+
  annotate("text", x=145.3832, y=21, label= ~underline(Mariana~Islands),size=8,color="#a1dab4")+
  annotate("text", x=180, y=32, label= ~underline(NW~Hawaiian~Islands),size=8,color="#41b6c4")+
  annotate("text", x=189, y=-6, label= ~underline(American~Samoa),size=8,color="#2c7fb8")+
  annotate("text", x=184, y=12, label= ~underline(Pacific~Remote~Island~Areas),size=8,color="#253494")+
  annotate("text", x=141, y=13.5, label= "Guam",size=6,color="#a1dab4")+
  annotate("text", x=142, y=18.1, label= "Pagan",size=6,color="#a1dab4")+
  annotate("text", x=150, y=15.1, label= "Saipan",size=6,color="#a1dab4")+
  annotate("text", x=176, y=28.5, label= "Kure Atoll",size=6,color="#41b6c4")+
  annotate("text", x=191, y=26.5, label= "Lisianski",size=6,color="#41b6c4")+
  annotate("text", x=183.5256, y=-1.3050296, label= "Baker",size=6,color="#253494")+
  annotate("text", x=183.3894, y=2.8093400, label= "Howland",size=6,color="#253494")+
  annotate("text", x=199.9971, y=1.030700, label= "Jarvis",size=6,color="#253494")+
  annotate("text", x=197.6121, y=8.4388300, label= "Kingman",size=6,color="#253494")+
  annotate("text", x=197.9694, y=4.8637900, label= "Palmyra",size=6,color="#253494")+
  annotate("text", x=166.6272, y=18.3160500, label= "Wake",size=6,color="#253494")+
  annotate("text", x=188.9229, y=-9.0456887, label= "Swains",size=6,color="#2c7fb8")+
  annotate("text", x=195, y=-15, label= "Rose",size=6,color="#2c7fb8")+
  annotate("text", x=190.5821, y=-12.511600, label= "Tau",size=6,color="#2c7fb8")+
  annotate("text", x=185.3399, y=-15, label= "Tutuila",size=6,color="#2c7fb8")+
  ylab("Latitude (°N)")+
  xlab("Longitude  (°E)")+
  #labs(title = "(a)") +
  theme_classic()+
  theme(text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())

Figure2=map_plot
ggsave("Figure2.pdf",Figure2,width=9.5,height=9.5)

###########Figure 3 Display island level DHW, G, and dnTA########

DHWplot=merge(CRAMP_OCC,maxDHWannual,by="OCC") %>%
  group_by(region,lat,island,year.y)%>%
  summarize(islandDHW=max(maxDHW)) %>%
  mutate(island = fct_reorder(island,lat)) %>%
  ggplot()+  
  geom_bar(stat='identity',position='dodge',aes(x=island,y=islandDHW,group=year.y,fill=year.y,width=0.6),size=1)+
  geom_bracket(xmin = "Guam", xmax = "Pagan", y.position = 31,label = "Mariana Islands",color="#a1dab4",size=1,label.size=5,tip.length=0.01)+
  geom_bracket(xmin = "Lisianski", xmax = "Kure", y.position = 31,label = "NW Hawaiian Islands",color="#41b6c4",size=1,label.size=5,tip.length=0.01)+
  geom_bracket(xmin = "Jarvis", xmax = "Wake", y.position = 31,label = "Pacific Remote Island Areas",color="#253494",size=1,label.size=5,tip.length=0.01)+
  geom_bracket(xmin = "Rose", xmax = "Swains", y.position = 31,label = "American Samoa",color="#2c7fb8",size=1,label.size=5,tip.length=0.01)+
  geom_hline(yintercept=4,linetype="dashed",size=1,color="#fdae61")+
  geom_hline(yintercept=8,linetype="dashed",size=1,color="#d73027")+
  geom_hline(yintercept=0)+
  scale_y_continuous(expand = c(0,0), limits=c(-2,35),breaks = seq(0,30,10))+
  scale_fill_manual(values = c("2014"="#d9d9d9","2015"="#969696","2016"="#525252","2017"="#252525"))+
  xlab("Island")+
  ylab("DHW (2014-2017)")+
  ylab(expression(DHW[max]~'(2014-2017)'))+
  labs(title = "(a)") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(vjust = -8),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="top",
        legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

G_island_plot=
  CRAMP %>%
  group_by(region,lat,island,year)%>%
  mutate(island = fct_reorder(island,lat)) %>%
  ggplot()+  geom_hline(yintercept=0)+
  geom_hline(yintercept=0,color='red',size=1)+
  geom_pointrange(mapping=aes(x=island,y=G,ymin=G_lwr,ymax=G_upr,color=region),size=1)+
  scale_color_manual(labels = c("Mariana Islands", "NW Hawaiian Islands", "Pacific Remote Island Areas", "American Samoa"),values = c("MARIAN"="#a1dab4","NWHI"="#41b6c4","PRIAs"="#253494","SAMOA"="#2c7fb8"))+
  scale_shape_manual(labels = c("Mariana Islands", "NW Hawaiian Islands", "Pacific Remote Island Areas", "American Samoa"),values=shapes)+
  scale_y_continuous(expand = c(0,0), breaks = seq(-24,24,2))+
  coord_cartesian(ylim=c(-6,12))+
  xlab("")+
  ylab(expression(G~(kg~CaCO[3]~m^-2~yr^-1)))+
  labs(title = "(b)") +
  theme_classic()+
  theme(text = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

dTA_island_plot=
  CRAMP %>%
  group_by(region,lat,island,year)%>%
  mutate(island = fct_reorder(island,lat)) %>%
  ggplot()+  geom_hline(yintercept=0)+
  geom_pointrange(mapping=aes(x=island,y=dnta,ymin=dnta_lwr,ymax=dnta_upr,color=region),size=1)+
  scale_color_manual(labels = c("Mariana Islands", "NW Hawaiian Islands", "Pacific Remote Island Areas", "American Samoa"),values = c("MARIAN"="#a1dab4","NWHI"="#41b6c4","PRIAs"="#253494","SAMOA"="#2c7fb8"))+
  scale_shape_manual(labels = c("Mariana Islands", "NW Hawaiian Islands", "Pacific Remote Island Areas", "American Samoa"),values=shapes)+
  scale_y_continuous(expand = c(0,0), limits=c(-150,302),breaks = seq(-150,300,50))+
  xlab("Island")+
  ylab(expression(dnTA~(µmol~kg^-1)))+
  labs(title = "(c)") +
  theme_classic()+
  theme(text = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

Figure3=DHWplot/G_island_plot/dTA_island_plot
ggsave("Figure3.pdf",Figure3,width=12,height=12)

###########Figure 4 Display summary of BCC, CDP, and CVI########
BCC_plot=CRAMP %>%
  group_by(region,lat, island,BCC) %>%
  summarize(reefs=length(BCC)) %>%
  mutate(island = fct_reorder(island,lat)) %>%
  ggplot()+  
  geom_bar(stat='identity',aes(x=island,y=reefs,group=BCC,fill=BCC,width=0.6),size=1)+
  geom_hline(yintercept=0)+
  geom_bracket(xmin = "Guam", xmax = "Pagan", y.position = 18,label = "Mariana Islands",color="#a1dab4",size=1,label.size=5,tip.length=0.01)+
  geom_bracket(xmin = "Lisianski", xmax = "Kure", y.position = 18,label = "NW Hawaiian Islands",color="#41b6c4",size=1,label.size=5,tip.length=0.01)+
  geom_bracket(xmin = "Jarvis", xmax = "Wake", y.position = 18,label = "Pacific Remote Island Areas",color="#253494",size=1,label.size=5,tip.length=0.01)+
  geom_bracket(xmin = "Rose", xmax = "Swains", y.position = 18,label = "American Samoa",color="#2c7fb8",size=1,label.size=5,tip.length=0.01)+
  scale_fill_manual(labels = c("Negative","Neutral","Positive"),values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"),guide = guide_legend(reverse = TRUE))+
  scale_y_continuous(expand = c(0,0), limits=c(0,20),breaks = seq(0,20,2))+
  ylab("Number of Reefs")+
  labs(title = "(a) Census: Net Carbonate Budget") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16, vjust = -8),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="top",
        legend.box = "vertical",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

CDP_plot=CRAMP %>%
  group_by(region,lat, island,CDP) %>%
  summarize(reefs=length(CDP)) %>%
  mutate(island = fct_reorder(island,lat)) %>%
  ggplot()+  
  geom_bar(stat='identity',aes(x=island,y=reefs,group=CDP,fill=CDP,width=0.6),size=1)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  scale_y_continuous(expand = c(0,0), limits=c(0,20),breaks = seq(0,20,2))+
  ylab("Number of Reefs")+
  labs(title = "(b) Chemistry: Net Calcification") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

CVI_plot=CRAMP %>%
  group_by(region,lat, island,CVI) %>%
  summarize(reefs=length(CVI)) %>%
  mutate(island = fct_reorder(island,lat)) %>%
  ggplot()+  
  geom_bar(stat='identity',aes(x=island,y=reefs,group=CVI,fill=CVI,width=0.6),size=1)+
  geom_hline(yintercept=0)+
  scale_fill_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  scale_y_continuous(expand = c(0,0), limits=c(0,20),breaks = seq(0,20,2))+
  ylab("Number of Reefs")+
  labs(title = "(c) Calcification Vulnerability Index") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

Figure4=BCC_plot/CDP_plot/CVI_plot
ggsave("Figure4.pdf",Figure4,width=11.75,height=10.5)

###########Figure 5 coral cover, G, dnTA, and maxDHW plots########
G_nta_plot=CRAMP %>%
  ggplot()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_pointrange(mapping=aes(x=G,y=dnta,ymin=dnta_lwr,ymax=dnta_upr),colour="#666666",size=0.25)+
  geom_errorbarh(mapping=aes(y=dnta,xmin=G_lwr,xmax=G_upr),colour="#666666",size=0.25)+
  scale_y_continuous(expand = c(0,0), limits=c(-150,302),breaks = seq(-150,300,50))+
  scale_x_continuous(expand = c(0,0), breaks = seq(-24,24,2))+
  coord_cartesian(xlim=c(-6,12))+
  ylab(expression(dnTA~(µmol~kg^-1)))+
  xlab(expression(G~(kg~CaCO[3]~m^-2~yr^-1)))+
  labs(title = "(a)") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

coral_G_plot=CRAMP %>%
  ggplot()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_pointrange(mapping=aes(x=coral,y=G,ymin=G_lwr,ymax=G_upr),colour="#666666",size=0.25)+
  geom_ribbon(data = newdat.m2, aes(x=CRAMP$coral,y=CRAMP$G,ymin = lower, ymax = upper, color = NULL),alpha = .15) +
  geom_line(data = newdat.m2, aes(x=CRAMP$coral,y = predlme), size = .75)+
  scale_y_continuous(expand = c(0,0), breaks = seq(-24,24,2))+
  coord_cartesian(ylim=c(-6,12))+
  scale_x_continuous(expand = c(0,0), limits=c(0,80),breaks = seq(0,80,10))+
  xlab("Percent Coral Cover")+
  ylab(expression(G~(kg~CaCO[3]~m^-2~yr^-1)))+
  labs(title = "(b)") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

coral_dnta_plot=CRAMP %>%
  ggplot()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_pointrange(mapping=aes(x=coral,y=dnta,ymin=dnta_lwr,ymax=dnta_upr),colour="#666666",size=0.25)+
  scale_y_continuous(expand = c(0,0), limits=c(-150,302),breaks = seq(-150,300,50))+
  scale_x_continuous(expand = c(0,0), limits=c(0,80),breaks = seq(0,80,10))+
  ylab(expression(dnTA~(µmol~kg^-1)))+
  xlab("Percent Coral Cover")+
  labs(title = "(c)") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

maxDHW_coral_plot=
  CRAMP %>%
  ggplot()+
  geom_vline(xintercept=4,size=1,color="#fdae61")+
  geom_vline(xintercept=8,size=1,color="#d73027")+
  geom_point(mapping=aes(x=maxDHW,y=coral),colour="#666666",size=1)+
  geom_ribbon(data = newdat.m5, aes(x=CRAMP$maxDHW,y=CRAMP$coral,ymin = lower, ymax = upper, color = NULL),alpha = .15) +
  geom_line(data = newdat.m5, aes(x=CRAMP$maxDHW,y = predlme), size = .75)+
  scale_x_continuous(expand = c(0,0), limits=c(0,32),breaks = seq(0,32,4))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,80,10))+
  coord_cartesian(ylim=c(0, 80))+
  xlab("Maximum Degree Heating Weeks")+
  ylab("Percent Coral Cover")+
  labs(title = "(d)") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

maxDHW_G_plot=
  CRAMP %>%
  ggplot()+
  geom_vline(xintercept=4,size=1,color="#fdae61")+
  geom_vline(xintercept=8,size=1,color="#d73027")+
  geom_hline(yintercept=0)+
  geom_pointrange(mapping=aes(x=maxDHW,y=G,ymin=G_lwr,ymax=G_upr),colour="#666666",size=0.25)+
  scale_x_continuous(expand = c(0,0), limits=c(0,32),breaks = seq(0,32,4))+
  scale_y_continuous(expand = c(0,0), breaks = seq(-24,24,2))+
  coord_cartesian(ylim=c(-6,12))+
  xlab("Maximum Degree Heating Weeks")+
  ylab(expression(G~(kg~CaCO[3]~m^-2~yr^-1)))+
  labs(title = "(e)") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

maxDHW_dnta_plot=
  CRAMP %>%
  ggplot()+
  geom_vline(xintercept=4,size=1,color="#fdae61")+
  geom_vline(xintercept=8,size=1,color="#d73027")+
  geom_hline(yintercept=0)+
  geom_pointrange(mapping=aes(x=maxDHW,y=dnta,ymin=dnta_lwr,ymax=dnta_upr),colour="#666666",size=0.25)+
  geom_ribbon(data = newdat.m7, aes(x=CRAMP$maxDHW,y=CRAMP$dnta,ymin = lower, ymax = upper, color = NULL),alpha = .15) +
  geom_line(data = newdat.m7, aes(x=CRAMP$maxDHW,y = predlme), size = .75)+
  scale_x_continuous(expand = c(0,0), limits=c(0,32),breaks = seq(0,32,4))+
  scale_y_continuous(expand = c(0,0), limits=c(-150,302),breaks = seq(-150,300,50))+
  xlab("Maximum Degree Heating Weeks")+
  ylab(expression(dnTA~(µmol~kg^-1)))+
  labs(title = "(f)") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

Figure5=(G_nta_plot + coral_G_plot + coral_dnta_plot)/
  (maxDHW_coral_plot + maxDHW_G_plot + maxDHW_dnta_plot)

ggsave("Figure5.pdf",Figure5,width=15,height=7)

###########Figure 6 Percent calcifier cover vs BCC, CDP, and CVI########
sigFunc = function(x){
  if(x < 0.01){"p < 0.05"}
  else{NA}}

BCC_DHW_plot=ggplot()+
  geom_point(data=CRAMP,mapping=aes(x=BCC,y=maxDHW,color=BCC),size=2,alpha=0.6)+
  geom_pointrange(data=BCC_coral_means_DHW,mapping=aes(x=BCC,y=y,ymin=ymin,ymax=ymax,color=BCC),size=1.5,shape=15)+
  geom_signif(data=CRAMP,mapping=aes(x=BCC,y=maxDHW),
              comparisons = list(c("neutral", "negative"),
                                 c("positive", "negative"),
                                 c("positive", "neutral")),
              test = "wilcox.test", step_increase = 0.075,map_signif_level = sigFunc, tip_length = 0.01)+
  scale_x_discrete(limits=c("positive","neutral","negative"),labels = c("Positive","Neutral","Negative"))+
  scale_y_continuous(expand = c(0,0), limits=c(-1,40),breaks = seq(0,100,10))+
  scale_color_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  #annotate("text",x=1,y=21,label="a",size=5)+
  #annotate("text",x=2,y=32,label="a",size=5)+
  #annotate("text",x=3,y=32,label="b",size=5)+
  ylab(expression(DHW[max]~'(2014-2017)'))+
  labs(title = "(a) Census: Net Carbonate Budget") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

CDP_DHW_plot=ggplot()+
  geom_point(data=CRAMP,mapping=aes(x=CDP,y=maxDHW,color=CDP),size=2,alpha=0.6)+
  geom_pointrange(data=CDP_coral_means_DHW,mapping=aes(x=CDP,y=y,ymin=ymin,ymax=ymax,color=CDP),size=1.5,shape=15)+
  geom_signif(data=CRAMP,mapping=aes(x=CDP,y=maxDHW),
              comparisons = list(c("positive", "neutral"),
                                 c("neutral", "negative"),
                                 c("positive", "negative")),
              test = "wilcox.test", step_increase = 0.075,map_signif_level = sigFunc, tip_length = 0.01)+
  scale_x_discrete(limits=c("positive","neutral","negative"),labels = c("Positive","Neutral","Negative"))+
  scale_y_continuous(expand = c(0,0), limits=c(-1,40),breaks = seq(0,100,10))+
  scale_color_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  ylab(expression(DHW[max]~'(2014-2017)'))+
  labs(title = "(b) Chemistry: Net Calcification") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

CVI_DHW_plot=ggplot()+
  geom_point(data=CRAMP,mapping=aes(x=CVI,y=maxDHW,color=CVI),size=2,alpha=0.6)+
  geom_pointrange(data=CVI_coral_means_DHW,mapping=aes(x=CVI,y=y,ymin=ymin,ymax=ymax,color=CVI),size=1.5,shape=15)+
  geom_signif(data=CRAMP,mapping=aes(x=CVI,y=maxDHW),
              comparisons = list(c("positive", "neutral")),
              test = "wilcox.test", step_increase = 0.075,map_signif_level = sigFunc, tip_length = 0.01)+
  scale_x_discrete(limits=c("positive","neutral","negative"),labels = c("Positive","Neutral","Negative"))+
  scale_y_continuous(expand = c(0,0), limits=c(-1,40),breaks = seq(0,100,10))+
  scale_color_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  ylab(expression(DHW[max]~'(2014-2017)'))+
  labs(title = "(c) Calcification Vulnerability Index") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_blank(),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

BCCcoralplot=ggplot()+
  geom_point(data=CRAMP,mapping=aes(x=BCC,y=coral,color=BCC),size=2,alpha=0.6)+
  geom_pointrange(data=BCC_coral_means,mapping=aes(x=BCC,y=y,ymin=ymin,ymax=ymax,color=BCC),size=1.5,shape=15)+
  geom_signif(data=CRAMP,mapping=aes(x=BCC,y=coral),
              comparisons = list(c("positive", "negative"),
                                 c("neutral", "negative"),
                                 c("positive", "neutral")),
              test = "wilcox.test", step_increase = 0.075,map_signif_level = sigFunc, tip_length = 0.01)+
  scale_x_discrete(limits=c("positive","neutral","negative"),labels = c("Positive","Neutral","Negative"))+
  scale_y_continuous(expand = c(0,0), limits=c(-2,90),breaks = seq(0,100,10))+
  scale_color_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  #annotate("text",x=1,y=60,label="a",size=5)+
  #annotate("text",x=2,y=77,label="a,b",size=5)+
  #annotate("text",x=3,y=8,label="b",size=5)+
  ylab("Percent Coral Cover")+
  labs(title = "(d) Census: Net Carbonate Budget") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

CDPcoralplot=ggplot()+
  geom_point(data=CRAMP,mapping=aes(x=CDP,y=coral,color=CDP),size=2,alpha=0.6)+
  geom_pointrange(data=CDP_coral_means,mapping=aes(x=CDP,y=y,ymin=ymin,ymax=ymax,color=CDP),size=1.5,shape=15)+
  geom_signif(data=CRAMP,mapping=aes(x=CDP,y=coral),
              comparisons = list(c("positive", "neutral"),
                                 c("neutral", "negative"),
                                 c("positive", "negative")),
              test = "wilcox.test", step_increase = 0.075,map_signif_level = sigFunc, tip_length = 0.01)+
  scale_x_discrete(limits=c("positive","neutral","negative"),labels = c("Positive","Neutral","Negative"))+
  scale_y_continuous(expand = c(0,0), limits=c(-2,90),breaks = seq(0,100,10))+
  scale_color_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  ylab("Percent Coral Cover")+
  labs(title = "(e) Chemistry: Net Calcification") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

CVIcoralplot=ggplot()+
  geom_point(data=CRAMP,mapping=aes(x=CVI,y=coral,color=CVI),size=2,alpha=0.6)+
  geom_pointrange(data=CVI_coral_means,mapping=aes(x=CVI,y=y,ymin=ymin,ymax=ymax,color=CVI),size=1.5,shape=15)+
  geom_signif(data=CRAMP,mapping=aes(x=CVI,y=coral),
              comparisons = list(c("positive", "neutral")),
              test = "wilcox.test", step_increase = 0.075,map_signif_level = sigFunc, tip_length = 0.01)+
  #annotate("text",x=1,y=60,label="a",size=5)+
  #annotate("text",x=2,y=78,label="b",size=5)+
  #annotate("text",x=3,y=8,label="b",size=5)+
  scale_x_discrete(limits=c("positive","neutral","negative"),labels = c("Positive","Neutral","Negative"))+
  scale_y_continuous(expand = c(0,0), limits=c(-2,90),breaks = seq(0,100,10))+
  scale_color_manual(values = c("negative"="#d7191c","neutral"="#fdae61","positive"="#2c7bb6"))+
  ylab("Percent Coral Cover")+
  labs(title = "(f) Calcification Vulnerability Index") +
  theme_classic()+
  theme(text = element_text(size=16),
        plot.title = element_text(size=16),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"))

Figure6=(BCC_DHW_plot|CDP_DHW_plot|CVI_DHW_plot)/(BCCcoralplot|CDPcoralplot|CVIcoralplot)
ggsave("Figure6.pdf",Figure6,width=12.5,height=8)

###########Electronic Supplementary Materials########

#Table S1: Parrotfish Bioerosion Scaling Data
write_csv(RB_parrotfish_a_b,"Table_S1-Parrotfish_Bioerosion_Scaling.csv")

#Table S2: Mixed Effects Model Summaries
tab_model(m1,m2,m3,m4,m5,m6,m7,
          pred.labels = c("Intercept", "G", "%CoralCover","maxDHW"),
          string.pred = "Coeffcient",
          string.ci = "Conf. Int (95%)",
          string.p = "p-Value",
          file="Table_S2-Model_Summaries.html")

#Table S3: Site Summary Data
write_csv(CRAMP,"Table_S3-Site_Summary_Data.csv")

#Figure S1: Time of Day of dnTA Samples
dnta_time_plot=ggplot()+
  geom_pointrange(data=OCC_sitelevel_dnta,mapping=aes(x=as.numeric(hour(local_time))+as.numeric(minute(local_time))/60+as.numeric(second(local_time))/3600,y=dnta,ymin=dnta_lwr,ymax=dnta_upr))+
  scale_x_continuous(limits=c(0,24),breaks=seq(0,24,1))+
  scale_y_continuous(limits=c(-150,302),breaks = seq(-150,300,50))+
  xlab("Hour of Day (local time)")+
  ylab(expression(dnTA~(µmol~kg^-1)))+
  geom_hline(yintercept=0)+
  theme_classic()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("FigureS1.pdf",dnta_time_plot,width=9,height=5)