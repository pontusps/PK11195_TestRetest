### === PK-11195 kinfit analysis === ###
### // Pontus P. Sigray, KI, April 2018

#install packages
#devtools::install_github("mathesong/kinfitr") #install kinfitr from github 

#load packages
library(tidyverse)
library(kinfitr) 

#load data
load(file = './DerivedData/TidyDataNested.Rdata')
load(file = './DerivedData/TidyBloodDataNestedLong.Rdata')

######### 
# 2TCM
#########

###Merge TAC, weights and blood long-formatted nested dataframe 
TidyDataNested <- left_join(TidyDataNested,TidyBloodDataNestedLong)

#Make 2TCM fun for fitting delay and vB using the graymatter ROI 
fitDelay_vB <- function(tacs, input) {
  twotcm(t_tac = tacs$time, tac = tacs$GM, input = input, weights = tacs$Weights)
}

#Fit delay and vB using 2TCM and GM ROI
TidyDataNested <- TidyDataNested %>%   
  mutate(delayFit = pmap(list(TACs,input), fitDelay_vB)) 

#Make TidyDataNested into even longer nested format by spreading ROIs
TidyDataNestedLong <- TidyDataNested %>%
  select(PETid, subject, PETno, TACs) %>%
  unnest() %>%
  gather(Region, TAC, -c(PETid, subject, PETno, CER,  IDblood, time, Weights)) %>%
  group_by(PETid, subject, PETno, Region) %>%
  nest(.key = 'TACs')

#Merge back blood data
TidyDataNestedLong <- TidyDataNested %>%
  select( PETid , subject, PETno, input) %>%
  left_join(TidyDataNestedLong, .)

#Merge back delayFit objects
TidyDataNestedLong <- TidyDataNested %>%
  select( PETid , subject, PETno, delayFit) %>%
  left_join(TidyDataNestedLong, .)

#Make function for fitting 2TCM
fit2tcm <- function(tacs, input, delayFit) {
  twotcm(t_tac = tacs$time, tac = tacs$TAC, input = input, weights = tacs$Weights,
         inpshift = delayFit$par$inpshift, vB_fixed = delayFit$par$vB) 
}

#Make function for fitting Logan VT
fitLogan <- function(tacs, input, delayFit) {
  Loganplot(t_tac = tacs$time, tstarIncludedFrames = 7,
            tac = tacs$TAC, input = input, weights = tacs$Weights,
            inpshift = delayFit$par$inpshift, vB = delayFit$par$vB) 
}

#Fit 2TCM and calculate VT, VS and BPND from rate constants
Tidy2TCM <- TidyDataNestedLong %>%  
  mutate(fit_2tcm = purrr::pmap(list(TACs, input, delayFit), fit2tcm)) %>%
  mutate(K1_2tcm = map_dbl(fit_2tcm, c('par', 'K1'))) %>%
  mutate(k2_2tcm = map_dbl(fit_2tcm, c('par', 'k2'))) %>%
  mutate(k3_2tcm = map_dbl(fit_2tcm, c('par', 'k3'))) %>%
  mutate(k4_2tcm = map_dbl(fit_2tcm, c('par', 'k4'))) %>%
  mutate(Vt.2tcm   = (K1_2tcm / k2_2tcm) * (1 + k3_2tcm / k4_2tcm)) %>%
  mutate(Vs.2tcm   = (K1_2tcm * k3_2tcm) / (k2_2tcm * k4_2tcm)) %>%
  mutate(BPND.2tcm = (k3_2tcm / k4_2tcm)) 

#Fit Logan and retrieve VT 
TidyLogan <- TidyDataNestedLong %>%  
  mutate(fit_logan = pmap(list(TACs, input, delayFit), fitLogan)) %>%
  mutate(Vt.logan = map_dbl(fit_logan, c('par', 'Vt'))) 

#Plot 2TCM fits
#plot_kinfit(Tidy2TCM$fit_2tcm[[26]]) + theme_bw()

#Place SVCA4refs and cerebellum back into parent df
TidyDataNestedLong <- TidyDataNested %>%
  select(PETid, subject, PETno, TACs) %>%
  unnest() %>%
  gather(Region, TAC, -c(PETid, subject, PETno, CER, SVCArefVUMC, SVCArefTurkheimer, SVCArefTPC, IDblood, time, Weights)) %>%
  group_by(PETid, subject, PETno, Region) %>%
  nest(.key = 'TACs')

#Merge back blood data
TidyDataNestedLong <- TidyDataNested %>%
  select( PETid , subject, PETno, input) %>%
  left_join(TidyDataNestedLong, .)


########################
# SRTM - Cerbellum
########################
set.seed(41)

#Make function for fitting SRTM with cerebellum as reference region
dofit_srtm_cer <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$CER, 
       roitac = tacs$TAC, weights = tacs$Weights,
       bp.lower = -1,
       bp.upper = 10,
       R1.lower = -1,
       k2.lower = -1,
       bp.start = 0,
       k2.start = 0.1,
       multstart_iter = 100)
}

#Fit srtm using cerebellum as reference region
TidySRTMCer <- TidyDataNestedLong %>%  
  mutate(fit_srtm_cer = pmap(list(TACs), dofit_srtm_cer)) %>%
  mutate(bp.srtm.cer = map_dbl(fit_srtm_cer, c('par', 'bp'))) 

#make function to refit SRTM with cerebellum as reference region for TACs that hit lower and upper bounds
dofit_srtm_cer <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$CER, 
       roitac = tacs$TAC, weights = tacs$Weights,
       bp.start = -0.1,
       k2.start = 0.9, 
       R1.start = 2.2)
}

#Refit srtm using cerebellum as reference region
TidySRTMCer_refit <- TidyDataNestedLong[c(5,42),] %>%  
  mutate(fit_srtm_cer = pmap(list(TACs), dofit_srtm_cer)) %>%
  mutate(bp.srtm.cer = map_dbl(fit_srtm_cer, c('par', 'bp'))) 

TidySRTMCer$fit_srtm_cer[c(5,42)] <- TidySRTMCer_refit$fit_srtm_cer
TidySRTMCer$bp.srtm.cer[c(5,42)] <- TidySRTMCer_refit$bp.srtm.cer

###################### 
# Logan - Cerebellum
######################

#make mrtm1 function to get k2prime
dofit_mrtm1_cer <- function(tacs) {
  mrtm1(t_tac = tacs$time, reftac = tacs$CER,roitac = tacs$TAC,weights = tacs$Weights) 
}

#Run mrtm1 to get k2prime
TidyLoganCer <- TidyDataNestedLong %>%  
  mutate(fit_mrtm1_cer = pmap(list(TACs), dofit_mrtm1_cer)) %>%
  mutate(bp.mrtm1.cer = map_dbl(fit_mrtm1_cer, c('par', 'bp')),
         k2prime.mrtm1.cer = map_dbl(fit_mrtm1_cer, c('par', 'k2prime')))

#make Logan fun, t* = 20min
dofit_logan_cer <- function(tacs, k2prime) {
  refLogan(t_tac = tacs$time,reftac = tacs$CER,
           roitac = tacs$TAC,k2prime = k2prime,
           tstarIncludedFrames = 7,weights =  tacs$Weights) 
}

#run logan 
TidyLoganCer <- TidyLoganCer %>%  
  mutate(fit_logan_cer = pmap(list(TACs,k2prime.mrtm1.cer), dofit_logan_cer)) %>%
  mutate(bp.logan.cer = map_dbl(fit_logan_cer, c('par', 'bp'))) 

plot(TidyLoganCer$bp.logan.cer,TidySRTMCer$bp.srtm.cer)
abline(0,1)

#################################################################################
# SRTMv - SVCA reference using Turkheimer 2006 population-based kinetic classes
#################################################################################
set.seed(42)

#Make function for fitting SRTMv with a SVCA derived reference TAC
dofit_srtmv_svca <- function(tacs) {
  srtm_v(t_tac = tacs$time, reftac = tacs$SVCArefTurkheimer,
         roitac = tacs$TAC, weights = tacs$Weights, 
         bloodtac = tacs$IDblood,
         bp.lower = -0.25,
         R1.lower = -0.1,
         k2.lower = -0.1,
         multstart_iter = 100
  )
}

#fit srtmv using SVCA reference from Turkheimer 2006 population-based kinetic classes
#NB: fitting takes a few minutes
TidySRTMvSCVATurkheimer <- TidyDataNestedLong %>%  
  mutate(fit_srtmv_svca_turkheimer = pmap(list(TACs), dofit_srtmv_svca)) %>%
  mutate(bp.srtmv.svca_turkheimer = map_dbl(fit_srtmv_svca_turkheimer, c('par', 'bp')))

#make function to refit SRTM with SVCA-Turkheimer for TACs that hit lower bounds
dofit_srtmv_svca <- function(tacs) {
  srtm_v(t_tac = tacs$time, reftac = tacs$SVCArefTurkheimer, 
         bloodtac = tacs$IDblood,
         roitac = tacs$TAC, 
         weights = tacs$Weights,
         bp.start = 0.1,
         k2.start = 0.4,
         k2.lower = -1,
         R1.start = 0.1,
         R1.lower = -1
  )
}

#Refit srtmv for subject that hit lower bound using SVCA reference from Turkheimer 2006 population-based kinetic classes
TidySRTMvSCVATurkheimer_refit <- TidySRTMvSCVATurkheimer[c(47),] %>%
  mutate(fit_srtmv_svca_turkheimer = pmap(list(TACs), dofit_srtmv_svca)) %>%
  mutate(bp.srtmv.svca_turkheimer = map_dbl(fit_srtmv_svca_turkheimer, c('par', 'bp')))  

TidySRTMvSCVATurkheimer$fit_srtmv_svca_turkheimer[c(47)] <- TidySRTMvSCVATurkheimer_refit$fit_srtmv_svca_turkheimer
TidySRTMvSCVATurkheimer$bp.srtmv.svca_turkheimer[c(47)] <- TidySRTMvSCVATurkheimer_refit$bp.srtmv.svca_turkheimer

###########################################################################################
# SRTMv - SVCA reference using VUMC (Yaqub et al, 2011) population-based kinetic classes
###########################################################################################
set.seed(43)

#Make function for fitting SRTMv with a SVCA derived reference TAC
dofit_srtmv_svca <- function(tacs) {
  srtm_v(t_tac = tacs$time, reftac = tacs$SVCArefVUMC,
         roitac = tacs$TAC, weights = tacs$Weights, 
         bloodtac = tacs$IDblood,
         bp.lower = -0.25,
         R1.lower = -0.1,
         k2.lower = -0.1,
         multstart_iter = 100
  )
}

#fit srtmv using SVCA reference from VUMC (Yaqub et al., 2011) population-based kinetic classes
#NB: fitting takes a few minutes
TidySRTMvSCVAVUMC <- TidyDataNestedLong %>%  
  mutate(fit_srtmv_svca_vumc = pmap(list(TACs), dofit_srtmv_svca)) %>%
  mutate(bp.srtmv.svca_vumc = map_dbl(fit_srtmv_svca_vumc, c('par', 'bp')))

#Make function for refitting SRTMv with a SVCA derived reference for TACs that hit lower bounds
dofit_srtmv_svca <- function(tacs) {
  srtm_v(t_tac = tacs$time, reftac = tacs$SVCArefVUMC,
         roitac = tacs$TAC, weights = tacs$Weights, 
         bloodtac = tacs$IDblood,
         bp.start = 0.1,
         k2.start = 0.4,
         k2.lower = -1,
         R1.start = 0.1,
         R1.lower = -1
  )
}

#Refit srtmv using SVCA reference from VUMC (Yaqub et al., 2011) population-based kinetic classes
#NB: fitting takes a few minutes
TidySRTMvSCVAVUMC_refit <- TidyDataNestedLong[c(38,42),] %>%  
  mutate(fit_srtmv_svca_vumc = pmap(list(TACs), dofit_srtmv_svca)) %>%
  mutate(bp.srtmv.svca_vumc = map_dbl(fit_srtmv_svca_vumc, c('par', 'bp')))

TidySRTMvSCVAVUMC$fit_srtmv_svca_vumc[c(38,42)] <- TidySRTMvSCVAVUMC_refit$fit_srtmv_svca_vumc
TidySRTMvSCVAVUMC$bp.srtmv.svca_vumc[c(38,42)] <- TidySRTMvSCVAVUMC_refit$bp.srtmv.svca_vumc

###########################################################################################
# SRTMv - SVCA reference using TPC population-based kinetic classes
###########################################################################################
set.seed(43)

#Make function for fitting SRTMv with a SVCA derived reference TAC
dofit_srtmv_svca <- function(tacs) {
  srtm_v(t_tac = tacs$time, reftac = tacs$SVCArefTPC,
         roitac = tacs$TAC, weights = tacs$Weights, 
         bloodtac = tacs$IDblood,
         bp.lower = -0.25,
         R1.lower = -0.1,
         k2.lower = -0.1,
         multstart_iter = 100
  )
}

#fit srtmv using SVCA reference from TPC population-based kinetic classes
#NB: fitting takes a few minutes
TidySRTMvSCVATPC <- TidyDataNestedLong %>%  
  mutate(fit_srtmv_svca_tpc = pmap(list(TACs), dofit_srtmv_svca)) %>%
  mutate(bp.srtmv.svca_tpc = map_dbl(fit_srtmv_svca_tpc, c('par', 'bp')))

#Make function for refitting SRTMv with a SVCA derived reference for TACs that hit upper bound
dofit_srtmv_svca <- function(tacs) {
  srtm_v(t_tac = tacs$time, reftac = tacs$SVCArefTPC,
         roitac = tacs$TAC, weights = tacs$Weights, 
         bloodtac = tacs$IDblood,
         bp.start = 0.1,
         k2.start = 0.4,
         k2.lower = -1,
         R1.start = 0.1,
         R1.lower = -1
  )
}

#Refit srtmv using SVCA reference from TPC population-based kinetic classes
#NB: fitting takes a few minutes
TidySRTMvSCVATPC_refit <- TidyDataNestedLong[c(27),] %>%  
  mutate(fit_srtmv_svca_tpc = pmap(list(TACs), dofit_srtmv_svca)) %>%
  mutate(bp.srtmv.svca_tpc = map_dbl(fit_srtmv_svca_tpc, c('par', 'bp')))

TidySRTMvSCVATPC$fit_srtmv_svca_tpc[c(27)] <- TidySRTMvSCVATPC_refit$fit_srtmv_svca_tpc
TidySRTMvSCVATPC$bp.srtmv.svca_tpc[c(27)] <- TidySRTMvSCVATPC_refit$bp.srtmv.svca_tpc

###############################################################################
# SRTM - SCVA reference using Turkheimer 2006 population-based kinetic classes
###############################################################################
set.seed(44)

#Make function for fitting srtm with a SVCA derived reference TAC
dofit_srtm_svca <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$SVCArefTurkheimer,
       roitac = tacs$TAC, weights = tacs$Weights, 
       bp.lower = -0.25,
       R1.lower = -0.1,
       k2.lower = -0.1,
       multstart_iter = 100
  )
}

#fit srtm using SVCA reference from Turkheimer 2006 population-based kinetic classes
#NB: fitting takes a few minutes
TidySRTMSCVATurkheimer <- TidyDataNestedLong %>%  
  mutate(fit_srtm_svca_turkheimer = pmap(list(TACs), dofit_srtm_svca)) %>%
  mutate(bp.srtm.svca_turkheimer = map_dbl(fit_srtm_svca_turkheimer, c('par', 'bp')))

#Make fun to refit value that hitter upper bound
dofit_srtm_svca <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$SVCArefTurkheimer,
       roitac = tacs$TAC, weights = tacs$Weights, 
       R1.upper = 0.99
  )
}

#Refit srtmv for subject that hit lower bound using SVCA reference 
TidySRTMSCVATurkheimer_refit <- TidySRTMSCVATurkheimer[c(39,45,48),] %>%
  mutate(fit_srtm_svca_turkheimer = pmap(list(TACs), dofit_srtm_svca)) %>%
  mutate(bp.srtm.svca_turkheimer = map_dbl(fit_srtm_svca_turkheimer, c('par', 'bp'))) 

TidySRTMSCVATurkheimer$fit_srtm_svca_turkheimer[c(39,45,48)] <- TidySRTMSCVATurkheimer_refit$fit_srtm_svca_turkheimer
TidySRTMSCVATurkheimer$bp.srtm.svca_turkheimer[c(39,45,48)] <- TidySRTMSCVATurkheimer_refit$bp.srtm.svca_turkheimer

#########################################################################################
# SRTM - SCVA reference using VUMC (Yaqub et al, 2011) population-based kinetic classes
#########################################################################################
set.seed(45)

#Make function for fitting srtm with a SVCA derived reference TAC. 
dofit_srtm_svca <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$SVCArefVUMC,
       roitac = tacs$TAC, weights = tacs$Weights, 
       bp.lower = -0.25,
       R1.lower = -0.1,
       k2.lower = -0.1,
       multstart_iter = 100
  )
}

#fit srtm using SVCA reference from VUMC (Yaqub, 2011)
#NB: fitting takes a few minutes
TidySRTMSCVAVUMC <- TidyDataNestedLong %>%  
  mutate(fit_srtm_svca_vumc = pmap(list(TACs), dofit_srtm_svca)) %>%
  mutate(bp.srtm.svca_vumc = map_dbl(fit_srtm_svca_vumc, c('par', 'bp')))

#Make fun to refit value that hit lower bound
dofit_srtm_svca <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$SVCArefVUMC,
       roitac = tacs$TAC, weights = tacs$Weights, 
       R1.upper = 1
  )
}

#Refit srtmv for subject that hit lower bound 
TidySRTMSCVAVUMC_refit <- TidySRTMSCVAVUMC[c(38),] %>%
  mutate(fit_srtm_svca_vumc = pmap(list(TACs), dofit_srtm_svca)) %>%
  mutate(bp.srtm.svca_vumc = map_dbl(fit_srtm_svca_vumc, c('par', 'bp'))) 

TidySRTMSCVAVUMC$fit_srtm_svca_vumc[c(38)] <- TidySRTMSCVAVUMC_refit$fit_srtm_svca_vumc
TidySRTMSCVAVUMC$bp.srtm.svca_vumc[c(38)] <- TidySRTMSCVAVUMC_refit$bp.srtm.svca_vumc

#########################################################################################
# SRTM - SCVA reference using TPC population-based kinetic classes
#########################################################################################
set.seed(45)

#Make function for fitting srtm with a SVCA derived reference TAC. 
dofit_srtm_svca <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$SVCArefTPC,
       roitac = tacs$TAC, weights = tacs$Weights, 
       bp.lower = -0.25,
       R1.lower = -0.1,
       k2.lower = -0.1,
       multstart_iter = 100
  )
}

#fit srtm using SVCA reference from TPC 
#NB: fitting takes a few minutes
TidySRTMSCVATPC <- TidyDataNestedLong %>%  
  mutate(fit_srtm_svca_tpc = pmap(list(TACs), dofit_srtm_svca)) %>%
  mutate(bp.srtm.svca_tpc = map_dbl(fit_srtm_svca_tpc, c('par', 'bp')))

#Make fun to refit value that hit upper bound
dofit_srtm_svca <- function(tacs) {
  srtm(t_tac = tacs$time, reftac = tacs$SVCArefTPC,
       roitac = tacs$TAC, weights = tacs$Weights, 
       bp.start = 0.1,
       k2.start = 0.4,
       k2.lower = -1,
       R1.start = 0.1,
       R1.lower = -1
  )
}

#Refit srtmv for subject that hit lower bound using SVCA reference 
TidySRTMSCVATPC_refit <- TidySRTMSCVATPC[c(31),] %>%
  mutate(fit_srtm_svca_tpc = pmap(list(TACs), dofit_srtm_svca)) %>%
  mutate(bp.srtm.svca_tpc = map_dbl(fit_srtm_svca_tpc, c('par', 'bp'))) 

TidySRTMSCVATPC$fit_srtm_svca_tpc[c(31)] <- TidySRTMSCVATPC_refit$fit_srtm_svca_tpc
TidySRTMSCVATPC$bp.srtm.svca_tpc[c(31)] <- TidySRTMSCVATPC_refit$bp.srtm.svca_tpc

######################
# SUVs - 40 to 60 min 
######################
library(pracma)
DemogRadioChem <- read.csv('./RawData/Demog/DemogRadioChem.csv')
TidyDataNestedLong <- left_join(TidyDataNestedLong,DemogRadioChem)

# Total SUV
calcSUV <- function(tacs, InjRadMbq, BodyWeightKg) {
  SUV(t_tac = tacs$time,
      tac = tacs$TAC*0.037,    # To kBq - because of kg
      injRad = InjRadMbq,
      bodymass = BodyWeightKg)
}

# 40-60 Minute SUV
calcSUV_4060 <- function(SUVout) {
  
  interptime = seq(SUVout$tacs$Time[1], rev(SUVout$tacs$Time)[1], by=1/60)
  interptac  = interp1(x = SUVout$tacs$Time, SUVout$tacs$SUV, 
                       xi = interptime, method = 'linear')
  step=interptime[2] - interptime[1]
  
  SUVdf <- data.frame(interptime, interptac) %>%
    filter(interptime > 40 & interptime <= 60)
  
  out <- mean(SUVdf$interptac)
  return(out)
}

#Calculate SUVs
TidySUV <- TidyDataNestedLong %>%
  mutate(suvout_tot = pmap(list(TACs, BodyWeightKg, InjRadMbq), calcSUV)) %>%
  mutate(SUV_tot = map_dbl(suvout_tot, c('par', 'intSUV'))) %>%
  mutate(SUV_4060 = map_dbl(suvout_tot, calcSUV_4060)) 


#################################
# Merge all kinetic model output
#################################

Tidy2TCM <- Tidy2TCM %>%
  select(PETid, subject, PETno, Region, Vt.2tcm,Vs.2tcm, BPND.2tcm, K1_2tcm)

TidyLogan<- TidyLogan %>%
  select(PETid, subject, PETno, Region, Vt.logan)

TidySRTMCer <- TidySRTMCer %>%
  select(PETid, subject, PETno, Region, bp.srtm.cer)

TidyLoganCer <- TidyLoganCer %>%
  select(PETid, subject, PETno, Region, bp.logan.cer)

TidySRTMSCVATurkheimer <- TidySRTMSCVATurkheimer %>%
  select(PETid, subject, PETno, Region,bp.srtm.svca_turkheimer)

TidySRTMSCVAVUMC <- TidySRTMSCVAVUMC %>%
  select(PETid, subject, PETno, Region,bp.srtm.svca_vumc)

TidySRTMSCVATPC <- TidySRTMSCVATPC %>%
  select(PETid, subject, PETno, Region,bp.srtm.svca_tpc)

TidySRTMvSCVATurkheimer <- TidySRTMvSCVATurkheimer %>%
  select(PETid, subject, PETno, Region,bp.srtmv.svca_turkheimer)

TidySRTMvSCVAVUMC <- TidySRTMvSCVAVUMC %>%
  select(PETid, subject, PETno, Region,bp.srtmv.svca_vumc)

TidySRTMvSCVATPC <- TidySRTMvSCVATPC %>%
  select(PETid, subject, PETno, Region,bp.srtmv.svca_tpc)

TidySUV <- TidySUV %>%
  select(PETid, subject, PETno, Region,SUV_4060)

TidyAll <- left_join(Tidy2TCM,TidySRTMCer)
TidyAll <- left_join(TidyAll,TidyLogan)
TidyAll <- left_join(TidyAll,TidyLoganCer)
TidyAll <- left_join(TidyAll,TidySRTMSCVATurkheimer)
TidyAll <- left_join(TidyAll,TidySRTMSCVAVUMC)
TidyAll <- left_join(TidyAll,TidySRTMSCVATPC)
TidyAll <- left_join(TidyAll,TidySRTMvSCVATurkheimer)
TidyAll <- left_join(TidyAll,TidySRTMvSCVAVUMC)
TidyAll <- left_join(TidyAll,TidySRTMvSCVATPC)
TidyAll <- left_join(TidyAll,TidySUV)                     

save(TidyAll,file = './DerivedData/TrTData.Rdata')

