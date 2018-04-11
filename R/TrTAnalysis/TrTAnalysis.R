### === PK-11195 TrT Analysis === ###
### // Pontus P. Sigray, KI, April 2018

#load packges
library(tidyverse)
library(granviller)

#Load data
load(file = './DerivedData/TrTData.Rdata')

trtdata <- TidyAll %>%
  select(subject, PETno, Region, Vt.2tcm, Vs.2tcm, BPND.2tcm, 
         bp.srtm.cer, bp.srtm.svca_turkheimer, bp.srtm.svca_vumc, bp.srtm.svca_tpc,
         bp.srtmv.svca_turkheimer, bp.srtmv.svca_vumc,bp.srtmv.svca_tpc,SUV_4060) %>%
  gather(Measure, Value, -subject, -PETno, -Region) %>%
  spread(Region, Value) 

trtout <- trtdata %>%
  gather(Region, Binding, -(subject:Measure)) %>%
  spread(PETno, Binding, drop = T) %>%
  group_by(Region, Measure) %>%
  do(trt = granviller::trt(.$`1`, .$`2`)$tidy) %>%
  ungroup() %>%
  unnest() %>%
  select(-cov, -se, -skew, -kurtosis, -md, -avgpercchange) %>%
  mutate('MD' = ((sem*1.96*sqrt(2))/abs(mean))*100) %>%
  arrange(Measure, Region)

trtoutSummary <- trtout %>%
  group_by(Measure) %>%
  select(mean, icc, aapd, MD ) %>%
  summarise( meanOutcome = mean(mean),
             sdOutcome = sd(mean),
            meanICC = mean(icc),
            sdICC = sd(icc),
            meanAAPD = mean (aapd),
            sdAAPD = sd(aapd),
            meanMD = mean(MD),
            sdMD = sd(MD) )

xlsx::write.xlsx(as.data.frame(trtout),file = './Results/Tables/TRTmetrics.xlsx',row.names = F)
xlsx::write.xlsx(as.data.frame(trtoutSummary),file = './Results/Tables/TRTmetricsSummary.xlsx',row.names = F)

