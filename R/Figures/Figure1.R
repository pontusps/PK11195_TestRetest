### === PK-11195 Figure1 === ###
### // Pontus P. Sigray, KI, Decemeber 2017
library(tidyverse)
load(file = './DerivedData/TrTData.Rdata')

TidyAllCorMat <- TidyAll %>% 
  select( -c(PETid, subject, PETno, Region ) ) %>%
  select(Vt.2tcm,Vs.2tcm,BPND.2tcm,
         bp.srtm.cer, bp.srtm.svca_turkheimer,bp.srtmv.svca_turkheimer,SUV_4060)

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r = ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # r2-value calculation
  r2 <- round(cor(x, y)*cor(x, y),2)
  txt <- bquote("R" ^2 ~ "=" ~ .(r2) )
  if(r2<0.01) txt <- bquote("R" ^2 ~ "< 0.01")
  text(0.5, 0.4, txt)
}

pdf(file = './Results/Figures/Figure1.pdf',width = 7,height = 7)
pairs(TidyAllCorMat, upper.panel = panel.cor, col = 'lightblue', cex = 1, pch = 19,
      labels = c (expression( paste('2TCM ',V[T]) ) ,expression( paste('2TCM ',V[S])),
                  expression( paste('2TCM ',BP[ND]) ), bquote( atop('SRTM-CER', 'BP' [ND])),
                  bquote( atop('SRTM-SVCA4', 'BP' [ND])), bquote( atop('SRTMv-SVCA4', 'BP' [ND])),
                  'SUV 40-60min'),
      cex.labels = 0.9
      )

dev.off()

