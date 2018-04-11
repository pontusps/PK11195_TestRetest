# PK11195_TestRetest
Re-analysis of the Karolinska Institutet's (R)-[11C]PK11195 Test-Retest dataset to assess the reliability and convergent validity of outcome measures and models with and without arterial input function.

This repository contains data and code to reproduce the results from the article "Test-retest reliability and convergent validity of (R)-[11C]PK11195 outcome measures without arterial input function" (LINK TO bioRxiv). 

Folder "DerivedData" contains three different .Rdata files: 

 - TidyDataNested.Rdata: contains all ROI TACs, reference TACs, TAC-time, image-derived blood curve ('IDblood') and weights
 - TidyBloodDataNestedLong.Rdata: contains the arterial sampled whole blood curve, the metabolite-corrected plasma curve and the plasma parent fraction curve
 - TrTData.Rdata: Contains derived outcome measures (e.g. VT and BPND values) from all models and methods evaluated in the study. 

Folder "R" contains R-scripts to do all kinetic modelling (requires the R-package 'kinfitr'), the test-retest analyses and to reproduce the Table and Figure in the article. 

Folder "Results" contains Figure 1 in pdf format and Table 1 in excel format. 

Folder "Preprint" contains the article preprint also uploaded to BioRxiv. The current version might differ from that on BioRxiv. 

Issues are turned off since this repo is meant for storage and is not actively maintained. 
