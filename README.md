# Codes for dBiRS
Codes for BiRS

## Data and software availability
- Previous data published in UK Biobank were used for this work and this research has been conducted using UK Biobank Resource under project 79237. One can refer to 
https://biobank.ctsu.ox.ac.uk/crystal/index.cgi for accessing and enabling data download or detailed description of UKB data.
- All WES studies are conducted in the UK Biobank Research Analysis Platform (RAP), see https://ukbiobank.dnanexus.com/landing for details.
- PLINK2.0 was used to pre-treat UKB data, see https://www.cog-genomics.org/plink/2.0/ for detailed manuls of using PLINK2.0.
- cosi2 was used to generated simulation data, see https://github.com/broadinstitute/cosi2 for detaild manuls of using cosi2.

## Code for dBiRS
1. Package *BiRS* contains the main functions of original BiRS algorithm and Maximum Marginal Score Test. One can install the package *BiRS* through the file ***BiRS_1.0.tar.gz***.  


2. Directory *simulations* contains the main functions for the simulations. Specifically, 
- The R scripts ***DistributedBiRS*** and ***BinaryGenerator*** contain the code for applying sBiRS method to blocks.
- The R script ***dBiRSAfter*** contains the code for applying sBiRS to detect significant blocks in central machine.
- The R script ***UniBlockQSCAN*** contains the code to apply QSCAN method to blocks.
- The R script ***QSCAN_After*** contains the code to summarize the detection results of QSCAN in each block.
- The R script ***KSAfter*** contains the code to summarize the detection results of KnockoffScreen in each block.
- The R script ***Impute_Func*** contains the code for imputation.

- The R script ***Simulation_Size*** contains the code for conducting simulation for calculating size of dBiRS and QSCAN.
- The R script ***Simulation-dBiRS*** contains the code for conducting simulation for dBiRS under alternative hypothesis.
- The R script ***Simulation-QSCAN*** contains the code for conducting simulation for QSCAN under alternative hypothesis.
- The R script ***Simulation-KS*** contains the code for conducting simulation for KnockoffScreen under alternative hypothesis.

- The R script ***Summary_Simulation*** contains the code for summarizing the simulation results.
- The R script ***Plot_Simulation*** contains the code for plotting the selection probability.


3. Directory *Application_Code* contains the main functions for WES studies. Since all the WES studies are conducted in RAP, we only provide the codes for analyzing Fuild Intelligence here and one can change the parameters to analyze Propective Memory. Specifically, 
- Please do quality controls and split blocks for these phenotypes before analysis.
- The R script ***distributed_BiRS*** and ***dBiRS_WES*** are codes for running sBiRS in blocks in RAP.
- The R script ***Summary_WES*** contains the code for applying sBiRS to detect significant blocks in central machine for WES studies.
