# Actin_Length_Simulations
Here we provide the main codes used in the work **Emergence and maintenance of variable-length actin filaments in a limiting pool of building blocks** ([BioRxiv link](https://www.biorxiv.org/content/10.1101/2021.11.07.467615v2.abstract)). These codes are written in FORTRAN90 language and tested in Linux. Commands to generate the ensembles using the codes were passed using shell scripting. Two MATLAB codes are provided for basic analysis and visuliation. Commands for running the codes related to each figures in the work are named likewise, e.g. commands for running the codes to generate Figure 2 of the manuscript is named as "Fig_2.sh". *Please note that we do not provide a full automated way to generate the figures of the work but provide the respective main codes with some examples*.
## growthcode_1.f90
This code was used to simulate the effective actin growth descriptions where explicit profilin and formin mediated nucleation was ignored. This code was used for Figures 2, 3, 4 and 6 in the manuscript.

## growthcode_2.f90
This code was used to simulate the actin growth descriptions with explicit profilin and formin mediated nucleation. This code was used in Figure 5 and many SI figures where the full model was explored.

## Figure 2: Spontaneous nucleation promotes F-actin length heterogeneity
- Data for Figure.2 can be generated by running Fig_2.sh which calls growthcode_1.f90 to simulate pure actin growth without any other actin binding proteins.
- The MATLAB codes Average_timeseries.m and Length_distribution.m can be run to generate the time evolution of average quantities and length distribution.

## Figure 3: F-actin capping increases filament length heterogeneity
- Data for Figure.3 can be generated by running Fig_3.sh which calls growthcode_1.f90 to simulate actin growth with capping proteins.
- The MATLAB codes Average_timeseries.m and Length_distribution.m can be run to generate the time evolution of average quantities and length distribution.

## Figure 4: Role of Formin on F-actin length control
- Data for Figure.4 can be generated by running Fig_4.sh which calls growthcode_1.f90 to simulate actin growth with formin proteins.
- The MATLAB codes Average_timeseries.m and Length_distribution.m can be run to generate the time evolution of average quantities and length distribution.

## Figure 5: Role of Formin on F-actin length control
- Data for Figure.5 can be generated by running Fig_5.sh which calls growthcode_2.f90 to simulate actin growth with profilin and formin proteins.
- The MATLAB codes Average_timeseries.m and Length_distribution.m can be run to generate the time evolution of average quantities and length distribution.

## Figure 6: Competition between Formin and Capping proteins results in bimodal F-actin length distribution
- Data for Figure.6 can be generated by running Fig_6.sh which calls growthcode_1.f90 to simulate actin growth with capping and formin proteins.
- The MATLAB codes Average_timeseries.m and Length_distribution.m can be run to generate the time evolution of average quantities and length distribution.

## Output data structure
The main codes produce two output files "m.txt" and "pl.txt" which recored the time evolution for the whole time duration and length of filaments during the last one minute of the simulation. 
### Format of m.txt:
9 column data. Each row represents - time | number of filaments | mean filament length in subunits | number of capping bound filaments | number of profilin | number of profilin-actin | number of formin-bound filaments | reaction count | number of actin monomers
### Format of pl.txt:
Single column data. Each row represents length of a filament in subunits.

These two files, "m.txt" and "pl.txt" are appended for each run and the combined data for the whole ensemble is saved as "P\_fr\_*formin concentration*\_cp\_*capping concentration*.txt" and "L\_fr\_*formin concentration*\_cp\_*capping concentration*.txt" respectively, e.g., for [FR]=5 nM and [CP]=0 nM the files will be named "P_fr_5_cp_0.txt" and "L_fr_5_cp_0.txt".  

## Codes for analysis
The two MATLAB codes "Average_timeseries.m" and "Length_distribution.m" read "P\_fr\_*formin concentration*\_cp\_*capping concentration*.txt" and "L\_fr\_*formin concentration*\_cp\_*capping concentration*.txt" to generate timeseries of averaged quantities (such as mean length) and length distribution respectively.




