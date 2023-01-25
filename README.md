# HOBMI
This is the Matlab code for the HOBMI algorithm proposed in the paper "Identification of Power System Oscillation Modes using Blind Source Separation based on Copula Statistic."
Instructions for running the code:
1. The code includes the dataset for case a A and case B of the 11-Bus 4-machine system in the Kundoor textbook, and toy example considered in the paper.
2. Please have PSAT installed on your system and keep it in your Matlab folder. 
3. Please avoid using the given data and step 1 if you have access to the actual measurements of relative angular frequency of your system.

Files:
HOBMI_2area_4machine_system: HOBMI i.e. HOBI with mode identification applied to 11-Bus 4-machine system and toy example where the data are generated from Archimedian copula family, namely, Frank, Clayton, Gumbel, and Gaussian copula.
HOBI: High-order blind identification algorithm based on copula statistics.
cosdv: calculates the copula statistics between two random variables. 
d_kundur1_mdl_05: data file for case A.
d_kundur1_mdl_1_3_20Hz: data file for case B.
