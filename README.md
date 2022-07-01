# RBD Intact Peak Analysis

An R script implemented to perform analyses accompanying the manuscript:
"Online Hydrophilic Interaction Chromatography (HILIC) Enhanced Top-Down Mass Spectrometry Characterization of SARS-Cov-2 Spike Receptor Binding Domain".
The code reads a list of deconvoluted peaks from intact protein MS analyses (matrix of mass, abundance, and elution time slice using Protein Metrics Intact Mass software version 4.2) and performs comparisons, filtering, statistical calculations and generates PDF figures and CSV reports.

### USAGE

1. Clone/download the repository.
2. Install the 2 required packages.
3. Update the dataBasePath variable according to the location of the data in your computer.
4. Variables dataFolder and inputReconstruction: only uncomment the 2 lines with the folder to be analyzed.
5. Execute all code in the script.

### Data

The 4 folders contain data for 4 intact recombinant RBDs: 
* Ray_WT
* Sino_N501Y
* Sino_WT
* Ray_N331Q

Each folder contains a list of deconvoluted peaks for each of the online separation techniques evaluated: 
* C2 reverse-phase liquid chromatography (C2)
* Capillary zone electrophoresis (CE)
* Acrylamide-based monolithic hydrophilic interaction chromatography (HILIC) 

Output:
* figures.pdf
* outputPeaks.csv

### Reference

If you use this code please cite:
Wilson et al. "Online Hydrophilic Interaction Chromatography (HILIC) Enhanced Top-Down Mass Spectrometry Characterization of SARS-Cov-2 Spike Receptor Binding Domain". Anal. Chem. 2022, 94, 15, 5909–5917. https://doi.org/10.1021/acs.analchem.2c00139.
