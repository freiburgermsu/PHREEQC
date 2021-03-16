# PHREEQC desalination simulations

Desalination is an increasingly critical technology for human sustainability. The geochemical phenomena within reverse osmosis (RO) desalination modules like brine formation and mineral scaling are researched mitigate fouling of the RO modules. The geochemical software PHREEQC is an established tool for evaluating gechemical phenonema, however, the software has been applied to depict scaling over the distance of an arbitrarily defined RO module. This work couples PHREEQC software with interactive Jupyter Notebooks that create PHREEQC input files and that process the PHREEQC output data into figures and data tables. The following sections stepwise detail the procedures from downloading PHREEQC software through processing simulation outputs.

## PHREEQC installation and execution

 The "PHREEQC execution instructions" PDF file provides explicit directions for installing and the executing PHREEQC on both Windows and Macintosh operating systems.


## Notebook usage

The two Jupyter notebooks are command-line interfaces. The "2021-03-09_APF_PHREEQC RO input file generation_01" notebook guides the user through creating an input file for PHREEQC  that will simulate RO desalination. The "2021-03-06_APF_excel output calculations_06" notebook guides the user through processing the PHREEQC SELECTED_OUTPUT data file into graphs and data tables.

The input file notebook presumes that the user has installed iPHREEQC software. The path directory of the iPHREEQC software must be known for the input file code to successfully function. The default path for a Windows 10 operating system is suggested in the code, however, the path directories for Macintosh and Linux systems must be identified by the user.     

The generated input file from the "2021-03-09_APF_PHREEQC RO input file generation_01" notebook must be imported and executed in iPHREEQC. The resultant SELECTED_OUTPUT output file from the iPHREEQC simulation must be placed within the same directory as the "2021-03-06_APF_excel output calculations_06" notebook for the output data to be processed through the code. The data processing notebook requests the title of the SELECTED_OUTPUT file and scans the file name for "brine" or "scaling" and "pitzer" or "phreeqc", which are used to classify the nature of the iPHREEQC simulation data. File names that lack the simulation details prompt user entries of the information, which directs the code to correctly process the data. The user will finally have the option to export the graph to the same directory as the code in a picture format.    
