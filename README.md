# PHREEQC desalination simulations

Desalination is an increasingly critical technology for human sustainability. The geochemical phenomena within reverse osmosis (RO) desalination modules like brine formation and mineral scaling are researched mitigate fouling of the RO modules. The geochemical software PHREEQC is an established tool for evaluating gechemical phenonema, however, the software has been applied to depict scaling over the distance of an arbitrarily defined RO module. This work couples PHREEQC software with interactive Jupyter Notebooks that create PHREEQC input files and that process the PHREEQC output data into figures and data tables. The following sections stepwise detail the procedures from downloading PHREEQC software through processing simulation outputs.

## PHREEQC installation and execution

 The "PHREEQC execution instructions" PDF file provides explicit directions for installing and the executing PHREEQC on both Windows and Macintosh operating systems.


## PHREEQC Notebooks 

The use of PHREEQC is facilitated by two Jupyter Notebooks. The Notebooks utilize command-line interfaces that circumvent the need of coding literacy. The "Input file generation" notebook guides the user through creating a PHREEQC input file that will simulate the geochemistry of RO desalination. The "Output processing" notebook guides the user through processing the PHREEQC SELECTED_OUTPUT data file into graphs and data tables.


### Input file generation

The "Input file generation" Notebook is organized into functions that each parameterize a different PHREEQC code block. The make_general_conditions function of the input file notebook defines initial details of the simulation and creates a title for the simulation. The other functions are more elaborate and are separately discussed in the following sections.


#### make_solutions


The make_solutions function generates the SOLUTION blocks of the PHREEQC simulation. The SOLUTION blocks of PHREEQC parameterize geochemical solution characteristics like elemental concentrations, alkalinity, and pH. The make_solutions function allows the user to define the feed geochemistry from one of the predefined water bodies, or to customize feed geochemistry through a series of command-line prompts. The predefined options for water bodies reflect experimental literature for the Mediterranean and Red Seas, which are archetypes of highly populated inland seas, and for groundwaters around the continential U.S.A. that represent the produced waters from fracking and are an intriguing water resource. The initial solution of the RO module is parameterized as pure H2O since the solution is rapidly flushed from the RO module during the simulation.


#### make_equilibrium_phases

The make_equilibrium_phases function generates the EQUILIBRIUM_PHASES block. The EQUILIBRIUM_PHASES block of PHREEQC defines the mineral scaling equilibrium during RO simulations. The function displays the alphabetized list of minerals and their chemical formulas that are defined by the user selected database, from which the user selects the minerals, and their preexisting quantities, that will be simulated. 


#### make_reactive_transport

The make_reactive_transport function generates both REACTION code blocks and the TRANSPORT code block. The function enables users to either simulate desalination through a  Dow-DuPont FILMTEC BW30-400 RO module, which is an archetype, or customize a module parameters through prompts in the notebook interface. The function contains a multitude of backend calculations that interpret permeate flux properties from the module parameters. The function currently only parameterizes single-domain model reactive transport.



#### make_selected_output

The make_selected_output function defines the data and the name of the output SELECTED_OUTPUT file. The SELECTED_OUTPUT file contains simulation data in a tab-delimited text file that is imported in the "Output processing" file and interpreted for simulation results. 


#### export

The export function exports the generated PHREEQC input file. The function simultaneously prints the complete input PHREEQC code in the notebook interface for the user's review and exports the input file with a predefined naming structure of the day’s date, the simulated water body selection, brine or scaling simulation perspective, and a whole number that increases incrementally until a unique name is achieved. The input file is saved as PQI file which denotes an iPHREEQC input file, although, the file extension is also aceptable for batch software versions of PHREEQC.












The input file notebook presumes that the user has installed iPHREEQC software. The path directory of the iPHREEQC software must be known for the input file code to successfully function. The default path for a Windows 10 operating system is suggested in the code, however, the path directories for Macintosh and Linux systems must be identified by the user.     

The generated input file from the "Input file generation" notebook must be imported and executed in iPHREEQC. The resultant SELECTED_OUTPUT output file from the iPHREEQC simulation must be placed within the same directory as the "Output processing" notebook for the output data to be processed through the code. The data processing notebook requests the title of the SELECTED_OUTPUT file and scans the file name for "brine" or "scaling" and "pitzer" or "phreeqc", which are used to classify the nature of the iPHREEQC simulation data. File names that lack the simulation details prompt user entries of the information, which directs the code to correctly process the data. The user will finally have the option to export the graph to the same directory as the code in a picture format.    
