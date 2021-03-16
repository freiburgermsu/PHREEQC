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

The make_reactive_transport function of the input file notebook parameterizes both code blocks. The function parameterizes either the Dow-DuPont FILMTEC BW30-400 RO module, which was selected as an archetypical industrial standard, or a customized module with guided prompts through the notebook interface. The input parameters are then calculated into module characteristics that govern the transport properties of the code. 

The REACTION blocks critically define the permeate fluxes of H¬2O from each module cell. The TRANSPORT block defines the spatiotemporal simulation variables like module distance and simulation timesteps, which are either directly or indirectly parameterized by user inputs. Simulations between seconds and a month in Figures S21-S22 were conducted through permutations of the TRANSPORT block. The removal of the TRANSPORT block while retaining the REACTION blocks simulates evaporation, where mineral precipitation occurs in the absence of feed flow or transport phenomena. Evaporation of Figure S7 was a crucial internal validation of transport scaling in Figure S8 at an equivalent CF.

The reactive transport of desalination can be modeled in either the single domain or the dual domain. A single-domain model considers the CP and bulk solution as a single blended solution while the dual-domain model differentiates the CP and the bulk solutions through parameterizing one set of cells [1,n],n∈W for the bulk solution and a second set of cells [n+2,m],m∈W>n+2 for the CP. The dual domain model contains additional parameters of the volume fractions between the bulk and CP solutions and of the rate of solvent exchange, termed the exchange  factor ε_F, between the bulk and CP solutions. 





	SELECTED_OUTPUT block

The SELECTED_OUTPUT block prints simulation data in a tab-delimited text file. The block importantly differentiates simulations that examine scaling throughout the module distance versus simulations that examine effluent brine. The make_selected_output function of the input file notebook importantly prompts whether brine or scaling will be simulated and to name the output file. 


	Export

The export function of the input file notebook exports the PHREEQC input file. The function both prints the complete input PHREEQC code in the notebook interface and exports the input code with a predefined naming structure of the day’s date, the simulated water body selection, brine or scaling simulation perspective, and a whole number that increases incrementally for each exported file of the otherwise same name. The input file is saved as PQI file which denotes an iPHREEQC input file.












The input file notebook presumes that the user has installed iPHREEQC software. The path directory of the iPHREEQC software must be known for the input file code to successfully function. The default path for a Windows 10 operating system is suggested in the code, however, the path directories for Macintosh and Linux systems must be identified by the user.     

The generated input file from the "Input file generation" notebook must be imported and executed in iPHREEQC. The resultant SELECTED_OUTPUT output file from the iPHREEQC simulation must be placed within the same directory as the "Output processing" notebook for the output data to be processed through the code. The data processing notebook requests the title of the SELECTED_OUTPUT file and scans the file name for "brine" or "scaling" and "pitzer" or "phreeqc", which are used to classify the nature of the iPHREEQC simulation data. File names that lack the simulation details prompt user entries of the information, which directs the code to correctly process the data. The user will finally have the option to export the graph to the same directory as the code in a picture format.    
