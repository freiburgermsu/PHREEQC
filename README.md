# PHREEQC desalination simulations

Desalination is an increasingly critical technology for human sustainability. The geochemical phenomena within reverse osmosis (RO) desalination modules like brine formation and mineral scaling are researched mitigate fouling of the RO modules. The geochemical software PHREEQC is an established tool for evaluating gechemical phenonema, however, the software has been applied to depict scaling over the distance of an arbitrarily defined RO module. This work couples PHREEQC software with interactive Jupyter Notebooks that create PHREEQC input files and that process the PHREEQC output data into figures and data tables. The following sections stepwise detail the procedures from downloading PHREEQC software through processing simulation outputs.

## PHREEQC installation and execution

 The "PHREEQC execution instructions" PDF file provides explicit directions for installing and the executing PHREEQC on both Windows and Macintosh operating systems.


## PHREEQC Notebooks 

The use of PHREEQC is facilitated by two Jupyter Notebooks. The Notebooks utilize command-line interfaces that circumvent the need of coding literacy. The "Input file generation" notebook guides the user through creating a PHREEQC input file that will simulate the geochemistry of RO desalination. The "Output processing" notebook guides the user through processing the PHREEQC SELECTED_OUTPUT data file into graphs and data tables.


### Input file generation


The make_general_conditions function of the input file notebook defines initial details. The function defines the database that the simulation will use and establishes the directory path on the user’s computer. The backend Python code ensures that the entered path actually exists in the user’s system and that the database is defined by the code. The function finally creates a title for the simulation.


	SOLUTION blocks


The SOLUTION blocks of PHREEQC parameterize geochemical solution characteristics. These block parameterize geochemical qualities like elemental concentrations, alkalinity, and pH that are accepted by the selected database. Charge imbalances can occur in the solution from either incompatible geochemical data sets32 or from neglected ionic species in the geochemical solution data. The charge imbalances can be automatically balanced by adjusting the pH, although, our sensitivity analysis revealed that maintaining a charge imbalance was no different than automatically correcting a charge imbalance in Tables S3-S4. Manual correction of the charge imbalance, however, through increasing the concentration of (-)  ions for (+) %-error or through increasing the concentration of (+)  ions for (-)%-error, augmented the charge imbalance and disrupted scaling results in Figures S37-S41.

The make_solutions function of the input file notebook generates SOLUTION blocks. The function populates the feed geochemistry from one of the predefined water bodies, or the user can elect to customize feed geochemistry through a series of command-line prompts. The predefined options reflect experimental geochemistry literature for the Mediterranean and Red Seas, which are archetypes of highly populated inland seas, and for groundwaters around the continential U.S.A. that represent produced fracking waters34,35,36,37, which are unutilized38 and hazardous39 waters of oil wells40–46 that are intriguing as a water resource. The initial solution within the RO module is parameterized as pure H2O, since the initial solution is arbitrary and is rapidly flushed from the RO module during the simulation.


	EQUILIBRIUM_PHASES blocks

The EQUILIBRIUM_PHASES block defines scaling equilibrium for RO simulations. This block parameterizes the minerals that will be examined for scaling in the simulation. The code block substantially influences the reactive transport geochemistry and changes the brine and scaling formation with different mineral parameters in Figures S42-S53.  The removal of the block entirely permits elemental oversaturation and can thereby emulate the addition of anti-scalants to a feed solution, which is demonstrated in Figures S9-S11.

The make_equilibrium_phases function of the input file notebook creates the code block. The notebook interface displays the alphabetized list of minerals and their chemical formulas that are defined by the selected database, from which list the user selects the minerals, and their preexisting quantities, that will be simulated. 


	REACTION and TRANSPORT blocks

The REACTION blocks critically define the permeate fluxes of H¬2O from each module cell. The TRANSPORT block defines the spatiotemporal simulation variables like module distance and simulation timesteps, which are either directly or indirectly parameterized by user inputs. Simulations between seconds and a month in Figures S21-S22 were conducted through permutations of the TRANSPORT block. The removal of the TRANSPORT block while retaining the REACTION blocks simulates evaporation, where mineral precipitation occurs in the absence of feed flow or transport phenomena. Evaporation of Figure S7 was a crucial internal validation of transport scaling in Figure S8 at an equivalent CF.

The reactive transport of desalination can be modeled in either the single domain or the dual domain. A single-domain model considers the CP and bulk solution as a single blended solution while the dual-domain model differentiates the CP and the bulk solutions through parameterizing one set of cells [1,n],n∈W for the bulk solution and a second set of cells [n+2,m],m∈W>n+2 for the CP. The dual domain model contains additional parameters of the volume fractions between the bulk and CP solutions and of the rate of solvent exchange, termed the exchange  factor ε_F, between the bulk and CP solutions. 

The make_reactive_transport function of the input file notebook parameterizes both code blocks. The function parameterizes either the Dow-DuPont FILMTEC BW30-400 RO module, which was selected as an archetypical industrial standard, or a customized module with guided prompts through the notebook interface. The input parameters are then calculated into module characteristics that govern the transport properties of the code. 



	SELECTED_OUTPUT block

The SELECTED_OUTPUT block prints simulation data in a tab-delimited text file. The block importantly differentiates simulations that examine scaling throughout the module distance versus simulations that examine effluent brine. The make_selected_output function of the input file notebook importantly prompts whether brine or scaling will be simulated and to name the output file. 


	Export

The export function of the input file notebook exports the PHREEQC input file. The function both prints the complete input PHREEQC code in the notebook interface and exports the input code with a predefined naming structure of the day’s date, the simulated water body selection, brine or scaling simulation perspective, and a whole number that increases incrementally for each exported file of the otherwise same name. The input file is saved as PQI file which denotes an iPHREEQC input file.












The input file notebook presumes that the user has installed iPHREEQC software. The path directory of the iPHREEQC software must be known for the input file code to successfully function. The default path for a Windows 10 operating system is suggested in the code, however, the path directories for Macintosh and Linux systems must be identified by the user.     

The generated input file from the "Input file generation" notebook must be imported and executed in iPHREEQC. The resultant SELECTED_OUTPUT output file from the iPHREEQC simulation must be placed within the same directory as the "Output processing" notebook for the output data to be processed through the code. The data processing notebook requests the title of the SELECTED_OUTPUT file and scans the file name for "brine" or "scaling" and "pitzer" or "phreeqc", which are used to classify the nature of the iPHREEQC simulation data. File names that lack the simulation details prompt user entries of the information, which directs the code to correctly process the data. The user will finally have the option to export the graph to the same directory as the code in a picture format.    
