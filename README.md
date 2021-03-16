# PHREEQC desalination simulations

Desalination is an increasingly critical technology for human sustainability. The geochemical phenomena within reverse osmosis (RO) desalination modules like brine formation and mineral scaling are researched mitigate fouling of the RO modules. The geochemical software PHREEQC is an established tool for evaluating gechemical phenonema, however, the software has been applied to depict scaling over the distance of an arbitrarily defined RO module. This work couples PHREEQC software with interactive Jupyter Notebooks that create PHREEQC input files and that process the PHREEQC output data into figures and data tables. The following sections stepwise detail the procedures from downloading PHREEQC software through processing simulation outputs.


# PHREEQC installation and execution

PHREEQC is freely downloaded from the US Geological Survey (https://www.usgs.gov/software/phreeqc-version-3). The "PHREEQC execution instructions" PDF file provides explicit directions for installing and the executing PHREEQC on both Windows and Macintosh operating systems. 


# PHREEQC Notebooks 

The use of PHREEQC is facilitated by two Jupyter Notebooks. The Notebooks utilize command-line interfaces that circumvent the need of coding literacy. The "Input file generation" notebook guides the user through creating a PHREEQC input file that will simulate the geochemistry of RO desalination. The "Output processing" notebook guides the user through processing the PHREEQC SELECTED_OUTPUT data file into graphs and data tables.


## Input file generation

The "Input file generation" Notebook is organized into functions that each parameterize a different PHREEQC code block. The make_general_conditions function of the input file notebook defines initial details of the simulation and creates a title for the simulation. The other functions are more elaborate and are separately discussed in the following sections.


### make_solutions


The make_solutions function generates the SOLUTION blocks of the PHREEQC simulation. The SOLUTION blocks of PHREEQC parameterize geochemical solution characteristics like elemental concentrations, alkalinity, and pH. The make_solutions function allows the user to define the feed geochemistry from one of the predefined water bodies, or to customize feed geochemistry through a series of command-line prompts. The predefined options for water bodies reflect experimental literature for the Mediterranean and Red Seas, which are archetypes of highly populated inland seas, and for groundwaters around the continential U.S.A. that represent the produced waters from fracking and are an intriguing water resource. The initial solution of the RO module is parameterized as pure H2O since the solution is rapidly flushed from the RO module during the simulation.


### make_equilibrium_phases

The make_equilibrium_phases function generates the EQUILIBRIUM_PHASES block. The EQUILIBRIUM_PHASES block of PHREEQC defines the mineral scaling equilibrium during RO simulations. The function displays the alphabetized list of minerals and their chemical formulas that are defined by the user selected database, from which the user selects the minerals, and their preexisting quantities, that will be simulated. 


### make_reactive_transport

The make_reactive_transport function generates both REACTION code blocks and the TRANSPORT code block. The function enables users to either simulate desalination through a  Dow-DuPont FILMTEC BW30-400 RO module, which is an archetype, or customize a module parameters through prompts in the notebook interface. The function contains a multitude of backend calculations that interpret permeate flux properties from the module parameters. The function currently only parameterizes single-domain model reactive transport.



### make_selected_output

The make_selected_output function defines the data and the name of the output SELECTED_OUTPUT file. The SELECTED_OUTPUT file contains simulation data in a tab-delimited text file that is imported in the "Output processing" file and interpreted for simulation results. 


### export

The export function exports the generated PHREEQC input file. The function simultaneously prints the complete input PHREEQC code in the notebook interface for the user's review and exports the input file with a predefined naming structure of the day’s date, the simulated water body selection, brine or scaling simulation perspective, and a whole number that increases incrementally until a unique name is achieved. The input file is saved as PQI file which denotes an iPHREEQC input file, although, the file extension is also aceptable for batch software versions of PHREEQC.



## Output processing

The "Output processing" Notebook is likewise organized into discrete functions. The functions each provide complementary ability to interpret the imported SELECTED_OUTPUT file of simulation output data through a few command-line requests and interactions. The SELECTED_OUTPUT file is imported through providing the Notebook with the file name of the output file in the current working directory. The tab-delimited output data is then imported into a Pandas DataFrame for subsequent data manipulations. The following functions of the Notebook generate figure(s) and\or a data table from the DataFrame data, and optionally exports the figure(s).  


### process_selected_output

The process_selected_output function automatically directs the user to the appropriate function for processing the data. The file name of the SELECTED_OUTPUT data file is scanned for the words “brine”, “scaling”, “pitzer”, and “phreeqc” that will specify the nature of the simulation as either brine or scaling and that will specify the database that was used in the simulation. The function then automatically executes the appropriate function for processing the data depending upon the identified words in the output file name or the user parameterizations where the aforementioned words are absent from the SELECTED_OUTPUT file name.


### make_brine_plot

The make_brine_plot function generates a figure for all simulated elements over time. The function first scans the DataFrame elemental data, which data is then plotted via the matplotlib Python library. The time domain for the plot is selected to be all times after the initial solution is flushed from the module, which is represented by the concentration of the most concentrated ion, chloride, first becoming non-zero. The function differentiates between simulations of scaling or brine perspectives and differentially processes the data to match the different data structures. The brine plots of both simulation perspectives, however, ultimately yield a figure, a title and a caption in the notebook interface, and are exported at the user selection. The function additionally yields a data table of the average elemental concentrations, over the same timeframe as the figure, for the brine simulation perspective, which facilitates identifying and quantifying the elemental concentrations that may be difficult to interpret from the generated figure.


### make_scaling_plot

The make_scaling_plot function generates a figure(s) of scaling over module distance. The function commences by quantifying the number of minerals that precipitated during the simulation. The user then selects, based upon the quantity of minerals, whether all of the minerals should be displayed on a single figure, or whether the each mineral should be separately expressed on a separate plots. The figure legends depict both the mineral names and their corresponding chemical formulas to remove barriers of geochemical knowledge from interpreting the figures. Each time series of a mineral is a separate plot, which conveys the accumulation of the mineral over time. The function concludes by displaying the figure(s) and by inquiring whether the figure(s) should be exported, which will direct the user to the export_plot function for the affirmative. 


### export_plot

The export_plot function exports the brine and scaling figure(s) to the working directory. The user is prompted to name the export figures, which are defaulted to the SELECTED_OUTPUT name. The user further defines the image format as either PNG, JPG, and SVG, which is defaulted to JPG. The file name finally consists of a whole number that incrementally increases until a unique file name + extension is determined for the current working directory. 
