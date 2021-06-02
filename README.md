# Reverse Osmosis desalination

Desalination is an increasingly critical technology for human sustainability. Reverse Osmosis (RO) is a prominent desalination technology, yet, some of the geochemical phenomena like brine formation and mineral scaling hinder RO efficacy. We developed a RO scaling simulation (ROSS) software that effectively evaluates brine formation and mineral scaling from desalination. The source code and executable files for this software, which incorporates batch PHREEQC into the executable files, is provided to support community research and user customization.


# Input file generation 

The software first facilitates user parameterization of simulation conditions through intuitive command-line prompts.
The functions and parameterized code blocks for PHREEQC file are detailed in the fillowing sub-sections. 


## make_general_conditions

The "make_general_conditions" function notebook defines initial simulation conditions


## make_solutions

The "make_solutions" function defines the feed solution for the SOLUTION block of PHREEQC. The geochemical characteristics like elemental concentrations, alkalinity, and pH are parameterized with either user-defined quantities or from predefined  options like the Mediterranean and Red Seas or American produced groundwaters that are sourced from literature. 


## make_equilibrium_phases

The "make_equilibrium_phases" function defines mineral scaling for the EQUILIBRIUM_PHASES block. The user selects minerals from an alphabetized list of mineral names and chemical formulae that are defined by the user-selected database. The default set of displayed minerals are those that may be assembled from the defined set of elements in the "make_solutions" function. 


## make_reactive_transport

The "make_reactive_transport" function defines the reactive transport conditions for the REACTION and TRANSPORT blocks. Users may selected predefined parameters of the Dow-DuPont FILMTEC BW30-400 RO module, or users may assign custom spatiotemporal variables for the module characteristics like permeate flux, feed flow rate, membrane thickness, and module length, et cetera. The function includes numerous intermediary values and calculates final PHREEQC parameters from the user inputs.


## make_selected_output

The "make_selected_output" function defines the output file and data for the SELECTED_OUTPUT block. The selected set of simulation data is output to the working directory as a tab-delimited file that subsequently processed for simulation results. 


## export

The "export" function prints and exports the generated PHREEQC input file to the working directory. The export file is established with a predefined naming structure that incorporates simulation details, which prevents overwriting files, facilitates user estimation of the contained simulation, and is used by subsequent data processing code to generate the data figures. 


# PHREEQC file execution
The 


# Output data processing

The ROSS simulations may generate either SELECTED_OUTPUT file of simulation output data through a few command-line requests and interactions. The SELECTED_OUTPUT file is imported through providing the Notebook with the file name of the output file in the current working directory. The tab-delimited output data is then imported into a Pandas DataFrame for subsequent data manipulations. The following functions of the Notebook generate figure(s) and\or a data table from the DataFrame data, and optionally exports the figure(s).  


### process_selected_output

The process_selected_output function automatically directs the user to the appropriate function for processing the data. The file name of the SELECTED_OUTPUT data file is scanned for the words “brine”, “scaling”, “pitzer”, and “phreeqc” that will specify the nature of the simulation as either brine or scaling and that will specify the database that was used in the simulation. The function then automatically executes the appropriate function for processing the data depending upon the identified words in the output file name or the user parameterizations where the aforementioned words are absent from the SELECTED_OUTPUT file name.


### make_brine_plot

The make_brine_plot function generates a figure for all simulated elements over time. The function first scans the DataFrame elemental data, which data is then plotted via the matplotlib Python library. The time domain for the plot is selected to be all times after the initial solution is flushed from the module, which is represented by the concentration of the most concentrated ion, chloride, first becoming non-zero. The function differentiates between simulations of scaling or brine perspectives and differentially processes the data to match the different data structures. The brine plots of both simulation perspectives, however, ultimately yield a figure, a title and a caption in the notebook interface, and are exported at the user selection. The function additionally yields a data table of the average elemental concentrations, over the same timeframe as the figure, for the brine simulation perspective, which facilitates identifying and quantifying the elemental concentrations that may be difficult to interpret from the generated figure.


### make_scaling_plot

The make_scaling_plot function generates a figure(s) of scaling over module distance. The function commences by quantifying the number of minerals that precipitated during the simulation. The user then selects, based upon the quantity of minerals, whether all of the minerals should be displayed on a single figure, or whether the each mineral should be separately expressed on a separate plots. The figure legends depict both the mineral names and their corresponding chemical formulas to remove barriers of geochemical knowledge from interpreting the figures. Each time series of a mineral is a separate plot, which conveys the accumulation of the mineral over time. The function concludes by displaying the figure(s) and by inquiring whether the figure(s) should be exported, which will direct the user to the export_plot function for the affirmative. 


### export_plot

The export_plot function exports the brine and scaling figure(s) to the working directory. The user is prompted to name the export figures, which are defaulted to the SELECTED_OUTPUT name. The user further defines the image format as either PNG, JPG, and SVG, which is defaulted to JPG. The file name finally consists of a whole number that incrementally increases until a unique file name + extension is determined for the current working directory. 
