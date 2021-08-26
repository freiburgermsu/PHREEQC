# import libraries
from matplotlib import pyplot 
from scipy.constants import nano, kilo, liter, day
from itertools import chain
from chempy.properties.water_density_tanaka_2001 import water_density
from pubchempy import get_compounds 
import subprocess
import datetime
import pandas
from math import pi, exp
import time
import os
import re


# calculation constants
grams_over_liters_h2o = 997.07 # @25 degrees celcius
grams_over_moles_h2o = 18.015 
kinematic_flow_velocity = 9.33E-7    #square meters / second
simulated_time_over_computational_time = 9.29    


# simulation constants
timesteps_over_cell = 1

# useful conditions and parameters
possible_answers = ['y', 'n']


class ROSSPkg():

    def __init__(self):
        # print the software introductory box
        print('\n\n')
        message = ('''* ROSS 1.1.1 *
        Reverse Osmosis Scaling Simulation
        by Andrew P. Freiburger, Ethan S. Chan, and Heather L. Buckley
        Summer 2021, Green Safe Water Lab, University of Victoria''')
        
        lines = message.split('\n')
        space = " " * indent
        width = max(map(len, lines))
        upper = f'╔{"═" * (width + indent * 2)}╗\n'
        middles = ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
        lower = f'╚{"═" * (width + indent * 2)}╝' 
        box = chain(upper, middles, lower)
        box_print = ''.join(box)
        print(box_print, '\n')

        
        # establish the general organization structures
        self.parameters = {}
        self.variables = {}
        self.results = {}

        
    def define_general(self, os, phreeqc_path, database_selection, simulation_type, simulation_title, indent = 1)
        '''
        Establish general conditions
        '''
        self.parameters['water_mw'] = get_compounds('water', 'name')[0].molecular_weight
        self.parameters['water_density'] = water_density()
        
        # parameterize the input file
        self.parameters['os'] =  os
        self.parameters['phreeqc_path'] = phreeqc_path
        self.parameters['database_selection'] = database_selection 
        self.parameters['simulation_type'] = simulation_type
        self.parameters['title'] = simulation_title
        
        title_line = 'TITLE\t %s' %(simulation_title)
        if os == 'Windows':
            self.parameters['database_path'] = phreeqc_path + '\\database\\%s.dat' %(database_selection)
            database_line = 'DATABASE %s' %(database_path)
            self.results['general_conditions'] = [database_line, title_line]
        else:
            self.results['general_conditions'] = [title_line]
            
            
    def database_parsing():
        db = self.parameters['database_selection']
        database = pandas.read_table(f'./databases/{db}.dat', sep='\n')
        database.columns = ['content']

        elements_rows = []
        minerals_rows = []
        for index, row in database.iterrows():
            if row['content'] == 'SOLUTION_MASTER_SPECIES':
                while database.at[index, 'content'] != 'SOLUTION_SPECIES':
                    split_row = database.at[index, 'content'].split('\t')
                    elements_rows.append(split_row)
                    index += 1
                continue

            if row['content'] == 'PHASES':
                loop = False
                while phreeqc_database.at[index, 'content'] != '# Gases from LLNL...':
                    print(phreeqc_database.at[index, 'content'])
                    if not loop:
                        minerals_rows.append(['phases', 'formula'])
                        loop = True
                        index += 1
                        continue

                    if re.search('(^\w+$)',phreeqc_database.at[index, 'content']):
                        reactants = phreeqc_database.at[index+1, 'content'].split(' = ')[0]
                        formula = reactants.split(' + ')[0].strip()
                        print(formula)
                        minerals_rows.append([phreeqc_database.at[index, 'content'], formula])
                    index += 1

        # define the elements content for the database
        elements = pandas.DataFrame(elements_rows)
        elements.columns = elements.iloc[2]
        elements.rename(columns = {'#element':'elements'}, inplace = True)
        elements = elements.iloc[pandas.RangeIndex(len(elements)).drop([x for x in range(4)])]
        
        element_list = list(elements['elements'])
        self.parameters['element_list'] = []
        for element in element_list:
            if not re.search('\#', element):
                self.parameters['elements_list'].append(element)

        # define the minerals content for the database
        minerals = pandas.DataFrame(minerals_rows)
        minerals.columns = minerals.iloc[0]
        minerals = minerals.drop(0)
        self.parameters['mineral_list'] = list(minerals['phases'])
        self.parameters['formula_list'] = list(minerals['formula'])
        

    def transport(self, module_selection, module_characteristics = {}, quantity_of_modules, self.parameters['cells_per_module'], domain, output_perspective):
        '''
        Define the TRANSPORT block
        '''
        # parameterize the module 
        if module_selection == 'BW30-400':
            module_diameter =  201                   #mm
            permeate_tube_diameter =  29             #mm
            module_length =  1.016                   #m
            permeate_flow = 40                       #cubic meters / day
            max_feed_flow = 15.9                     #cubic meters / hour
            membrane_thickness = 250 * (constants.milli / constants.nano)   #mm
            feed_thickness = 0.7112                  #mm
            permeate_thickness = 0.3                 #mm
            polysulfonic_layer_thickness = 0.05      #mm 
            support_layer_thickness = 0.15           #mm
            repeated_membrane_thickness = 2 * membrane_thickness + feed_thickness + permeate_thickness + 2 * polysulfonic_layer_thickness + 2 * support_layer_thickness      #mm
            print('\nMembrane thickness:', '%s (mm)' %(repeated_membrane_thickness))

        elif module_selection == 'Custom':
            module_diameter = module_characteristics['diameter']                             #mm
            permeate_tube_diameter = module_characteristics['permeate_diameter']             #mm
            module_length = module_characteristics['length']                                 #m
            permeate_flow = module_characteristics['permeate_flow']                          #cubic meters / day
            max_feed_flow = module_characteristics['feed_flow']                              #cubic meters / day
            membrane_thickness = module_characteristics['membrane_thickness']                #mm
            feed_thickness = module_characteristics['feed_spacer_thickness']
            permeate_thickness = input('- What is the permeate spacer thickness? (mm) __ ')
            polysulfonic_layer_thickness = module_characteristics['polysulfonic_thickness']
            support_layer_thickness = module_characteristics['support_thickness']
            repeated_membrane_thickness = 2 * membrane_thickness + feed_thickness + permeate_thickness + 2 * polysulfonic_layer_thickness + 2 * support_layer_thickness
            print('Thickness of the repeated membrane unit (mm): %s' %(repeated_membrane_thickness))

        self.parameters['quantity_of_modules'] = quantity_of_modules
        self.parameters['cells_per_module'] = cells_per_module
        cell_length = module_length / self.parameters['cells_per_module']     #meters

        # calculate module properties
        module_cross_sectional_area = module_diameter**2 * math.pi / 4        #squared millimeters
        permeate_tube_cross_sectional_area = permeate_tube_diameter**2 * math.pi / 4     #squared millimeters
#             filtering_layer_thicknes = (module_diameter - permeate_thickness) / 2         #millimeters
        filtration_cross_sectional_area = module_cross_sectional_area - permeate_tube_cross_sectional_area        #squared millimeters
        feed_cross_sectional_area = (feed_thickness / repeated_membrane_thickness) * filtration_cross_sectional_area     #squared millimeters
        feed_volume = feed_cross_sectional_area * module_length * kilo**2   #cubic meters
        feed_mass = feed_volume * self.parameters['water_density'] * (1 / constants.liter) / constants.kilo    #kilograms, which assumes pure water for mass
        feed_moles = feed_mass * constants.kilo / self.parameters['water_mw'] 

        # calculate fluid flow characteristics
        feed_velocity = max_feed_flow / (feed_cross_sectional_area * kilo**2) / day     #meters / second
        reynolds_number = feed_velocity * module_length / kinematic_flow_velocity
        print('Reynold\'s number: ', reynolds_number)

        # calculate module cell characteristics
        self.parameters['feed_mass_cell'] = feed_mass / self.parameters['cells_per_module']      #kg
        feed_moles_cell = feed_moles / self.parameters['cells_per_module']    #moles

        # calculate simulation parameters that will populate the PHREEQC simulation 
        maximal_timestep = cell_length / feed_velocity * timesteps_over_cell         #seconds, from the Courant condition
        self.parameters['permeate_removal_per_cell'] = maximal_timestep * permeate_flow / day * liter * self.parameters['water_density'] / self.parameters['water_mw'] / self.parameters['cells_per_module']      #moles / cell

        # define the transport black
        transport_line = '\nTRANSPORT'
        cells_line = '-cells\t\t\t%s' %(self.parameters['cells_per_module'])
        
        simulation_shifts = (2*self.parameters['cells_per_module']*self.parameters['quantity_of_modules'])
        shifts_line = '-shifts\t\t\t%s' %(simulation_shifts)
        lengths_line = '-lengths\t\t%s' %(cell_length)
        timestep_line = '-time_step\t\t%s\t# the Courant condition is satisfied with the cell_length of %s m and the feed velocity of %s m/s' %(maximal_timestep, cell_length, feed_velocity)
        initial_time_line = '-initial_time\t\t0'    
        boundary_conditions_line = '-boundary_conditions\tconstant\tconstant \t # Dirichlet boundary condition'
        
        if domain == 'single':
            domain_line = '-stagnant\t\t0\t0\t0\t0 \t # single domain\n#^stagnant cells\t^exchange factor\t^CP volume\t^bulk volume'
        elif domain == 'dual':
            domain_line = '-stagnant\t\t1\t1\t0.1\t0.9 \t # dual domain\n#^stagnant cells\t^exchange factor\t^CP volume\t^bulk volume'

        self.parameters['output_perspective'] = output_perspective
        if self.parameters['output_perspective'] == 'Scaling':
            punch_cells_line = '-punch_cells\t\t1-%s' %(self.parameters['cells_per_module'])
            punch_frequency_line = '-punch_frequency\t%s' %(self.parameters['cells_per_module'])
        elif self.parameters['output_perspective'] == 'Brine':
            punch_cells_line = '-punch_cells\t\t%s' %(self.parameters['cells_per_module'])
            punch_frequency_line = '-punch_frequency\t1'       

        # create the transport block
        self.results['transport_block'] = []
        self.results['transport_block'].extend((transport_line, cells_line, shifts_line, lengths_line, timestep_line, initial_time_line, boundary_conditions_line, domain_line, punch_cells_line, punch_frequency_line))


    def reaction(self, permeate_approach, permeate_efficiency, head_loss, final_cf):
        '''
        Define the REACTION block
        '''
        cfs = []
        cell_moles = []
        reaction_parameters = []
        iteration = 0
        cumulative_cf = 1 
        self.results['reaction_block'] = []
        for module in range(quantity_of_modules):
            module_previous_moles_removed = 0
            if permeate_approach == 'Linear permeate':
                if iteration == 0:
                    initial_moles_cell = feed_moles_cell

                initial_moles_removed = self.parameters['permeate_removal_per_cell'] * 2 / (1 + math.exp(head_loss))
                final_moles_removed = initial_moles_removed * math.exp(head_loss)
                try:
                    removed_moles_slope = ((final_moles_removed - initial_moles_removed) / (self.parameters['cells_per_module'])) / permeate_efficiency
                except:
                    removed_moles_slope = 0
                average_moles_removed = (final_moles_removed + initial_moles_removed) / 2
                print('(Removed moles / cell) slope: ', removed_moles_slope)

                for cell in range(self.parameters['cells_per_module']):
                    removed_moles_in_cell = (cell * removed_moles_slope + initial_moles_removed)
                    reaction_parameters.append(removed_moles_in_cell)
                    module_previous_moles_removed += removed_moles_in_cell
                    
                cf = initial_moles_cell / (initial_moles_cell - module_previous_moles_removed)
                cumulative_cf *= cf

                initial_moles_cell -= module_previous_moles_removed
                iteration += 1

            if permeate_approach == 'Linear CF':
                module_iteration = 0
                initial_cf = 1

                cf_slope = (final_cf - initial_cf) / self.parameters['cells_per_module']
                for cell in range(self.parameters['cells_per_module']):
                    cell_cf = (cell+1) * cf_slope + initial_cf
                    cfs.append(cell_cf)    

                for cf in cfs:
                    moles_to_be_removed =  feed_moles_cell - (feed_moles_cell / cf)
                    if module_iteration == 0:
                        initial_moles_cell = feed_moles_cell
                        reaction_parameters.append(moles_to_be_removed)
                    if module_iteration > 0:
                        module_previous_moles_removed += reaction_parameters[-1] 
                        reaction_parameter = moles_to_be_removed - module_previous_moles_removed
                        reaction_parameters.append(reaction_parameter)
                        initial_moles_cell -= moles_to_be_removed

                    module_iteration += 1

                cf = cfs[-1]
                cumulative_cf *= cf
                initial_moles_cell = feed_moles_cell - moles_to_be_removed   # moles_to_be_removed = the final quantity of moles is calculated with naming from the linear permeate flux method to be consistent  

            if self.parameters['simulation_type'] == 'transport':
                reaction_line.append('\n')
                for cell in range(self.parameters['cells_per_module']):
                    cell_number = (cell+1) + self.parameters['cells_per_module'] * module
                    reaction_line = 'REACTION %s' %(cell_number)
                    
                    if (cell+1) < self.parameters['cells_per_module']:
                        reaction_line += f'\n\tH2O -1; {reaction_parameters[cell_number-1]}' 
                    elif (cell+1) == self.parameters['cells_per_module']:
                        reaction_line += f'''\n\tH2O -1; {reaction_parameters[cell_number-1]}
        INCREMENTAL_REACTIONS \ttrue''' %(reaction_parameters[cell_number-1])     

                    self.results['reaction_block'].append(reaction_line)

            elif self.parameters['simulation_type'] == 'evaporation':
                parameter_quantity = 15                          
                recursive_assymtote_multiplier = 1.335449219     # ??? arbitrary assignment of kg of water in the simulation?
                total_moles_removed = sum(reaction_parameters)
                initial_evaporation_parameter = total_moles_removed / recursive_assymptote_multiplier
                evaporation_reaction_parameters = ['0', initial_evaporation_parameter]  # ???
                for parameter in range(1, parameter_quantity):
                    evaporation_reaction_parameter = evaporation_reaction_parameters[parameter] * 1/4
                    evaporation_reaction_parameters.append(evaporation_reaction_parameter)

                # define the reaction block
                reaction_line = 'REACTION 1'
                reaction_line += '\n\tH2O -1; '
                reaction_line += ' '.join(evaporation_reaction_parameters) 
                self.parameters['reaction_block'] = [reaction_line]

            # the calculated reaction parameters will be added and printed to a generated PHREEQC input file
            final_solution_mass = initial_moles_cell * grams_over_moles_h2o / constants.kilo  #kg water mass
            final_cf_cell = self.parameters['feed_mass_cell'] / final_solution_mass
            print('Effluent module %s CF:' %(module + 1), final_cf_cell)

            if self.parameters['os'] == 'Windows':
                self.results['reaction_block'].append('#%s' %(permeate_approach))
                if permeate_approach == 'Linear permeate':
                    self.results['reaction_block'].append('''
        #Permeate efficiency parameter: %s
        #Head loss parameter: %s''' %(permeate_efficiency, head_loss))
                self.results['reaction_block'].append('''    #Effluent module %s:
        #Estimated CF: %s
        #Estimated solution mass: %s\n\n''' %(module + 1, cumulative_cf, final_solution_mass))


    def make_solutions(self, water_selection, custom_water_parameters = {}):
        """
        Specify the SOLUTION block of the simulation.
        """
        import json

        # create the solution line of the input file
        self.results['solution_block'] = []
        self.parameters['water_selection'] = water_selection

        if self.parameters['simulation_type'] == 'transport':
            initial_solution_line = '\nSOLUTION 0\t%s' % (solution_description)
        elif self.parameters['simulation_type'] == 'evaporation':
            initial_solution_line = '\nSOLUTION 1\t%s' % (solution_description)
        self.results['solution_block'].append(initial_solution_line)

        #=============================================================================
        # determine which predefined solution should be simulated

        if water_selection is not None:       
            # import the predefined water body
            water_file_path = f'./water_bodies/{water_selection}.json'
            water_body = json.load_json(water_file_path)
            
            self.parameters['elements'] = []
            elements_lines = []
            for content, information in water_body.items():
                if content == 'element':
                    for element, information2 in information.items():
                        if element in self.parameters['elements_list']:
                            self.parameters['elements'].append(element)
                            conc = information2['concentration (ppm)']
                            ref = information2['reference']
                            if len(conc) <= 3:
                                elements_lines.append(f'{element}\t\t{conc}\t#{ref}')
                            else:
                                elements_lines.append(f'{element}\t{conc}\t#{ref}')
                        else:
                            print('\n--> ERROR: The {} element is not accepted by the {} database')
                                
                elif content == 'temperature':
                    temperature = information['celcius']
                    temperature_reference = information['reference']
                elif content == 'pe':
                    pe = information['value']
                    pe_reference = information['reference']
                elif content == 'Alkalinity':
                    alkalinity = information['value']
                    alkalinity_reference = information['reference'] 
                elif content == 'pH':
                    ph = information['value']
                    ph_reference = information['reference']

        if custom_water_parameters != {}:
            for content, information in custom_water_parameters.items():
                
                if content == 'element':
                    for element, information2 in information.items():
                        conc = information2['concentration (ppm)']
                        ref = information2['reference']
                        if element in self.parameters['elements_list']:
                            if element in self.parameters['elements']:
                                element_index = self.parameters['elements'].index(element)                           
                                if len(conc) <= 3:
                                    elements_lines[element_index] = (f'{element}\t\t{conc}\t#{ref}')
                                else:
                                    elements_lines[element_index] = (f'{element}\t{conc}\t#{ref}')
                            else:
                                self.parameters['elements'].append(element)
                                if len(conc) <= 3:
                                    elements_lines.append(f'{element}\t\t{conc}\t#{ref}')
                                else:
                                    elements_lines.append(f'{element}\t{conc}\t#{ref}')
                        else:
                            print('\n--> ERROR: The {} element is not accepted by the {} database')
                                    
                            
                # create the temperature line of the input file
                elif content == 'temperature':                    
                    temperature = custom_water_parameters['temperature']['value']
                    temperature_reference = custom_water_parameters['temperature']['reference']

                elif content == 'pe':       
                    pe = custom_water_parameters['pe']['value']
                    pe_reference = custom_water_parameters['pe']['reference']

                elif content == 'Alkalinity':
                    alkalinity = custom_water_parameters['Alkalinity']['value']
                    alkalinity_reference = custom_water_parameters['Alkalinity']['reference'] 

                elif content == 'pH':
                    ph = custom_water_parameters['ph']['value']
                    ph_reference = custom_water_parameters['ph']['reference']
                    
                    
        # parameterize the lines of the SOLUTIONS block
        temperature_line = f'temp \t {temperature} \t #{temperature_reference}.'
        ph_line = f'pH \t\t {ph} charge #{ph_reference}'
        pe_line = f'pe \t\t {pe} \t   #{pe_reference} // 4.00 is the default (?)'
        alkalinity_line = f'Alkalinity \t {alkalinity} #{alkalinity_reference}'
        unit_line = 'units \t ppm' 
        elements_line = '\n'.join(elements_lines)
        water_mass = self.parameters['feed_mass_cell']
        if water_selection == 'Bakken formation':
            water_line = f'-water \t{water_mass}\t#TDS=300 ppthousand [before fudging]' 
        elif water_selection == 'German Basin':
            water_line = f'-water \t{water_mass}\t#TDS=314 ppthousand [before fudging]'
        else:
            water_line = f'-water \t{water_mass}'

        self.results['solution_block'].extend([temperature_line, ph_line, pe_line, alkalinity_line, unit_line, elements_line, water_line])

        #parameterize the initial module solution
        if self.parameters['simulation_type'] == 'transport':
            total_cells = self.parameters['cells_per_module'] * quantity_of_modules
            feed_solution_line = f'\nSOLUTION 1-{total_cells}\tInitial solution in the RO module'
            self.results['solution_block'].extend([feed_solution_line,'temp \t 25','units \t ppm'])

            for element in elements:
                element_concentration = 0
                element_line = f'{element}\t{element_concentration}'    
                self.results['solution_block'].append(element_line)

            water_line = '-water \t %s' %(self.parameters['feed_mass_cell'])
            self.results['solution_block'].append(water_line)


    def make_equilibrium_phases(self, block_comment, ignored_minerals = [], existing_parameters = {}):
        """
        Specify the EQUILIBRIUM_PHASES block of the simulation.
        """
        # define mineral sizes for later spacing
        short_mineral_names = ['Barite', 'Gypsum','Halite', 'Natron', 'Quartz', 'Talc','Trona', 'Borax', 'Albite', 'K-mica','Illite', 'Pyrite', 'Sulfur',]
        long_mineral_names = ['Anthophyllite', 'Hexahydrite', 'Leonhardite', 'Nesquehonite', 'Pentahydrite', 'Portlandite','Sepiolite(d)', 'Boric_acid,s', 'K2B4O7:4H2O', 'NaB5O8:5H2O', 'Rhodochrosite', 'Strontianite','Hydroxyapatite', 'Chlorite(14A)', 'Mackinawite', 'Hausmannite', 'Pyrochroite']

        # determine the set of possible minerals 
        remaining_characters = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '(', ')', 'I', '•']
        permitted_auxillary_reagents = ['OH', 'H2O','O']

        self.variables['described_minerals'] = []
        self.variables['described_mineral_formulas'] = []
        for mineral in self.parameters['mineral_list']:
            mineral_index = self.parameters['mineral_list'].index(mineral)
            mineral_formula = self.parameters['formula_list'][mineral_index]

            for element in self.parameters['elements']:
                if re.search(element, mineral_formula):
                    mineral_formula = re.sub(element, '', mineral_formula)

            for reagent in permitted_auxillary_reagents:
                if re.search(reagent, mineral_formula):
                    mineral_formula = re.sub(reagent, '', mineral_formula)

            if all(character in remaining_characters for character in mineral_formula):
                self.variables['described_minerals'].append(mineral)
                self.variables['described_mineral_formulas'].append(formulas[-1])
            else:
                print(f'The {mineral} mineral contained {mineral_formula} remaining characters after substitution with the parameterized elements.')


        # define the equilibrium_phases block
        self.results['equilibrium_phases_block'] = []
        if self.parameters['simulation_type'] == 'transport':
            equilibrium_phases_number = '1-{}'.format(self.parameters['cells_per_module']*quantity_of_modules)
        elif self.parameters['simulation_type'] == 'evaporation':
            equilibrium_phases_number = '1'
            
        equilibrium_phases_line = f'\nEQUILIBRIUM_PHASES {equilibrium_phases_number}\t{block_comment}'
        self.results['equilibrium_phases_block'].append(equilibrium_phases_line)

        # define the equilibrium_phases lines for the code block
        for possible_mineral in self.variables['described_minerals']:
            if possible_mineral not in ignored_minerals:
                formula = {self.variables['described_mineral_formulas'][described_minerals.index(possible_mineral)]
                if possible_mineral in short_mineral_names:
                    mineral_line = f'{possible_mineral}\t\t' 
                elif possible_mineral == 'Ca-Montmorillonite':
                    mineral_line = f'{possible_mineral}'
                else:
                    mineral_line = f'{possible_mineral}\t'

                if possible_mineral in existing_parameters:
                    for key, value in existing_parameters[possible_mineral].items():
                        if key == 'saturation':
                            mineral_saturation = value['saturation']
                            mineral_line += f'\t{mineral_saturation}'
                        if key == 'initial_moles':
                            initial_moles = value['initial_moles']
                            mineral_line += f'\t{initial_moles}'
                else:
                    mineral_line += f'\t0\t0'
                    
                self.results['equilibrium_phases_block'].append(mineral_line)        


    def make_selected_output(self, output_filename = None):
        '''
        Specify the output file after a PHREEQC simulation
        '''
        # create parameter lines 
        if output_filename is not None:
            count = 0
            selected_output_file_name = '{}_{}_{}_{}_{}_{}'.format(datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['output_perspective'], count) 
            while os.path.exists(f'{selected_output_file_name}.txt':
                count += 1
                selected_output_file_name = '{}_{}_{}_{}_{}_{}'.format(datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['output_perspective'], count) 
        else:
            selected_output_file_name = output_filename

        selected_output_file_name += '.txt'
        self.parameters['selected_output_file_name'] = selected_output_file_name

        minerals_line = ''
        for mineral in self.parameters['elements']:
             minerals_line += f' {mineral}'

        elements_line = ''
        for element in self.parameters['elements']:
             elements_line += f' {element}'

        # define parameter lines
        first_line = '\nSELECTED_OUTPUT'
        file_name_line = f'-file\t\t\t{selected_output_file_name}'
        reaction_line = '-reaction\t\ttrue'
        temperature_line = '-temperature\t\ttrue'
        total_elements_line = '-totals\t\t\t' + elements_line    
        saturation_indices_line = f'-saturation_indices\t{minerals_line}'
        equilibrium_phases_line = f'-equilibrium_phases\t{minerals_line}'
        ph_line = '-pH\t\t\ttrue'
        solution_line = '-solution'
        time_line = '-time\t\t\ttrue'
        distance_line = '-distance\t\ttrue'
        simulation_line = '-simulation\t\ttrue'
        high_precision_line = '-high_precision\ttrue'
        step_line = '-step'
        water_line = '-water'

        # establish the selected_output_block
        self.results['selected_output_block'] = []
        self.results['selected_output_block'].extend((first_line, file_name_line, reaction_line, temperature_line, total_elements_line, saturation_indices_line, equilibrium_phases_line, ph_line, time_line, distance_line, simulation_line, high_precision_line, solution_line, step_line,water_line))


    def export(self, print_block = True):
        """
        View and export the PHREEQC input file
        """
        global complete_lines
        global input_file_name

        complete_lines = chain(general_conditions, self.results['solution_block'], self.results['equilibrium_phases_block'], self.results['reaction_block'], self.results['selected_output_block'], self.results['transport_block']) 

        # create the export input file
        file_number = 0
        proposed_file_name = '%s_ROSS_%s_%s_%s.pqi' %(datetime.date.today(), self.parameters['water_selection'], self.parameters['output_perspective'], file_number)
        if not os.path.exists(proposed_file_name):
            self.parameters['input_file_name'] = proposed_file_name
        elif os.path.exists(proposed_file_name):
            while os.path.exists('%s_ROSS_%s_%s_%s_%s.pqi' %(datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['output_perspective'], file_number)):
                file_number += 1
                self.parameters['input_file_name'] = '%s_ROSS_%s_%s_%s_%s.pqi' %(datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['output_perspective'], file_number)

        # printing and exporting the input file
        input_file = open(input_file_name,'a')
        for line in complete_lines:
            if printing_block:
                print(line)
            input_file.write(line + '\n')
        input_file.close()


    def execute():
        '''
        Execute the input file through the batch PHREEQC software 
        '''
        global database_selection
        global input_file_name

        announcement = '\nExecute the input file:'
        print(announcement, '\n', '='*len(announcement))
        #subprocess.call("start powershell", shell=True)
        working_directory = os.getcwd()
        phreeqc_path = 'C:\\Program Files\\USGS\\phreeqc-3.6.2-15100-x64'
        if input_selection == 'n':
            database_selection = input('''- What database do you select?
            < pitzer > or < phreeqc >  __ ''')
            input_file_name = input('- What is the input file name?')
            input_path = os.path.join(working_directory, input_file_name)
            while not os.path.exists(input_path):
                print('''ERROR: The input file path {} does not exist. Provide a valid input file path.'''.format(input_path))
                input_file_name = input('What is the input file name?')
                input_path = os.path.join(working_directory, input_file_name)

            output_file_name = re.sub('(?<=\.)(.+)' , 'pqo', input_file_name)
            database_path = ''

        else:
            input_path = os.path.join(working_directory, input_file_name)
            output_file_name = re.sub('pqi', 'pqo', input_file_name)
            database_path = os.path.join(phreeqc_path, 'database\\%s.dat' %(database_selection))        

        output_path = os.path.join(working_directory, output_file_name)
        bat_path = os.path.join(phreeqc_path, 'bin\\phreeqc.bat')
        print('\ninput file path: {}'.format(input_path))

        proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
        command = str.encode("\"" + bat_path + "\" \"" + input_path + "\" \"" + output_path + "\"  \"" + database_path + "\"\n") 
        proc.stdin.write(command)
        proc.stdin.close()  
        proc.wait()

        if os.path.exists(output_path):
            if input_selection == 'y':
                if os.path.exists(selected_output_file_name):
                    print('The execution is complete.')
                else:
                    print('ERROR: The SELECTED_OUTPUT file is missing. The simulation failed to simulate.')
            else:
                print('The execution is complete.')
        else:
            print('\nERROR: The simulation failed to execute.')


    def process_selected_output():
        """
        Interpreting the PHREEQC SELECTED_OUTPUT file and conducting the plotting functions
        """
        global selected_output_file_name
        global self.parameters['simulation_type']
        global simulation_perspective
        global initial_solution_mass
        global selected_output_file
        global graphical_selection
        global self.parameters['output_perspective']
        global simulation_cf
        global csv_data   
        global database


        announcement = '\nProcess the selected_output file:'
        print(announcement, '\n', '='*len(announcement))

        # determining the appropriate variables
        if input_selection == 'y':
            selected_output_file = selected_output_file_name
            simulation_perspective = self.parameters['output_perspective']
            database = database_selection 

        else:       
            if execute_selection == 'y':
                input_file = open(input_file_name, 'r')
                for line in input_file:
                    if re.search('(-file\s+)', line):
                        selected_output_file = re.sub('(-file\s+)', '', line)
                        selected_output_file = re.sub('(\n)', '', selected_output_file)

                    if re.search('(TRANSPORT)', line):
                        self.parameters['simulation_type'] = 'transport'
                    else:
                        self.parameters['simulation_type'] = 'evaporation'

            else:  
                selected_output_file = input('''- What is the name and\or path of the simulation file?
                Include the extension, like < .txt > __ ''')
                selected_output_file_name = re.sub('(?<=\.)(.+)', '', selected_output_file)
                while not os.path.exists(selected_output_file):
                    print('ERROR: The simulation file is missing from the current directory.')  
                    selected_output_file = input('What is the name and\or path of the simulation file?')  

                for line in selected_output_file:
                    if re.search('(TRANSPORT)', line):
                        self.parameters['simulation_type'] = 'transport'
                    else:
                        self.parameters['simulation_type'] = 'evaporation'

            #determining the scope of the simulation data
            if re.search('(Scaling)', selected_output_file, flags=re.IGNORECASE):
                simulation_perspective = 'Scaling'

            elif re.search('(Brine)', selected_output_file, flags=re.IGNORECASE):
                simulation_perspective = 'Brine'

            else:     
                simulation_perspective = input('''- Is the output file representative of a simulation for scaling or brine?
                < Scaling > or < Brine >''')
                while simulation_perspective != 'Brine' and simulation_perspective != 'Scaling':
                    print('''ERROR: The value is not one of the options. Select one of the choices to proceed.''')  
                    simulation_perspective = input('- Is the output file representative of a simulation for scaling or brine?')

            self.parameters['output_perspective'] = ''

        graphical_selection = input('''- Would you like to view the effluent brine or the module scaling from your {} simulation?
        < Brine > or < Scaling >'''.format(self.parameters['output_perspective']))
        while graphical_selection != 'Brine' and graphical_selection != 'Scaling':
            print('''ERROR: The value is not one of the options. Select one of the choices to proceed.''')  
            graphical_selection = input('''- Would you like to view brine over time in the module < Brine > , 
            or would your like to view scaling over distance in the module < Scaling > ? __ ''')

        # preparing the SELECTED_OUTPUT file into a dataframe
        selected_output = open(selected_output_file, 'r')
        original_data = pandas.read_table(selected_output, sep = '\t')
        csv_data = pandas.DataFrame(original_data)
        for column in csv_data.columns:
            new_column = column.strip()
            csv_data.rename(columns={column:new_column}, inplace = True)

        initial_solution_mass = csv_data.at[0, 'mass_H2O']
        final_solution_mass = csv_data['mass_H2O'].iloc[-1]
        simulation_cf = initial_solution_mass / final_solution_mass

        # conducting the appropriate visualization function
        if graphical_selection == 'Brine':
            make_brine_plot()
        elif graphical_selection == 'Scaling':
            make_scaling_plot()
        else:
            print('''ERROR: Reconfigure the SELECTED_OUTPUT simulation file name and\or the above visualizaiton parameters.''')


    def make_brine_plot():
        """
        Generate plots of the elemental concentrations from effluent brine in the PHREEQC SELECTED_OUTPUT file  
        """
        global export_name
        global plot_title 
        global figure

        elements = []
        for column in csv_data.columns:
            if re.search('([A-Z][a-z]?(\(\d\))?){1}$', column) and not re.search('(_|H2O|pH)', column):
                elements.append(column)

        csv_data.drop(csv_data.index[:3], inplace=True)

        # plot the brine concentrations figure
        unit = 'mol/kgw'
        concentration_array = []
        pyplot.figure(figsize = (17,10))
        plot_title = (input('''- What is the title of the plot?
        Default = Effluent brine elemental concentrations ____ ''')) or 'Effluent brine elemental concentrations'
        pyplot.title(plot_title, fontsize = 'xx-large')
        pyplot.xlabel('time (s)', fontsize = 'x-large')
        pyplot.ylabel('concentration (%s)' %(unit), fontsize = 'x-large')
        pyplot.grid(True)

        for element in elements:  
            concentration_serie = []
            time_serie = []
            initial_solution_time = 0
            for index, row in csv_data.iterrows():
                if csv_data.at[index, 'Cl'] == 0:
                    initial_solution_time += 1
                    #print('yes')
                else:
                    concentration_serie.append(csv_data.at[index, element])
                    time_serie.append(csv_data.at[index, 'time'] - initial_solution_time * (csv_data.at[index, 'time'] - csv_data.at[index-1, 'time']))

            pyplot.plot(time_serie,concentration_serie)

        plot_caption = '''\n\nBrine Figure:\n%s 
        The effluent concentrations of each existing element in the brine. Brine plots from brine data rapidly reach a steady state elemental concentration. Brine plots from scaling data consist of vertical concentrations that represent the concentrations at each distance throughout the RO module at the specified time, where the low end represents the influent concentration while the high end represents the effluent concentration.''' %('='*len('Brine Figure'))

        pyplot.legend(elements, loc='best', title = 'non-zero elements', fontsize = 'x-large')
        pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
        pyplot.yscale('log')
        figure = pyplot.gcf()
        print('\nClose the figure to proceed.')
        print(plot_caption)
        pyplot.show()
        export_figure = input('''Will you export the figure? < n > or < y >  __  ''')
        while export_figure not in possible_answers:
            print('ERROR: The entered mineral is not among the options.')  
            export_figure = input('''\nWill you export the figure? < n > or < y >  __  ''')

        # create a complementary concentration data table for the scaling figure 
        loop_iteration = 1
        table_view = {}
        average_concentrations_table = pandas.DataFrame()
        index_elements = []
        if simulation_perspective == 'Scaling':
            for element in elements:
                quantity_of_steps_index = 0
                average_iteration = 0
                time_serie = []            
                time_averages = {}
                for index, row in csv_data.iterrows():
                    if csv_data.at[index, 'time'] == 0:
                        time_serie.append(csv_data.at[index,element])                 
                        quantity_of_steps_index += 1                    

                    elif csv_data.at[index-1,'soln'] == quantity_of_steps_index:       
                        #process the complete time serie
                        try:
                            average_concentration = sum(time_serie) / len(time_serie)
                        except:
                            average_concentration = 0

                        #print(average_concentration)
                        table_view['Time (s): %s' %(average_iteration * quantity_of_steps_index)] = average_concentration
                        average_iteration += 1
                        #print('mid-End')   

                        #begin the new time serie
                        time_serie = []
                        time_averages = {}
                        time_serie.append(csv_data.at[index,element])
                        #print('middle')    

                    elif index == len(csv_data[element]) + 2:       
                        time_serie.append(csv_data.at[index,element])            
                        try:
                            average_concentration = sum(time_serie) / len(time_serie)
                        except:
                            average_concentration = 0
                        table_view['Time (s): %s' %(round(average_iteration * quantity_of_steps_index), 1)] = average_concentration
                        average_iteration += 1
                        #print('end-End')  
                        index_elements.append(element)
                        average_concentrations_table = average_concentrations_table.append(table_view, ignore_index=True)
                        loop_iteration += 1

                    else:
                        #print('middle')      
                        time_serie.append(csv_data.at[index,element])

            average_concentrations_table.index = index_elements
            average_concentrations_table.index.name = 'Elements'
            dataframe_title = 'Average elemental molal concentrations of the feed water in the RO module for each %s seconds of simulation:' %(quantity_of_steps_index)


        # create a complementary concentration data table for the brine figure 
        elif simulation_perspective == 'Brine':
            total_time = csv_data['time'].iloc[-1]
            for element in elements:  
                concentration_serie = []
                time_serie = []
                for index, row in csv_data.iterrows():
                    if csv_data.at[index, 'Cl'] != 0:
                        concentration_serie.append(csv_data.at[index,element])      
                average_concentration = sum(concentration_serie) / len(concentration_serie)
                table_view['%s' %(element)] = average_concentration

            average_concentrations_table = average_concentrations_table.append(table_view, ignore_index=True)
            average_concentrations_table.rename(index = {0:'Concentrations (molal)'}, inplace = True)
            dataframe_title = 'Average elemental molal concentrations of the feed water in the RO module over %s seconds of simulation:' %(total_time)

        print('\n\n\n',dataframe_title,'\n%s'%('='*len(dataframe_title)))
        print(average_concentrations_table)

        # export the output graphic
        if export_figure == 'y':
            export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file_name)
            if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
            else:
                export_name = export_filename_progenitor

            export_option = input('''- Would you like to export the figure? < y > or < n > __ ''')
            if export_option == 'y':
                export_plot()
            while export_option not in possible_answers:
                print('''ERROR: The value is not one of the options.''')  
                export_option = input('''- Would you like to export the figure?
                < y > or < n >''')

    def make_scaling_plot():
        """
        Generate plots of scaling along the module distance in the PHREEQC SELECTED_OUTPUT file  
        """
        global initial_solution_mass
        global individual_plots
        global mineral_formulas
        global export_name
        global plot_title
        global mineral
        global minerals
        global figure

        if input_selection == 'n':
            minerals = ['Akermanite', 'Al(OH)3(a)', 'Albite', 'Anhydrite', 'Anorthite', 'Anthophyllite', 'Antigorite', 'Aragonite', 'Arcanite', 'Artinite', 'Barite', 'Bischofite', 'Bloedite', 'Borax', 'Boric_acid,s', 'Brucite', 'Burkeite', 'Ca-Montmorillonite', 'Calcite', 'Carnallite', 'Celestite', 'Chalcedony', 'Chlorite(14A)', 'Chrysotile', 'Diopside', 'Dolomite', 'Enstatite', 'Epsomite', 'Fe(OH)3(a)', 'FeS(ppt)', 'Fluorite', 'Forsterite', 'Gaylussite', 'Gibbsite', 'Glaserite', 'Glauberite', 'Goergeyite', 'Goethite', 'Gypsum', 'Halite', 'Hausmannite', 'Hematite', 'Hexahydrite', 'Huntite', 'Hydroxyapatite', 'Illite', 'K-feldspar', 'K-mica', 'K2B4O7:4H2O', 'KB5O8:4H2O', 'Kainite', 'Kalicinite', 'Kaolinite', 'Kieserite', 'Labile_S', 'Leonhardite', 'Leonite', 'Mackinawite', 'Magnesite', 'Manganite', 'MgCl2_2H2O', 'MgCl2_4H2O', 'Mirabilite', 'Misenite', 'NaB5O8:5H2O', 'NaBO2:4H2O', 'Nahcolite', 'Natron', 'Nesquehonite', 'Pentahydrite', 'Pirssonite', 'Polyhalite', 'Portlandite', 'Pyrite', 'Pyrochroite', 'Pyrolusite', 'Quartz', 'Rhodochrosite', 'Schoenite', 'Sepiolite', 'Sepiolite(d)', 'SiO2(a)', 'Siderite', 'Strontianite', 'Sulfur', 'Sylvite', 'Syngenite', 'Talc', 'Teepleite', 'Thenardite', 'Trona', 'Vivianite', 'Witherite']

            mineral_formulas = ['Ca2Mg[Si2O7]', 'Al(OH)3', 'NaAlSi3O8', 'CaSO4', 'CaAl2Si2O8', '☐Mg2Mg5Si8O22(OH)2', 'Mg48Si34O85(OH)62', 'CaCO3', 'K2SO4', 'Mg2(CO3)(OH)2·3H2O', 'BaSO4', 'MgCl2·6H2O', 'Na2Mg(SO4)2·4H2O', 'Na2B4O5(OH)4•8(H2O)', 'B(OH)3', 'Mg(OH)2', 'Na6(CO3)(SO4)2', 'Ca0.165Al2.33Si3.67O10(OH)2', 'CaCO3', 'KMgCl3•6(H2O)', 'SrSO4', 'SiO2', 'Mg5Al2Si3O10(OH)8', 'Mg3Si2O5(OH)4', 'CaMgSi2O6', 'CaMg(CO3)2', 'MgSiO3', 'MgSO4•7(H2O)', 'Fe(OH)3', 'FeS', 'CaF2', 'Mg2SiO4', 'Na2Ca(CO3)2•5(H2O)', 'Al(OH)3', 'NaK3(SO4)2', 'Na2Ca(SO4)2', 'K2Ca5(SO4)6•(H2O)', 'FeO(OH)', 'CaSO4•2(H2O)', 'NaCl', 'Mn(II)Mn(III)2O4', 'Fe2O3', 'MgSO4•6(H2O)', 'CaMg3(CO3)4', 'Ca5(PO4)3(OH)', 'K0.6Mg0.25Al2.3Si3.5O10(OH)2', 'KAlSi3O8', 'KAl3Si3O10(OH)2', 'K2B4O7•4H2O', 'KB5O8•4H2O', 'MgSO4•KCl•3(H2O)', 'KHCO3', 'Al2Si2O5(OH)4', 'MgSO4•(H2O)', 'Na4Ca(SO4)3•2H2O', 'MgSO4•4(H2O)', 'K2Mg(SO4)2•4(H2O)', 'FeS', 'MgCO3', 'MnO(OH)', 'MgCl2:2H2O', 'MgCl2•4H2O', 'Na2SO4•10(H2O)', 'K8H6(SO4)7', 'NaB5O8•5H2O', 'NaBO2•4H2O', 'NaHCO3', 'Na2CO3•10(H2O)', 'Mg(HCO3)(OH)•2(H2O)', 'MgSO4•5(H2O)', 'Na2Ca(CO3)2•2(H2O)', 'K2Ca2Mg(SO4)4•2(H2O)', 'Ca(OH)2', 'FeS2', 'Mn(OH)2', 'MnO2', 'SiO2', 'MnCO3', 'K2Mg(SO4)2•6(H2O)', 'Mg4Si6O15(OH)2•6(H2O)', 'Mg4Si6O15(OH)2•6(H2O)', 'SiO2', 'FeCO3', 'SrCO3', 'S8', 'KCl', 'K2Ca(SO4)2•(H2O)', 'Mg3Si4O10(OH)2', 'Na2B(OH)4Cl', 'Na2SO4', 'Na3(CO3)(HCO3)•2(H2O)', 'Fe3(PO4)2•8(H2O)', 'BaCO3']      

        # the complete list of all minerals is created
        csv_minerals = []
        for column in csv_data.columns:
            if re.search('([A-Z].{3,})', column) and not re.search('(\(|\_|\:)', column):
                csv_minerals.append(column)

        if self.parameters['simulation_type'] == 'transport':
            csv_data.drop(csv_data.index[:3], inplace=True)

        plot_title = (input('''- What is the title of the plot? 
            Default = Scaling throughout the RO module  ___ ''')) or 'Scaling throughout the RO module'



        # all of the non-zero minerals are identified and the chemical formulas are sorted into a list
        non_zero_minerals = []
        for mineral in csv_minerals:
            for row in csv_data[mineral]:
                if row != 0 and mineral not in non_zero_minerals:
                    non_zero_minerals.append(mineral)

        quantity_nonzero_minerals = len(non_zero_minerals)      

        non_zero_mineral_formulas = []
        for mineral in non_zero_minerals:
            mineral_index = minerals.index(mineral)
            mineral_formula = mineral_formulas[mineral_index]
            non_zero_mineral_formulas.append(mineral_formula)      

        # plot the simulation depending upon the simulation perspective
        unit = 'moles'
        if simulation_perspective == "Brine":
            individual_plots = 'n'
            pyplot.figure(figsize = (17,10))
            pyplot.title(plot_title, fontsize = 'xx-large')
            pyplot.xlabel('Simulation time (s)', fontsize = 'x-large')
            pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  
            experimental_loop = []
            formula_index = 0
            for mineral in non_zero_minerals:
                mineral_serie = []
                time_serie = []
                for index, row in csv_data.iterrows():
                    mineral_serie.append(csv_data.at[index, mineral]) 
                    time = csv_data.at[index, 'time']
                    time_serie.append(time)

                pyplot.plot(time_serie,mineral_serie)
                pyplot.scatter(time_serie,mineral_serie)

                experimental_loop.append('%s [%s]' %(mineral,non_zero_mineral_formulas[formula_index]))   

            pyplot.legend(experimental_loop, loc='best', fontsize = 'x-large')
            pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
            figure = pyplot.gcf()
            print('\nClose the figure to proceed.\n')
            pyplot.show()

            export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file_name)
            if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
            else:
                export_name = export_filename_progenitor

            export_option = input('''- Would you like to export the figure? < y > or < n > __ ''')
            while export_option not in possible_options:
                print('''ERROR: The value is not one of the options.''')  
                export_option = input('''Would you like to export the figure? < y > or < n > __ ''')  

            if export_option == 'y':
                export_plot()

        elif simulation_perspective == 'Scaling':
            if quantity_nonzero_minerals < 2:
                individual_plots = 'n'
            elif quantity_nonzero_minerals >= 2:
                individual_plots = 'y'
            print('\nQuantity of precipitated minerals: %s' %(quantity_nonzero_minerals))
            #print('Number of timesteps per mineral: %s')
            individual_plots = input('''- Would you like to plot each mineral on a separate figure?
        Suggestion = < %s >  ;  < y > or < n >  __ '''  %(individual_plots)) or individual_plots
            while individual_plots not in possible_answers:
                print('ERROR: The entered value is not accepted.')  
                individual_plots = input('- Would you like to plot each mineral on a separate figure?')

            if individual_plots == 'y':
                viewing_plots = input('''How many mineral figures would you like to view?
        < All >, < A few >, or < None >''')
                while viewing_plots != 'All' and exporting_plots != 'A few' and exporting_plots != 'None':
                    print('ERROR: The entered value is not accepted.')  
                    viewing_plots = input('- How many mineral figures would you like to view?')

                if viewing_plots == ('All' or 'A few'):
                    viewing_minerals = []
                    if viewing_plots == 'All':
                        viewing_minerals = non_zero_minerals

                    elif exporting_plots == 'A few':
                        for mineral in non_zero_minerals:
                            print('< %s >' %(mineral))

                        viewing_mineral = input('''- Which minerals will you view?
            Type < done > when you are finished.''')
                        while viewing_mineral not in non_zero_minerals:
                            print('ERROR: The entered mineral is not among the options.')  
                            viewing_mineral = input('- Which minerals will you view?')
                        while viewing_mineral != 'done':
                            viewing_minerals.append(exporting_mineral)
                            viewing_mineral = input('- Which minerals will you view?')
                            while viewing_mineral not in non_zero_minerals:
                                print('ERROR: The entered mineral is not among the options.')  
                                viewing_mineral = input('- Which minerals will you view?')

                    for mineral in viewing_minerals:
                        pyplot.figure(figsize = (17,10))
                        pyplot.title(plot_title, fontsize = 'xx-large')

                        if self.parameters['simulation_type'] == 'transport':
                            pyplot.xlabel('Midpoint module distance (m)', fontsize = 'x-large')
                            pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  
                            experimental_loop = []

                            iteration = 0
                            distance_serie = []
                            time_serie = []
                            quantity_of_steps_index = 0   
                            for index, row in csv_data.iterrows():
                                if csv_data.at[index, 'time'] == 0:
                                    time_serie.append(csv_data.at[index, mineral]) 
                                    distance_serie.append(csv_data.at[index, 'dist_x'])
                                    quantity_of_steps_index += 1   
                                    time = 0

                                elif csv_data.at[index-1, 'soln'] == quantity_of_steps_index:
                                    experimental_loop.append('%s [%s] ; time (s): %s' 
                                                             %(mineral, 
                                                                mineral_formulas[minerals.index(mineral)], 
                                                                round(time, 2)))
                                    pyplot.plot(distance_serie,time_serie)
                                    distance_serie = []
                                    time_serie = []
                                    time_serie.append(csv_data.at[index, mineral])
                                    distance_serie.append(csv_data.at[index, 'dist_x'])
                                    time = csv_data.at[index, 'time']

                                elif index == len(csv_data[mineral]) + 2:   
                                    experimental_loop.append('%s [%s] ; time (s): %s' 
                                                             %(mineral, 
                                                                mineral_formulas[minerals.index(mineral)], 
                                                                round(time, 2)))
                                    pyplot.plot(distance_serie,time_serie)

                                else:
                                    time_serie.append(csv_data.at[index, mineral])
                                    distance_serie.append(csv_data.at[index, 'dist_x'])
                                    iteration += 1

                        elif self.parameters['simulation_type'] == 'evaporation':
                            pyplot.xlabel('Concentration Factor (CF)', fontsize = 'x-large')
                            pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  

                            experimental_loop = []
                            cf_series = []
                            concentration_series = []
                            data_length = len(csv_data['mass_H2O'])
                            for index, row in csv_data.iterrows():
                                if index < data_length:
                                    if csv_data.at[index, 'step'] >= 1:
                                        concentration_series.append(csv_data.at[index, mineral]) 
                                        solution_mass = csv_data.at[index, 'mass_H2O']
                                        cf_series.append(initial_solution_mass / solution_mass)  

                                    elif index > 1:
                                        print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')

                            concentration_series.append(csv_data.at[index, mineral]) 
                            solution_mass = csv_data.at[index, 'mass_H2O']
                            cf_series.append(initial_solution_mass / solution_mass)  

                            experimental_loop.append('%s [%s]' %(mineral, mineral_formulas[minerals.index(mineral)]))
                            pyplot.plot(cf_series,concentration_series)                    


                        pyplot.legend(experimental_loop, loc='best', fontsize = 'x-large')
                        pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
                        figure = pyplot.gcf()
                        print('\nClose the figure to proceed.')
                        pyplot.show()
                        export_figure = input('''Will you export the figure? < n > or < y >  __  ''')
                        while export_figure not in possible_answers:
                            print('ERROR: The entered mineral is not among the options.')  
                            export_figure = input('''\nWill you export the figure? < n > or < y >  __  ''')

                        # export the direct figures
                        if export_figure == 'y':
                            export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file)
                            if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                                export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
                            else:
                                export_name = export_filename_progenitor

                            export_plot()


            elif individual_plots == 'n':
                pyplot.figure(figsize = (17,10))
                pyplot.title(plot_title, fontsize = 'xx-large')
                pyplot.xlabel('Midpoint module distance (m)', fontsize = 'x-large')
                pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  
                experimental_loop = []
                for mineral in non_zero_minerals:
                    if self.parameters['simulation_type'] == 'transport':
                        iteration = 0
                        distance_serie = []
                        time_serie = []
                        quantity_of_steps_index = 0   
                        for index, row in csv_data.iterrows():
                            if csv_data.at[index, 'time'] == 0:
                                time_serie.append(csv_data.at[index, mineral]) 
                                distance_serie.append(csv_data.at[index, 'dist_x'])
                                quantity_of_steps_index += 1   
                                time = 0

                            elif csv_data.at[index-1, 'soln'] == quantity_of_steps_index:
                                experimental_loop.append('%s [%s] ; time (s): %s' %(mineral,                                                       mineral_formulas[minerals.index(mineral)], round(time, 2)))
                                pyplot.plot(distance_serie,time_serie)
                                distance_serie = []
                                time_serie = []
                                time_serie.append(csv_data.at[index, mineral])
                                distance_serie.append(csv_data.at[index, 'dist_x'])
                                time = csv_data.at[index, 'time']

                            elif index == len(csv_data[mineral]) + 2:   
                                experimental_loop.append('%s [%s] ; time (s): %s' %(mineral, mineral_formulas[minerals.index(mineral)], round(time, 2)))
                                pyplot.plot(distance_serie,time_serie)

                            else:
                                time_serie.append(csv_data.at[index, mineral])
                                distance_serie.append(csv_data.at[index, 'dist_x'])
                                iteration += 1


                    elif self.parameters['simulation_type'] == 'evaporation':
                        pyplot.xlabel('Concentration Factor (CF)', fontsize = 'x-large')
                        pyplot.ylabel('Quantity (%s)' %(unit), fontsize = 'x-large')  

                        experimental_loop = []
                        iteration = 0
                        mass_series = []
                        concentration_series = []
                        quantity_of_steps_index = 0  
                        initial_solution_mass = csv_data.at[0, 'mass_H2O']
                        for index, row in csv_data.iterrows():
                            try:
                                if csv_data.at[index+1, 'mass_H2O']:
                                    if csv_data.at[index, 'step'] >= 1:
                                        concentration_series.append(csv_data.at[index, mineral]) 
                                        solution_mass = csv_data.at[index, 'mass_H2O']
                                        mass_series.append(initial_solution_mass / solution_mass)   

                                    else:
                                        print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')


                            except:
                                concentration_series.append(csv_data.at[index, mineral]) 
                                solution_mass = csv_data.at[index, 'mass_H2O']
                                mass_series.append(initial_solution_mass / solution_mass)  

                                experimental_loop.append('%s [%s]' %(mineral, mineral_formulas[minerals.index(mineral)]))
                                pyplot.plot(mass_series,concentration_series)

                                mass_series = []
                                concentration_series = []
                                concentration_series.append(csv_data.at[index, mineral]) 
                                solution_mass = csv_data.at[index, 'mass_H2O']
                                mass_series.append(initial_solution_mass / solution_mass)         


                pyplot.legend(experimental_loop, loc='best', fontsize = 'x-large')
                pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(simulation_cf), wrap=True, horizontalalignment='left', fontsize=12)
                figure = pyplot.gcf()
                print('\nClose the figure to proceed.')
                pyplot.show()
                export_figure = input('''Will you export the figure? < n > or < y >  __  ''')
                while export_figure not in possible_answers:
                    print('ERROR: The entered mineral is not among the options.')  
                    export_figure = input('''\nWill you export the figure? < n > or < y >  __  ''')

                # export the output graphic
                if export_figure == 'y':
                    export_filename_progenitor = re.sub('(\.\w+)', '', selected_output_file)
                    if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
                        export_name = '{}, {}'.format(export_filename_progenitor, simulation_perspective)
                    else:
                        export_name = export_filename_progenitor

                    export_option = input('''Would you like to export the figure? < y > or < n > __ ''')
                    while export_option not in possible_answers:
                        print('''ERROR: The value is not one of the options.''')  
                        export_option = input('''Would you like to export the figure? < y > or < n > __ ''')       

                    if export_option == 'y':
                        export_plot()


    def export_plot():
        """
        Export the plots to the current working directory  
        """

        if graphical_selection == 'Brine' or (graphical_selection == 'Scaling' and individual_plots == 'n'):
            export_file_name = input('''- What is the name of your export figure?
            Omit < . > and < \ > in the name.
            Default = %s''' %(export_name)) or export_name
            while re.search('(\.|\\))', export_name):
                print('''ERROR: Remove < . > and < \ > from the figure name.''')  
                export_file_name = input('- What will be the name of your export figure?') 
        else:
            export_file_name = input('''- What is the name of your export figure?
            Omit < . > and < \ > in the name.
            Default = %s, %s''' %(export_name, mineral)) or '%s, %s' %(export_name, mineral)
            while re.search('(\.|\\))',export_name):
                print('''ERROR: Remove < . > and < \ > from the figure name.''')  
                export_file_name = input('- What will be the name of your export figure?')

        available_formats = ['jpg', 'png', 'svg']
        export_format = input('''- What will be the format of your export figure?
        Select from < %s >, < %s >, and < %s >.
        Default = < jpg >''' %('jpg', 'png', 'svg')) or 'jpg'
        while export_format not in available_formats:
            print('''ERROR: Select from < %s >, < %s >, and < %s >.''' %('jpg', 'png', 'svg'))  
            export_format = input('- What will be the format of your export figure?')       

        file_number = 0
        if not os.path.exists('%s.%s' %(export_file_name, export_format)):
            figure.savefig('%s.%s' %(export_file_name, export_format))
        elif os.path.exists('%s.%s' %(export_file_name, export_format)):
            while os.path.exists('%s_%s.%s' %(export_file_name, file_number, export_format)):
                file_number += 1
            figure.savefig('%s_%s.%s' %(export_file_name, file_number, export_format))


    def choose_your_path():
        """
        The software user directs the code functions based upon the use case  
        """
        global visualize_selection
        global input_selection
        global execute_selection

        welcome()

        input_selection = input('''- Will you create an input file? < y > or < n > ___ ''')
        while input_selection not in possible_answers:
            print('ERROR: provide < y > or < n > as your answer.')
            input_selection = input('''- Will you create an input file? < y > or < n > ___ ''')

        execute_selection = input('''- Will you execute a simulation?
        < y > or < n >''')
        while execute_selection not in possible_answers:
            print('ERROR: provide < y > or < n > as your answer.')
            execute_selection = input('''- Will you execute a simulation? < y > or < n > ___ ''')

        visualize_selection = input('''- Will you visualize simulation results? < y > or < n > ___ ''')
        while visualize_selection not in possible_answers:
            print('ERROR: provide < y > or < n > as your answer.')
            visualize_selection = input('''- Will you visualize simulation results? < y > or < n > ___ ''')

        while input_selection == execute_selection == visualize_selection == 'n':
            print('ERROR: one of the software options equate < y >.')
            simulation_purpose = input('''- Will you create an input file, execute a simulation, or visualize simulation results?
            < input >, < execute >, or < visualize > __ ''')
            while simulation_purpose not in ['input', 'execute', 'visualize']:
                print('ERROR: provide the purpose of the simulation, or exit the software.')
                simulation_purpose = input('''- Will you create an input file, execute a simulation, or visualize simulation results?
                < input >, < execute >, or < visualize > __ ''')

            if simulation_purpose == 'input':
                input_selection = 'y'

            elif simulation_purpose == 'execute':
                execute_selection = 'y'

            elif simulation_purpose == ' visualize':
                visualize_selection = 'y'

        if visualize_selection == input_selection == execute_selection == 'y':  
            print('\nSimulate and visualize a created input\n', '+'*len('Simulate and visualize a created input'))
            make_general_conditions()
            make_reactive_transport()
            make_solutions()
            make_equilibrium_phases()
            make_selected_output()
            export()
            execute()
            process_selected_output()

        elif execute_selection == visualize_selection == 'y' and input_selection == 'n':
            print('\nSimulation and visualization\n', '+'*len('Simulation and visualization'))
            execute()
            process_selected_output()

        elif input_selection == execute_selection == 'y' and visualize_selection == 'n':
            print('\nSimulate a created input\n', '+'*len('Simulate a created input'))
            make_general_conditions()
            make_reactive_transport()
            make_solutions()
            make_equilibrium_phases()
            make_selected_output()
            export()
            execute()
            process_selected_output()

        elif visualize_selection == 'y' and input_selection == execute_selection == 'n':
            print('\nVisualize\n', '+'*len('Visualize'))
            process_selected_output()

        elif input_selection == 'y' and visualize_selection == execute_selection == 'n':
            print('\nInput creation\n', '+'*len('Input creation'))
            make_general_conditions()
            make_reactive_transport()
            make_solutions()
            make_equilibrium_phases()
            make_selected_output()
            export()

        elif execute_selection == 'y' and input_selection == visualize_selection == 'n':
            print('\nSimulation\n', '+'*len('Simulation'))
            execute()

        final_message = '\n\nThe simulation is complete.'
        print('%s\n%s' %(final_message, '='*len(final_message)))
        
        
        
    def define(self):
        pass
    
    
    def calculate(self):
        if self.parameters['simulation_type'] == 'transport':
            transport()
    
    def execute(self):
        pass
    
    def simulate(self):
        pass