# import libraries
from matplotlib import pyplot 
from scipy.constants import nano, kilo, milli, liter, day
from itertools import chain
from chempy.properties.water_density_tanaka_2001 import water_density
from pubchempy import get_compounds 
import subprocess
import datetime
import pandas
from math import pi, exp
from glob import glob
import time
import json
import os
import re

# calculation constants
simulated_time_over_computational_time = 9.29    


class ROSSPkg():
    def __init__(self):       
        # establish the general organization structures
        self.parameters = {}
        self.variables = {}
        self.results = {}
        self.results['figures'] = {}

        
    def define_general(self, operating_system, phreeqc_path, database_selection, simulation_type, simulation_title):
        '''
        Establish general conditions
        '''
        self.parameters['water_mw'] = float(get_compounds('water', 'name')[0].molecular_weight)
        self.parameters['water_density'] = water_density()
        
        # parameterize the input file
        self.parameters['os'] =  operating_system
        self.parameters['phreeqc_path'] = phreeqc_path
        self.parameters['simulation_type'] = simulation_type
        database_path = f'databases/{database_selection}.json'
        
        title_line = 'TITLE\t %s' %(simulation_title)
        if operating_system == 'Windows':
            database_line = 'DATABASE %s' %(database_path)
            self.results['general_conditions'] = [database_line, title_line]
        else:
            self.results['general_conditions'] = [title_line]
            
        # establish the database content
        self.parameters['database_selection'] = database_selection 
        database = json.load(open(database_path, 'r'))
        self.parameters['elements'] = database['elements']
        self.parameters['minerals'] = database['minerals']

    def transport(self, module_selection = 'BW30-400', module_characteristics = {}, quantity_of_modules = 1, cells_per_module = 12, domain = 'dual', output_perspective = 'scaling'):
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
            membrane_thickness = 250 * (nano / milli)   #mm
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
        module_cross_sectional_area = module_diameter**2 * pi / 4        #squared millimeters
        permeate_tube_cross_sectional_area = permeate_tube_diameter**2 * pi / 4     #squared millimeters
#             filtering_layer_thicknes = (module_diameter - permeate_thickness) / 2         #millimeters
        filtration_cross_sectional_area = module_cross_sectional_area - permeate_tube_cross_sectional_area        #squared millimeters
        feed_cross_sectional_area = (feed_thickness / repeated_membrane_thickness) * filtration_cross_sectional_area     #squared millimeters
        feed_volume = feed_cross_sectional_area * module_length * kilo**2   #cubic meters
        feed_mass = feed_volume * self.parameters['water_density'] * (1 / liter) / kilo    #kilograms, which assumes pure water for mass
        feed_moles = feed_mass * kilo / self.parameters['water_mw'] 

        # calculate fluid flow characteristics
        kinematic_flow_velocity = 9.33E-7    #square meters / second
        feed_velocity = max_feed_flow / (feed_cross_sectional_area * kilo**2) / day     #meters / second
        reynolds_number = feed_velocity * module_length / kinematic_flow_velocity
        print('Reynold\'s number: ', reynolds_number)

        # calculate module cell characteristics
        self.parameters['feed_mass_cell'] = feed_mass / self.parameters['cells_per_module']      #kg
        self.parameters['feed_moles_cell'] = feed_moles / self.parameters['cells_per_module']    #moles

        # calculate simulation parameters that will populate the PHREEQC simulation 
        timesteps_over_cell = 1
        print('cell length: ', cell_length)
        print('feed velocity: ', feed_velocity)
        print('max_feed_flow: ', max_feed_flow)        
        print('feed_cross_sectional_area: ', feed_cross_sectional_area)
        
        
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
        if self.parameters['output_perspective'] == 'scaling':
            punch_cells_line = '-punch_cells\t\t1-%s' %(self.parameters['cells_per_module'])
            punch_frequency_line = '-punch_frequency\t%s' %(self.parameters['cells_per_module'])
        elif self.parameters['output_perspective'] == 'brine':
            punch_cells_line = '-punch_cells\t\t%s' %(self.parameters['cells_per_module'])
            punch_frequency_line = '-punch_frequency\t1'       

        # create the transport block
        self.results['transport_block'] = []
        self.results['transport_block'].extend((transport_line, cells_line, shifts_line, lengths_line, timestep_line, initial_time_line, boundary_conditions_line, domain_line, punch_cells_line, punch_frequency_line))


    def reaction(self, permeate_approach = 'linear permeate', permeate_efficiency = 1, head_loss = -0.15, final_cf = 2):
        '''
        Define the REACTION block
        '''
        self.parameters['permeate_approach'] = permeate_approach
        cfs = []
        cell_moles = []
        reaction_parameters = []
        iteration = 0
        cumulative_cf = 1 
        self.results['reaction_block'] = []
        for module in range(quantity_of_modules):
            module_previous_moles_removed = 0
            if permeate_approach == 'linear permeate':
                if iteration == 0:
                    initial_moles_cell = self.parameters['feed_moles_cell']

                initial_moles_removed = self.parameters['permeate_removal_per_cell'] * 2 / (1 + exp(head_loss))
                final_moles_removed = initial_moles_removed * exp(head_loss)
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

            if permeate_approach == 'linear CF':
                module_iteration = 0
                initial_cf = 1

                cf_slope = (final_cf - initial_cf) / self.parameters['cells_per_module']
                for cell in range(self.parameters['cells_per_module']):
                    cell_cf = (cell+1) * cf_slope + initial_cf
                    cfs.append(cell_cf)    

                for cf in cfs:
                    moles_to_be_removed =  self.parameters['feed_moles_cell'] - (self.parameters['feed_moles_cell'] / cf)
                    if module_iteration == 0:
                        initial_moles_cell = self.parameters['feed_moles_cell']
                        reaction_parameters.append(moles_to_be_removed)
                    if module_iteration > 0:
                        module_previous_moles_removed += reaction_parameters[-1] 
                        reaction_parameter = moles_to_be_removed - module_previous_moles_removed
                        reaction_parameters.append(reaction_parameter)
                        initial_moles_cell -= moles_to_be_removed

                    module_iteration += 1

                cf = cfs[-1]
                cumulative_cf *= cf
                initial_moles_cell = self.parameters['feed_moles_cell'] - moles_to_be_removed   # moles_to_be_removed = the final quantity of moles is calculated with naming from the linear permeate flux method to be consistent  

            if self.parameters['simulation_type'] == 'transport':
                self.results['reaction_block'].append('\n')
                for cell in range(self.parameters['cells_per_module']):
                    cell_number = (cell) + self.parameters['cells_per_module'] * module
                    reaction_line = 'REACTION %s' %(cell_number)
                    
                    if (cell+1) < self.parameters['cells_per_module']:
                        reaction_line += f'\n\tH2O -1; {reaction_parameters[cell_number]}' 
                    elif (cell+1) == self.parameters['cells_per_module']:
                        reaction_line += f'''\n\tH2O -1; {reaction_parameters[cell_number]}
        INCREMENTAL_REACTIONS \ttrue'''     

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
            final_solution_mass = initial_moles_cell * self.parameters['water_mw'] / kilo  #kg water mass
            final_cf_cell = self.parameters['feed_mass_cell'] / final_solution_mass
            print('Effluent module %s CF:' %(module + 1), final_cf_cell)

            if self.parameters['os'] == 'windows':
                self.results['reaction_block'].append('#%s' %(permeate_approach))
                if permeate_approach == 'linear permeate':
                    self.results['reaction_block'].append('''
        #Permeate efficiency parameter: %s
        #Head loss parameter: %s''' %(permeate_efficiency, head_loss))
                self.results['reaction_block'].append('''    #Effluent module %s:
        #Estimated CF: %s
        #Estimated solution mass: %s\n\n''' %(module + 1, cumulative_cf, final_solution_mass))


    def solutions(self, water_selection, custom_water_parameters = {}, solution_description = ''):
        """
        Specify the SOLUTION block of the simulation.
        """
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
            water_body = json.load(open(water_file_path))
            
            elements_lines = []
            for content, information in water_body.items():
                if content == 'element':
                    for element, information2 in information.items():
                        if element in self.parameters['elements']:
                            conc = information2['concentration (ppm)']
                            ref = information2['reference']
                            if len(str(conc)) <= 3:
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
                        if element in self.parameters['elements']:                          
                            if len(str(conc)) <= 3:
                                elements_lines[element] = (f'{element}\t\t{conc}\t#{ref}')
                            else:
                                elements_lines[element] = (f'{element}\t{conc}\t#{ref}')
                        else:
                            print('\n--> ERROR: The {} element is not accepted by the {} database'.format(element, self.parameters['database_selection']))
                                    
                            
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

            for element in self.parameters['elements']:
                element_concentration = 0
                element_line = f'{element}\t{element_concentration}'    
                self.results['solution_block'].append(element_line)

            water_line = '-water \t %s' %(self.parameters['feed_mass_cell'])
            self.results['solution_block'].append(water_line)


    def equilibrium_phases(self, block_comment = '', ignored_minerals = [], existing_parameters = {}, verbose = True):
        """
        Specify the EQUILIBRIUM_PHASES block of the simulation.
        """
        # define mineral sizes for later spacing
        short_mineral_names = ['Barite', 'Gypsum','Halite', 'Natron', 'Quartz', 'Talc','Trona', 'Borax', 'Albite', 'K-mica','Illite', 'Pyrite', 'Sulfur',]
        long_mineral_names = ['Anthophyllite', 'Hexahydrite', 'Leonhardite', 'Nesquehonite', 'Pentahydrite', 'Portlandite','Sepiolite(d)', 'Boric_acid,s', 'K2B4O7:4H2O', 'NaB5O8:5H2O', 'Rhodochrosite', 'Strontianite','Hydroxyapatite', 'Chlorite(14A)', 'Mackinawite', 'Hausmannite', 'Pyrochroite']

        # determine the set of possible minerals 
        elements_list = list(self.parameters['elements'])
        compounds_list = [self.parameters['elements'][x]['gfw_formula'] for x in elements_list]
        total_accepted_characters = list(chain(compounds_list, elements_list))

        self.variables['described_minerals'] = {}
        for mineral in self.parameters['minerals']:
            mineral_formula = self.parameters['minerals'][mineral]
            original_formula = mineral_formula
            
            # remove entities is an ordered fashion
            mineral_formula = re.sub('(SO4|H2O|OH|CO3)', '', mineral_formula)
            mineral_formula = re.sub('([0-9]|\(|\)|•|:|\.)', '', mineral_formula)
            mineral_elements = re.findall('[A-Z][a-z]?', mineral_formula)
            
            for element in mineral_elements:
                if element in elements_list:
                    mineral_formula = re.sub(element, '', mineral_formula)
            if mineral_formula == '':
                self.variables['described_minerals'][mineral] = original_formula
            else:
                if verbose:
                    print(f'--> < {mineral_formula} > characters of < {mineral} >, formula = < {original_formula} > , are undescribed by the parameterized elements.')
                    
        print(self.variables['described_minerals'])


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


    def selected_output(self, output_filename = None):
        '''
        Specify the output file after a PHREEQC simulation
        '''
        # create parameter lines 
        if output_filename is None:
            count = 0
            selected_output_file_name = '_'.join([str(x) for x in [datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['output_perspective'], count]]) 
            while os.path.exists(f'{selected_output_file_name}.txt'):
                count += 1
                selected_output_file_name = '_'.join([str(x) for x in [datetime.date.today(), self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'], self.parameters['output_perspective'], count]]) 
        else:
            selected_output_file_name = output_filename

        selected_output_file_name += '.txt'
        self.parameters['selected_output_file_name'] = selected_output_file_name

        minerals_line = ''
        for mineral in self.parameters['minerals']:
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
        self.results['complete_lines'] = chain(self.results['general_conditions'], self.results['solution_block'], self.results['equilibrium_phases_block'], self.results['reaction_block'], self.results['selected_output_block'], self.results['transport_block']) 

        # create the export input file
        file_number = 0
        if self.parameters['permeate_approach'] == 'linear permeate':
            permeate_approach_name = 'LinPerm'
        elif self.parameters['permeate_approach'] == 'linear CF':
            permeate_approach_name = 'LinCF'
        
        self.parameters['input_file_name'] = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'],self.parameters['output_perspective'], permeate_approach_name, file_number]]) + '.pqi'
        while os.path.exists(self.parameters['input_file_name']):
            file_number += 1
            self.parameters['input_file_name'] = '_'.join([str(x) for x in [datetime.date.today(), 'ROSS', self.parameters['water_selection'], self.parameters['simulation_type'], self.parameters['database_selection'],self.parameters['output_perspective'], permeate_approach_name, file_number]]) + '.pqi'

        # printing and exporting the input file
        with open(self.parameters['input_file_name'],'w') as input_file:
            for line in self.results['complete_lines']:
                if print_block:
                    print(line)
                input_file.write(line + '\n')

    def execute(self, input_path = None, output_path = None):
        '''
        Execute the input file through the batch PHREEQC software 
        '''
        # define the file paths for the simulation input and output
        if input_path is None:
            working_directory = os.getcwd()
            input_path = os.path.join(working_directory, self.parameters['input_file_name'])
        if output_path is None:
            output_file_name = re.sub('(?<=\.)(.+)' , 'pqo', self.parameters['input_file_name'])        
            output_path = os.path.join(working_directory, output_file_name)

        # locate the batch PHREEQC software
        bat_path = os.path.join(self.parameters['phreeqc_path'], 'bin\\phreeqc.bat')
        database_path = os.path.join(self.parameters['phreeqc_path'], 'database\\{}.dat'.format(self.parameters['database_selection']))

        # execute command for the PHREEQC calculations
        proc = subprocess.Popen('cmd.exe', stdin=subprocess.PIPE)
        command = str.encode("\"" + bat_path + "\" \"" + input_path + "\" \"" + output_path + "\"  \"" + database_path + "\"\n") 
        proc.stdin.write(command)
        proc.stdin.close()  
        proc.wait()

        # verify that the PHREEQC executed and generated the appropriate files
        if os.path.exists(output_path):
            if os.path.exists(self.parameters['selected_output_file_name']):
                print('The execution is complete.')
            else:
                print('ERROR: The SELECTED_OUTPUT file is missing. The simulation failed to simulate.')
        else:
            print('\nERROR: The simulation failed to execute.')


    def process_selected_output(self, selected_output_path = None, graphical_selection = 'scaling', plot_title = '', title_font = 'xx-large', label_font = 'x-large', plot_caption = '', table_title = None, export_figure = True, export_format = 'svg', individual_plots = True):
        """
        Interpreting the PHREEQC SELECTED_OUTPUT file and conducting the plotting functions
        """
        self.parameters['graphical_selection'] = graphical_selection
        
        databases = [database for database in glob('./databases/*.json')]
        databases = [re.search('(\w+)(?=.json)', database).group() for database in databases]

        # determining the appropriate variables
        if selected_output_path is not None:
            # define the simulation perspective
            self.parameters['output_perspective'] = 'scaling'
            if re.search('(brine)', selected_output_path, flags=re.IGNORECASE):
                self.parameters['output_perspective'] = 'brine'
                      
            # define the simulation type
            self.parameters['simulation_type'] = 'transport'
            for line in selected_output_path:
                if re.search('(evaporation)', line, re.IGNORECASE):
                    self.parameters['simulation_type'] = 'evaporation'
            
            # define the database contents
            self.parameters['database_selection'] = 'llnl'
            for database in databases:
                names = database.split('_')
                print(names)
                print(selected_output_path)
                if all(re.search(name, selected_output_path, re.IGNORECASE) for name in names):
                    self.parameters['database_selection'] = database                   
                    break                    
                    
            database_json = json.load(open('./databases/{}.json'.format(self.parameters['database_selection'])))
            self.parameters['minerals'] = database_json['minerals']
            self.parameters['elements'] = database_json['elements']
            
            # define the simulation
            self.parameters['selected_output_file_name'] = re.search('(\w+)(?=\.)', selected_output_path).group()
                    
        else:
            working_directory = os.getcwd()
            selected_output_path = os.path.join(working_directory, self.parameters['selected_output_file_name'])                     

        # preparing the SELECTED_OUTPUT file into a dataframe
        selected_output = open(selected_output_path, 'r')
        original_data = pandas.read_table(selected_output, sep = '\t')
        self.results['csv_data'] = pandas.DataFrame(original_data)
        for column in self.results['csv_data'].columns:
            new_column = column.strip()
            self.results['csv_data'].rename(columns={column:new_column}, inplace = True)

        self.variables['initial_solution_mass'] = self.results['csv_data'].at[0, 'mass_H2O']
        final_solution_mass = self.results['csv_data']['mass_H2O'].iloc[-1]
        self.variables['simulation_cf'] = self.variables['initial_solution_mass'] / final_solution_mass
                                 
        # conducting the appropriate visualization function
        pyplot = self.plot_definition(plot_title, title_font, label_font,)
        if self.parameters['graphical_selection'] == 'brine':
            self.brine_plot(pyplot, plot_title, title_font, label_font, plot_caption, table_title, export_figure, export_format)
        elif self.parameters['graphical_selection'] == 'scaling':
            self.scaling_plot(pyplot, plot_title, title_font, label_font, plot_caption, table_title, individual_plots, export_figure, export_format)

                                 
    def plot_definition(self, plot_title, title_font, label_font):
        pyplot.figure(figsize = (17,10))
        pyplot.title(plot_title, fontsize = title_font)                
        pyplot.grid(True)
                                 
        if (self.parameters['graphical_selection'] or self.parameters['output_perspective']) == 'brine':
            if plot_title == '':
                 plot_title = 'Effluent brine elemental concentrations'
            pyplot.figure(figsize = (17,10))
            pyplot.xlabel('Time (s)', fontsize = label_font)
            pyplot.ylabel('Concentration (ppm)', fontsize = label_font)
            pyplot.yscale('log')
            pyplot.legend(elements, loc='best', title = 'non-zero elements', fontsize = 'x-large')
        else:
            if plot_title == '':
                 plot_title = 'scaling throughout the RO module'
            pyplot.xlabel('Midpoint module distance (m)', fontsize = label_font)
            pyplot.ylabel('Quantity (moles)', fontsize = label_font)  
                                 
        return pyplot
                                 
    def brine_plot(self, pyplot, plot_title, title_font, label_font, plot_caption, table_title, export_figure, export_format):
        """
        Generate plots of the elemental concentrations from effluent brine in the PHREEQC SELECTED_OUTPUT file  
        """
        elements = []
        for column in self.results['csv_data'].columns:
            if re.search('([A-Z][a-z]?(\(\d\))?){1}$', column) and not re.search('(_|H2O|pH)', column):
                elements.append(column)

        self.results['csv_data'].drop(self.results['csv_data'].index[:3], inplace=True)

        # plot the brine concentrations figure
        unit = 'mol/kgw'
        concentration_array = []

        for element in elements:  
            concentration_serie = []
            time_serie = []
            initial_solution_time = 0
            for index, row in self.results['csv_data'].iterrows():
                if self.results['csv_data'].at[index, 'Cl'] == 0:
                    initial_solution_time += 1
                else:
                    concentration_serie.append(self.results['csv_data'].at[index, element])
                    delta_time = self.results['csv_data'].at[index, 'time'] - self.results['csv_data'].at[index-1, 'time']
                    time_serie.append(self.results['csv_data'].at[index, 'time'] - initial_solution_time * delta_time)

            pyplot.plot(time_serie,concentration_serie)

        pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(self.variables['simulation_cf']), wrap=True, horizontalalignment='left', fontsize=12)
        figure = pyplot.gcf()
        if plot_caption == '':
            plot_caption = '''\n\nbrine Figure:\n%s 
                    The effluent concentrations of each existing element in the brine. brine plots from brine data rapidly reach a steady state elemental concentration. brine plots from scaling data consist of vertical concentrations that represent the concentrations at each distance throughout the RO module at the specified time, where the low end represents the influent concentration while the high end represents the effluent concentration.''' %('='*len('brine Figure'))
                                 
        print(plot_caption)
        pyplot.show()

        # create a complementary concentration data table
        loop_iteration = 1
        table_view = {}
        average_concentrations_table = pandas.DataFrame()
        index_elements = []
        if self.parameters['output_perspective'] == 'scaling':
            for element in self.parameters['elements']:
                quantity_of_steps_index = 0
                average_iteration = 0
                time_serie = []            
                time_averages = {}
                for index, row in self.results['csv_data'].iterrows():
                    if self.results['csv_data'].at[index, 'time'] == 0:
                        time_serie.append(self.results['csv_data'].at[index,element])                 
                        quantity_of_steps_index += 1                    

                    elif self.results['csv_data'].at[index-1,'soln'] == quantity_of_steps_index:       
                        #process the complete time serie
                        try:
                            average_concentration = sum(time_serie) / len(time_serie)
                        except:
                            average_concentration = 0

                        #print(average_concentration)
                        table_view['Time (s): %s' %(average_iteration * quantity_of_steps_index)] = average_concentration
                        average_iteration += 1 

                        #begin the new time serie
                        time_serie = []
                        time_averages = {}
                        time_serie.append(self.results['csv_data'].at[index,element])  

                    elif index == len(self.results['csv_data'][element]) + 2:       
                        time_serie.append(self.results['csv_data'].at[index,element])            
                        try:
                            average_concentration = sum(time_serie) / len(time_serie)
                        except:
                            average_concentration = 0
                        table_view['Time (s): %s' %(round(average_iteration * quantity_of_steps_index), 1)] = average_concentration
                        average_iteration += 1
                        index_elements.append(element)
                        average_concentrations_table = average_concentrations_table.append(table_view, ignore_index=True)
                        loop_iteration += 1

                    else:
                        time_serie.append(self.results['csv_data'].at[index,element])

            # defining the index column of the DataFrame
            average_concentrations_table.index = index_elements
            average_concentrations_table.index.name = 'Elements'
            if table_title is None:
                table_title = 'Average elemental molal concentrations of the feed water in the RO module for each %s seconds of simulation:' %(quantity_of_steps_index)

        # create a complementary concentration data table for the brine figure 
        elif self.parameters['output_perspective'] == 'brine':
            total_time = self.results['csv_data']['time'].iloc[-1]
            for element in elements:  
                concentration_serie = []
                time_serie = []
                for index, row in self.results['csv_data'].iterrows():
                    if self.results['csv_data'].at[index, 'Cl'] != 0:
                        concentration_serie.append(self.results['csv_data'].at[index,element])      
                average_concentration = sum(concentration_serie) / len(concentration_serie)
                table_view['%s' %(element)] = average_concentration

            average_concentrations_table = average_concentrations_table.append(table_view, ignore_index=True)
            average_concentrations_table.rename(index = {0:'Concentrations (molal)'}, inplace = True)
            if table_title is None:
                table_title = 'Average elemental molal concentrations of the feed water in the RO module over %s seconds of simulation:' %(total_time)

        print('\n\n\n',table_title,'\n%s'%('='*len(table_title)))
        print(average_concentrations_table)

        # export the graphic
        if export_figure:
            self.export_plot(figure, plot_title, export_format)
                                 

    def scaling_plot(self, pyplot, plot_title, title_font, label_font, plot_caption, table_title, individual_plots, export_figure, export_format = 'svg'):
        """
        Generate plots of scaling along the module distance in the PHREEQC SELECTED_OUTPUT file  
        """
        # define reused functions
        def series_creation():
            iteration = 0
            distance_serie = []
            time_serie = []
            quantity_of_steps_index = 0   
            for index, row in self.results['csv_data'].iterrows():
                if self.results['csv_data'].at[index, 'time'] == 0:
                    time_serie.append(self.results['csv_data'].at[index, mineral]) 
                    distance_serie.append(self.results['csv_data'].at[index, 'dist_x'])
                    quantity_of_steps_index += 1   
                    time = 0

                elif self.results['csv_data'].at[index-1, 'soln'] == quantity_of_steps_index:
                    experimental_loop.append(f'{mineral} [{mineral_formula}] ; time: {round(time)} s')
                    pyplot.plot(distance_serie,time_serie)
                    distance_serie = []
                    time_serie = []
                    time_serie.append(self.results['csv_data'].at[index, mineral])
                    distance_serie.append(self.results['csv_data'].at[index, 'dist_x'])
                    time = self.results['csv_data'].at[index, 'time']

                elif index == len(self.results['csv_data'][mineral]) + 2:   
                    experimental_loop.append(f'{mineral} [{mineral_formula}] ; time: {round(time)} s')
                    pyplot.plot(distance_serie,time_serie)

                else:
                    time_serie.append(self.results['csv_data'].at[index, mineral])
                    distance_serie.append(self.results['csv_data'].at[index, 'dist_x'])
                    iteration += 1
                                 
            return experimental_loop     
                                 
        def illustrate(pyplot):                     
            pyplot.legend(experimental_loop, loc='best', fontsize = 'x-large')
            pyplot.figtext(0.2, 0, 'Desalination CF: %s' %(self.variables['simulation_cf']), wrap=True, horizontalalignment='left', fontsize=12)
            figure = pyplot.gcf()
            pyplot.show()
                        
            return figure
                                 
        # the complete list of all minerals is created
        csv_minerals = []
        for column in self.results['csv_data'].columns:
            if re.search('([A-Z].{3,})', column) and not re.search('(\(|\_|\:)', column):
                csv_minerals.append(column)

        if self.parameters['simulation_type'] == 'transport':
            self.results['csv_data'].drop(self.results['csv_data'].index[:3], inplace=True)

        # all of the non-zero minerals are identified and the chemical formulas are sorted into a list
        non_zero_minerals = set()
        for mineral in csv_minerals:
            for value in self.results['csv_data'][mineral]:
                if value != 0:
                    non_zero_minerals.add(mineral)

        quantity_nonzero_minerals = len(non_zero_minerals)         

        # plot the simulation depending upon the simulation perspective
        unit = 'moles'
        if self.parameters['output_perspective'] == "brine":
            experimental_loop = []
            formula_index = 0
            for mineral in non_zero_minerals:
                mineral_formula = self.parameters['minerals'][mineral]
                                 
                mineral_serie = []
                time_serie = []
                for index, row in self.results['csv_data'].iterrows():
                    mineral_serie.append(self.results['csv_data'].at[index, mineral]) 
                    time = self.results['csv_data'].at[index, 'time']
                    time_serie.append(time)

                pyplot.plot(time_serie,mineral_serie)
                pyplot.scatter(time_serie,mineral_serie)

                experimental_loop.append(f'{mineral} [{mineral_formula}]')   

            # export the figure
            figure = illustrate(pyplot)
            if export_figure:
                self.results['figures'][self.parameters['selected_output_file_name']] = {'figure':figure, 'title':plot_title}
                self.export_plot(figure, plot_title, export_format)

        elif self.parameters['output_perspective'] == 'scaling':
            if quantity_nonzero_minerals < 2:
                individual_plots = 'n'
            elif quantity_nonzero_minerals >= 2:
                individual_plots = 'y'
            print(f'\nQuantity of precipitated minerals: {quantity_nonzero_minerals}')

            if individual_plots:
                for mineral in non_zero_minerals:
                    print(mineral)
                    mineral_formula = self.parameters['minerals'][mineral]
                                 
                    pyplot.figure(figsize = (17,10))
                    pyplot.title(plot_title, fontsize = title_font)

                    if self.parameters['simulation_type'] == 'transport':
                        pyplot.xlabel('Midpoint module distance (m)', fontsize = label_font)
                        pyplot.ylabel('Quantity (ppm)', fontsize = label_font)  
                        experimental_loop = []

                        experimental_loop = series_creation()
                        

                    elif self.parameters['simulation_type'] == 'evaporation':
                        pyplot.xlabel('Concentration Factor (CF)', fontsize = label_font)
                        pyplot.ylabel('Quantity (%s)' %(unit), fontsize = label_font)  

                        experimental_loop = []
                        cf_series = []
                        concentration_series = []
                        data_length = len(self.results['csv_data']['mass_H2O'])
                        for index, row in self.results['csv_data'].iterrows():
                            if index < data_length:
                                if self.results['csv_data'].at[index, 'step'] >= 1:
                                    concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                    solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                    cf_series.append(self.variables['initial_solution_mass'] / solution_mass)  
                                elif index > 1:
                                    print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')

                        concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                        solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                        cf_series.append(self.variables['initial_solution_mass'] / solution_mass)  

                        experimental_loop.append(f'{mineral} [{mineral_formula}]')
                        pyplot.plot(cf_series,concentration_series)                    

                    # export the figure
                    figure = illustrate(pyplot)
                    if export_figure:
                        self.export_plot(figure, plot_title, mineral, export_format)


            elif not individual_plots:
                experimental_loop = []
                for mineral in non_zero_minerals:
                    if self.parameters['simulation_type'] == 'transport':
                        experimental_loop = series_creation()

                    elif self.parameters['simulation_type'] == 'evaporation':
                        pyplot.xlabel('Concentration Factor (CF)', fontsize = label_font)
                        pyplot.ylabel('Quantity (%s)' %(unit), fontsize = label_font)  

                        experimental_loop = []
                        iteration = 0
                        cf_series = []
                        concentration_series = []
                        quantity_of_steps_index = 0  
                        self.variables['initial_solution_mass'] = self.results['csv_data'].at[0, 'mass_H2O']
                        for index, row in self.results['csv_data'].iterrows():
                            try:
                                if self.results['csv_data'].at[index+1, 'mass_H2O']:
                                    if self.results['csv_data'].at[index, 'step'] >= 1:
                                        concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                        solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                        cf_series.append(self.variables['initial_solution_mass'] / solution_mass)   
                                    else:
                                        print('ERROR: The SELECTED_OUTPUT file possesses an unexcepted data structure.')

                            except:
                                concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                cf_series.append(self.variables['initial_solution_mass'] / solution_mass)  

                                experimental_loop.append(f'{mineral} [{mineral_formula}]')
                                pyplot.plot(mass_series,concentration_series)

                                cf_series = []
                                concentration_series = []
                                concentration_series.append(self.results['csv_data'].at[index, mineral]) 
                                solution_mass = self.results['csv_data'].at[index, 'mass_H2O']
                                cf_series.append(self.variables['initial_solution_mass'] / solution_mass)         
                                                         
                    # export the figure
                    figure = illustrate(pyplot)
                    self.export_plot(figure, plot_title, export_format)


    def export_plot(self, figure, plot_title, mineral = None, individual_plots = False, export_format = 'svg'):
        """
        Export the plots to the current working directory  
        """
        # define the output name
        export_filename_progenitor = re.sub('(\.\w+)', '', self.parameters['selected_output_file_name'])
        if not re.search('(scaling|brine)', export_filename_progenitor, flags=re.IGNORECASE):
            export_name = '{}, {}'.format(export_filename_progenitor, self.parameters['output_perspective'])
        else:
            export_name = export_filename_progenitor
        if individual_plots:
            export_name += f'_{mineral}'
                                                         
        self.results['figures'][export_name] = {'figure':figure, 'title':plot_title}
                                
        # export the plot
        file_number = 0
        if not os.path.exists('%s.%s' %(export_name, export_format)):
            self.results['figures'][export_name]['figure'].savefig('%s.%s' %(export_name, export_format))
        elif os.path.exists('%s.%s' %(export_name, export_format)):
            while os.path.exists('%s_%s.%s' %(export_name, file_number, export_format)):
                file_number += 1
            figure.savefig('%s_%s.%s' %(export_name, file_number, export_format))


    def input_file(self, operating_system, phreeqc_path, database_selection, simulation_type, simulation_title, water_selection, quantity_of_modules = 1, module_characteristics = {}, output_perspective = 'scaling', domain = 'dual', permeate_approach = 'linear permeate', permeate_efficiency = 1, head_loss = -0.15, final_cf = 2, custom_water_parameters = {}, ignored_minerals = [], existing_parameters = {}, export_figure = True):
        """
        Concisely create an input file of the software  
        """
        self.define_general(operating_system, phreeqc_path, database_selection, simulation_type, simulation_title)
        self.transport(module_characteristics = module_characteristics, quantity_of_modules = quantity_of_modules, domain = domain, output_perspective = output_perspective)
        self.reaction(permeate_approach = permeate_approach, permeate_efficiency = permeate_efficiency, head_loss = head_loss, final_cf = final_cf)
        self.solutions(water_selection = water_selection, custom_water_parameters = custom_water_parameters)
        self.equilibrium_phases(ignored_minerals = ignored_minerals, existing_parameters = existing_parameters)
        self.selected_output()
        self.export()
        
            
    def complete_simulation():
        self.execute()
        self.process_selected_output(export_figure = export_figure)