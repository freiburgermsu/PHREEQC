from glob import glob
import json

def database_parsing(db):
    print(f'\n\n\n{db}')
    
    database = pandas.read_table(db, sep='\n')
    start_master = False
    if database.columns == ['SOLUTION_MASTER_SPECIES']:
        start_master = True
    database.columns = ['content']

    elements_rows = []
    minerals_rows = []
    elemental_parsing = False
    mineral_parsing = False
    for index, row in database.iterrows():
        if (re.search('SOLUTION_MASTER_SPECIES', row['content']) or start_master) and not elemental_parsing :
            while not re.search('SOLUTION_SPECIES', database.at[index, 'content']):
                split_row = database.at[index, 'content'].split()
                if all(not re.search('^#', entity) for entity in split_row):
                    elements_rows.append(split_row)
                index+= 1
            elemental_parsing = True

        if re.search('PHASES', row['content']) and not mineral_parsing:
            loop = False
            while not re.search('PITZER|EXCHANGE_MASTER_SPECIES|SURFACE_MASTER_SPECIES', database.at[index, 'content']):
                if not loop:
                    minerals_rows.append(['phases', 'formula'])
                    loop = True
                    index += 1
                    continue

                if re.search('(^\w+\s*\d*$)',database.at[index, 'content']):
                    reactants = database.at[index+1, 'content'].split(' = ')[0]
                    if all('#' not in entity for entity in reactants):
                        formula = reactants.split(' + ')[0].strip()
                        name = database.at[index, 'content']
                        name = re.sub('\s*\d*', '', name)
                        minerals_rows.append([name, formula])
                index+= 1
                
                if index == len(database):
                    break
            mineral_parsing = True

    # define the elements content for the database
    elements = pandas.DataFrame(elements_rows)
    elements.fillna(' ')
    elements.drop([0], inplace = True)
    for column in elements:
        nan_entries = 0
        alphanumeric_entries = 0
        for entry in elements[column]:
            if entry is not None:
                if re.search('[a-z]|[0-9]', entry, re.IGNORECASE):
                    alphanumeric_entries += 1
                else:
                    nan_entries += 1
            else:
                nan_entries += 1
        if nan_entries > alphanumeric_entries and len(elements.columns) > 5:
            print('deleted column: ', column)
            del elements[column]
    
    elements.columns = ['elements', 'species', 'alk', 'gfw_formula', 'element_gfw']
    print(elements)

    # define the minerals content for the database
    minerals = pandas.DataFrame(minerals_rows)
    minerals.columns = minerals.iloc[0]
    minerals = minerals.drop(0)
    print(minerals)

    return elements, minerals
    
    
def database_json_creation(database, elements, minerals):
    database_json = {'elements': {}, 'minerals': {}}
    
    # create the elements JSON
    for index, element in elements.iterrows():
        database_json['elements'][element['elements']] = {'charge_specie': element['alk'], 'gfw_formula':element['gfw_formula'], 'element_gfw':element['element_gfw']}
    
    # create the mienrals JSON
    for index, mineral in minerals.iterrows():
        database_json['minerals'][mineral['phases']] = mineral['formula']     
        
    # export the JSON files
    database_json_name = re.sub('.dat$', '.json', database)
    with open(database_json_name, 'w') as output:
        json.dump(database_json, output, indent = 4)

# execute the functions over all databases
erroneous_databases = []
for database in glob('./databases/*.dat'):
    try:
        elements, minerals = database_parsing(database)
        database_json_creation(database, elements, minerals)
    except:
        erroneous_databases.append(database)
        print(f'ERROR: {database}')

print(erroneous_databases)