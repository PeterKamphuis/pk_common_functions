# -*- coding: future_fstrings -*-
import astroquery
from astroquery.ipac.ned import Ned
import pickle
import os
import numpy as np
#This file allow the NED distance table to be pickled or to be read for 
# an individual galaxy  
# 
 
def obtain_table(force_original=False,csv_name='NED30.5.1-D-17.1.2-20200415.csv'
                 ,pickle_name='Pickled_NED_Table.pkl',directory = ''):
    #csv_name = '/home/peter/Galaxies/ESO_358-60/TF/NED-D_distances/NED30.5.1-D-17.1.2-20200415.csv'
   # pickle_name = '/home/peter/Galaxies/ESO_358-60/TF/NED-D_distances/Pickled_NED_Table.pkl'

    if not os.path.isfile(f'{directory}{pickle_name}') or force_original:
        distance_table = read_ned_table(f'{directory}{csv_name}')
        with open(pickle_name,'wb') as tmp:
            pickle.dump(f'{directory}{pickle_name}',tmp) 
    else:
        with open(f'{directory}{pickle_name}','rb') as tmp:
            distance_table = pickle.load(tmp) 
    return distance_table


   
def read_ned_table(name):
    #read the NED Redshift independent distance (https://ned.ipac.caltech.edu/Library/Distances/)
    #Into a python dictionary    
    with open(name) as file:
        lines = file.readlines()
        
    table = {}
    method = []
    trig =  False
    for line in lines: 
        values = [x.strip('"') for x in line.split('","')]
       
        if values[1] == 'Exclusion Code':
            columns = values[1:]+['Comment']
            trig = True
            continue
        if trig:
           
            line_dict = {}    
            for i,col in enumerate(columns):
                line_dict[f'{col}'] = values[i]
            if len(values) != len(columns):
                print(line_dict,values,columns,len(values),len(columns))
                exit()
            if line_dict['Galaxy ID'] in [x for x in table]:
                table[line_dict['Galaxy ID']].append({'Distance': line_dict['D (Mpc)'], 'Method': line_dict['Method']})
            else:
                print(f'processing galaxy {line_dict["Galaxy ID"]}')
                table[line_dict['Galaxy ID']] = [{'Distance': line_dict['D (Mpc)'], 'Method': line_dict['Method']}] 
           
            
            if line_dict['Method'] not in method:
                method.append(line_dict['Method'])  
    return table





def get_ned_distance(names, method = None,force_original=False,directory = ''):


    table = obtain_table(force_original=force_original,directory=directory)

    #If method is not set we use all non rulers
    # based on table 2 at  (https://ned.ipac.caltech.edu/Library/Distances/)
    standards = ['SNIa SDSS', 'SNIa', 'Cepheids', 'CMD', 'FGLR', 'PNLF', 'TRGB',\
            'Brightest Stars', 'Statistical', 'SBF', 'Quasar spectrum', 'BCG',\
            'GRB', 'AGN time lag', 'Black Hole',   'Carbon Stars', 'Red Clump',\
            'RR Lyrae', 'SGRB', 'SZ effect', 'Wolf-Rayet', 'Horizontal Branch',\
            'Miras', 'BL Lac Luminosity',  'PAGB Stars',  'Type II Cepheids',\
            'GCLF', 'GC SBF', 'Novae', 'S Doradus Stars', 'Blue Supergiant', \
            'Delta Scuti', 'M Stars', 'OB Stars', 'RV Stars', 'HII LF','RSV Stars',\
            'SNII radio', 'AGB', 'B Stars', 'Subdwarf fitting', 'White Dwarfs', \
            'SX Phe Stars' ]
    
  
    #Grav. Wave is not listed in table 2, added it two secondary 
    secondary = [ 'FP', 'Tully-Fisher','IRAS','Sosies','Tully est', 'D-Sigma',\
        'Faber-Jackson',  'Tertiary', 'Magnitude', 'GeV TeV ratio','Diameter',\
        'GC FP', 'Radio Brightness', 'L(H{beta})-{sigma}', 'Mass Model','Grav. Wave',\
        'Dwarf Ellipticals', 'LSB galaxies', 'H I + optical distribution']
    rulers = ['Ring Diameter', 'SNII optical','G Lens', 'Proper Motion',\
        'HII region diameter''Orbital Mech.', 'Eclipsing Binary', 'GC K vs. (J-K)',\
        'CO ring diameter', 'GC radius', 'Maser',  'Grav. Stability Gas. Disk',\
        'Dwarf Galaxy Diameter', 'Jet Proper Motion']
    if method is None:
        method = standards
    if 'standards' in method:
        method.pop(method.index('standards'))
        method = method+standards
    if 'rulers' in method:
        method.pop(method.index('rulers'))
        method = method+rulers
    if 'secondary' in method:
        method.pop(method.index('secondary'))
        method = method+secondary
 
    distances = {}
    for name in names:
        distances[name] = {'Distance': [None,None], 'Uni_Name': None}
        print(f'Start search for data of {name}')
        try:
            result_table = Ned.query_object(name,get_query_payload=False)
        except astroquery.exceptions.RemoteServiceError:   
            continue
    
        distances[name]['Uni_Name'] = result_table['Object Name'].value[0]
        if result_table['Object Name'].value[0] in table:
           all_distances = table[result_table['Object Name'].value[0]] 
        else:
            continue
    
        found_distances = []
        for ind_div in all_distances:
            if ind_div['Method'] in method:
                found_distances.append(float(ind_div['Distance']))
    
        if len(found_distances) == 0:
            final_distance = [None,None]
        elif len(found_distances) > 2:
            final_distance = [np.mean(found_distances), np.std(found_distances)]
        else:
            final_distance = [np.mean(found_distances), np.mean(found_distances)*0.1]
        distances[name]['Distance'] = final_distance
    return distances