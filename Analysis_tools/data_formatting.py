# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 12:05:16 2016

@author: imchugh
"""
# Python modules
import os
import numpy as np
import pdb

# My modules
import DataIO as io
import data_filtering as filt

def get_model_NEE_from_OzFluxQCncL6(f):
    """
    Builds a continuous time series (i.e. day and night) from OzFluxQC L6
    model data (specifically nocturnal ER and daytime NEE)
    Pass: valid file target (path / name string)
    Returns: data dictionary containing 2 key / value pairs:
                 - 'date_time': numpy array of Python datetimes for each datum
                                in original time series
                 - 'Fc_model': numpy array of corresponding model data (float)
    """
    # Grab data    
    data_dict = io.OzFluxQCnc_to_data_structure(f,
                                                var_list = ['Fsd',
                                                            'ER_SOLO_all', 
                                                            'Fc_SOLO'])

    # Create continuous model series from ER and Fc series
    temp_dict = {'Fc_model': np.concatenate([data_dict['ER_SOLO_all']
                                             [data_dict['Fsd'] < 10],
                                             data_dict['Fc_SOLO']
                                             [data_dict['Fsd'] >= 10]]),
                 'date_time': np.concatenate([data_dict['date_time']
                                             [data_dict['Fsd'] < 10],
                                             data_dict['date_time']
                                             [data_dict['Fsd'] >= 10]])}
    temp_dict = filt.sort_dict_on_index_variable(temp_dict, 'date_time')
    data_dict['Fc_model'] = temp_dict['Fc_model']
    
    return {'date_time': temp_dict['date_time'],
            'Fc_model': temp_dict['Fc_model']}

def rename_data_dict_vars(data_dict, names_dict, pass_through = True):
    """
    Renames variables in a data dictionary according to key / value pairs in 
    a names dictionary specifying the old (key) and new (value) names (note 
    that by default variables in the data dictionary that are not in the keys
    of the names dictionary are passed through - switching pass_through to 
    'False' will stop these variables from being passed to dictionary that is
    returned).
    Pass the following args: 'data_dict' - a dictionary containing key / value
                                           pairs of variable name and numpy
                                           array
                             'names_dict' - a dictionary containing key / value
                                            pairs of old variable name and new
                                            variable name
    Optional kwargs: 'pass_through' - boolean specifying whether to keep 
                                      variables that are not specified in the 
                                      keys of the names dictionary (this can 
                                      be used to pass through variables that 
                                      do not require name changes)
    """
    wrong_keys = [i for i in names_dict.keys() if not i in data_dict.keys()]
    if not len(wrong_keys) == 0:
        wrong_key_str = ', '.join(wrong_keys)
        raise Exception('The following keys specified in names dictionary ' \
                        'not present in data dictionary: ' + wrong_key_str)

    new_dict = {}
    for key in names_dict.keys():
        new_dict[names_dict[key]] = data_dict.pop(key)

    if pass_through:
        for key in data_dict.keys():
            new_dict[key] = data_dict[key]
        
    return new_dict

def standard_names_dictionary():
    
    return {'carbon_flux':'NEE_series',
            'carbon_storage': 'Sc',
            'temperature': 'TempC',
            'solar_radiation': 'Fsd',
            'vapour_pressure_deficit': 'VPD',
            'friction_velocity': 'ustar',
            'wind_speed': 'ws',
            'modelled_carbon_flux': 'NEE_model',
            'soil_moisture': 'Sws',
            'generic': 'filter_var'}

# Get the data and format appropriately
def get_data(configs_dict):

    # Get file extension and target
    paths_dict = configs_dict['files']
    ext = os.path.splitext(paths_dict['input_file'])[1]
    data_input_target = os.path.join(paths_dict['input_path'],
                                     paths_dict['input_file'])

    # Initialise name change dictionary with new names via common keys
    oldNames_dict = configs_dict['variables']
    newNames_dict = {'carbon_flux':'Fc',
                     'temperature': 'TempC',
                     'solar_radiation': 'Fsd',
                     'vapour_pressure_deficit': 'VPD',
                     'friction_velocity': 'ustar',
                     'wind_speed': 'ws'}
    names_dict = {oldNames_dict[key]: newNames_dict[key] for key in oldNames_dict}                     

    # get data (screen only the Fc data to obs only)
    if ext == '.nc':
        Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                  var_list = [oldNames_dict
                                                              ['carbon_flux']],
                                                  QC_accept_codes = [0])
        Fc_dict.pop('date_time')
        ancillary_vars = [oldNames_dict[var] for var in oldNames_dict.keys() 
                          if not var == 'carbon_flux']
        ancillary_dict, global_attr = io.OzFluxQCnc_to_data_structure(
                                          data_input_target,
                                          var_list = ancillary_vars,
                                          return_global_attr = True)
        data_dict = dict(Fc_dict, **ancillary_dict)
    elif ext == '.df':
        data_dict, global_attr = io.DINGO_df_to_data_structure(
                                     data_input_target,
                                     var_list = oldNames_dict.values(),
                                     return_global_attr = True)

    data_dict['Fc_model'] = get_model_NEE_from_OzFluxQCncL6(data_input_target)
    data_dict = rename_data_dict_vars(data_dict, names_dict)

    return data_dict, global_attr