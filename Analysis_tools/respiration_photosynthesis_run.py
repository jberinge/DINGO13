# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 16:08:41 2016

@author: imchugh
"""

# Python standard modules
import os
import copy as cp
import numpy as np
import pdb

# My modules
import DataIO as io
import respiration as re
import photosynthesis as ps
import data_formatting as dt_fm
import data_filtering as data_filter

reload(dt_fm)
reload(ps)
reload(re)
reload(data_filter)

#------------------------------------------------------------------------------
# Fetch data from configurations
def get_data(configs_dict):

    # Get data (screen Fc data to obs only - keep gap-filled drivers etc)
    data_input_target = os.path.join(configs_dict['files']['input_path'],
                                     configs_dict['files']['input_file'])
    
    ext = os.path.splitext(data_input_target)[1]
    if ext == '.nc':
        Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                  var_list = [configs_dict['variables']
                                                                          ['carbon_flux']],
                                                  QC_accept_codes = [0])
        Fc_dict.pop('date_time')
        ancillary_vars = [configs_dict['variables'][var] for var in 
                          configs_dict['variables'] if not var == 'carbon_flux']
        ancillary_dict, attr = io.OzFluxQCnc_to_data_structure(
                                   data_input_target,
                                   var_list = ancillary_vars,
                                   return_global_attr = True)
        data_dict = dict(Fc_dict, **ancillary_dict)
    elif ext == '.df':
        data_dict, attr = io.DINGO_df_to_data_structure(data_input_target,
                              var_list = configs_dict['variables'].values(),
                              return_global_attr = True)

    # Rename to generic names used by scripts
    old_names_dict = configs_dict['variables']
    std_names_dict = dt_fm.standard_names_dictionary()
    try:
        new_dict = {std_names_dict[key]: data_dict[old_names_dict[key]] 
                    for key in old_names_dict.keys()}
        new_dict['date_time'] = data_dict.pop('date_time')
    except:
        raise Exception()
    
    return new_dict, attr
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Rebuild the master configuration file for passing to respiration and light 
# response (if requested) functions
def build_config_file(configs_master_dict, do_light_response):
    
    # Build a specific configuration file
    configs_dict = {'files': configs_master_dict['global_configs']['files'],
                    'global_options': (configs_master_dict['global_configs']
                                                          ['options'])}                                                          
    if do_light_response:
        configs_dict['variables'] = dict(configs_master_dict
                                         ['respiration_configs']['variables'],
                                         ** configs_master_dict
                                            ['photosynthesis_configs']
                                            ['variables'])
    else:
        configs_dict['variables'] = (configs_master_dict['respiration_configs']
                                                        ['variables'])
    return configs_dict                                                         
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------    
def main(use_storage = 'from_config', storage_var = 'from_config',
         ustar_threshold = 'from_config', 
         config_file = False, 
         do_light_response = False):
    """
    No positional arguments - prompts for a configuration file
    Kwargs: use_storage - if True then algorithm looks for a variable called 
                          Fc_storage and then sums that with Fc; if False uses
                          Fc alone; if set to 'from_config', looks in 
                          config file
            ustar_threshold - set to a particular value to override the ustar
                              threshold set in the global configs root item of
                              the configuration file (this is done so that can
                              be set in function call from another script)
    """
    
    # Do the respiration fit

    # Get master configuration file
    if not config_file:
        configs_master_dict = io.config_to_dict(io.file_select_dialog())
    else:
        configs_master_dict = io.config_to_dict(config_file)
    
    # Build custom configuration file for this script
    configs_dict = build_config_file(configs_master_dict, do_light_response)

    # Override default configuration storage_var if not set to from_config
    if not storage_var == 'from_config':
        if isinstance(storage_var, str):
            configs_dict['variables']['carbon_storage'] = storage_var
        else:
            raise Exception('storage_var kwarg must be a string if not '\
                            'set to from_config! Quitting...')

    # Get data
    data_dict, attr = get_data(configs_dict)

    # Override default configuration file use_storage if not set to from_config
    if not use_storage == 'from_config':
        if isinstance(use_storage, (bool, str)):
            configs_dict['global_options']['use_storage'] = use_storage
        else:
            raise Exception('use_storage kwarg must be a boolean if not '\
                            'set to from_config! Quitting...')

    # Override default configuration file ustar_threshold if requested by user
    if not ustar_threshold == 'from_config':
        if isinstance(ustar_threshold, (int, float, dict)):
            configs_dict['global_options']['ustar_threshold'] = ustar_threshold
        else:
            raise Exception('ustar_threshold kwarg must be set to a number or '\
                            'dictionary of numbers (year[int] as key, ' \
                            'threshold [float] as value)! Quitting...')

    # Sum Fc and Sc if storage is to be included, otherwise if requested, 
    # remove all Fc where Sc is missing
    if configs_dict['global_options']['use_storage']:
        data_dict['NEE_series'] = (data_dict['NEE_series'] + 
                                   data_dict['Sc'])
    elif configs_dict['global_options']['unify_flux_storage_cases']:
#        data_dict['NEE_series'][np.isnan(data_dict['Sc'])] = np.nan
        data_dict['NEE_series'][np.isnan(data_dict['filter_var'])] = np.nan

    # Remove low ustar data
    data_filter.screen_low_ustar(data_dict, 
                                 configs_dict['global_options']['ustar_threshold'],
                                 configs_dict['global_options']['noct_threshold'])     
    
    # Set up respiration configs and add measurement interval and output path
    re_configs_dict = configs_master_dict['respiration_configs']['options']
    re_configs_dict['measurement_interval'] = int(attr['time_step'])
    re_full_path = os.path.join(configs_dict['files']['output_path'],
                                re_configs_dict['output_folder'])
    if not os.path.isdir(re_full_path): os.makedirs(re_full_path)
    re_configs_dict['output_path'] = re_full_path
    
    # Calculate Re                            
    re_rslt_dict, re_params_dict, re_error_dict = re.main(cp.copy(data_dict), 
                                                          re_configs_dict)
    
    # Add original time series                                                              
    re_rslt_dict['NEE_series'] = data_dict['NEE_series']

    

    # Do the light response parameters (or not)
    if do_light_response:
    
        # Set up light response configs and add measurement interval and output path
        li_configs_dict = configs_master_dict['photosynthesis_configs']['options']
        li_configs_dict['measurement_interval'] = int(attr['time_step'])
        li_full_path = os.path.join(configs_dict['files']['output_path'],
                                    li_configs_dict['output_folder'])
        if not os.path.isdir(li_full_path): os.makedirs(li_full_path)
        li_configs_dict['output_path'] = li_full_path
        
        # Convert insolation to PPFD
        data_dict['PAR'] = data_dict['Fsd'] * 0.46 * 4.6
        
        # Call light response function
        li_rslt_dict, li_params_dict, li_error_dict = ps.main(data_dict, 
                                                              li_configs_dict, 
                                                              re_params_dict)

        # Add original time series                                                              
        li_rslt_dict['NEE_series'] = data_dict['NEE_series']
        
        # Write data to file
        var_order_list = ['date', 'Eo', 'Eo_error_code', 'rb', 'alpha', 'beta', 'k', 
                          'light_response_error_code']
        if li_configs_dict['use_nocturnal_rb']:
            var_order_list.insert(4, 'rb_error_code')
        io.array_dict_to_csv(li_params_dict, 
                             os.path.join(li_full_path, 'params.csv'), 
                             var_order_list)
        io.array_dict_to_csv(li_rslt_dict, os.path.join(li_full_path, 'Re_GPP.csv'), 
                             ['date_time', 'Re', 'GPP'])
        io.text_dict_to_text_file(li_error_dict, os.path.join(li_full_path, 'error_codes.txt'))

        return li_rslt_dict, li_params_dict
        
    else:
        
        # Write data to file
        io.array_dict_to_csv(re_params_dict, 
                             os.path.join(re_full_path, 'params.csv'), 
                             ['date', 'Eo', 'Eo_error_code', 'rb', 'rb_error_code'])
        io.array_dict_to_csv(re_rslt_dict, os.path.join(re_full_path, 'Re.csv'), 
                             ['date_time', 'Re'])
        io.text_dict_to_text_file(re_error_dict, os.path.join(re_full_path, 'error_codes.txt'))        
        
        return re_rslt_dict, re_params_dict

if __name__ == "__main__":

   main()