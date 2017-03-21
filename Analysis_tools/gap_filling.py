# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:36:03 2015

Generic imputation routines

@author: ian_mchugh
"""
import numpy as np
import pdb
from scipy.interpolate import griddata as griddata_sc
import ffnet as ffnet_class
from ffnet import ffnet, tmlgraph
import warnings
import os

def generic_2d_linear(data_2d):

    """
    Takes a 2d array as input and;
     1) tiles this into a 3 x 3 space (9 repeats of the original 2d array in 3 
        columns and 3 rows)
     2) removes the missing data (c.missing_value) from the tiled array
     3) does a bi-linear interpolation to replace the the missing data
     4) returns the central tile
     Note: the effect is to replace missing data in the original 2d array with 
     data from a bi-linear interpolation, the tiling repeats the original array 
     along its boundaries to avoid problems at the array edges.
    """
    
    data_2d_tiled = np.tile(data_2d, (3,3))   
    
    num_x = np.shape(data_2d_tiled)[1]
    flat_x = np.arange(0, num_x)
    num_y = np.shape(data_2d_tiled)[0]
    flat_y = np.arange(0, num_y)

    # Make the regular grid to project the data onto
    coords_x, coords_y = np.meshgrid(flat_x, flat_y)
    
    # Make a flat array of the tiled data
    data_flat = data_2d_tiled.flatten()

    # Define an index that will return all valid data for the array
    index = np.where(~np.isnan(data_flat))
    
    # Generate a 2d array with existing 2d coordinates of the tiled data
    data_coords = np.column_stack([coords_x.flatten(), 
                                   coords_y.flatten()])

    # Do the interpolation
    grid_z = griddata_sc(data_coords[index], data_flat[index], 
                         (coords_x, coords_y), method = 'linear')
    
    # Return the central tile
    return grid_z[num_y / 3: num_y / 3 * 2, num_x / 3: num_x / 3 * 2].T[:, 0]

# Simple linear interpolation
def interp_params(param_rslt_array):

    def do_interp(array_1D):
        xp = np.arange(len(arr))
        fp = array_1D[:]
        nan_index = np.isnan(fp)
        fp[nan_index] = np.interp(xp[nan_index], xp[~nan_index], fp[~nan_index])
        return fp   
    
    arr = param_rslt_array.copy()    
    num_vars = np.shape(arr)
    if len(num_vars) == 1:
        arr = do_interp(arr)
    else:
        num_vars = num_vars[1]
        for i in range(num_vars):
            arr[:, i] = do_interp(arr[:, i])

    return arr            
    
def train_ANN(inputs_array, target_array, iterations, node_architecture, 
              **configs_dict):

    # Same first dimension?
    if not inputs_array.shape[0] == target_array.shape[0]:
        raise Exception('Input and target arrays must have same first ' \
                        'dimension!')

    # Specified number of input nodes matches second dim of input array?
    n_input_nodes = node_architecture[0]
    if len(inputs_array.shape) == 1:
        sec_dim_inputs = 1
    else: 
        sec_dim_inputs = inputs_array.shape[1]
    if not n_input_nodes == sec_dim_inputs:
        raise Exception('Specified input node architecture (n = %s) ' \
                        'incompatible with passed input arrays... Returning!'
                        %str(n_input_nodes))

    # Specified number of target nodes matches second dim of target array?
    n_target_nodes = node_architecture[-1]
    if len(target_array.shape) == 1:
        sec_dim_target = 1
    else: 
        sec_dim_target = target_array.shape[1]
    if not n_target_nodes == sec_dim_target:
        raise Exception('Specified target node architecture (n = %s) ' \
                        'incompatible with passed input arrays... Returning!'
                        %str(n_target_nodes))        

    # Missing data in inputs array? (Warning only)
    if np.isnan(inputs_array).any():
        missing_inputs_flag = True
        warnings.warn('Specified ANN training input variables contain missing ' \
                      'data. NaNs will be inserted into prediction series!')
    else:
        missing_inputs_flag = False

    # Missing data in target array? (Warning only)
    if np.isnan(target_array).any():
        missing_target_flag = True
        warnings.warn('Specified ANN training target variables contain missing ' \
                      'data. These will be removed for training!')
    else:
        missing_target_flag = False

    # Check if saving trained network
    save_flag = False
    if 'save_network' in configs_dict.keys():
        if configs_dict['save_network']: 
            save_flag = True
        if not 'network_filepath' in configs_dict.keys():
            raise Exception('You must specify a file path if you wish to ' \
                            'save a new network!')
        else:
            split_pathname_list = os.path.split(configs_dict['network_filepath'])
            if not os.path.isdir(split_pathname_list[0]):
                raise Exception('The specified file path is not valid!')
            if split_pathname_list[1] == '':
                print 'Filename not supplied - using this_net.ann!'
                configs_dict['network_filepath'] = os.path.join(split_pathname_list[0],
                                                                'this_net.ann')
                      
    # Check if doing testing
    test_flag = False
    if 'test' in configs_dict:
        if configs_dict['test']:
            test_flag = True

    # Create a second series with nans dropped
    if missing_inputs_flag or missing_target_flag:
        new_array = np.empty([inputs_array.shape[0], 
                              sec_dim_inputs + sec_dim_target])
        new_array[:, :sec_dim_target] = target_array
        new_array[:, sec_dim_target:] = inputs_array
        new_array = new_array[~np.isnan(new_array).any(axis = 1)]
        clean_target_array = new_array[:, :sec_dim_target]
        clean_inputs_array = new_array[:, sec_dim_target:]

    # Generate network and train
    conec = tmlgraph(node_architecture)
    net = ffnet(conec)
    net.train_tnc(clean_inputs_array, clean_target_array, 
                  maxfun = iterations, messages = 1)

    # Save network if requested
    if save_flag:
        ffnet_class.savenet(net, configs_dict['network_filepath'])

    # Generate full series from inputs
    predict_array = net.call(inputs_array)

    # Do testing if requested
    if test_flag:
        vars_list = ['slope', 'intercept', 'r-value', 'p-value', 
                     'slope stderr', 'estim. stderr']
        valid_predict_array, stats_list = net.test(clean_inputs_array, 
                                                   clean_target_array)
        stats_dict = {var: stats_list[0][i] for i, var in enumerate(vars_list)}
        return predict_array, stats_dict
    else:
        return predict_array