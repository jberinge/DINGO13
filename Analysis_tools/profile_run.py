# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 14:34:44 2015

@author: imchugh
"""
import pandas as pd

import profile_scripts as ps
import DataIO as io

reload (ps)

configs_dict = {'met_file_in': '/home/imchugh/Ozflux/Sites/Whroo/Data/' \
                               'Processed/all/Whroo_2011_to_2015_L4.nc'}

outer_dict = {'C_names_dict': {'Cc_1': 'Cc_LI840_0.5m',
                               'Cc_2': 'Cc_LI840_2m',
                               'Cc_3': 'Cc_LI840_4m',
                               'Cc_4': 'Cc_LI840_8m',
                               'Cc_5': 'Cc_LI840_16m',
                               'Cc_6': 'Cc_LI840_36m'},
              'met_names_dict': {'Ta': 'Ta',
                                 'Press': 'ps'},
              'C_heights_dict': {'Cc_1': 0.5,
                                 'Cc_2': 2,
                                 'Cc_3': 4,
                                 'Cc_4': 8,
                                 'Cc_5': 16,
                                 'Cc_6': 36}}
                
              

full_C_df = ps.get_CO2_data()

C_df = ps.truncate_data(full_C_df)

met_df = io.OzFluxQCnc_to_pandasDF(configs_dict['met_file_in'], ['Ta', 'ps'])[0]

# Check the dates and use only those which overlap, then reindex and join
start = met_df.index[0] if C_df.index[0] < met_df.index[0] else C_df.index[0]
end = met_df.index[-1] if C_df.index[-1] > met_df.index[-1] else C_df.index[-1]
new_index = pd.date_range(start, end, freq = '30T')
C_df = C_df.reindex(new_index)
met_df = met_df.reindex(new_index)
df = C_df.join(met_df)

storage_df = ps.process_data(df, outer_dict)

storage_df.to_csv('/home/imchugh/Analysis/Whroo/Data/Flux_and_met/temp.csv',
                  index_label = 'date_time')