#!/usr/bin/env python3

import numpy as np
import math
import pandas as pd
import os
import time
from glob import glob
import shutil
from tqdm import tqdm
import pickle

from nest_values import *
from funciones   import *

############################################################ Move Files #######################################################################

create_folder(sd_path); remove_contents(sd_path)
create_folder(df_folder); remove_contents(df_folder)
create_folder(positions_path); remove_contents(positions_path)

files = glob('*spike_detector*')
for file in files:
    shutil.move(file, sd_path + '/' + file)

files = glob('*positions-*')
for file in files:
    shutil.move(file, positions_path + '/' + file)
    
########################################################### Read files ###################################################################
files = glob(positions_path + '/*')
pos = []
for file in files:
    positions = pd.read_table(file,names = ['Number','x_pos','y_pos'], index_col=False, sep = ' ')
    pos.append(positions)

positions = pd.concat(pos)
del(pos)

layers_to_record = load_dict('to_record_layer')
spike_detectors = load_dict('to_record_sd')

######################################################### To dataframe ###################################################################

print('Fixing the data...')
t = time.time()
spike_detectors = list(spike_detectors.values())
for layer, spk in zip(layers_to_record,spike_detectors):

    data = []
    num_layer = layers_to_record[layer][0]
    files = glob('spk_detectors_folder/spike_detector-' + str(spk[0]) + '-*')
    if files == []:
        files = glob('spk_detectors_folder/spike_detector-' + '0' + str(spk[0]) + '-*')
    for spk_detector in files:

        df = pd.read_table(spk_detector,names = ['Number','Time'], index_col=False)
        data.append(df)
    
    data = pd.concat(data)
    data = data.set_index(([pd.Index([i for i in range(0,len(data))])]))
    data['Number'] = data.Number.astype(float)
    data = pd.merge(data,positions,how = 'left',on = 'Number' )
    
    #Save dataframe
    data.to_pickle(df_folder + '/data_' + str(layer) +'.pkl')
    print('Numero de spikes ' +str(layer)+': '+str(len(data)))
    
del(data)
del(positions)
print('Tiempo tratamiento de datos: '+str(time.time() - t))










