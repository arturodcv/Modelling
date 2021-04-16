#####!/usr/bin/env python3
#!/bin/env python

import nest
import pylab
import nest.topology as tp
nest.Install('mymodule')

import cv2
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
from PIL import Image
import time
from glob import glob
import shutil
from tqdm import tqdm
from termcolor import cprint
import pickle

from nest_values import *
from funciones   import *


######################################################### Imagen ###################################################################
input_images_path = 'input_images_folder' 
#create_folder(input_images_path)

freq_sin = 4
#sin = create_sin2d(500,freq_sin)
#plt.imshow(sin)
#plt.colorbar()
#plt.show()
#plt.imsave(input_images_path + '/sinusoide.png',sin)

#img = cv2.imread("circle.png")
#img = cv2.imread("Captura.PNG")
img = cv2.imread(input_images_path + "/sinusoide.png")
gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

########################################################### Nest ###################################################################

output_gabor_0 = gabor(gray_img,0,'on')
input_spike = output_gabor_0 
flat_list_0_on = input_treatment(input_spike,x_cortex_size,y_cortex_size,max_normalized_value)


nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads':local_num_threads})
nest.SetKernelStatus({"total_num_virtual_procs": total_num_virtual_procs})
nest.SetKernelStatus({'print_time':True})
nest.SetKernelStatus({'overwrite_files':True})
nest.SetKernelStatus({"resolution": resolution})


################################################################# Modelos ################################################################

RS_dict =  {'a':0.02, 'b':0.2, 'c':-65.,'d':8.0, 'V_th':30.}
FS_dict =  {'a':0.1, 'b':0.2, 'c':-65., 'd':2.0, 'V_th':30.}
nest.CopyModel("izhikevich","exc", RS_dict)
nest.CopyModel("izhikevich","inh", FS_dict) 

################################################################# Layers ################################################################# 

l_exc_0_on  = create_layer(x_cortex_size,y_cortex_size,extent,'exc',neurons_per_column_exc)
l_inh_0_on  = create_layer(x_cortex_size,y_cortex_size,extent,'inh',neurons_per_column_inh)

########################################################### Poisson generator ##############################################################

l_poiss_0_on = create_layer(x_cortex_size,y_cortex_size,extent,'poisson_generator',1)
fixed_list = [i*factor + poisson_bias for i in flat_list_0_on]
nest.SetStatus(nest.GetNodes(l_poiss_0_on)[0],'rate', fixed_list)

########################################################### Spike detectors ################################################################

#sd_parrot = nest.Create('spike_detector', params = {'to_file': True})
sd_exc_0_on = nest.Create('spike_detector', params = {'to_file': True, 'to_memory': False})
#sd_inh_0_on = nest.Create('spike_detector', params = {'to_file': True})

########################################################### Parrot layers #################################################################

#l_parrots = create_layer(x_cortex_size,y_cortex_size,extent,'parrot_neuron',1)

########################################################### Dictionaries ##################################################################

dict_poiss_to_exc  = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},'weights': input_weight_exc, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_poiss_to_inh  = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},'weights': input_weight_inh, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_divergent = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

def create_lat_dict(conn_type,weight_type,min_delay_type,slowness_type, synapse_model):
    return  {'connection_type': 'convergent',                   
             'weights': weight_type,
             'delays': {'linear':{'c':min_delay_type,'a':slowness_type}},
             'synapse_model': synapse_model,
             'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses
            }

dict_lat_exc = create_lat_dict('convergent',0.2 * weight_exc,delay_exc_min,slowness_exc, syn_model_exc)

dict_lat_inh =  {'connection_type': 'convergent',
                 'mask': {'grid':{'rows':n_rows_latconn_inh_4cbeta,'columns':n_cols_latconn_inh_4cbeta}, 
                          'anchor':{'row':(n_rows_latconn_inh_4cbeta-1)//2,'column':(n_cols_latconn_inh_4cbeta-1)//2}},
                 'delays': {'linear':{'c':delay_inh_min,'a':slowness_inh}},
                 'kernel': {'gaussian2D':{'p_center':1., 'sigma_x':stddev_lat_conn_inh, 'sigma_y':stddev_lat_conn_inh}},
                 'weights': 0.2 * weight_inh , 
                 'synapse_model':syn_model_inh,
                 'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}


########################################################### Connetivity ##################################################################
#Poisson to layer
#tp.ConnectLayers(l_poiss_0_on,l_parrots, dict_divergent)
tp.ConnectLayers(l_poiss_0_on, l_exc_0_on, dict_poiss_to_exc)
tp.ConnectLayers(l_poiss_0_on, l_inh_0_on, dict_poiss_to_inh)

tp.ConnectLayers(l_exc_0_on, l_exc_0_on, update_dict(dict_lat_exc,'PlosOne_J',0.0,0.0))
tp.ConnectLayers(l_exc_0_on, l_inh_0_on, update_dict(dict_lat_exc,'PlosOne_W',0.0,0.0))
tp.ConnectLayers(l_inh_0_on, l_exc_0_on, dict_lat_inh)
tp.ConnectLayers(l_inh_0_on, l_inh_0_on, dict_lat_inh)

#Record connections
#leaves_parrot = nest.GetLeaves(l_parrots, local_only=True)[0]
#nest.Connect(leaves_parrot, sd_parrot)
#spike_detectors = [sd_parrot]
#layers_to_record = {'l_parrots': l_parrots} 


#leaves_inh_0_on = nest.GetLeaves(l_inh_0_on, local_only=True)[0]
#nest.Connect(leaves_inh_0_on, sd_inh_0_on)
#spike_detectors = [sd_exc_0_on, sd_inh_0_on]
#layers_to_record = {'l_exc_0_on': l_exc_0_on, 'l_inh_0_on': l_inh_0_on}

leaves_exc_0_on = nest.GetLeaves(l_exc_0_on, local_only=True)[0]
nest.Connect(leaves_exc_0_on, sd_exc_0_on)
layers_to_record = {'l_exc_0_on': l_exc_0_on}
spike_detectors = [sd_exc_0_on]

to_record_layer = open("to_record_layer.pkl", "wb")
pickle.dump(layers_to_record, to_record_layer)
to_record_layer.close()

with open('to_record_sd.pkl', 'wb') as f:
      pickle.dump(spike_detectors, f)

############################################################ Simulation ##################################################################

t = time.time()
nest.Simulate(simulation_time)
print('tiempo de simulacion: ', time.time() - t)

######################################################### Data Treatment #################################################################

for layer,j in zip(layers_to_record,range(0,len(layers_to_record))):
    tp.DumpLayerNodes(layers_to_record[layer],'positions-'+str(layer))































