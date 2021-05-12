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
import pickle

from nest_values import *
from funciones   import *


########################################################### Nest ###################################################################

nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads':local_num_threads})
nest.SetKernelStatus({'print_time':True})
nest.SetKernelStatus({'overwrite_files':True})
nest.SetKernelStatus({"resolution": resolution})

nest.CopyModel("izhikevich","exc", RS_dict)
nest.CopyModel("izhikevich","inh", FS_dict) 

################################################################# Main ################################################################## 

layers, poiss_layers = main_all_orientations(num_orientations)
print("All layers succesfully connected!")

############################################################## Simulation ################################################################

#input_files = [input_images_path + '/sinusoide_' + str(int(i * 180 / num_images_to_simulate) ) + '.png' for i in range(0,num_images_to_simulate)]
#input_files = [input_images_path + '/gaussian.png' ]
input_files = [input_images_path + '/sinusoide_0.png' ]


create_folder(gabor_folder)
remove_contents(gabor_folder)
 
t = time.time()
for i in range(0,num_images_to_simulate):
    set_poisson_values(input_files[i], poiss_layers, num_orientations)
    nest.Simulate(ms_per_stimuli)
print('tiempo de simulacion: ', time.time() - t) 

######################################################### Data Treatment #################################################################

layers_to_record = {}
spike_detectors = {}
for i in layers:
    layers_to_record.update(dict(list(layers[i].items())[:2]))
    spike_detectors.update(dict(list(layers[i].items())[2:]))
    
save_dict(layers_to_record,'to_record_layer')
save_dict(spike_detectors,'to_record_sd')

for layer,j in zip(layers_to_record,range(0,len(layers_to_record))):
    tp.DumpLayerNodes(layers_to_record[layer],'positions-'+str(layer))

