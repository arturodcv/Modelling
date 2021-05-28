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


########################################################### Image processing ####################################################################

create_folder(gabor_folder); remove_contents(gabor_folder)
create_folder(sd_path); remove_contents(sd_path)
create_folder(df_folder); remove_contents(df_folder)
create_folder(positions_path); remove_contents(positions_path)

t = time.time()
gabors_to_nest = full_img_filtering(images_to_simulate,num_orientations)
gabors_time = time.time() - t
#save_gabors(gabors_to_nest, images_to_simulate,num_orientations);

############################################################  Connectivity ######################################################################

t = time.time()
layers, poiss_layers = main_all_orientations(num_orientations)
conn_time = time.time() - t
print("All layers succesfully connected!")


############################################################## Simulation #######################################################################
 
t = time.time()
#print("Lets simulate!")
for i in range(0,num_images_to_simulate):
    set_poisson_values(gabors_to_nest['image_' + str(i)], poiss_layers, num_orientations)
    nest.Simulate(ms_per_stimuli)
sim_time = time.time() - t

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
    
print("Times: \n Building architecture: " + str(conn_time) + "\n Image processing: " + str(gabors_time) + "\n Simulation: " + str(sim_time))

