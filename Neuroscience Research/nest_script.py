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

freq_sin = 4
#sin = create_sin2d(500,freq_sin)
#plt.imshow(sin)
#plt.colorbar()
#plt.show()
#plt.imsave(input_images_path + '/sinusoide.png',sin)


########################################################### Nest ###################################################################

nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads':local_num_threads})
nest.SetKernelStatus({"total_num_virtual_procs": total_num_virtual_procs})
nest.SetKernelStatus({'print_time':True})
nest.SetKernelStatus({'overwrite_files':True})
nest.SetKernelStatus({"resolution": resolution})

nest.CopyModel("izhikevich","exc", RS_dict)
nest.CopyModel("izhikevich","inh", FS_dict) 

################################################################# Main ################################################################## 

num_orientations = 2
result = main_all_orientations(input_images_path + '/sinusoide.png',num_orientations)

############################################################## Simulation ################################################################

t = time.time()
nest.Simulate(simulation_time)
print('tiempo de simulacion: ', time.time() - t)

######################################################### Data Treatment #################################################################

layers_to_record = {}
spike_detectors = {}
for i in result:
    layers_to_record.update(dict(list(result[i].items())[:2]))
    spike_detectors.update(dict(list(result[i].items())[2:]))
    
save_dict(layers_to_record,'to_record_layer')
save_dict(spike_detectors,'to_record_sd')

for layer,j in zip(layers_to_record,range(0,len(layers_to_record))):
    tp.DumpLayerNodes(layers_to_record[layer],'positions-'+str(layer))

