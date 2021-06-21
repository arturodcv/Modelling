#!/usr/bin/env python3
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
from collections import Counter
from collections import OrderedDict
from scipy.fft import rfft
import pickle
import scipy.signal

from nest_values import *
from funciones   import *

import datatable as dt

print("Results for total data: ")
path = results_path + '/data_total'; create_folder(path) ; remove_contents(path)   

data = read_and_fix_dataframe('','total')
times = data['Time'].tolist() ; times = (np.around(data['Time'],window_time)).tolist()
unique_times = np.unique(np.array(times)).tolist()
full_times = [i for i in range(0,int(simulation_time) + 1)]
complementary_time_list = list(set(full_times) - set(unique_times))
eeg = get_eeg(times, complementary_time_list, 'total', '_', path)
get_frequencies(eeg,'total','_', path)



orientations = [0.0,45.0,90.0,135.0]
neuron_types = ['l_exc', 'l_inh']

for orientation_to_read in orientations:
    for exc_or_inh in neuron_types:

        path = results_path + '/results_' + str(orientation_to_read) + '_' + str(exc_or_inh)
        create_folder(path); remove_contents(path)
        
        # Read dataframe
        data = read_and_fix_dataframe(orientation_to_read,exc_or_inh)

        if len(data) < 2:
            #print("Not enough results for orientation ",orientation_to_read, " ",exc_or_inh ); 
            continue
        #print("\nResults for layer " + exc_or_inh + " and orientation " + str(orientation_to_read) + ': ')
        #Generate frames
        times = generate_frames(data)
        frames, complementary_time_list = generate_empty_frames(times)
        # Read frames and create avg image and video
        img_array = read_frames(frames)
        create_video(img_array,orientation_to_read ,exc_or_inh, path)
        create_avg_img(img_array,orientation_to_read ,exc_or_inh , path)
        # Create eeg plot
        #eeg = get_eeg(times, complementary_time_list, orientation_to_read, exc_or_inh, path)
        # Get eeg frequencies plot
        #get_frequencies(eeg,orientation_to_read,exc_or_inh, path)
        

print("\n")




