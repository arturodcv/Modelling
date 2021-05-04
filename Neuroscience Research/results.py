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


orientations = [-45.0,0.0,45.0,90.0]
neuron_types = ['exc', 'inh']

for orientation_to_read in orientations:
    for exc_or_inh in neuron_types:
        print("_______________________________________________________________________________________________ \n")
        print("_____________________ Results for layer " + exc_or_inh + " and orientation " + str(orientation_to_read) + ' ______________________')
        print("\n _______________________________________________________________________________________________")
        
        path = results_path + '/results_' + str(orientation_to_read) + '_' + str(exc_or_inh)
        create_folder(path)
        remove_contents(path)
        # Read dataframe
        data = read__and_fix_dataframe(orientation_to_read,exc_or_inh)
        #Generate frames
        times = generate_frames(data)
        frames, complementary_time_list = generate_empty_frames(times)
        # Read frames and create avg image and video
        img_array = read_frames(frames)
        create_video(img_array,orientation_to_read ,exc_or_inh, path)
        create_avg_img(img_array,orientation_to_read ,exc_or_inh , path)
        # Create eeg plot
        eeg = get_eeg(times, complementary_time_list, orientation_to_read, exc_or_inh, path)
        # Get eeg frequencies plot
        get_frequencies(eeg,orientation_to_read,exc_or_inh, path)






