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
import nest
import pylab
import nest.topology as tp
from nest_values import *

#################################################### Folders ################################################################

def create_folder(path_name):
    if not os.path.exists(path_name):
        os.makedirs(path_name)

def remove_contents(path_name):
    for filename in os.listdir(path_name):
        file_path = os.path.join(path_name, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))
            


##################################################### Create estimuli #########################################################

def create_sin2d(size,freq):
    array = np.zeros((size,size))
    for j in range(0,size):
        for i in range(0,size):
            array[i][j] = (math.sin( (i / (size/(4*freq)) ) * math.pi/2  ))
    return array   


##################################################### Gabor ####################################################################

def gabor_filter(K_size,Lambda, Theta, Sigma, Gamma, Psi, mode):

    gabor = np.zeros((K_size, K_size), dtype=np.float32)

    for x in range(K_size):
        for y in range(K_size):
            
            # distancias
            dx = x - K_size // 2
            dy = y - K_size // 2

            x_ = np.cos(Theta) * dx + np.sin(Theta) * dy
            y_ = -np.sin(Theta) * dx + np.cos(Theta) * dy

            gabor[x, y] = np.exp(-(x_**2 + Gamma**2 * y_**2) / (2 * Sigma**2)) * np.cos(2*np.pi*x_/Lambda + Psi)
    
    if mode == 'on':
        gabor = gabor

    elif mode == 'off':
        gabor = -gabor

    return gabor

def apply_filter(gray_img, K_size, Lambda, Theta, Sigma, Gamma, Psi, mode):

    gray = np.pad(gray_img, (K_size//2, K_size//2), 'edge')

    output = np.zeros((gray_img.shape[0],gray_img.shape[1]), dtype=np.float32)

    gabor = gabor_filter(K_size = K_size, Lambda = Lambda, Theta = Theta, Sigma = Sigma, Gamma = Gamma, Psi = 0, mode = mode)
    for x in range(gray_img.shape[0]):
        for y in range(gray_img.shape[1]):
            output[x, y] = np.sum(gray[x : x + K_size, y : y + K_size] * gabor)
            
    return output

def gabor(gray_img,orientation_in_radians, mode):

    output = np.zeros((gray_img.shape[0],gray_img.shape[1]), dtype=np.float32) 
    
    orientation = orientation_in_radians*math.pi/180
    
    output = apply_filter(gray_img, K_size=K_size, Lambda=Lambda, Theta=orientation, Sigma=Sigma, Gamma=Sigma,Psi = Psi, mode = mode )
    output = np.clip(output, 0, np.max(output))
    
    return output

########################################################### Nest ###################################################################

def input_treatment(input_spike,x_cortex_size,y_cortex_size,max_value):
    input_as_img = Image.fromarray(input_spike)
    input_resized = input_as_img.resize((x_cortex_size,y_cortex_size))
    input_norm = np.divide(input_resized,np.max(input_resized)) * max_value
    
    input_transposed = input_norm.transpose()
    input_as_list = input_transposed.tolist()
    flat_list = [item for sublist in input_as_list for item in sublist]
    return flat_list

################################################################# Layers #################################################################

def create_layer(rows,columns,extent,elements,neurons_per_column):
    return tp.CreateLayer({'rows': rows,'columns': columns,'extent': extent,'elements': [elements,neurons_per_column],'edge_wrap': False}) 

########################################################### Connetivity ##################################################################

def update_dict(dictionary, kernel_type, orientation_i, orientation_j):
    new_dict = {'kernel': {kernel_type: {'rows': float(x_cortex_size), 'kappa': 1.,
                                         'orientation_i': orientation_i, 'orientation_j': orientation_j }}}
    new_dict.update(dictionary)
    return new_dict
    








