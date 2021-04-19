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
import pickle
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
            dx = x - K_size // 2; dy = y - K_size // 2
            x_ = np.cos(Theta) * dx + np.sin(Theta) * dy; y_ = -np.sin(Theta) * dx + np.cos(Theta) * dy
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
    
def main_img(img,orientation):
    img = cv2.imread(img)
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    output_gabor = gabor(gray_img,orientation,'on')
    flat_list = input_treatment(output_gabor,x_cortex_size,y_cortex_size,max_normalized_value)
    return flat_list
    
def main_create(flat_list):
    l_exc  = create_layer(x_cortex_size,y_cortex_size,extent,'exc',neurons_per_column_exc)
    l_inh  = create_layer(x_cortex_size,y_cortex_size,extent,'inh',neurons_per_column_inh)
    sd_exc = nest.Create('spike_detector', params = {'to_file': True, 'to_memory': False})
    sd_inh = nest.Create('spike_detector', params = {'to_file': True, 'to_memory': False})
    l_poiss = create_layer(x_cortex_size,y_cortex_size,extent,'poisson_generator',1)
    fixed_list = [i*factor + poisson_bias for i in flat_list]
    nest.SetStatus(nest.GetNodes(l_poiss)[0],'rate', fixed_list)
    return l_exc,l_inh,sd_exc,sd_inh,l_poiss
    
def main_one_orientation(img,orientation):
    orientation = orientation*math.pi/180
    flat_list = main_img(img,orientation)
    l_exc,l_inh,sd_exc,sd_inh,l_poiss = main_create(flat_list)
    main_self_connections(l_exc,l_inh,sd_exc,sd_inh,l_poiss, orientation)
    orientation = orientation*180 / math.pi
    out_dict = {'l_exc_'+str(orientation): l_exc,'l_inh_'+str(orientation): l_inh,'sd_exc_'+str(orientation): sd_exc,'sd_inh_'+str(orientation): sd_inh}
    return out_dict
    
def main_all_orientations(img,num_orientations):
    d = {}
    for i in range(0,num_orientations):
        orientation = i*180/num_orientations
        d["orientation_" + str(orientation)] = main_one_orientation(img,orientation)
    for i in range(0,num_orientations):
        for j in range(0,num_orientations):
            if i!=j:
                or_i = i*180/num_orientations;
                or_j = j*180/num_orientations
                l_exc_i = d['orientation_'+str(or_i)]['l_exc_'+str(or_i)]
                l_exc_j = d['orientation_'+str(or_j)]['l_exc_'+str(or_j)]
                l_inh_i = d['orientation_'+str(or_i)]['l_inh_'+str(or_i)]
                l_inh_j = d['orientation_'+str(or_j)]['l_inh_'+str(or_j)]
                main_lat_connections(l_exc_i,l_exc_j,l_inh_i,l_inh_j,i*math.pi/num_orientations,j*math.pi/num_orientations)                  
    return d

################################################################# Layers #################################################################

def create_layer(rows,columns,extent,elements,neurons_per_column):
    return tp.CreateLayer({'rows': rows,'columns': columns,'extent': extent,'elements': [elements,neurons_per_column],'edge_wrap': False}) 

########################################################### Connetivity ##################################################################

def create_lat_exc(kernel_type,kappa,orientation_i,orientation_j):
    return  {'connection_type': 'convergent',    
             'kernel': {kernel_type: {'rows': float(x_cortex_size), 'kappa': kappa,
                                      'orientation_i': orientation_i, 'orientation_j': orientation_j }},
             'weights': 0.2 * weight_exc,
             'delays': {'linear':{'c':delay_exc_min,'a':slowness_exc}},
             'synapse_model': syn_model_exc,
             'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses
            }    

def main_self_connections(l_exc,l_inh,sd_exc,sd_inh,l_poiss, self_orientation):
    tp.ConnectLayers(l_poiss, l_exc, dict_poiss_to_exc)
    tp.ConnectLayers(l_poiss, l_inh, dict_poiss_to_inh)
    tp.ConnectLayers(l_exc, l_exc, create_lat_exc('PlosOne_J',kappa_j,self_orientation,self_orientation))
    tp.ConnectLayers(l_exc, l_inh, create_lat_exc('PlosOne_W',kappa_w,self_orientation,self_orientation))
    tp.ConnectLayers(l_inh, l_exc, dict_lat_inh)
    tp.ConnectLayers(l_inh, l_inh, dict_lat_inh)
    leaves_exc = nest.GetLeaves(l_exc, local_only=True)[0]
    nest.Connect(leaves_exc, sd_exc)
    leaves_inh = nest.GetLeaves(l_inh, local_only=True)[0]
    nest.Connect(leaves_inh, sd_inh)
    
def main_lat_connections(l_exc_i,l_exc_j,l_inh_i,l_inh_j,orientation_i,orientation_j): 
    tp.ConnectLayers(l_exc_i, l_exc_j, create_lat_exc('PlosOne_J',kappa_j,orientation_i,orientation_j))
    tp.ConnectLayers(l_exc_i, l_inh_j, create_lat_exc('PlosOne_W',kappa_w,orientation_i,orientation_j))
    tp.ConnectLayers(l_inh_i, l_exc_j, dict_lat_inh)
    tp.ConnectLayers(l_inh_i, l_inh_j, dict_lat_inh)  
    
############################################################## Data ###############################################################

def save_dict(to_save,name_to_save):
    to_record = open( name_to_save + ".pkl", "wb")
    pickle.dump(to_save, to_record)
    to_record.close()

def load_dict(name):
    to_load = open(name +str('.pkl'), "rb")
    data = pickle.load(to_load)
    to_load.close()
    return data








