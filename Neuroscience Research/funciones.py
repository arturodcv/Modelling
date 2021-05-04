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
from collections import Counter
from collections import OrderedDict
from scipy.fft import rfft
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

def create_sin2d(size,freq,a,b):
    array = np.zeros((size,size))
    for j in range(0,size):
        for i in range(0,size):
            array[i][j] =  math.sin( 2*math.pi*freq*(a*i+b*j) )
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
    gabor = gabor_filter(K_size = K_size, Lambda = Lambda, Theta = Theta, Sigma = Sigma, Gamma = Gamma, Psi = Psi, mode = mode)
    for x in range(gray_img.shape[0]):
        for y in range(gray_img.shape[1]):
            output[x, y] = np.sum(gray[x : x + K_size, y : y + K_size] * gabor) 
    return output

    
def gabor(gray_img,orientation_in_radians, mode):
    output = np.zeros((gray_img.shape[0],gray_img.shape[1]), dtype=np.float32) 
    #orientation = (90 - orientation_in_radians)*math.pi/180
    orientation = orientation_in_radians*math.pi/180
    output = apply_filter(gray_img, K_size=K_size, Lambda=Lambda, Theta=orientation, Sigma=Sigma, Gamma=Gamma,Psi = Psi, mode = mode )
    output = np.clip(output, 0, np.max(output))
    return output

########################################################### Nest ###################################################################

def input_treatment(input_spike,x_cortex_size,y_cortex_size):
    input_as_img = Image.fromarray(input_spike)
    input_resized = input_as_img.resize((x_cortex_size,y_cortex_size), resample = Image.NEAREST  )
    input_norm = np.divide(input_resized,1) * 1 ## cambiar
    input_transposed = input_norm.transpose()
    input_as_list = input_transposed.tolist()
    flat_list = [item for sublist in input_as_list for item in sublist]
    return flat_list
    
def main_img(img,orientation):
    img = cv2.imread(img)
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) 
    output_gabor = gabor(gray_img,orientation,'on')
    output_gabor = output_gabor[cut_pixels:-cut_pixels,cut_pixels:-cut_pixels]
    flat_list = input_treatment(output_gabor,x_cortex_size,y_cortex_size)
    return flat_list
    
def main_create():
    l_exc  = create_layer(x_cortex_size,y_cortex_size,extent,'exc',neurons_per_column_exc)
    l_inh  = create_layer(x_cortex_size,y_cortex_size,extent,'inh',neurons_per_column_inh)
    sd_exc = nest.Create('spike_detector', params = {'to_file': True, 'to_memory': False})
    sd_inh = nest.Create('spike_detector', params = {'to_file': True, 'to_memory': False})
    l_poiss = create_layer(x_cortex_size,y_cortex_size,extent,'poisson_generator',1)
    return l_exc,l_inh,sd_exc,sd_inh,l_poiss
    
def main_one_orientation(orientation):
    orientation = orientation*math.pi/180
    l_exc,l_inh,sd_exc,sd_inh,l_poiss = main_create()
    main_self_connections(l_exc,l_inh,sd_exc,sd_inh,l_poiss ,orientation)
    orientation = orientation*180 / math.pi
    out_dict = {'l_exc_'+str(orientation): l_exc,'l_inh_'+str(orientation): l_inh,'sd_exc_'+str(orientation): sd_exc,'sd_inh_'+str(orientation): sd_inh}
    poiss_dict = {'l_poiss_'+str(orientation): l_poiss}
    return out_dict, poiss_dict



def main_all_orientations(num_orientations):
    lyrs = {}
    poiss = {}
    for i in range(0,num_orientations):
        #orientation = 90 - i*180/num_orientations ###
        orientation =  i*180/num_orientations 
        lyrs["orientation_" + str(orientation)], poiss["orientation"+str(orientation)] = main_one_orientation(orientation)
    for i in range(0,num_orientations):
        for j in range(0,num_orientations):
            if i!=j:
                #or_i = 90 - i*180/num_orientations; ##
                or_i = i*180/num_orientations;
                #or_j = 90 - j*180/num_orientations;  ##
                or_j = j*180/num_orientations;
                l_exc_i = lyrs['orientation_'+str(or_i)]['l_exc_'+str(or_i)]
                l_exc_j = lyrs['orientation_'+str(or_j)]['l_exc_'+str(or_j)]
                l_inh_i = lyrs['orientation_'+str(or_i)]['l_inh_'+str(or_i)]
                l_inh_j = lyrs['orientation_'+str(or_j)]['l_inh_'+str(or_j)]
                

                #main_lat_connections(l_exc_i,l_exc_j,l_inh_i,l_inh_j,math.pi/2 - i*math.pi/num_orientations,math.pi/2 - j*math.pi/num_orientations) 
                main_lat_connections(l_exc_i,l_exc_j,l_inh_i,l_inh_j,i*math.pi/num_orientations,j*math.pi/num_orientations)
    return lyrs, poiss



def set_poisson_values(img, poiss_layers,num_orientations):
    flat_lists = []
    for i in range(0,num_orientations):
        #orientation = 90 - i*180/num_orientations; ##
        orientation = i*180/num_orientations;
        flat_list = main_img(img,orientation)
        l_poiss = list(poiss_layers['orientation'+str(orientation)].values())[0]
        fixed_list = [k*factor / 41 + poisson_bias for k in flat_list] ################## cuidado con el 41
        nest.SetStatus(nest.GetNodes(l_poiss)[0],'rate', fixed_list)
    


################################################################# Layers #################################################################

def create_layer(rows,columns,extent,elements,neurons_per_column):
    return tp.CreateLayer({'rows': rows,'columns': columns,'extent': extent,'elements': [elements,neurons_per_column],'edge_wrap': edge_wrap}) 

########################################################### Connetivity ##################################################################

def create_lat_exc(kernel_type,kappa,orientation_i,orientation_j):
    return  {'connection_type': 'convergent',    
             'kernel': {kernel_type: {'max_dist': float(max_dist_plosone), 'kappa': kappa,
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
    #tp.ConnectLayers(l_inh_i, l_inh_j, dict_lat_inh) 
    
    
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
    
########################################################## Results ##################################################################

def read__and_fix_dataframe(orientation_to_read,exc_or_inh):
    data = pd.read_pickle(df_folder + '/data_l_' + exc_or_inh + '_' + str(orientation_to_read) + '.pkl')
    print("Dataframe loaded!")
    data = data.sort_values(by=['Time'])
    data = data.set_index([pd.Index([i for i in range(len(data))])])
    return data

def generate_frames(data):
    create_folder(plots_path)
    remove_contents(plots_path)

    times = data['Time'].tolist() 
    times = (np.around(data['Time'],window_time)).tolist()
    actual_time = times[0]
    x_pos = ( np.around( (data['x_pos'] + extent[0] / 2) * 10 - 0.5)).astype(int).tolist()
    y_pos = ( np.around( (data['y_pos'] + extent[0] / 2) * 10 - 0.5)).astype(int).tolist()

    num_spikes = len(data)
    array = [0] * x_cortex_size * y_cortex_size
    
    print('Generating frames ')
    for i in tqdm(range(0,num_spikes)): 
        if times[i] - actual_time != 0 or i == num_spikes:
            name = plots_path+'/plot_time_'+str(actual_time)+'.tiff'
            array = np.reshape(array,(x_cortex_size,y_cortex_size))
            img = Image.fromarray(np.uint8(array),'L').resize(re_size).transpose(Image.FLIP_TOP_BOTTOM)
            img.save(name,compress_level = 1)
            actual_time = times[i]
            array = [0] * x_cortex_size * y_cortex_size
        array[ x_pos[i] + y_cortex_size * y_pos[i] ] += 1

    name = plots_path+'/plot_time_'+str(actual_time)+'.tiff'
    array = np.reshape(array,(x_cortex_size,y_cortex_size))
    img = Image.fromarray(np.uint8(array),'L').resize(re_size).transpose(Image.FLIP_TOP_BOTTOM)
    img.save(name,compress_level = 1)
    
    return times
    
def generate_empty_frames(times):
    unique_times = np.unique(np.array(times)).tolist()
    full_times = [i for i in range(0,int(simulation_time) + 1)]
    complementary_time_list = list(set(full_times) - set(unique_times))

    array = [0] * x_cortex_size * y_cortex_size
    array = np.reshape(array,(x_cortex_size,y_cortex_size))
    img = Image.fromarray(np.uint8(array),'L').resize(re_size).transpose(Image.FLIP_TOP_BOTTOM)
    for i in tqdm(complementary_time_list):
        name = plots_path+'/plot_time_'+str(i)+'.tiff'
        img.save(name,compress_level = 1)
    frames = unique_times + complementary_time_list
    frames.sort()
    return frames, complementary_time_list
    
def read_frames(frames):
    img_array = []
    print('Reading frames')
    for milisecond in tqdm(frames[:-1]):
        filename = plots_path + '/plot_time_' + str(milisecond)+'.tiff'
        img = cv2.imread(filename,cv2.IMREAD_GRAYSCALE)
        img_array.append(img)
    return img_array
    
def create_video(img_array,orientation_to_read ,exc_or_inh, path):
    size = (img_array[0].shape[0],img_array[0].shape[1])
    video_out = cv2.VideoWriter(path + '/neurons_video_' + str(orientation_to_read) + '_' + str(exc_or_inh)+'_.avi',0, frames_per_second, size)

    max_array = np.max(img_array)
    print('Writing frames in video')
    for img in tqdm(img_array):
        img = (np.multiply(img, 255 / max_array)).astype(np.uint8)
        img = cv2.applyColorMap(img,cv2.COLORMAP_JET)
        video_out.write(img)
    video_out.release()
    print("Video created succesfully!")
    
def create_avg_img(img_array,orientation_to_read ,exc_or_inh, path ):
    print('Getting average image')
    img_sum = np.zeros(re_size)
    for img in tqdm(img_array):
        img_sum = img_sum + np.divide(img,len(img_array)*neurons_per_column_exc)

    plt.imshow(img_sum)
    plt.title('Average image')
    plt.colorbar()
    plt.imsave(path + '/Average_img_' + str(orientation_to_read) + '_' + str(exc_or_inh)+'_.png',img_sum)
    plt.close('all')
    print("Average image created succesfully!")
    
def get_eeg(times, complementary_time_list, orientation_to_read, exc_or_inh, path):
    eeg = Counter(times)
    eeg.update({i:0 for i in complementary_time_list})
    eeg = list(OrderedDict(sorted(eeg.items())).values())

    eeg_cut = eeg[100:1000]
    plt.plot(eeg_cut); plt.title('eeg_cut 100:500');
    plt.savefig(path + '/eeg_cut_' + str(orientation_to_read) + '_' + str(exc_or_inh)+'_.png')
    plt.close('all')

    plt.plot(eeg); plt.title('eeg');
    plt.savefig(path + '/eeg_' + str(orientation_to_read) + '_' + str(exc_or_inh)+'_.png')
    plt.close('all')

    print("Eeg created succesfully!")
    return eeg
    
def get_frequencies(eeg,orientation_to_read,exc_or_inh, path):
    frequencies = np.abs(rfft(np.asarray(eeg)))
    plt.plot(frequencies); plt.xlabel("Frequency (Hz)"); plt.ylabel("Frequency Domain (Spectrum) Magnitude")
    plt.grid(); plt.savefig(path + '/frequencies_' + str(orientation_to_read) + '_' + str(exc_or_inh)+'_.png')
    plt.close('all')

    print("Frequencies plot created succesfully!")
  









