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
from scipy.fft import fft, fftfreq
import pickle

from nest_values import *
from funciones   import *

################################################## Data Treatment ################################################################
create_folder(sd_path); #remove_contents(sd_path)
create_folder(df_folder); #remove_contents(df_folder)
create_folder(positions_path); #remove_contents(positions_path)

files = glob('*spike_detector*')
for file in files:
    shutil.move(file, sd_path + '/' + file)

files = glob('*positions-*')
for file in files:
    shutil.move(file, positions_path + '/' + file)

files = glob(positions_path + '/*')
pos = []
for file in files:
    positions = pd.read_table(file,names = ['Number','x_pos','y_pos'], index_col=False, sep = ' ')
    pos.append(positions)
  
positions = pd.concat(pos)
del(pos)

to_record_layer = open("to_record_layer.pkl", "rb")
layers_to_record = pickle.load(to_record_layer)
to_record_layer.close()

with open('to_record_sd.pkl', 'rb') as f:
      spike_detectors = pickle.load(f)

print('Fixing the data...')
t = time.time()
for layer,j in zip(layers_to_record,range(0,len(layers_to_record))):
    data = []
    for i in range(0,total_num_virtual_procs): 
        spk_number = str(i)
        if total_num_virtual_procs > 9 and i < 10:
            spk_number = str(0) + spk_number
        name = sd_path + '/spike_detector-' + str(spike_detectors[j][0]) + '-' + spk_number + '.gdf'
        df = pd.read_table(name,names = ['Number','Time'], index_col=False)
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


#################################################### Read data ###################################################################
#Read dataframe
#data = pd.read_pickle(df_folder + "/data_l_parrots.pkl")
data = pd.read_pickle(df_folder + "/data_l_exc_0_on.pkl")
#data = pd.read_pickle(df_folder + "/data_l_inh_0_on.pkl")

data = data.sort_values(by=['Time'])
data = data.set_index([pd.Index([i for i in range(len(data))])])

################################################### Generate frames ##############################################################

create_folder(plots_path)
remove_contents(plots_path)

array = [0] * x_cortex_size * y_cortex_size

times = data['Time'].tolist() 
times = (np.around(data['Time'],window_time)).tolist()
actual_time = times[0]
x_pos = (np.around(data['x_pos'] * x_cortex_size + x_cortex_size/2 - 0.5)).astype(int).tolist()
y_pos = (np.around(data['y_pos'] * y_cortex_size + x_cortex_size/2 - 0.5)).astype(int).tolist()
num_spikes = len(data)

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

del(data)

##################################################### Create video and images ######################################################
create_folder(results_path)

frames = np.unique(np.array(times)).tolist()
img_array = []

print('Reading frames')
for milisecond in tqdm(frames[:-1]):
    filename = plots_path + '/plot_time_' + str(milisecond)+'.tiff'
    img = cv2.imread(filename,cv2.IMREAD_GRAYSCALE)
    img_array.append(img)

size = (img.shape[0],img.shape[1])
video_out = cv2.VideoWriter(results_path + '/neurons_video.avi',0, frames_per_second, size)

max_array = np.max(img_array)

print('Writing frames in video')
for img in tqdm(img_array):
    img = (np.multiply(img, 255 / max_array)).astype(np.uint8)
    img = cv2.applyColorMap(img,cv2.COLORMAP_JET)
    video_out.write(img)
video_out.release()

#Average Image
print('Getting average image')
img_sum = np.zeros(re_size)
for img in tqdm(img_array):
    img_sum = img_sum + np.divide(img,len(img_array)*neurons_per_column_exc)

plt.imshow(img_sum)
plt.title('Average image')
plt.colorbar()
#plt.show()
#plt.close(img_sum)
plt.imsave(results_path + '/Average_img.png',img_sum)
plt.close('all')

############################################################## EEG ################################################################
counts = Counter(times)
i = 0; eeg = {}
while counts[i] < 1 and i < simulation_time:
    eeg.update({float(i):0})
    i += 1
eeg.update(counts)
eeg = list(eeg.values())

eeg_cut = eeg[0:400]
plt.plot(eeg_cut);
plt.title('eeg_cut')
#plt.close(eeg_cut)
plt.savefig(results_path + '/eeg_cut')

plt.title('eeg')
plt.plot(eeg);
plt.savefig(results_path + '/eeg')
plt.close('all')

############################################### Fourier transform and frequencies #####################################################

fourier_transform = fft(np.asarray(eeg))
N = len(fourier_transform)
T = 1.0 / simulation_time
x_points = fftfreq(N, T)[:N//2]

plt.plot(x_points[1:], np.abs(fourier_transform[1:N//2]) )
plt.xlabel("Frequency (Hz?)")
plt.ylabel("Frequency Domain (Spectrum) Magnitude")
plt.grid()
#plt.close(eeg)
plt.savefig(results_path + '/Frequencies')
plt.close('all')


