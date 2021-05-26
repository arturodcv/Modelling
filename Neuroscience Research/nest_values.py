#Paths
sd_path = 'spk_detectors_folder' 
df_folder = 'dataframes_folder'
plots_path = 'plots_folder' 
results_path = 'results_folder'
positions_path = 'positions_folder'
input_images_path = 'input_images_folder' 
gabor_folder = 'gabor_outputs'

#Gabor
#K_size = 500; Lambda = 120;Gamma = 0.1 ;Psi = 0 ;Sigma = 18 # Para el vertical
K_size = 500; Lambda = 100;Gamma = 1.6 ;Psi = 0 ;Sigma = 65.5 # Para las gaussianas en linea horizontal
#K_size = 100; Lambda = 50;Gamma = 0.9 ;Psi = 3.14159265/2 ;Sigma = 25.25832176208496 # Para la gabor sola

cut_pixels = 0
max_normalized_value = 100

#Simulation 
#images_to_simulate = [input_images_path + '/sinusoide_' + str(int(i * 180 / 2) ) + '.png' for i in range(0,2)]
#images_to_simulate = [input_images_path + '/gaussian.png' ]
#images_to_simulate = [input_images_path + '/sinusoide_0_freq.png' ]
#images_to_simulate = [input_images_path + '/double_vertical_gauss.png' ]
#images_to_simulate = [input_images_path + '/sinusoid_vertical_freq_half.png', input_images_path + '/sinusoide_0_freq.png' ]
images_to_simulate = [input_images_path + '/pixil-frame-0.png' ] 

num_images_to_simulate = len(images_to_simulate)
ms_per_stimuli = 1000.0
simulation_time = ms_per_stimuli * num_images_to_simulate 

#Size
num_hipercolumns = 5
columns_in_hipercolumns = 10
x_cortex_size = num_hipercolumns * columns_in_hipercolumns 
y_cortex_size = num_hipercolumns * columns_in_hipercolumns
cortex_size = x_cortex_size * y_cortex_size

#Nest
local_num_threads = 2
resolution = 0.1

#Number of orientations
num_orientations = 4

#Layers
extent = [float(num_hipercolumns), float(num_hipercolumns)]
ratio_exc_inh = 4
neurons_per_column_inh = 6
neurons_per_column_exc = ratio_exc_inh * neurons_per_column_inh
poisson_bias = 5.0


#Poisson
factor = 430


############ Dictionaries

#Layers
edge_wrap = True 
allow_autapses = True
allow_multapses = False
input_weight_exc = 1.0
input_weight_inh = input_weight_exc * 1

#Weights
ratio_inh_exc_w = 4.0
default_synapse_weight_exc = 1.0 
default_synapse_weight_inh = - ratio_inh_exc_w * default_synapse_weight_exc

#Delays
delay_exc_min = 2.0
delay_inh_min = 0.5
slowness_exc = 5.0
slowness_inh = 0.01 * slowness_exc

#Kernel
stddev_c_rf = 0.08
factor_exc_inh_2 = 1.0

p_center_inh = 1.0
mean_lat_conn_inh_4cbeta = 3 * stddev_c_rf
stddev_lat_conn_inh_4cbeta = 1 * stddev_c_rf 

p_center_exc = p_center_inh / factor_exc_inh_2
mean_lat_conn_exc_4cbeta = 0.0
stddev_lat_conn_exc_4cbeta = stddev_lat_conn_inh_4cbeta / 2.0

n_sigmas_inside_mask_4cbeta = 2.25
n_microcolumn_height = 10
n_microcolumn_width = 10
n_rows_latconn_inh_4cbeta = int(stddev_lat_conn_inh_4cbeta * n_sigmas_inside_mask_4cbeta * 2 * n_microcolumn_height) + 1
n_cols_latconn_inh_4cbeta = int(stddev_lat_conn_inh_4cbeta * n_sigmas_inside_mask_4cbeta * 2 * n_microcolumn_width) + 1
n_rows_latconn_exc_4cbeta = int(stddev_lat_conn_exc_4cbeta * n_sigmas_inside_mask_4cbeta * 2 * n_microcolumn_height) + 1
n_cols_latconn_exc_4cbeta = int(stddev_lat_conn_exc_4cbeta * n_sigmas_inside_mask_4cbeta * 2 * n_microcolumn_height) + 1

#Synapse model
syn_model_inh = 'static_synapse_hpc'
syn_model_exc = 'static_synapse_hpc'

# PlosOne
kappa_j = 0.126 * 7.14
kappa_w = 0.14  * 7.14
rescale = 5.0

dict_poiss_to_exc  = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},'weights': input_weight_exc, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_poiss_to_inh  = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},'weights': input_weight_inh, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_divergent     = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

short_range_inh       =  {'connection_type': 'convergent',
                       'mask': {'grid':{'rows':n_rows_latconn_inh_4cbeta,'columns':n_cols_latconn_inh_4cbeta}, 
                                'anchor':{'row':(n_rows_latconn_inh_4cbeta-1)//2,'column':(n_cols_latconn_inh_4cbeta-1)//2}},
                       'delays': {'linear':{'c':delay_inh_min,'a':slowness_inh}},
                       'kernel': {'gaussian':{'p_center': p_center_inh, 'sigma':stddev_lat_conn_inh_4cbeta , 'mean': mean_lat_conn_inh_4cbeta}},
                       'weights': 0.1 * 0.2 * default_synapse_weight_inh ,  
                       'synapse_model':syn_model_inh,
                       'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}
                       
short_range_exc    =  {'connection_type': 'convergent',
                       'mask': {'grid':{'rows':n_rows_latconn_exc_4cbeta,'columns':n_cols_latconn_exc_4cbeta}, 
                                'anchor':{'row':(n_rows_latconn_exc_4cbeta-1)//2,'column':(n_cols_latconn_exc_4cbeta-1)//2}},
                       'delays': {'linear':{'c':delay_exc_min,'a':slowness_exc}},
                       'kernel': {'gaussian':{'p_center': p_center_exc, 'sigma':stddev_lat_conn_exc_4cbeta , 'mean':mean_lat_conn_exc_4cbeta}},
                       'weights': 0.1 * 0.2 * default_synapse_weight_exc , ###!!!!!
                       'synapse_model':syn_model_exc,
                       'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}
                 

#Modelos 
RS_dict =  {'a':0.02, 'b':0.2, 'c':-65.,'d':8.0, 'V_th':30.}
FS_dict =  {'a':0.1, 'b':0.2, 'c':-65., 'd':2.0, 'V_th':30.}





#Save dictionary
def save_dict(to_save,name_to_save):
    to_record = open( name_to_save + ".pkl", "wb")
    pickle.dump(to_save, to_record)
    to_record.close()

#Video and image
window_time = 0
re_size = (x_cortex_size,y_cortex_size)
frames_per_second = 20
