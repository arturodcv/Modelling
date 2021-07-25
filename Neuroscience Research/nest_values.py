#Paths
sd_path = 'spk_detectors_folder' 
df_folder = 'dataframes_folder'
plots_path = 'plots_folder' 
results_path = 'results_folder'
positions_path = 'positions_folder'
input_images_path = 'input_images_folder' 
gabor_folder = 'gabor_outputs'

#Gabor
K_size = 300; Lambda = 75 ;Psi = 0 ;Sigma = 47; Gamma = 0.4
cut_pixels = 0
get_output_gabors = False

#Simulation 
images_to_simulate = [input_images_path + '/small.png' ] # big_5 small_5 big_9 small_9 big small test

 
num_images_to_simulate = len(images_to_simulate)
ms_per_stimuli = 1000.0
simulation_time = ms_per_stimuli * num_images_to_simulate 



#Size
num_hipercolumns = 7 
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
poisson_bias = 2.0 




#Poisson
factor_exc = 109.0
factor_inh = 20.4 
factor_inh_poiss = 4.78 * factor_inh
factor_exc_poiss = 97.5

############ Dictionaries

#Layers
edge_wrap = True 
allow_autapses = True
allow_multapses = False
  

#Kernel
stddev_c_rf = 0.08 

p_center_inh = 1.0
mean_lat_conn_inh = 3 * stddev_c_rf
stddev_lat_conn_inh = 1 * stddev_c_rf 
input_stddev_inh = stddev_c_rf ## 0.0

p_center_exc = p_center_inh * 1
mean_lat_conn_exc = 0.0
stddev_lat_conn_exc = stddev_lat_conn_inh *  2
input_stddev_exc = stddev_c_rf

n_sigmas_inside_mask = 4
n_microcolumn_height = 10 
n_microcolumn_width = 10
n_rows_latconn_inh = int(stddev_lat_conn_inh * n_sigmas_inside_mask * 2 * n_microcolumn_height) + 1 
n_cols_latconn_inh = int(stddev_lat_conn_inh * n_sigmas_inside_mask * 2 * n_microcolumn_width)  + 1 
n_rows_latconn_exc = int(stddev_lat_conn_exc * n_sigmas_inside_mask * 2 * n_microcolumn_height) + 1
n_cols_latconn_exc = int(stddev_lat_conn_exc * n_sigmas_inside_mask * 2 * n_microcolumn_width)  + 1  


#Synapse model
syn_model_inh = 'static_synapse_hpc'
syn_model_exc = 'static_synapse_hpc'

#Delays
delay_exc_min = 1.0
delay_exc_min_large = 1.0
delay_inh_min = 0.5 

slowness_exc = 0.0
slowness_exc_large = 0.0
slowness_inh = 0.01 * slowness_exc # 0.01


# PlosOne
lateral_connections = False

kappa_j = 0.15
kappa_w = 1.0
rescale = 2.0 
radius_lat = 3.5 
weight_large_range_exc_exc = 0.01 
weight_large_range_exc_inh = 0.01 

input_p_center_inh = 1.0 * 1.0 ;  input_weight_inh = 1.0
input_p_center_exc = 1.0 * 1.0 ;  input_weight_exc = 1.0


dict_poiss_to_inh  = {'connection_type': 'divergent','weights': input_weight_inh ,
                      'kernel': {'gaussian':{'p_center': input_p_center_inh, 'sigma':input_stddev_inh , 'mean': 0.0}}, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_poiss_to_exc = {'connection_type': 'divergent','weights': input_weight_exc ,
                      'kernel': {'gaussian':{'p_center': input_p_center_exc , 'sigma':input_stddev_exc , 'mean': 0.0}}, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}



dict_divergent    = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}


ratio_inh_exc_w = 4.0 
p_center_inh_exc = 1.0 / 1 ; weight_inh_exc = - ratio_inh_exc_w * .7
p_center_inh_inh = 1.0 / 6 ; weight_inh_inh = - ratio_inh_exc_w * 0.7
p_center_exc_exc = 1.0 / 4 ; weight_exc_exc = 1.0 * 0.7
p_center_exc_inh = 1.0 / 2 ; weight_exc_inh = 1.0 * 0.7

short_range_inh_exc =  {'connection_type': 'convergent',
                        'mask': {'grid':{'rows':n_rows_latconn_inh ,'columns':n_cols_latconn_inh }, 
                                 'anchor':{'row':(n_rows_latconn_inh - 1)//2,'column':(n_cols_latconn_inh - 1)//2}},
                        'delays': {'linear':{'c':delay_inh_min,'a':slowness_inh}},
                        'kernel': {'gaussian':{'p_center': p_center_inh_exc , 'sigma':stddev_lat_conn_inh , 'mean': mean_lat_conn_inh}},
                        'weights': weight_inh_exc ,
                        'synapse_model':syn_model_inh,
                        'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}
                       
short_range_inh_inh =  {'connection_type': 'convergent',
                        'mask': {'grid':{'rows':n_rows_latconn_inh,'columns':n_cols_latconn_inh}, 
                                 'anchor':{'row':(n_rows_latconn_inh - 1)//2,'column':(n_cols_latconn_inh - 1)//2}},
                        'delays': {'linear':{'c':delay_inh_min,'a':slowness_inh}},
                        'kernel': {'gaussian':{'p_center': p_center_inh_inh  , 'sigma':stddev_lat_conn_inh , 'mean': mean_lat_conn_inh}}, 
                        'weights': weight_inh_inh ,
                        'synapse_model':syn_model_inh,
                        'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}
                       
short_range_exc_exc =  {'connection_type': 'convergent',
                        'mask': {'grid':{'rows':n_rows_latconn_exc,'columns':n_cols_latconn_exc}, 
                                 'anchor':{'row':(n_rows_latconn_exc - 1)//2,'column':(n_cols_latconn_exc - 1)//2}},
                        'delays': {'linear':{'c':delay_exc_min,'a':slowness_exc}},
                        'kernel': {'gaussian':{'p_center': p_center_exc_exc, 'sigma':stddev_lat_conn_exc , 'mean':mean_lat_conn_exc}}, 
                        'weights':  weight_exc_exc  , 
                        'synapse_model':syn_model_exc,
                        'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}
                       
short_range_exc_inh =  {'connection_type': 'convergent',
                        'mask': {'grid':{'rows':n_rows_latconn_exc,'columns':n_cols_latconn_exc}, 
                                 'anchor':{'row':(n_rows_latconn_exc - 1)//2,'column':(n_cols_latconn_exc - 1)//2}},
                        'delays': {'linear':{'c':delay_exc_min,'a':slowness_exc}},
                        'kernel': {'gaussian':{'p_center': p_center_exc_inh  , 'sigma':stddev_lat_conn_exc, 'mean':mean_lat_conn_exc}}, 
                        'weights':  weight_exc_inh , 
                        'synapse_model':syn_model_exc,
                        'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}
                       

#Modelos 
RS_dict =  {'a':0.02, 'b':0.2, 'c':-65.,'d':8.0, 'V_th':30.}
FS_dict =  {'a':0.1, 'b':0.2, 'c':-65., 'd':2.0, 'V_th':30.}



#Results (image, video, EEG, frequencies)
window_time = 0
re_size = (x_cortex_size,y_cortex_size)
frames_per_second = 20
num_max_frequencies = 2
broadband_initial = 30
image_from = 100
eeg_from = 100
radius = columns_in_hipercolumns / 2 * 1

