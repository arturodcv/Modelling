#Paths
sd_path = 'spk_detectors_folder' 
df_folder = 'dataframes_folder'
plots_path = 'plots_folder' 
results_path = 'results_folder'
positions_path = 'positions_folder'
input_images_path = 'input_images_folder' 

#Gabor
max_normalized_value = 100.0
K_size=11
Lambda=3
Sigma=1
Gamma=1.2
Psi = 0

#Size
x_cortex_size = 15
y_cortex_size = 15
cortex_size = x_cortex_size * y_cortex_size

#Nest
local_num_threads = 2
total_num_virtual_procs = 2
resolution = 0.1

#Layers
extent = [1.0,1.0]
ratio_exc_inh = 4
neurons_per_column_inh = 6
neurons_per_column_exc = ratio_exc_inh * neurons_per_column_inh
poisson_bias = 5.0
#poisson_bias = 0.0

#Poisson
#factor = 1
factor = 430

#Dictionaries

input_weight_exc = 1.0
input_weight_inh = input_weight_exc * 1
ratio_inh_exc_w = 4
weight_exc = 1.0 
weight_inh = - ratio_inh_exc_w * weight_exc
delay_exc_min = 2.0
slowness_exc = 5.0
delay_inh_min = 0.5
slowness_inh = 0.01 * slowness_exc
crf_conn_factor = 2.5
stddev_c_rf = 0.08
stddev_lat_conn_inh = 3 * stddev_c_rf
n_sigmas_inside_mask_4cbeta = 2.25
stddev_lat_conn_inh_4cbeta = 3 * stddev_c_rf
n_microcolumn_height = 10
n_microcolumn_width = 10
n_rows_latconn_inh_4cbeta = int(stddev_lat_conn_inh_4cbeta * n_sigmas_inside_mask_4cbeta * 2 * n_microcolumn_height) + 1
n_cols_latconn_inh_4cbeta = int(stddev_lat_conn_inh_4cbeta * n_sigmas_inside_mask_4cbeta * 2 * n_microcolumn_width) + 1
syn_model_inh = 'static_synapse_hpc'
syn_model_exc = 'static_synapse_hpc'
allow_autapses = True
allow_multapses = False
kappa_j = 1.0
kappa_w = 1.0

dict_poiss_to_exc  = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},'weights': input_weight_exc, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_poiss_to_inh = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},'weights': input_weight_inh, 
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_divergent = {'connection_type': 'divergent','mask': {'grid': {'rows': 1, 'columns': 1}},
                      'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}

dict_lat_inh =  {'connection_type': 'convergent',
                 'mask': {'grid':{'rows':n_rows_latconn_inh_4cbeta,'columns':n_cols_latconn_inh_4cbeta}, 
                          'anchor':{'row':(n_rows_latconn_inh_4cbeta-1)//2,'column':(n_cols_latconn_inh_4cbeta-1)//2}},
                 'delays': {'linear':{'c':delay_inh_min,'a':slowness_inh}},
                 'kernel': {'gaussian2D':{'p_center':1., 'sigma_x':stddev_lat_conn_inh, 'sigma_y':stddev_lat_conn_inh}},
                 'weights': 0.2 * weight_inh , 
                 'synapse_model':syn_model_inh,
                 'allow_autapses': allow_autapses, 'allow_multapses': allow_multapses}
                 

#Modelos 
RS_dict =  {'a':0.02, 'b':0.2, 'c':-65.,'d':8.0, 'V_th':30.}
FS_dict =  {'a':0.1, 'b':0.2, 'c':-65., 'd':2.0, 'V_th':30.}

#Simulation
simulation_time = 1000.0

#Save dictionary
def save_dict(to_save,name_to_save):
    to_record = open( name_to_save + ".pkl", "wb")
    pickle.dump(to_save, to_record)
    to_record.close()

#Video and image
window_time = 0
re_size = (200,200)
frames_per_second = 15
#COLORMAP_HOT #COLORMAP_COOL #COLORMAP_JET
#colormap = cv2.COLORMAP_JET