#Paths
sd_path = 'spk_detectors_folder' 
df_folder = 'dataframes_folder'
plots_path = 'plots_folder' 
results_path = 'results_folder'
positions_path = 'positions_folder'

#Gabor
max_normalized_value = 100.0
K_size=11
Lambda=3
Sigma=1
Gamma=1.2
Psi = 0

#Size
x_cortex_size = 50
y_cortex_size = 50
cortex_size = x_cortex_size * y_cortex_size

#Nest
local_num_threads = 32
total_num_virtual_procs = 32
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

#Simulation
simulation_time = 1000.0

#Video and image
window_time = 0
re_size = (200,200)
frames_per_second = 15
#COLORMAP_HOT #COLORMAP_COOL #COLORMAP_JET
#colormap = cv2.COLORMAP_JET