import numpy as np
import pandas as pd
import os
import time
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

#Some folder functions for our purposes
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

#This is the directory where the frames will be saved. If you want to change the directory, rewrite the path variable.
path = 'plots_folder'
create_folder(path)


#Read and fix the data
data = pd.read_table('output.txt', sep = ';', header = None)
data = data.dropna(axis = 1)
data = data.to_numpy()

#Generate and save the frames for the video
x_points = [k for k in range(0,len(data[0]))] 

for i in range(0,len(data)):
    plt.xlim(x_points[0] - 1, x_points[-1] + 1)
    plt.ylim(np.min(data) -1 , np.max(data) + 1)
    plt.scatter(x_points,data[i])
    plt.savefig(path + '/plot_' + str(i) + '.png')
    plt.clf()

#Read and release the video
img = [] 
frames = [] 
fig = plt.figure()
plt.axis('off')
for i in range(len(data)):
    img.append(plt.imread(path + '/plot_' + str(i) + '.png'))
    frames.append([plt.imshow(img[i], cmap=cm.Greys_r,animated=True)])
ani = animation.ArtistAnimation(fig, frames, interval= 150, blit=True,
                                repeat_delay=1000)
ani.save('results.mp4')


#Remove the created files and the directory
remove_contents(path)
os.rmdir(path)

print("Program finished! Reproduce the 'results.mp4' video to see the results\n")





