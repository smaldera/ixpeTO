import pandas as pd
import numpy as np
import h5py as h5

from tqdm import tqdm

import matplotlib.pyplot as plt

s1 = list()
s2 = list()

with open("simIxpe_5Kev.txt", "r") as f:
    lines = f.readlines()
    print(len(lines))
    for i in tqdm(range(0,len(lines))):
        if (len(lines[i]) > 100):
            s1.append(i)

for i in tqdm(range(0,len(lines))):
    s2.append(i)
    
for i in tqdm((s1)):
    s2.remove(i)
    
pixels = pd.read_csv("simIxpe_5Kev.txt", skiprows=s1, header=None, delimiter = "  ", names = ["x","y","int"])
labels = pd.read_csv("simIxpe_5Kev.txt", skiprows=s2, header=None, delimiter = " ", names = ["garbage1","ID","garbage2","nPixels","garbage3","Energy","garbage4","x","garbage5","y","garbage6","phi","garbage7","theta"])

number_of_all_pixels = int(pixels.size/3)
number_of_all_events = int(labels.size/14)

pixels.x *= 40
pixels.y *= 40/3.464

x_min = 0
y_min = 0

x_max = 0
y_max = 0

for i in tqdm(range(0,number_of_all_pixels)):
    if pixels.x.values[i] < x_min:
        x_min = pixels.x[i]
    if pixels.y.values[i] < y_min:
        y_min = pixels.y[i]

pixels.x -= x_min
pixels.y -= y_min

for i in tqdm(range(0,number_of_all_pixels)):
    if pixels.x.values[i] > x_max:
        x_max = pixels.x[i]
    if pixels.y.values[i] + 0.5*pixels.x.values[i] > y_max:
        y_max = pixels.y[i] + 0.5*pixels.x.values[i]

sum_pixels = 0
num_pixels = 0

max_dif_x = 0
max_dif_y = 0

picture_size_list = np.zeros((number_of_all_events, 4), dtype = np.float16)

for i in tqdm(range(0,int(number_of_all_events))):
    num_pixels = labels.nPixels[i]
    sum_pixels += num_pixels
    temp_pixels = pixels[sum_pixels-num_pixels:sum_pixels]
    temp_x_min = 500  
    temp_x_max = 0    
    temp_y_min = 500  
    temp_y_max = 0    
    for k in range(0,int(temp_pixels.size/3)):
        x = int(temp_pixels.x[sum_pixels-num_pixels+k])
        y = int(temp_pixels.y[sum_pixels-num_pixels+k]+0.5*temp_pixels.x[sum_pixels-num_pixels+k])
        if (x < temp_x_min):
            temp_x_min = x
        if (x > temp_x_max):
            temp_x_max = x
        if (y < temp_y_min):
            temp_y_min = y
        if (y > temp_y_max):
            temp_y_max = y
    size_x = temp_x_max - temp_x_min + 1
    size_y = temp_y_max - temp_y_min + 1
    picture_size_list[i,0] = size_x
    picture_size_list[i,1] = size_y
    picture_size_list[i,2] = temp_x_min
    picture_size_list[i,3] = temp_y_min
    if size_x > max_dif_x:
        max_dif_x = size_x
    if size_y > max_dif_y:
        max_dif_y = size_y
        
sum_pixels = 0
num_pixels = 0

event_list = np.zeros((number_of_all_events, max_dif_x, max_dif_x, 1), dtype = np.float16)
label_list = np.zeros((number_of_all_events, 7))

for i in tqdm(range(0,int(number_of_all_events))):
    num_pixels = labels.nPixels[i]
    sum_pixels += num_pixels
    temp_pixels = pixels[sum_pixels-num_pixels:sum_pixels]
    for k in range(0,int(temp_pixels.size/3)):
        event_list[i, int(temp_pixels.x[sum_pixels-num_pixels+k]-picture_size_list[i,2] + (max_dif_x-picture_size_list[i,0])/2),int(temp_pixels.y[sum_pixels-num_pixels+k]+0.5*temp_pixels.x[sum_pixels-num_pixels+k]-picture_size_list[i,3] + (max_dif_x-picture_size_list[i,1])/2), 0] = temp_pixels.int[sum_pixels-num_pixels+k]
    
for i in tqdm(range(0,int(number_of_all_events))):
    label_list[i,0] = labels.ID[i]
    label_list[i,1] = labels.nPixels[i]
    label_list[i,2] = labels.Energy[i]
    label_list[i,3] = labels.x[i]
    label_list[i,4] = labels.y[i]
    label_list[i,5] = labels.phi[i]
    label_list[i,6] = labels.theta[i]

outfile = h5.File("pixels.hdf5", "w")
outfile.create_dataset("data", data=event_list)
outfile.close()

outfile = h5.File("labels.hdf5", "w")
outfile.create_dataset("data", data=label_list)
outfile.close()

for i in range(0,10):
    plt.imshow(event_list[100*i,:,:,0], interpolation="nearest")
    plt.colorbar()
    plt.show()