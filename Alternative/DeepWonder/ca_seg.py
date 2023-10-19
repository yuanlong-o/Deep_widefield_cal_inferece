# %%

import numpy as np
import time
import os
import tifffile as tiff
from skimage import io
import matplotlib.pyplot as plt

from scipy.io import savemat, loadmat

# load modules
from PostProcessing.par3 import fastthreshold
from PostProcessing.combine import segs_results, unique_neurons2_simp, \
    group_neurons, piece_neurons_IOU, piece_neurons_consume
from PostProcessing.complete_post import complete_segment

from caiman_plots.visualization import nb_view_patches

import argparse


experiment_path = 'Experiments/'

def convert_mask(mask):
    row_sum = np.sum(mask, axis = 0)
    col_sum = np.sum(mask, axis = 1)
    row_mask = np.where(row_sum > 0)
    col_mask = np.where(col_sum > 0)
    roi = [np.min(col_mask), np.max(col_mask)+1, np.min(row_mask), np.max(row_mask)+1] # notice +1!
    mask_valid = mask[roi[0]:roi[1],roi[2]:roi[3]]
    return mask_valid, roi

# %%
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--path','-p','-path of the recording')
    args = parser.parse_args()
    path = args.path

    display = True
    start = time.time()
    # load prob_map
    im_folder = path
    input_dir = im_folder+ '/process/RMBG/RMBG.tif'
    output_folder = im_folder + '/process/SEG/'

    prob_map = tiff.imread(input_dir)
    prob_map = prob_map.astype(np.float32)
    prob_map = (prob_map-prob_map.min()).astype(np.float32) / prob_map.max()
    nframes = prob_map.shape[0]
    Lx =  prob_map.shape[1]
    Ly =  prob_map.shape[2]

    pixel_size = 2

    # set parameters for post-processing
    # Optimization_Info = loadmat(os.path.join('.', 'Optimization_Info_{}.mat'.format(0)))
    # Params_post_mat = Optimization_Info['Params'][0]
    Params_post={
            # minimum area of a neuron (unit: pixels).
            'minArea': 15, 
            # average area of a typical neuron (unit: pixels) 
            'avgArea': 150,
            # uint8 threshould of probablity map (uint8 variable, = float probablity * 256 - 1)
            'thresh_pmap':0.1 * 256, 
            # values higher than "thresh_mask" times the maximum value of the mask are set to one.
            'thresh_mask': 0.5, 
            # maximum COM distance of two masks to be considered the same neuron in the initial merging (unit: pixels)
            'thresh_COM0': 4 / pixel_size, 
            # maximum COM distance of two masks to be considered the same neuron (unit: pixels)
            'thresh_COM': 8 / pixel_size, 
            # minimum IoU of two masks to be considered the same neuron
            'thresh_IOU': 0.5, 
            # minimum consume ratio of two masks to be considered the same neuron
            'thresh_consume': 0.75, 
            # minimum consecutive number of frames of active neurons
            'cons':3}

    # post-processing
    print(Params_post)
    Params_post_copy = Params_post.copy()
    Params_post_copy['thresh_pmap'] = None # Avoid repeated thresholding in postprocessing

    pmaps_b = np.zeros(prob_map.shape, dtype='uint8')
    thresh_pmap_float = (Params_post['thresh_pmap']+1)/256 # for published version

    # threshold the probability map to binary activity
    fastthreshold(prob_map, pmaps_b, thresh_pmap_float)
    # io.imsave(im_folder + '\\' + im_name + '_seg' + 'threshold.tiff', pmaps_b, check_contrast=False)

    # the rest of post-processing. The result is a 2D sparse matrix of the segmented neurons
    p = None
    useWT = False
    useMP = False
    Masks_2 = complete_segment(pmaps_b, Params_post_copy, useMP=useMP, p=p, useWT=useWT, display=display)
    if display:
        finish = time.time()
        time_post = finish-start
        time_frame_post = time_post/nframes*1000
        print('Processing Done: {:6f} s, {:6f} ms/frame'.format(time_post, time_frame_post))

    # convert to a 3D array of the segmented neurons
    Masks = np.reshape(Masks_2.toarray(), (Masks_2.shape[0], Lx, Ly)).astype('bool')

    # save the mask
    output_img = Masks
    # output_img = output_img1[0:raw_noise_img.shape[0],0:raw_noise_img.shape[1],0:raw_noise_img.shape[2]]


    output_img = output_img.astype('uint16')
    result_name = output_folder + '/SEG.tiff'
    io.imsave(result_name, output_img, check_contrast=False)

    result_name = output_folder + '/SEG_SUM.png'
    mask_sum = np.sum(output_img, axis = 0).astype('uint8')
    mask_sum = mask_sum * int(255/np.max(mask_sum))
    io.imsave(result_name, mask_sum)

    if display:
        finish = time.time()
        time_post = finish-start
        print('Mask saved: {:6f} s'.format(time_post))
    


# %% further show the temporal information
    N = output_img.shape[0] # component number
    print('Final Neuron number: ', N)

    C = np.zeros((N, nframes), dtype = 'float32')
    for i in range(N):
        # read the target positions
        mask_valid, roi = convert_mask(output_img[i])
        curr_temporal_signal = np.squeeze(np.sum(prob_map[:,roi[0]:roi[1],roi[2]:roi[3]] * mask_valid.astype(np.float32), axis=(1, 2))) # broadcast here.
        C[i] = curr_temporal_signal

    if display:
        finish = time.time()
        time_post = finish-start
        print('Temporal signal extracted: {:6f} s'.format(time_post))

    # save the time 
    savemat(output_folder + '/infer_results.mat', {'C': C})
    
    if display:
        finish = time.time()
        time_post = finish-start
        print('Temporal signal saved: {:6f} s'.format(time_post))
    
    exit()
    


#%% plot the temporal signals
    YrA = np.reshape(prob_map, (nframes, Lx * Ly))
    YrA = np.transpose(YrA)
    A = np.reshape(Masks, (N, Lx* Ly)).astype('float32')
    A = np.transpose(A)
    b = np.zeros((Lx * Ly, 1), dtype = 'float32')
    f = np.zeros((1, nframes), dtype = 'float32')

    print('Plotting Traces.....................')
    output_file_path = output_folder + '/Traces.html'
    # fig = plt.ion() # pop up a window
    nb_view_patches(YrA, A, C, b, f,  # values from CNMF
                    Lx, Ly, thr=0.9,output_file_path = output_file_path) # dimensions    # how can I save that?
    print(f'Done. Trace saved to {output_file_path}')
    if display:
        finish = time.time()
        time_post = finish-start
        print('All task completed: {:6f} s'.format(time_post))
