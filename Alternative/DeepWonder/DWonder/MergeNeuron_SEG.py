import tifffile as tiff
import cv2
import numpy as np
from skimage import io

import time
import datetime
import scipy.io as scio
import multiprocessing as mp
import math
import os

from DWonder.MN.utils import Split_Neuron, listAddtrace
from DWonder.MN.utils import Neuron_List_Initial, list2mask, listAddtrace, centroid_distance, list_union, initial_mask_list
from DWonder.MN.utils import list2contours, Joint_Mask_List_Simple, Joint_Mask_List_Mul
from DWonder.MN.utils import listAddcontours_Laplacian_pytorch, list2contours
from DWonder.MN.utils import group_mask, calculate_trace, Mining_rest_neuron, clear_neuron

def merge_neuron_SEG_mul(mask_stack, 
                    raw_image, 
                    quit_round_rate = 0.5, 
                    good_round_rate = 0.85,
                    good_round_size_rate = 0.4, 
                    corr_mark=0.9, 
                    area_mark=0.9, 
                    max_value=1000,
                    smallest_neuron_area = 60,
                    if_nmf = True):
    '''
    save_folder = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    if not os.path.exists(save_folder): 
        os.mkdir(save_folder)
    '''

    threshold = max_value*0.3
    mask_stack[mask_stack>threshold] = max_value
    mask_stack[mask_stack<threshold] = 0

    w_good_neuron_list = []
    w_bad_neuron_list = []
    g_f_contours = np.zeros(mask_stack.shape)
    b_f_contours = np.zeros(mask_stack.shape)

    prev_time = time.time()
    time_start = time.time()

    if_mul=1
    if if_mul:
        num_cores = int(mp.cpu_count())
        pool = mp.Pool(num_cores)
        mask_dict = {}
        for i in range(0, mask_stack.shape[0]):
            mask_dict[str(i)] = mask_stack[i,:,:].squeeze()
        # results = [pool.apply_async(initial_mask_list, args=(mask, 0.5, 0.85, 60)) for name, mask in mask_dict.items()]
        results = [pool.apply_async(initial_mask_list, args=(mask, quit_round_rate, good_round_rate, smallest_neuron_area)) for name, mask in mask_dict.items()]
        for p in results:
            good_neuron_list, bad_neuron_list = p.get()
            w_good_neuron_list = good_neuron_list + w_good_neuron_list
            w_bad_neuron_list = w_bad_neuron_list + bad_neuron_list
            # print('good_neuron_list ---> ',len(good_neuron_list),' bad_neuron_list ---> ',len(bad_neuron_list))
        time_end = time.time()
        time_cost = time_end - time_start
        print('[Neuron Segmentation Mask Initialize Time Cost: %.0d s] \n'%(time_cost), end=' ')
        pool.close()
        pool.join()
    if not if_mul:
        for i in range(0, mask_stack.shape[0]):
            now_mask = mask_stack[i,:,:].squeeze()
            good_neuron_list, bad_neuron_list = initial_mask_list(now_mask, 
                                                            quit_round_rate = quit_round_rate, 
                                                            good_round_rate = good_round_rate, 
                                                            smallest_neuron_area = smallest_neuron_area)

            w_good_neuron_list = good_neuron_list + w_good_neuron_list
            w_bad_neuron_list = w_bad_neuron_list + bad_neuron_list
            ################################################################################################################
            # Determine approximate time left
            batches_done = i+1
            batches_left = 1 * mask_stack.shape[0] - batches_done
            time_left_seconds = int(batches_left * (time.time() - prev_time))
            time_left = datetime.timedelta(seconds=time_left_seconds)
            prev_time = time.time()
            ################################################################################################################
            if i % 1 == 0:
                time_end = time.time()
                time_cost = time_end - time_start  # datetime.timedelta(seconds= (time_end - time_start))
                print(
                    '\r\033[1;33m[MERGE]\033[0m [Patch %d/%d] [Time Cost: %.0d s] [ETA: %s s]     '
                    % ( 
                    i ,
                    mask_stack.shape[0],
                    time_cost,
                    time_left_seconds
                    ), end=' ')
            if (i + 1) %  mask_stack.shape[0] == 0:
                print('\n', end=' ')
            ################################################################################################################

    w_good_neuron_list = listAddtrace(w_good_neuron_list, raw_image)

    if_time = True
    if if_time:
        time_cost = time.time() - time_start
        print('[1 Time Cost: %.0d s] \n'%(time_cost), end=' ')
    
    if_split = True
    # corr_mark=0.9, area_mark=0.9,
    if not if_split:
        w_good_neuron_list = Joint_Mask_List_Simple(w_good_neuron_list, [], corr_mark=corr_mark, area_mark=area_mark, active_rate=0, if_coor=True, if_area=True, if_merge=False)
    if if_split:
        w_good_neuron_list = Joint_Mask_List_Mul(w_good_neuron_list, corr_mark=corr_mark, area_mark=area_mark, active_rate=0, if_coor=True, if_area=True, if_merge=False)

    if if_time:
        time_cost = time.time() - time_start
        print('[2 Time Cost: %.0d s] \n'%(time_cost), end=' ')
    # w_good_neuron_list = listAddtrace(w_good_neuron_list, raw_image)
    # print('w_good_neuron_list ---> ',len(w_good_neuron_list))
    w_good_neuron_list = clear_neuron(w_good_neuron_list, [], area_mark=0.2, area_size=100)

    if if_time:
        time_cost = time.time() - time_start
        print('[3 Time Cost: %.0d s] \n'%(time_cost), end=' ')

    '''
    if not if_split:
        w_good_neuron_list = Joint_Mask_List_Simple(w_good_neuron_list, [], corr_mark=0.7, area_mark=0.6, active_rate=0, if_coor=False, if_area=True, if_merge=False)
    if if_split:
        w_good_neuron_list = Joint_Mask_List_Mul(w_good_neuron_list, corr_mark=0.7, area_mark=0.6, active_rate=0, if_coor=False, if_area=True, if_merge=False)
    '''
    # print('w_good_neuron_list ---> ',len(w_good_neuron_list))
    # w_good_neuron_list = Joint_Mask_List_Simple(w_good_neuron_list, [], corr_mark=0.7, area_mark=0.6, active_rate=0, if_coor=False, if_area=True, if_merge=False)
    
    
    if if_time:
        time_cost = time.time() - time_start
        print('[4 Time Cost: %.0d s] \n'%(time_cost), end=' ')
    # ,  save_folder=save_folder
    w_bad_neuron_list = listAddtrace(w_bad_neuron_list, raw_image)
    if not if_split:
        w_bad_neuron_list = Joint_Mask_List_Simple(w_bad_neuron_list, [], corr_mark=corr_mark, area_mark=area_mark, active_rate=0, if_coor=True, if_area=True, if_merge=True)
    if if_split:
        w_bad_neuron_list = Joint_Mask_List_Mul(w_bad_neuron_list, corr_mark=corr_mark, area_mark=area_mark, active_rate=0, if_coor=True, if_area=True, if_merge=True)

    if if_time:
        time_cost = time.time() - time_start
        print('[5 Time Cost: %.0d s] \n'%(time_cost), end=' ')

    add_neuron_list = Mining_rest_neuron(w_good_neuron_list, w_bad_neuron_list, raw_image, quit_round_rate = 0.5, smallest_neuron_area = 200)
    add_neuron_list = listAddtrace(add_neuron_list, raw_image)
    if not if_split:
        add_neuron_list = Joint_Mask_List_Simple(add_neuron_list, [], corr_mark=corr_mark, area_mark=area_mark, active_rate=0, if_coor=True, if_area=True, if_merge=True)
    if if_split:
        add_neuron_list = Joint_Mask_List_Mul(add_neuron_list, corr_mark=corr_mark, area_mark=area_mark, active_rate=0, if_coor=True, if_area=True, if_merge=True)
    
    if if_time:
        time_cost = time.time() - time_start
        print('[6 Time Cost: %.0d s] \n'%(time_cost), end=' ')

    w_g_neuron_list = w_good_neuron_list+add_neuron_list
    # corr_mark, area_mark=0.9, active_rate=0.02, if_coor=True, if_area=True, if_merge=True)

    # w_g_mask_list = listAddcontours_Laplacian_pytorch(w_g_neuron_list, mask_stack.shape[1], mask_stack.shape[2])
    # w_g_f_contours, w_g_w_contours = list2contours(w_g_mask_list, mask_stack.shape[1], mask_stack.shape[2])
    # io.imsave(save_folder+'//'+'good.tif', w_g_f_contours)
    # io.imsave(save_folder+'//'+'good_w.tif', w_g_w_contours)


    # print('w_good_neuron_list -----> ',len(w_good_neuron_list))
    if if_nmf:
        group_neuron_list, arranged_index, cc_mask_j = group_mask(w_g_neuron_list, raw_image)
        '''
        group_neuron_list1 = listAddcontours_Laplacian_pytorch(group_neuron_list, mask_stack.shape[1], mask_stack.shape[2])
        g_n_f_contours, g_n_w_contours = list2contours(group_neuron_list1, mask_stack.shape[1], mask_stack.shape[2])
        io.imsave(save_folder+'//'+'group.tif', g_n_f_contours)
        io.imsave(save_folder+'//'+'cc_mask_j.tif', cc_mask_j)
        '''
        w_g_neuron_list = calculate_trace(raw_image, group_neuron_list, w_g_neuron_list, arranged_index)

    # print('w_good_neuron_list -----> ',len(w_good_neuron_list))
    # return final_mask_list_ct
    '''
    mat_save_name = save_folder+'//'+'neuron_list.mat'
    scio.savemat(mat_save_name, {'w_good_neuron_list':w_good_neuron_list})
    '''
    w_g_neuron_list = Joint_Mask_List_Simple(w_g_neuron_list, [], corr_mark=0.7, area_mark=0.6, active_rate=0, if_coor=False, if_area=True, if_merge=False)

    return w_g_neuron_list



def merge_neuron_SEG(mask_stack, 
                    raw_image, 
                    quit_round_rate = 0.5, 
                    good_round_rate = 0.8,
                    good_round_size_rate = 0.4, 
                    corr_mark=0.9, 
                    max_value=1000,
                    if_nmf = True):

    save_folder = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    if not os.path.exists(save_folder): 
        os.mkdir(save_folder)


    threshold = max_value*0.3
    mask_stack[mask_stack>threshold] = max_value
    mask_stack[mask_stack<threshold] = 0

    w_good_neuron_list = []
    w_bad_neuron_list = []
    g_f_contours = np.zeros(mask_stack.shape)
    b_f_contours = np.zeros(mask_stack.shape)

    prev_time = time.time()
    time_start = time.time()
    for i in range(0, mask_stack.shape[0]):
        now_mask = mask_stack[i,:,:].squeeze()
        good_neuron_list, bad_neuron_list = initial_mask_list(now_mask, 
                                                        quit_round_rate = 0.5, 
                                                        good_round_rate = 0.85, 
                                                        smallest_neuron_area = 60)

        # good_mask_list = listAddcontours_Laplacian_pytorch(good_neuron_list, mask_stack.shape[1], mask_stack.shape[2])
        # good_final_contours, good_whole_contours = list2contours(good_mask_list, mask_stack.shape[1], mask_stack.shape[2])
        # g_f_contours[i,:,:] = good_final_contours
        # io.imsave(save_folder+'//'+str(i)+'_good.tif', good_final_contours)

        # bad_mask_list = listAddcontours_Laplacian_pytorch(bad_neuron_list, mask_stack.shape[1], mask_stack.shape[2])
        # bad_final_contours, bad_whole_contours = list2contours(bad_mask_list, mask_stack.shape[1], mask_stack.shape[2])
        # b_f_contours[i,:,:] = bad_final_contours
        # io.imsave(save_folder+'//'+str(i)+'_bad.tif', bad_final_contours)

        w_good_neuron_list = good_neuron_list + w_good_neuron_list
        w_bad_neuron_list = w_bad_neuron_list + bad_neuron_list
        ################################################################################################################
        # Determine approximate time left
        batches_done = i+1
        batches_left = 1 * mask_stack.shape[0] - batches_done
        time_left_seconds = int(batches_left * (time.time() - prev_time))
        time_left = datetime.timedelta(seconds=time_left_seconds)
        prev_time = time.time()
        ################################################################################################################
        if i % 1 == 0:
            time_end = time.time()
            time_cost = time_end - time_start  # datetime.timedelta(seconds= (time_end - time_start))
            print(
                '\r\033[1;33m[MERGE]\033[0m [Patch %d/%d] [Time Cost: %.0d s] [ETA: %s s]     '
                % ( 
                i ,
                mask_stack.shape[0],
                time_cost,
                time_left_seconds
                ), end=' ')
        if (i + 1) %  mask_stack.shape[0] == 0:
            print('\n', end=' ')
        ################################################################################################################

    # io.imsave(save_folder+'//'+'g_f_contours.tif', g_f_contours)
    # io.imsave(save_folder+'//'+'b_f_contours.tif', b_f_contours)

    w_good_neuron_list = listAddtrace(w_good_neuron_list, raw_image)
    w_good_neuron_list = Joint_Mask_List_Simple(w_good_neuron_list, [], corr_mark=0.9, area_mark=0.9, active_rate=0, if_coor=True, if_area=True, if_merge=False)
    # w_good_neuron_list = listAddtrace(w_good_neuron_list, raw_image)
    print('w_good_neuron_list ---> ',len(w_good_neuron_list))
    w_good_neuron_list = clear_neuron(w_good_neuron_list, [], area_mark=0.2, area_size=100)
    w_good_neuron_list = Joint_Mask_List_Simple(w_good_neuron_list, [], corr_mark=0.7, area_mark=0.7, active_rate=0, if_coor=False, if_area=True, if_merge=False)
    print('w_good_neuron_list ---> ',len(w_good_neuron_list))
    # w_good_neuron_list = Joint_Mask_List_Simple(w_good_neuron_list, [], corr_mark=0.7, area_mark=0.6, active_rate=0, if_coor=False, if_area=True, if_merge=False)

    # ,  save_folder=save_folder
    w_bad_neuron_list = listAddtrace(w_bad_neuron_list, raw_image)
    w_bad_neuron_list = Joint_Mask_List_Simple(w_bad_neuron_list, [], corr_mark=0.9, area_mark=0.9, active_rate=0, if_coor=True, if_area=True, if_merge=True)
    add_neuron_list = Mining_rest_neuron(w_good_neuron_list, w_bad_neuron_list, raw_image, quit_round_rate = 0.5, smallest_neuron_area = 200)
    add_neuron_list = listAddtrace(add_neuron_list, raw_image)
    add_neuron_list = Joint_Mask_List_Simple(add_neuron_list, [], corr_mark=0.9, area_mark=0.9, active_rate=0, if_coor=True, if_area=True, if_merge=True)
    
    w_g_neuron_list = w_good_neuron_list+add_neuron_list
    # corr_mark, area_mark=0.9, active_rate=0.02, if_coor=True, if_area=True, if_merge=True)

    w_g_mask_list = listAddcontours_Laplacian_pytorch(w_g_neuron_list, mask_stack.shape[1], mask_stack.shape[2])
    w_g_f_contours, w_g_w_contours = list2contours(w_g_mask_list, mask_stack.shape[1], mask_stack.shape[2])
    io.imsave(save_folder+'//'+'good.tif', w_g_f_contours)
    io.imsave(save_folder+'//'+'good_w.tif', w_g_w_contours)


    # print('w_good_neuron_list -----> ',len(w_good_neuron_list))
    if if_nmf:
        group_neuron_list, arranged_index, cc_mask_j = group_mask(w_good_neuron_list, raw_image)
        '''
        group_neuron_list1 = listAddcontours_Laplacian_pytorch(group_neuron_list, mask_stack.shape[1], mask_stack.shape[2])
        g_n_f_contours, g_n_w_contours = list2contours(group_neuron_list1, mask_stack.shape[1], mask_stack.shape[2])
        io.imsave(save_folder+'//'+'group.tif', g_n_f_contours)
        io.imsave(save_folder+'//'+'cc_mask_j.tif', cc_mask_j)
        '''
        w_good_neuron_list = calculate_trace(raw_image, group_neuron_list, w_good_neuron_list, arranged_index)

    # print('w_good_neuron_list -----> ',len(w_good_neuron_list))
    # return final_mask_list_ct
    '''
    mat_save_name = save_folder+'//'+'neuron_list.mat'
    scio.savemat(mat_save_name, {'w_good_neuron_list':w_good_neuron_list})
    '''
    return w_good_neuron_list




