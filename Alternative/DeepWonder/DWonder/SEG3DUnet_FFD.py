import os
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.utils.data import DataLoader
import argparse
import time
import datetime
import sys
import math
import scipy.io as scio
from skimage import io
import numpy as np
import math

from DWonder.SEG.network import SEG_Network_3D_Unet
from DWonder.SEG.data_process import test_preprocess_lessMemoryNoTail_SubImgSEG, testset, singlebatch_test_save, multibatch_test_save
from DWonder.SEG.utils import FFDrealign4, inv_FFDrealign4

def seg_3dunet_ffd(net,
                sub_img,
                SEG_ffd = True,
                if_use_GPU = True,
                SEG_GPU='0', 
                SEG_batch_size=1, 
                SEG_img_w=256, 
                SEG_img_h=256, 
                SEG_img_s=64, 
                SEG_gap_w=224, 
                SEG_gap_h=224, 
                SEG_gap_s=32, 
                SEG_normalize_factor=1000):
    #############################################################################################################################################
    opt = {}
    #############################################################################################################################################
    opt['GPU'] = SEG_GPU
    opt['batch_size'] = SEG_batch_size
    opt['img_w'] = SEG_img_w
    opt['img_h'] = SEG_img_h
    opt['img_s'] = SEG_img_s

    opt['gap_w'] = SEG_gap_w
    opt['gap_h'] = SEG_gap_h
    opt['gap_s'] = SEG_gap_s

    opt['normalize_factor'] = SEG_normalize_factor
    #############################################################################################################################################
    os.environ["CUDA_VISIBLE_DEVICES"] = str(opt['GPU'])
    ##############################################################################################################################################################
    if if_use_GPU:
        if torch.cuda.is_available():
            net = net.cuda()
    ##############################################################################################################################################################
    name_list, noise_img, coordinate_list = test_preprocess_lessMemoryNoTail_SubImgSEG(opt, sub_img)

        #print(len(name_list))
    prev_time = time.time()
    time_start = time.time()

    num_s = math.ceil((noise_img.shape[0]-SEG_img_s+SEG_gap_s)/SEG_gap_s)
    denoise_img = np.zeros((num_s, noise_img.shape[1], noise_img.shape[2]))

    test_data = testset(name_list, coordinate_list, noise_img)
    testloader = DataLoader(test_data, batch_size=opt['batch_size'], shuffle=False, num_workers=4)
    for iteration, (noise_patch,single_coordinate) in enumerate(testloader):
        noise_patch=noise_patch.float()

        if SEG_ffd:
            real_A = FFDrealign4(noise_patch)
        if not SEG_ffd:
            real_A = noise_patch
        real_A = Variable(real_A)
        if if_use_GPU:
            real_A = real_A.cuda()
        
        fake_B = net(real_A)
        ################################################################################################################
        # Determine approximate time left
        batches_done = iteration
        batches_left = 1 * len(testloader) - batches_done
        time_left_seconds = int(batches_left * (time.time() - prev_time))
        time_left = datetime.timedelta(seconds=time_left_seconds)
        prev_time = time.time()
        ################################################################################################################
        if iteration % 1 == 0:
            time_end = time.time()
            time_cost = time_end - time_start  # datetime.timedelta(seconds= (time_end - time_start))
            # printf('\033[1;40;32m color!!! \033[0m hello\n')
            # print('\033[1;31mUsing GPU for training -----> \033[0m')
            print(
                '\r\033[1;31m[SEG]\033[0m [Patch %d/%d] [Time Cost: %.0d s] [ETA: %s s]    '
                % ( iteration + 1,
                    len(testloader),
                    time_cost,
                    time_left_seconds
                ), end=' ')
            # print(noise_patch.cpu().detach().numpy().max(),' ---> ', noise_patch.cpu().detach().numpy().min())

        if (iteration + 1) % len(testloader) == 0:
            print('\n', end=' ')
        ################################################################################################################
        if SEG_ffd:
            fake_B = fake_B.cpu()
            fake_B_realign = inv_FFDrealign4(fake_B)
        if not SEG_ffd:
            fake_B = fake_B.cpu()
            fake_B_realign = fake_B

        # fake_B_realign = fake_B
        output_image = np.squeeze(fake_B_realign.cpu().detach().numpy())
        # real_A_realign = inv_FFDrealign4(real_A)
        # real_A_realign = real_A
        # raw_image = np.squeeze(real_A_realign.cpu().detach().numpy())
        if(output_image.ndim==2):
            turn=1
        else:
            turn=output_image.shape[0]
        #print(turn)
        if(turn>1):
            for id in range(turn):
                #print('shape of output_image -----> ',output_image.shape)
                aaaa, stack_start_w, stack_end_w, stack_start_h, stack_end_h, stack_start_s= \
                multibatch_test_save(single_coordinate,id,output_image)
                denoise_img[stack_start_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
                = aaaa 

        else:
            aaaa, stack_start_w, stack_end_w, stack_start_h, stack_end_h, stack_start_s = \
            singlebatch_test_save(single_coordinate, output_image)
            # print(stack_start_s)
            denoise_img[stack_start_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
                            = aaaa 
    return denoise_img

