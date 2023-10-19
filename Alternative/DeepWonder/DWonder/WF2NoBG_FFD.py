import os
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.utils.data import DataLoader
import argparse
import time
import datetime
from skimage import io
import numpy as np

from DWonder.RMBG.data_process import test_preprocess_lessMemoryNoTail_SubImg, testset, multibatch_test_save, singlebatch_test_save
from DWonder.RMBG.network import Network_3D_Unet
from DWonder.RMBG.utils import FFDrealign4, inv_FFDrealign4
#############################################################################################################################################
def wf2nobg_ffd(net,
                sub_img,
                if_use_GPU,
                RMBG_GPU='0', 
                RMBG_batch_size=1, 
                RMBG_img_w=256, 
                RMBG_img_h=256, 
                RMBG_img_s=256, 
                RMBG_gap_w=224, 
                RMBG_gap_h=224, 
                RMBG_gap_s=224, 
                RMBG_normalize_factor=1):
    #############################################################################################################################################
    opt = {}
    #############################################################################################################################################
    opt['GPU'] = RMBG_GPU
    opt['batch_size'] = RMBG_batch_size
    opt['img_w'] = RMBG_img_w
    opt['img_h'] = RMBG_img_h
    opt['img_s'] = RMBG_img_s

    opt['gap_w'] = RMBG_gap_w
    opt['gap_h'] = RMBG_gap_h
    opt['gap_s'] = RMBG_gap_s

    opt['normalize_factor'] = RMBG_normalize_factor
    #############################################################################################################################################
    if if_use_GPU:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(opt['GPU'])

    ##############################################################################################################################################################
    if if_use_GPU:
        if torch.cuda.is_available():
            net = net.cuda()
    ##############################################################################################################################################################
    name_list, noise_img, coordinate_list = test_preprocess_lessMemoryNoTail_SubImg(opt, sub_img)

    #print(len(name_list))
    prev_time = time.time()
    time_start = time.time()
    denoise_img = np.zeros(noise_img.shape)
    input_img = np.zeros(noise_img.shape)

    test_data = testset(name_list, coordinate_list, noise_img)
    testloader = DataLoader(test_data, batch_size=opt['batch_size'], shuffle=False, num_workers=0)
    for iteration, (noise_patch,single_coordinate) in enumerate(testloader):
        '''
        if if_use_GPU:
            noise_patch=noise_patch.cuda().float()
        if not if_use_GPU:
        '''
        noise_patch=noise_patch.float()
        # print('noise_patch -----> ',noise_patch.shape)
        real_A = FFDrealign4(noise_patch)
        real_A = Variable(real_A)
        if if_use_GPU:
            real_A = real_A.cuda()
        
        # print('real_A -----> ',real_A.shape)
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
            print(
                '\r\033[1;32m[RMBG]\033[0m [Patch %d/%d] [Time Cost: %.0d s] [ETA: %s s]     '
                % ( iteration + 1,
                    len(testloader),
                    time_cost,
                    time_left_seconds
                ), end=' ')

        if (iteration + 1) % len(testloader) == 0:
            print('\n', end=' ')
        ################################################################################################################
        fake_B = fake_B.cpu()
        fake_B_realign = inv_FFDrealign4(fake_B)
        output_image = np.squeeze(fake_B_realign.cpu().detach().numpy())
        real_A = real_A.cpu()
        real_A_realign = inv_FFDrealign4(real_A)
        raw_image = np.squeeze(real_A_realign.cpu().detach().numpy())
        if(output_image.ndim==3):
            turn=1
        else:
            turn=output_image.shape[0]
        #print(turn)
        if(turn>1):
            for id in range(turn):
                #print('shape of output_image -----> ',output_image.shape)
                aaaa,bbbb,stack_start_w,stack_end_w,stack_start_h,stack_end_h,stack_start_s,stack_end_s=multibatch_test_save(single_coordinate,id,output_image,raw_image)
                denoise_img[stack_start_s:stack_end_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
                = aaaa 
                # input_img[stack_start_s:stack_end_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
                # = bbbb
        else:
            aaaa, bbbb, stack_start_w, stack_end_w, stack_start_h, stack_end_h, stack_start_s, stack_end_s = singlebatch_test_save(
                single_coordinate, output_image, raw_image)
            denoise_img[stack_start_s:stack_end_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
                            = aaaa 
            # input_img[stack_start_s:stack_end_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
            #     = bbbb
            

    # del noise_img
    output_img = denoise_img.squeeze().astype(np.float32) * opt['normalize_factor']
    del denoise_img
    # print(np.min(output_img),' -----> ',np.max(output_img))
    output_img = np.clip(output_img, 0, 65535).astype('uint16')
    # output_img = output_img - output_img.min()
    # output_img = output_img.astype('uint16')

    return output_img

