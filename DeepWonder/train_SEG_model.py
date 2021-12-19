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
import numpy as np
from skimage import io

import random

from deepwonder.SEG.data_process import train_preprocess_lessMemory_seg, shuffle_datasets_lessMemory
from deepwonder.SEG.network import SEG_Network_3D_Unet
from deepwonder.SEG.utils import FFDrealign4, inv_FFDrealign4
#############################################################################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("--epoch", type=int, default=0, help="epoch to start training from")
parser.add_argument("--n_epochs", type=int, default=100, help="number of training epochs")
parser.add_argument('--GPU', type=int, default=3, help="the index of GPU you will use for computation")
parser.add_argument('--cuda', action='store_true', help='use GPU computation')
parser.add_argument('--output_dir', type=str, default='results', help="the output folder")

parser.add_argument('--img_w', type=int, default=512, help="")
parser.add_argument('--img_h', type=int, default=512, help="")
parser.add_argument('--input_nc', type=int, default=1, help="")
parser.add_argument('--output_nc', type=int, default=1, help="")
parser.add_argument('--f_maps', type=int, default=1, help="")
parser.add_argument('--frame_num', type=int, default=1, help="")

parser.add_argument('--lr', type=float, default=0.0001, help='initial learning rate')
parser.add_argument("--b1", type=float, default=0.5, help="adam: decay of first order momentum of gradient")
parser.add_argument("--b2", type=float, default=0.999, help="adam: decay of first order momentum of gradient")
parser.add_argument('--normalize_factor', type=int, default=10000, help='actions: train or predict')

parser.add_argument('--datasets_folder', type=str, default='rawdata', help="the name of your project")
parser.add_argument('--datasets_path', type=str, default='datasets', help="the name of your project")
parser.add_argument('--input_folder', type=str, default='image', help="")
parser.add_argument('--GT_folder', type=str, default='mask', help="")
parser.add_argument('--pth_path', type=str, default='SEG_pth', help="the name of your project")
parser.add_argument('--train_datasets_size', type=int, default=1000, help='actions: train or predict')

opt = parser.parse_args()
print('the parameter of your training ----->')
print(opt)
########################################################################################################################
if not os.path.exists(opt.output_dir): 
    os.mkdir(opt.output_dir)
current_time = 'TS3DUnetFFD_'+opt.datasets_path+'_'+opt.datasets_folder \
    +'_ic'+str(opt.input_nc)+'_oc'+str(opt.output_nc)+'_lr'+str(opt.lr)+'_fm'+str(opt.f_maps)+'_'+datetime.datetime.now().strftime("%Y%m%d-%H%M")

output_path = opt.output_dir + '//' + current_time
pth_folder = opt.pth_path+'//'+ current_time

print('output_path ---> ',output_path)
print('pth_folder ---> ',pth_folder)

if not os.path.exists(opt.output_dir): 
    os.mkdir(opt.output_dir)
if not os.path.exists(output_path): 
    os.mkdir(output_path)
if not os.path.exists(opt.pth_path): 
    os.mkdir(opt.pth_path)
if not os.path.exists(pth_folder): 
    os.mkdir(pth_folder)


os.environ["CUDA_VISIBLE_DEVICES"] = str(opt.GPU)
lr = opt.lr

coor_list, GT_list, input_list = train_preprocess_lessMemory_seg(opt)

net_seg = SEG_Network_3D_Unet(UNet_type= 'TS_UNet3D',
                            in_channels = 4,
                            out_channels = 4,
                            frame_num = opt.frame_num,
                            final_sigmoid = True,
                            f_maps = opt.f_maps)

L1loss_function = torch.nn.L1Loss()
L2loss_function = torch.nn.MSELoss()
BCEloss_function = torch.nn.BCELoss()
if torch.cuda.is_available():
    print('Using GPU.')
    net_seg.cuda()
    L1loss_function.cuda()
    L2loss_function.cuda()
    BCEloss_function.cuda()


optimizer_seg = torch.optim.Adam(net_seg.parameters(), lr=opt.lr, betas=(opt.b1, 0.999))
########################################################################################################################

time_start=time.time()
prev_time = time.time()
# print('train_im_list key -----> ',train_im_list.keys())
iteration_num = 0
for epoch in range(opt.epoch, opt.n_epochs):
    coor_list = shuffle_datasets_lessMemory(coor_list)
    # print('name list -----> ',name_list)   
    for index in range(len(coor_list)):
        per_coor = coor_list[index]
        train_im_name = per_coor['name']
        init_w = per_coor['init_w']
        init_h = per_coor['init_h']
        GT_im = GT_list[train_im_name]
        input_im = input_list[train_im_name]

        # print('GT_im -----> ',GT_im.shape)
        # print('input_im -----> ',input_im.shape)
        GT_patch = GT_im[:, init_w:init_w+opt.img_w, init_h:init_h+opt.img_h]
        input_patch = input_im[:, init_w:init_w+opt.img_w, init_h:init_h+opt.img_h].copy()

        rand_bg = np.random.randint(0, 100)
        rand_gama = np.random.randint(1000, 5000)/1000
        input_patch = (input_patch+rand_bg)/rand_gama

        GT_patch = torch.from_numpy(GT_patch.astype(np.float32))
        GT_patch = torch.from_numpy(np.expand_dims(np.expand_dims(GT_patch,0),0)).cuda()
        input_patch = torch.from_numpy(input_patch.astype(np.float32))
        input_patch = torch.from_numpy(np.expand_dims(np.expand_dims(input_patch,0),0)).cuda()
        # print('input_patch -----> ',input_patch.shape)
        # print('GT_patch -----> ',GT_patch.shape)
        GT_patch = FFDrealign4(GT_patch).cuda()
        input_patch = FFDrealign4(input_patch).cuda()

        # print('input_patch -----> ',input_patch.shape)
        pred_patch = net_seg(input_patch)
        # print('pred_patch -----> ',pred_patch.shape)

        optimizer_seg.zero_grad()
        # L1Loss_A2B = L1loss_function(train_imB, pred_imA)
        BCELoss_B2A_SDNN = BCEloss_function(pred_patch, GT_patch)
        L2Loss_B2A_SDNN = L2loss_function(pred_patch, GT_patch)
        loss_SDNN = BCELoss_B2A_SDNN+L2Loss_B2A_SDNN
        loss_SDNN.backward()
        optimizer_seg.step()

        iteration_num = iteration_num +1
        ################################################################################################################
        batches_done = epoch * opt.train_datasets_size + index
        batches_left = opt.n_epochs * opt.train_datasets_size - batches_done
        time_left = datetime.timedelta(seconds=batches_left * (time.time() - prev_time))
        prev_time = time.time()
        ################################################################################################################ 8_HCC20
        if (index%1000 == 0):
            time_end=time.time()
            print('time cost',time_end-time_start,'s')
            sys.stdout.write("\r[Epoch %d/%d] [Batch %d/%d]  ETA: %s"
            % (epoch, opt.n_epochs, index, opt.train_datasets_size, time_left ))
            print('    loss_SDNN ',loss_SDNN.cpu().detach().numpy())
            print('iteration_num ',iteration_num)
            print(GT_patch.cpu().detach().numpy().max(),' ---> ', GT_patch.cpu().detach().numpy().min(),' ---> ',\
            pred_patch.cpu().detach().numpy().max(),' ---> ',  pred_patch.cpu().detach().numpy().min(),' ---> ',\
            input_patch.cpu().detach().numpy().max(),' ---> ',  input_patch.cpu().detach().numpy().min())

        # if (index%50 == 0): # or ((epoch+1)%1 == 0):
        if (iteration_num+1)%1000 == 0:
            print('save image')
            train_im_name = per_coor['name']
            input_output_path = output_path + '/input'
            pred_output_path = output_path + '/pred'
            GT_output_path = output_path + '/GT'

            if not os.path.exists(input_output_path): 
                os.mkdir(input_output_path)
            if not os.path.exists(pred_output_path): 
                os.mkdir(pred_output_path)
            if not os.path.exists(GT_output_path): 
                os.mkdir(GT_output_path)

            input_patch_realign = inv_FFDrealign4(input_patch)
            input_patch_realign = input_patch_realign.cpu().detach().numpy()
            pred_patch_realign = inv_FFDrealign4(pred_patch)
            pred_patch_realign = pred_patch_realign.cpu().detach().numpy()
            GT_patch_realign = inv_FFDrealign4(GT_patch)
            GT_patch_realign = GT_patch_realign.cpu().detach().numpy()

            input_patch_realign = input_patch_realign.squeeze().astype(np.float32)*opt.normalize_factor
            pred_patch_realign = pred_patch_realign.squeeze().astype(np.float32)*opt.normalize_factor
            GT_patch_realign = GT_patch_realign.squeeze().astype(np.float32)*opt.normalize_factor

            input_patch_realign = np.clip(input_patch_realign, 0, 65535).astype('uint16')
            pred_patch_realign = np.clip(pred_patch_realign, 0, 65535).astype('uint16')
            GT_patch_realign = np.clip(GT_patch_realign, 0, 65535).astype('uint16')

            input_name = input_output_path + '/' + str(epoch) + '_' + str(index) + '_' + train_im_name+'_input.tif'
            pred_name = pred_output_path + '/' + str(epoch) + '_' + str(index) + '_' + train_im_name+'_pred.tif'
            GT_name = GT_output_path + '/' + str(epoch) + '_' + str(index) + '_' + train_im_name+'_GT.tif'

            io.imsave(input_name, input_patch_realign)
            io.imsave(pred_name, pred_patch_realign)
            io.imsave(GT_name, GT_patch_realign)

    if (epoch+1)%1 == 0:
        torch.save(net_seg.state_dict(), pth_folder + '//seg_' + str(epoch) + '.pth')





