import os
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.utils.data import DataLoader
import argparse
import time
import datetime
import numpy as np
from skimage import io
import math

from deepwonder.data_process import Img2Subimg

from deepwonder.RMBG.network import Network_3D_Unet
from deepwonder.SEG.network import SEG_Network_3D_Unet

from deepwonder.WF2NoBG_FFD import wf2nobg_ffd
from deepwonder.SEG3DUnet import seg_3dunet
from deepwonder.SEG3DUnet_FFD import seg_3dunet_ffd

from deepwonder.MergeNeuron_SEG import merge_neuron_SEG, merge_neuron_SEG_mul
from deepwonder.MN.utils import joint_neuron
from deepwonder.MN.utils import listAddcontours_Laplacian_pytorch, list2contours
from deepwonder.MN.utils import listAddtrace
import scipy.io as scio

import warnings
warnings.filterwarnings('ignore')
############################################################################################################################################# 
parser = argparse.ArgumentParser()
parser.add_argument('--GPU', type=str, default='0', help="")

parser.add_argument('--sub_img_w', type=int, default=304, help="")
parser.add_argument('--sub_img_h', type=int, default=304, help="")
parser.add_argument('--sub_img_s', type=int, default=2000, help="")

parser.add_argument('--sub_gap_w', type=int, default=272, help='')
parser.add_argument('--sub_gap_h', type=int, default=272, help='')
parser.add_argument('--sub_gap_s', type=int, default=1960, help='')

parser.add_argument('--datasets_path', type=str, default='datasets', help="")
parser.add_argument('--datasets_folder', type=str, default='test', help="")
parser.add_argument('--test_datasize', type=int, default=2000, help='dataset size to be tested')
###############################################################################
parser.add_argument('--RMBG_GPU', type=str, default='0', help="")
parser.add_argument('--RMBG_batch_size', type=int, default=1, help="")
parser.add_argument('--RMBG_normalize_factor', type=int, default=1, help='normalize factor')
parser.add_argument('--RMBG_fmap', type=int, default=32, help='number of feature maps')

parser.add_argument('--RMBG_img_w', type=int, default=256, help="")
parser.add_argument('--RMBG_img_h', type=int, default=256, help="")
parser.add_argument('--RMBG_img_s', type=int, default=128, help="")

parser.add_argument('--RMBG_gap_w', type=int, default=224, help='')
parser.add_argument('--RMBG_gap_h', type=int, default=224, help='')
parser.add_argument('--RMBG_gap_s', type=int, default=96, help='')

parser.add_argument('--RMBG_model_path', type=str, default='RMBG_pth', help="")
parser.add_argument('--RMBG_model_folder', type=str, default='pth', help="")
parser.add_argument('--RMBG_model_name', type=str, default='train_20210401_1712', help='')
###############################################################################
parser.add_argument('--SEG_GPU', type=str, default='0', help="")
parser.add_argument('--SEG_batch_size', type=int, default=1, help="")
parser.add_argument('--SEG_normalize_factor', type=int, default=1000, help='normalize factor')
parser.add_argument('--SEG_fmap', type=int, default=16, help='number of feature maps')

parser.add_argument('--SEG_img_w', type=int, default=256, help="")
parser.add_argument('--SEG_img_h', type=int, default=256, help="")
parser.add_argument('--SEG_img_s', type=int, default=64, help="")

parser.add_argument('--SEG_gap_w', type=int, default=224, help='')
parser.add_argument('--SEG_gap_h', type=int, default=224, help='')
parser.add_argument('--SEG_gap_s', type=int, default=32, help='')

parser.add_argument('--SEG_model_path', type=str, default='SEG_pth', help="")
parser.add_argument('--SEG_model_folder', type=str, default='pth', help="")
parser.add_argument('--SEG_model_name', type=str, default='train_20210401_1712', help='')
###############################################################################
parser.add_argument('--output_dir', type=str, default='./results', help="output directory")
opt = parser.parse_args()
print(opt.datasets_folder)
############################################################################################################################################# 
opt.RMBG_GPU = opt.GPU
opt.SEG_GPU = opt.GPU
os.environ["CUDA_VISIBLE_DEVICES"] = str(opt.GPU)
###############################################################################
opt.sub_gap_s = opt.sub_img_s-40


opt.RMBG_img_w = opt.sub_img_w
opt.RMBG_img_h = opt.sub_img_h

opt.RMBG_gap_w = opt.sub_gap_w
opt.RMBG_gap_h = opt.sub_gap_h

opt.SEG_img_w = opt.sub_img_w
opt.SEG_img_h = opt.sub_img_h

opt.SEG_gap_w = opt.sub_gap_w
opt.SEG_gap_h = opt.sub_gap_h
###############################################################################
opt_dic={}

opt_dic['GPU'] = opt.GPU
opt_dic['sub_img_w'] = opt.sub_img_w
opt_dic['sub_img_h'] = opt.sub_img_h
opt_dic['sub_img_s'] = opt.sub_img_s

opt_dic['sub_gap_w'] = opt.sub_gap_w
opt_dic['sub_gap_h'] = opt.sub_gap_h
opt_dic['sub_gap_s'] = opt.sub_gap_s

opt_dic['datasets_path'] = opt.datasets_path
opt_dic['datasets_folder'] = opt.datasets_folder
opt_dic['test_datasize'] = opt.test_datasize

opt_dic['RMBG_GPU'] = opt.RMBG_GPU
opt_dic['RMBG_batch_size'] = opt.RMBG_batch_size
opt_dic['RMBG_normalize_factor'] = opt.RMBG_normalize_factor
opt_dic['RMBG_fmap'] = opt.RMBG_fmap

opt_dic['RMBG_img_w'] = opt.RMBG_img_w
opt_dic['RMBG_img_h'] = opt.RMBG_img_h
opt_dic['RMBG_img_s'] = opt.RMBG_img_s

opt_dic['RMBG_gap_w'] = opt.RMBG_gap_w
opt_dic['RMBG_gap_h'] = opt.RMBG_gap_h
opt_dic['RMBG_gap_s'] = opt.RMBG_gap_s

opt_dic['RMBG_model_path'] = opt.RMBG_model_path
opt_dic['RMBG_model_folder'] = opt.RMBG_model_folder
opt_dic['RMBG_model_name'] = opt.RMBG_model_name

opt_dic['SEG_GPU'] = opt.SEG_GPU
opt_dic['SEG_batch_size'] = opt.SEG_batch_size
opt_dic['SEG_normalize_factor'] = opt.SEG_normalize_factor
opt_dic['SEG_fmap'] = opt.SEG_fmap

opt_dic['SEG_img_w'] = opt.SEG_img_w
opt_dic['SEG_img_h'] = opt.SEG_img_h
opt_dic['SEG_img_s'] = opt.SEG_img_s

opt_dic['SEG_gap_w'] = opt.SEG_gap_w
opt_dic['SEG_gap_h'] = opt.SEG_gap_h
opt_dic['SEG_gap_s'] = opt.SEG_gap_s

opt_dic['SEG_model_path'] = opt.SEG_model_path
opt_dic['SEG_model_folder'] = opt.SEG_model_folder
opt_dic['SEG_model_name'] = opt.SEG_model_name

opt_dic['output_dir'] = opt.output_dir
###############################################################################
name_list, img_list, coordinate_list = Img2Subimg(opt_dic)



if_use_GPU = False
###############################################################################
RMBG_net = Network_3D_Unet(in_channels=4,
                            out_channels=4,
                            f_maps=opt_dic['RMBG_fmap'],
                            final_sigmoid=True)
RMBG_w_model_name = opt_dic['RMBG_model_path']+'//'+opt_dic['RMBG_model_folder']+'//'+opt_dic['RMBG_model_name']+'.pth'
if if_use_GPU:
    if torch.cuda.is_available():
        RMBG_net.load_state_dict(torch.load(RMBG_w_model_name))
if not if_use_GPU:
    RMBG_net.load_state_dict(torch.load(RMBG_w_model_name, map_location='cpu'))
###############################################################################
SEG_net = SEG_Network_3D_Unet(UNet_type = 'TS_UNet3D',
                            in_channels = 4,
                            out_channels = 4,
                            frame_num = opt_dic['SEG_img_s'],
                            final_sigmoid = True,
                            f_maps = opt_dic['SEG_fmap'])
SEG_w_model_name = opt_dic['SEG_model_path']+'//'+opt_dic['SEG_model_folder']+'//'+opt_dic['SEG_model_name']+'.pth'
if if_use_GPU:
    if torch.cuda.is_available():
        SEG_net.load_state_dict(torch.load(SEG_w_model_name))
if not if_use_GPU:
    SEG_net.load_state_dict(torch.load(SEG_w_model_name, map_location='cpu'))
###############################################################################

if not os.path.exists(opt_dic['output_dir']): 
    os.mkdir(opt_dic['output_dir'])

ccm = 0.95
cgrr = 0.75
current_time = datetime.datetime.now().strftime("%Y%m%d-%H%M")
all_output_path = opt_dic['output_dir'] + '//' + 'RSM_'+opt_dic['datasets_folder']+'_'+current_time
if not os.path.exists(all_output_path): 
    os.mkdir(all_output_path)

SEG_RMBG_output_path = all_output_path+ '//' + 'RMBG'
if not os.path.exists(SEG_RMBG_output_path): 
    os.mkdir(SEG_RMBG_output_path)
SEG_SEG_output_path = all_output_path+ '//' + 'SEG'
if not os.path.exists(SEG_SEG_output_path): 
    os.mkdir(SEG_SEG_output_path)
SEG_w_con_output_path =all_output_path+ '//' + 'w_con'
if not os.path.exists(SEG_w_con_output_path): 
    os.mkdir(SEG_w_con_output_path)
SEG_f_con_output_path = all_output_path+ '//' + 'f_con'
if not os.path.exists(SEG_f_con_output_path): 
    os.mkdir(SEG_f_con_output_path)
SEG_mat_con_output_path = all_output_path+ '//' + 'mat'
if not os.path.exists(SEG_mat_con_output_path): 
    os.mkdir(SEG_mat_con_output_path)

w_time_start = time.time()
prev_time = time.time()
time_start = time.time()
# print('len(name_list) -----> ',len(name_list))
for im_index in range(0, len(name_list)):
    im_name = name_list[im_index]
    per_coor_list = coordinate_list[im_name]
    img = img_list[im_name]

    img_RMBG = np.zeros(img.shape,dtype='uint16')
    w_num_s = math.ceil((img.shape[0]-opt_dic['sub_img_s']+opt_dic['sub_gap_s'])/opt_dic['sub_gap_s'])
    num_s = math.ceil((opt_dic['sub_img_s']-opt_dic['SEG_img_s']+opt_dic['SEG_gap_s'])/opt_dic['SEG_gap_s'])
    # print('num_s --- ',num_s)
    img_SEG = np.zeros((num_s*w_num_s, img.shape[1], img.shape[2])).astype('uint16')
    whole_mask_list = []
    for coor_index in range (0, len(per_coor_list)):
        per_coor = per_coor_list[coor_index]
        
        init_h = per_coor['init_h']
        end_h = per_coor['end_h']
        init_w = per_coor['init_w']
        end_w = per_coor['end_w']
        init_s = per_coor['init_s']
        end_s = per_coor['end_s']

        sub_img = img[init_s:end_s,init_h:end_h,init_w:end_w]
        #############################################################################
        sub_img_RMBG = wf2nobg_ffd(RMBG_net,
                sub_img,
                if_use_GPU = if_use_GPU,
                RMBG_GPU = opt_dic['RMBG_GPU'], 
                RMBG_batch_size = opt_dic['RMBG_batch_size'], 
                RMBG_img_w=opt_dic['RMBG_img_w'], 
                RMBG_img_h=opt_dic['RMBG_img_h'], 
                RMBG_img_s = opt_dic['RMBG_img_s'], 
                RMBG_gap_w = opt_dic['RMBG_gap_w'], 
                RMBG_gap_h = opt_dic['RMBG_gap_h'], 
                RMBG_gap_s = opt_dic['RMBG_gap_s'], 
                RMBG_normalize_factor = opt_dic['RMBG_normalize_factor'])
        #############################################################################
        sub_img_SEG = seg_3dunet_ffd(SEG_net,
                sub_img_RMBG,
                if_use_GPU = if_use_GPU,
                SEG_GPU = opt_dic['SEG_GPU'],
                SEG_batch_size = opt_dic['SEG_batch_size'], 
                SEG_img_w = opt_dic['SEG_img_w'], 
                SEG_img_h = opt_dic['SEG_img_h'], 
                SEG_img_s = opt_dic['SEG_img_s'], 
                SEG_gap_w = opt_dic['SEG_gap_w'], 
                SEG_gap_h = opt_dic['SEG_gap_h'], 
                SEG_gap_s = opt_dic['SEG_gap_s'], 
                SEG_normalize_factor = opt_dic['SEG_normalize_factor'])
        sub_img_SEG = sub_img_SEG*opt_dic['SEG_normalize_factor']
        #############################################################################

        stack_start_w = int(per_coor['stack_start_w'])
        stack_end_w = int(per_coor['stack_end_w'])
        patch_start_w = int(per_coor['patch_start_w'])
        patch_end_w = int(per_coor['patch_end_w'])

        stack_start_h = int(per_coor['stack_start_h'])
        stack_end_h = int(per_coor['stack_end_h'])
        patch_start_h = int(per_coor['patch_start_h'])
        patch_end_h = int(per_coor['patch_end_h'])

        stack_start_s = int(per_coor['stack_start_s'])
        stack_end_s = int(per_coor['stack_end_s'])
        patch_start_s = int(per_coor['patch_start_s'])
        patch_end_s = int(per_coor['patch_end_s'])
        # print('stack_start_w -----> ',stack_start_w, ' patch_start_w -----> ',patch_start_w,'stack_end_w -----> ',stack_end_w, ' patch_end_w -----> ',patch_end_w)
        # print('stack_start_h -----> ',stack_start_h, ' patch_start_h -----> ',patch_start_h, 'stack_end_h -----> ',stack_end_h, ' patch_end_h -----> ',patch_end_h)
        img_RMBG[stack_start_s:stack_end_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
        = sub_img_RMBG[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
        #############################################################################
        z = int(per_coor['z'])
        # print('sub_img_SEG ----- ',sub_img_SEG.shape)
        # print('z ----- ',z, z*num_s, z*num_s+num_s)
        img_SEG[z*num_s:z*num_s+num_s, stack_start_h:stack_end_h, stack_start_w:stack_end_w] \
        = sub_img_SEG[:, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
        aa = sub_img_SEG[:, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
        ################################################################################################################
        # Determine approximate time left
        batches_done = coor_index+im_index*len(per_coor_list)+1
        batches_left = len(name_list)*len(per_coor_list) - batches_done
        time_left_seconds = int(batches_left * (time.time() - prev_time))
        time_left = datetime.timedelta(seconds=time_left_seconds)
        prev_time = time.time()
        ################################################################################################################
        # if coor_index % 1 == 0:
        time_now = time.time()
        time_cost = time_now - time_start
        print(
            '\r\033[1;34m[WHOLE]\033[0m [Patch %d/%d] [Time Cost: %.0d s] [ETA: %s s]    '
            % ( batches_done,
                len(name_list)*len(per_coor_list),
                time_cost,
                time_left_seconds
            ), end=' ')
        print(opt_dic['datasets_folder'],'\n', end=' ')


    img_SEG1 = img_SEG.copy()
    img_RMBG1 = img_RMBG.copy()
    #############################################################################

    whole_mask_list = merge_neuron_SEG_mul(img_SEG1, 
                                    img_RMBG1, 
                                    quit_round_rate = 0.5, 
                                    good_round_rate = 0.8,
                                    good_round_size_rate = 0.6, 
                                    corr_mark=0.9, 
                                    max_value=1000 ,
                                    if_nmf = False)
    print('whole_mask_list -----> ',len(whole_mask_list))

    img_RMBG_name = SEG_RMBG_output_path+'//'+opt_dic['RMBG_model_folder']+'_'+opt_dic['RMBG_model_name']+'_'+im_name+'_RMBG.tif'
    img_RMBG = img_RMBG.astype('uint16')
    io.imsave(img_RMBG_name, img_RMBG)

    img_SEG_name = SEG_SEG_output_path+'//'+opt_dic['SEG_model_name']+'_'+im_name+'_SEG.tif'
    img_SEG = img_SEG.astype('uint16')
    io.imsave(img_SEG_name, img_SEG)


    # whole_mask_list = listAddcontours_Laplacian(whole_mask_list, img.shape[1], img.shape[2])
    whole_mask_list = listAddtrace(whole_mask_list, img_RMBG, mode='update', trace_mode='all')
    whole_mask_list = listAddcontours_Laplacian_pytorch(whole_mask_list, img.shape[1], img.shape[2])
    final_contours,whole_contours = list2contours(whole_mask_list, img.shape[1], img.shape[2])

    img_f_contours_name = SEG_f_con_output_path+'//'+opt_dic['SEG_model_name']+'_'+im_name+'_f_con.tif'
    final_contours = final_contours.astype('uint16')
    io.imsave(img_f_contours_name, final_contours)

    mat_save_name = SEG_mat_con_output_path+'//'+opt_dic['SEG_model_name']+'_'+im_name+'.mat'
    data = {'a':whole_mask_list, 'final_contours':final_contours}
    scio.savemat(mat_save_name, {'final_mask_list':data['a'], 'final_contours':data['final_contours']})

    img_w_contours_name = SEG_w_con_output_path+'//'+opt_dic['SEG_model_name']+'_'+im_name+'_w_con.tif'
    whole_contours = whole_contours.astype('uint16')

time_end = time.time()
time_left_seconds = int( (time_end - time_start))
print('Time Cost: ',time_left_seconds)
