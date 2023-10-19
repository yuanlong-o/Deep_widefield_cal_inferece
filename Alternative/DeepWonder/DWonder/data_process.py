import numpy as np
import os
import tifffile as tiff
import random
import math
import torch
from torch.utils.data import Dataset


def Img2Subimg(args):

    im_folder = args['datasets_path']+'//'+args['datasets_folder']

    name_list = []
    # train_raw = []
    image_list={}
    coordinate_list={}
    print('im_folder ---> ',im_folder)
    print('list(os.walk(im_folder, topdown=False))[-1][-1] ---> ',list(os.walk(im_folder, topdown=False))[-1][-1])
    for im_name in list(os.walk(im_folder, topdown=False))[-1][-1]:
        img_h = args['sub_img_h']
        img_w = args['sub_img_w']
        img_s2 = args['sub_img_s']
        gap_h = args['sub_gap_h']
        gap_w = args['sub_gap_w']
        gap_s2 = args['sub_gap_s']
        cut_w = (img_w - gap_w)/2
        cut_h = (img_h - gap_h)/2
        cut_s = (img_s2 - gap_s2)/2

        print('im_name -----> ',im_name)
        im_dir = im_folder+'//'+im_name
        im = tiff.imread(im_dir)
        print('im -----> ',im.shape)

        if '.tiff' in im_name:
            im_name = im_name.replace('.tiff','')
        if ('.tif' in im_name)&('.tiff' not in im_name):
            im_name = im_name.replace('.tif','')
        name_list.append(im_name)

        if im.shape[0]>args['test_datasize']:
            im = im[0:args['test_datasize'],:,:].copy()
        if im.shape[0]<img_s2:
            img_s2 = im.shape[0]
            args['sub_img_s'] = img_s2
            gap_s2 = img_s2
            cut_s = 0
        if im.shape[1]<img_h:
            img_h = im.shape[1]
            args['sub_img_h'] = img_h
            gap_h = img_h
            cut_h = 0
        if im.shape[2]<img_w:
            img_w = im.shape[2]
            args['sub_img_w'] = img_w
            gap_w = img_w
            cut_w = 0
        
        image_list[im_name.replace('.tif','')] = im

        whole_w = im.shape[2]
        whole_h = im.shape[1]
        whole_s = im.shape[0]

        num_w = math.ceil((whole_w-img_w+gap_w)/gap_w)
        num_h = math.ceil((whole_h-img_h+gap_h)/gap_h)
        num_s = math.ceil((whole_s-img_s2+gap_s2)/gap_s2)

        per_coor_list = []
        for x in range(0,num_h):
            for y in range(0,num_w):
                for z in range(0,num_s):
                    # print('x -----> ',x,num_h,'y -----> ',y,num_w,'z -----> ',z,num_s)
                    per_coor = {}
                    if x != (num_h-1):
                        init_h = gap_h*x
                        end_h = gap_h*x + img_h
                    elif x == (num_h-1):
                        init_h = whole_h - img_h
                        end_h = whole_h

                    if y != (num_w-1):
                        init_w = gap_w*y
                        end_w = gap_w*y + img_w
                    elif y == (num_w-1):
                        init_w = whole_w - img_w
                        end_w = whole_w

                    if z != (num_s-1):
                        init_s = gap_s2*z
                        end_s = gap_s2*z + img_s2
                    elif z == (num_s-1):
                        init_s = whole_s - img_s2
                        end_s = whole_s
                    per_coor['init_h'] = init_h
                    per_coor['end_h'] = end_h
                    per_coor['init_w'] = init_w
                    per_coor['end_w'] = end_w
                    per_coor['init_s'] = init_s
                    per_coor['end_s'] = end_s

                    if (num_w-1)==0:
                        per_coor['stack_start_w'] = 0
                        per_coor['stack_end_w'] = whole_w
                        per_coor['patch_start_w'] = 0
                        per_coor['patch_end_w'] = img_w
                    elif (num_w-1)>0:
                        if y == 0:
                            per_coor['stack_start_w'] = y*gap_w
                            per_coor['stack_end_w'] = y*gap_w+img_w-cut_w
                            per_coor['patch_start_w'] = 0
                            per_coor['patch_end_w'] = img_w-cut_w
                        elif y == num_w-1:
                            per_coor['stack_start_w'] = whole_w-img_w+cut_w
                            per_coor['stack_end_w'] = whole_w
                            per_coor['patch_start_w'] = cut_w
                            per_coor['patch_end_w'] = img_w
                        else:
                            per_coor['stack_start_w'] = y*gap_w+cut_w
                            per_coor['stack_end_w'] = y*gap_w+img_w-cut_w
                            per_coor['patch_start_w'] = cut_w
                            per_coor['patch_end_w'] = img_w-cut_w

                    if (num_h-1)==0:
                        per_coor['stack_start_h'] = 0
                        per_coor['stack_end_h'] = whole_h
                        per_coor['patch_start_h'] = 0
                        per_coor['patch_end_h'] = img_h
                    elif (num_h-1)>0:
                        if x == 0:
                            per_coor['stack_start_h'] = x*gap_h
                            per_coor['stack_end_h'] = x*gap_h+img_h-cut_h
                            per_coor['patch_start_h'] = 0
                            per_coor['patch_end_h'] = img_h-cut_h
                        elif x == num_h-1:
                            per_coor['stack_start_h'] = whole_h-img_h+cut_h
                            per_coor['stack_end_h'] = whole_h
                            per_coor['patch_start_h'] = cut_h
                            per_coor['patch_end_h'] = img_h
                        else:
                            per_coor['stack_start_h'] = x*gap_h+cut_h
                            per_coor['stack_end_h'] = x*gap_h+img_h-cut_h
                            per_coor['patch_start_h'] = cut_h
                            per_coor['patch_end_h'] = img_h-cut_h

                    if (num_s-1)==0:
                        per_coor['stack_start_s'] = 0
                        per_coor['stack_end_s'] = whole_s
                        per_coor['patch_start_s'] = 0
                        per_coor['patch_end_s'] = img_s2 
                    elif (num_s-1)>0:             
                        if z == 0:
                            per_coor['stack_start_s'] = z*gap_s2
                            per_coor['stack_end_s'] = z*gap_s2+img_s2-cut_s
                            per_coor['patch_start_s'] = 0
                            per_coor['patch_end_s'] = img_s2-cut_s
                        elif z == num_s-1:
                            per_coor['stack_start_s'] = whole_s-img_s2+cut_s
                            per_coor['stack_end_s'] = whole_s
                            per_coor['patch_start_s'] = cut_s
                            per_coor['patch_end_s'] = img_s2
                        else:
                            per_coor['stack_start_s'] = z*gap_s2+cut_s
                            per_coor['stack_end_s'] = z*gap_s2+img_s2-cut_s
                            per_coor['patch_start_s'] = cut_s
                            per_coor['patch_end_s'] = img_s2-cut_s
                    
                    per_coor['z'] = z

                    # noise_patch1 = noise_im[init_s:end_s,init_h:end_h,init_w:end_w]
                    patch_name = args['datasets_folder']+'_x'+str(x)+'_y'+str(y)+'_z'+str(z)
                    per_coor['name'] = im_name.replace('.tif','')
                    # print(' single_coordinate -----> ',single_coordinate)
                    per_coor_list.append(per_coor)
        # print(' per_coor_list -----> ',len(per_coor_list))
        '''
        while len(per_coor_list)%args.batch_size!=0:
            per_coor_list.append(per_coor)
            print(' per_coor_list -----> ',len(per_coor_list),' -----> ',len(per_coor_list)%args.batch_size,' -----> ',args.batch_size)
        '''
        coordinate_list[im_name.replace('.tif','')] = per_coor_list
    return  name_list, image_list, coordinate_list