import os
import sys
import io

import numpy as np
import torch
import cv2
from skimage import io

########################################################################################################################
def create_feature_maps(init_channel_number, number_of_fmaps):
    return [init_channel_number * 2 ** k for k in range(number_of_fmaps)]


def name2index(opt, input_name, num_h, num_w, num_s):
    # print(input_name)
    name_list = input_name.split('_')
    # print(name_list)
    z_part = name_list[-1]
    # print(z_part)
    y_part = name_list[-2]
    # print(y_part)
    x_part = name_list[-3]
    # print(x_part)
    z_index = int(z_part.replace('z',''))
    y_index = int(y_part.replace('y',''))
    x_index = int(x_part.replace('x',''))
    # print("x_index ---> ",x_index,"y_index ---> ", y_index,"z_index ---> ", z_index)

    cut_w = (opt.img_w - opt.gap_w)/2
    cut_h = (opt.img_h - opt.gap_h)/2
    cut_s = (opt.img_s - opt.gap_s)/2
    # print("z_index ---> ",cut_w, "cut_h ---> ",cut_h, "cut_s ---> ",cut_s)
    if x_index == 0:
        stack_start_w = x_index*opt.gap_w
        stack_end_w = x_index*opt.gap_w+opt.img_w-cut_w
        patch_start_w = 0
        patch_end_w = opt.img_w-cut_w
    elif x_index == num_w-1:
        stack_start_w = x_index*opt.gap_w+cut_w
        stack_end_w = x_index*opt.gap_w+opt.img_w
        patch_start_w = cut_w
        patch_end_w = opt.img_w
    else:
        stack_start_w = x_index*opt.gap_w+cut_w
        stack_end_w = x_index*opt.gap_w+opt.img_w-cut_w
        patch_start_w = cut_w
        patch_end_w = opt.img_w-cut_w

    if y_index == 0:
        stack_start_h = y_index*opt.gap_h
        stack_end_h = y_index*opt.gap_h+opt.img_h-cut_h
        patch_start_h = 0
        patch_end_h = opt.img_h-cut_h
    elif y_index == num_h-1:
        stack_start_h = y_index*opt.gap_h+cut_h
        stack_end_h = y_index*opt.gap_h+opt.img_h
        patch_start_h = cut_h
        patch_end_h = opt.img_h
    else:
        stack_start_h = y_index*opt.gap_h+cut_h
        stack_end_h = y_index*opt.gap_h+opt.img_h-cut_h
        patch_start_h = cut_h
        patch_end_h = opt.img_h-cut_h

    if z_index == 0:
        stack_start_s = z_index*opt.gap_s
        stack_end_s = z_index*opt.gap_s+opt.img_s-cut_s
        patch_start_s = 0
        patch_end_s = opt.img_s-cut_s
    elif z_index == num_s-1:
        stack_start_s = z_index*opt.gap_s+cut_s
        stack_end_s = z_index*opt.gap_s+opt.img_s
        patch_start_s = cut_s
        patch_end_s = opt.img_s
    else:
        stack_start_s = z_index*opt.gap_s+cut_s
        stack_end_s = z_index*opt.gap_s+opt.img_s-cut_s
        patch_start_s = cut_s
        patch_end_s = opt.img_s-cut_s
    return int(stack_start_w) ,int(stack_end_w) ,int(patch_start_w) ,int(patch_end_w) ,\
    int(stack_start_h) ,int(stack_end_h) ,int(patch_start_h) ,int(patch_end_h), \
    int(stack_start_s) ,int(stack_end_s) ,int(patch_start_s) ,int(patch_end_s)

def save_tiff_image(args, image_tensor, image_path):
    image = (image_tensor.cpu().detach().numpy()*args.normalize_factor) #.astype(np.uint16)
    # print('image max ',np.max(image),' image min ',np.min(image))
    '''
    image_tensor = torch.from_numpy(image_numpy)
    image_new_numpy = image_tensor.permute(2, 0, 1).numpy()
    # image_new_numpy = image_numpy
    if '_A' in image_path:
        image_new_numpy = image_new_numpy[...,np.newaxis] 
    '''
    save_tiff_path = image_path.replace('.png','.tif')
    io.imsave(save_tiff_path,image)

def save_numpy_image(args, image_tensor, image_path):
    image = (image_tensor*args.normalize_factor) #.astype(np.uint16)
    save_tiff_path = image_path.replace('.png','.tif')
    io.imsave(save_tiff_path,image)

def save_feature_tiff_image(image_tensor, image_path):
    image = image_tensor.cpu().detach().numpy()*255 #).astype(np.uint8)
    # image = np.clip(image, 0, 65535).astype('uint16')
    image = np.clip(image, 0, 255).astype('uint8')
    # print('image max ',np.max(image),' image min ',np.min(image))
    save_tiff_path = image_path.replace('.png','.tif')
    io.imsave(save_tiff_path,image)


def FFDrealign4(input):
    # batch channel time height width
    realign_input = torch.FloatTensor(input.shape[0], input.shape[1]*4, input.shape[2], int(input.shape[3]/2), int(input.shape[4]/2))
    # print('realign_input -----> ',realign_input.shape)
    # print('input -----> ',input.shape)
    # print('realign_input[:,0,:,:,:] -----> ',realign_input[:,0,:,:,:].shape)
    # print('input[:,:,:,0::2, 0::2] -----> ',input[:,:,:,0::2, 0::2].shape)
    realign_input[:,0,:,:,:] = torch.squeeze(input[:,:,:,0::2, 0::2])
    realign_input[:,1,:,:,:] = torch.squeeze(input[:,:,:,0::2, 1::2])
    realign_input[:,2,:,:,:] = torch.squeeze(input[:,:,:,1::2, 0::2])
    realign_input[:,3,:,:,:] = torch.squeeze(input[:,:,:,1::2, 1::2])
    return realign_input

def inv_FFDrealign4(input):
    # batch channel time height width
    realign_input = torch.FloatTensor(input.shape[0], int(input.shape[1]/4), input.shape[2], int(input.shape[3]*2), int(input.shape[4]*2))

    # print('realign_input[:,:,:,0::2, 0::2] -----> ',realign_input[:,:,:,0::2, 0::2].shape)
    # print('input[:,0,:,:,:] -----> ',input[:,0,:,:,:].shape)
    realign_input[:,:,:,0::2, 0::2] = torch.unsqueeze(input[:,0,:,:,:], 1)
    realign_input[:,:,:,0::2, 1::2] = torch.unsqueeze(input[:,1,:,:,:], 1)
    realign_input[:,:,:,1::2, 0::2] = torch.unsqueeze(input[:,2,:,:,:], 1)
    realign_input[:,:,:,1::2, 1::2] = torch.unsqueeze(input[:,3,:,:,:], 1)
    return realign_input