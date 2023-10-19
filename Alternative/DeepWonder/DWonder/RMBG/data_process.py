import numpy as np
import os
import tifffile as tiff
import random
import math
import torch
from torch.utils.data import Dataset

###
#  this file is for pre-processing images for training and test.
#  currently we use splitted dataset


###############################
# dataset class

class trainset(Dataset):
    def __init__(self,name_list,coordinate_list,noise_img, noise_img2):
        self.name_list = name_list
        self.coordinate_list=coordinate_list
        self.noise_img = noise_img
        self.noise_img2 = noise_img2

    def __getitem__(self, index):
        #fn = self.images[index]
        single_coordinate = self.coordinate_list[self.name_list[index]]
        init_h = single_coordinate['init_h']
        end_h = single_coordinate['end_h']
        init_w = single_coordinate['init_w']
        end_w = single_coordinate['end_w']
        init_s = single_coordinate['init_s']
        end_s = single_coordinate['end_s']
        input = self.noise_img[init_s: end_s : 1, init_h:end_h, init_w:end_w]
        target = self.noise_img2[init_s: end_s : 1, init_h:end_h, init_w:end_w]

        input=torch.from_numpy(np.expand_dims(input, 0))
        target=torch.from_numpy(np.expand_dims(target, 0))
        #target = self.target[index]
        return input, target

    def __len__(self):
        return len(self.name_list)


class trainsetMul(Dataset):
    def __init__(self,coor_list, GT_list, input_list, opt):
        self.coor_list = coor_list
        self.GT_list = GT_list
        self.input_list = input_list
        self.opt = opt

    def __getitem__(self, index):
        #fn = self.images[index]
        per_coor = self.coor_list[index]
        train_im_name = per_coor['name']
        init_w = per_coor['init_w']
        init_h = per_coor['init_h']
        init_s = per_coor['init_s']
        GT_im = self.GT_list[train_im_name]
        input_im = self.input_list[train_im_name]
        
        GT_patch =       GT_im[init_s:init_s+self.opt.img_s, init_w:init_w+self.opt.img_w, init_h:init_h+self.opt.img_h]
        input_patch = input_im[init_s:init_s+self.opt.img_s, init_w:init_w+self.opt.img_w, init_h:init_h+self.opt.img_h]

        # input = self.noise_img[init_s: end_s : 1, init_h:end_h, init_w:end_w]
        # target = self.noise_img2[init_s: end_s : 1, init_h:end_h, init_w:end_w]

        input=torch.from_numpy(np.expand_dims(input_patch, 0))
        target=torch.from_numpy(np.expand_dims(GT_patch, 0))
        #target = self.target[index]
        return input, target

    def __len__(self):
        return len(self.coor_list)


class testset(Dataset):
    def __init__(self,name_list,coordinate_list,noise_img):
        self.name_list = name_list
        self.coordinate_list=coordinate_list
        self.noise_img = noise_img

    def __getitem__(self, index):
        #fn = self.images[index]
        single_coordinate = self.coordinate_list[self.name_list[index]]
        init_h = single_coordinate['init_h']
        end_h = single_coordinate['end_h']
        init_w = single_coordinate['init_w']
        end_w = single_coordinate['end_w']
        init_s = single_coordinate['init_s']
        end_s = single_coordinate['end_s']
        noise_patch = self.noise_img[init_s:end_s, init_h:end_h, init_w:end_w]
        noise_patch=torch.from_numpy(np.expand_dims(noise_patch, 0))
        #target = self.target[index]
        return noise_patch,single_coordinate

    def __len__(self):
        return len(self.name_list)
###############################
# save
def singlebatch_test_save(single_coordinate,output_image,raw_image):
    stack_start_w = int(single_coordinate['stack_start_w'])
    stack_end_w = int(single_coordinate['stack_end_w'])
    patch_start_w = int(single_coordinate['patch_start_w'])
    patch_end_w = int(single_coordinate['patch_end_w'])

    stack_start_h = int(single_coordinate['stack_start_h'])
    stack_end_h = int(single_coordinate['stack_end_h'])
    patch_start_h = int(single_coordinate['patch_start_h'])
    patch_end_h = int(single_coordinate['patch_end_h'])

    stack_start_s = int(single_coordinate['stack_start_s'])
    stack_end_s = int(single_coordinate['stack_end_s'])
    patch_start_s = int(single_coordinate['patch_start_s'])
    patch_end_s = int(single_coordinate['patch_end_s'])
    '''
    print('RMBG stack_start_w -----> ',stack_start_w, ' patch_start_w -----> ', \
    patch_start_w,'stack_end_w -----> ',stack_end_w, ' patch_end_w -----> ',patch_end_w)
    print('RMBG stack_start_h -----> ',stack_start_h, ' patch_start_h -----> ', \
    patch_start_h, 'stack_end_h -----> ',stack_end_h, ' patch_end_h -----> ',patch_end_h)
    '''
    aaaa = output_image[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
    bbbb = raw_image[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
    return aaaa,bbbb,stack_start_w,stack_end_w,stack_start_h,stack_end_h,stack_start_s,stack_end_s


def multibatch_test_save(single_coordinate,id,output_image,raw_image):
    stack_start_w_id = single_coordinate['stack_start_w'].numpy()
    stack_start_w = int(stack_start_w_id[id])
    stack_end_w_id = single_coordinate['stack_end_w'].numpy()
    stack_end_w=int(stack_end_w_id[id])
    patch_start_w_id = single_coordinate['patch_start_w'].numpy()
    patch_start_w=int(patch_start_w_id[id])
    patch_end_w_id = single_coordinate['patch_end_w'].numpy()
    patch_end_w=int(patch_end_w_id[id])

    stack_start_h_id = single_coordinate['stack_start_h'].numpy()
    stack_start_h = int(stack_start_h_id[id])
    stack_end_h_id = single_coordinate['stack_end_h'].numpy()
    stack_end_h = int(stack_end_h_id[id])
    patch_start_h_id = single_coordinate['patch_start_h'].numpy()
    patch_start_h = int(patch_start_h_id[id])
    patch_end_h_id = single_coordinate['patch_end_h'].numpy()
    patch_end_h = int(patch_end_h_id[id])

    stack_start_s_id = single_coordinate['stack_start_s'].numpy()
    stack_start_s = int(stack_start_s_id[id])
    stack_end_s_id = single_coordinate['stack_end_s'].numpy()
    stack_end_s = int(stack_end_s_id[id])
    patch_start_s_id = single_coordinate['patch_start_s'].numpy()
    patch_start_s = int(patch_start_s_id[id])
    patch_end_s_id = single_coordinate['patch_end_s'].numpy()
    patch_end_s = int(patch_end_s_id[id])

    output_image_id=output_image[id]
    raw_image_id=raw_image[id]
    aaaa = output_image_id[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
    bbbb = raw_image_id[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]

    return aaaa,bbbb,stack_start_w,stack_end_w,stack_start_h,stack_end_h,stack_start_s,stack_end_s

###############################
# shuffle
def shuffle_datasets(train_raw, train_GT, name_list):
    # target: train_raw, train_GT, name_list
    index_list = list(range(0, len(name_list)))
    # print('index_list -----> ',index_list)
    random.shuffle(index_list)
    random_index_list = index_list
    # print('index_list -----> ',index_list)
    new_name_list = list(range(0, len(name_list)))
    train_raw = np.array(train_raw)
    # print('train_raw shape -----> ',train_raw.shape)
    train_GT = np.array(train_GT)
    # print('train_GT shape -----> ',train_GT.shape)
    new_train_raw = train_raw
    new_train_GT = train_GT
    for i in range(0,len(random_index_list)):
        # print('i -----> ',i)
        new_train_raw[i,:,:,:] = train_raw[random_index_list[i],:,:,:] 
        new_train_GT[i,:,:,:] = train_GT[random_index_list[i],:,:,:]
        new_name_list[i] = name_list[random_index_list[i]]
    # new_train_raw = np.expand_dims(new_train_raw, 4)
    # new_train_GT = np.expand_dims(new_train_GT, 4)
    return new_train_raw, new_train_GT, new_name_list

###############################
# pre-process
def get_gap_s(args, img, stack_num):
    whole_w = img.shape[2]
    whole_h = img.shape[1]
    whole_s = img.shape[0]
    # print('whole_w -----> ',whole_w)
    # print('whole_h -----> ',whole_h)
    # print('whole_s -----> ',whole_s)
    w_num = math.floor((whole_w-args.img_w)/args.gap_w)+1
    h_num = math.floor((whole_h-args.img_h)/args.gap_h)+1
    s_num = math.ceil(args.train_datasets_size/w_num/h_num/stack_num)
    # print('w_num -----> ',w_num)
    # print('h_num -----> ',h_num)
    # print('s_num -----> ',s_num)
    gap_s = math.floor((whole_s-args.img_s*2)/(s_num-1))
    # print('gap_s -----> ',gap_s)
    return gap_s


def test_preprocess_lessMemoryNoTail_chooseOne (args, N):
    img_h = args.img_h
    img_w = args.img_w
    img_s2 = args.img_s
    gap_h = args.gap_h
    gap_w = args.gap_w
    gap_s2 = args.gap_s
    cut_w = (img_w - gap_w)/2
    cut_h = (img_h - gap_h)/2
    cut_s = (img_s2 - gap_s2)/2
    im_folder = args.datasets_path+'//'+args.datasets_folder

    name_list = []
    # train_raw = []
    coordinate_list={}
    img_list = list(os.walk(im_folder, topdown=False))[-1][-1]
    img_list.sort()
    # print(img_list)

    im_name = img_list[N]

    im_dir = im_folder+'//'+im_name
    print('im_dir -----> ',im_dir)
    noise_im = tiff.imread(im_dir)
    # print('noise_im shape -----> ',noise_im.shape)
    # print('noise_im max -----> ',noise_im.max())
    # print('noise_im min -----> ',noise_im.min())
    if noise_im.shape[0]>args.test_datasize:
        noise_im = noise_im[0:args.test_datasize,:,:]
    noise_im = (noise_im).astype(np.float32)/args.normalize_factor

    noise_im_ave_single = np.mean(noise_im, axis=0)
    noise_im_ave = np.zeros(noise_im.shape)
    for i in range(0, noise_im.shape[0]):
        noise_im_ave[i,:,:] = noise_im_ave_single
    noise_im = noise_im-noise_im_ave


    whole_w = noise_im.shape[2]
    whole_h = noise_im.shape[1]
    whole_s = noise_im.shape[0]

    num_w = math.ceil((whole_w-img_w+gap_w)/gap_w)
    num_h = math.ceil((whole_h-img_h+gap_h)/gap_h)
    num_s = math.ceil((whole_s-img_s2+gap_s2)/gap_s2)
    # print('int((whole_h-img_h+gap_h)/gap_h) -----> ',int((whole_h-img_h+gap_h)/gap_h))
    # print('int((whole_w-img_w+gap_w)/gap_w) -----> ',int((whole_w-img_w+gap_w)/gap_w))
    # print('int((whole_s-img_s2+gap_s2)/gap_s2) -----> ',int((whole_s-img_s2+gap_s2)/gap_s2))
    for x in range(0,num_h):
        for y in range(0,num_w):
            for z in range(0,num_s):
                single_coordinate={'init_h':0, 'end_h':0, 'init_w':0, 'end_w':0, 'init_s':0, 'end_s':0}
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
                single_coordinate['init_h'] = init_h
                single_coordinate['end_h'] = end_h
                single_coordinate['init_w'] = init_w
                single_coordinate['end_w'] = end_w
                single_coordinate['init_s'] = init_s
                single_coordinate['end_s'] = end_s

                if (num_w-1)==0:
                    single_coordinate['stack_start_w'] = 0
                    single_coordinate['stack_end_w'] = whole_w
                    single_coordinate['patch_start_w'] = 0
                    single_coordinate['patch_end_w'] = img_w
                elif (num_w-1)>0:
                    if y == 0:
                        single_coordinate['stack_start_w'] = y*gap_w
                        single_coordinate['stack_end_w'] = y*gap_w+img_w-cut_w
                        single_coordinate['patch_start_w'] = 0
                        single_coordinate['patch_end_w'] = img_w-cut_w
                    elif y == num_w-1:
                        single_coordinate['stack_start_w'] = whole_w-img_w+cut_w
                        single_coordinate['stack_end_w'] = whole_w
                        single_coordinate['patch_start_w'] = cut_w
                        single_coordinate['patch_end_w'] = img_w
                    else:
                        single_coordinate['stack_start_w'] = y*gap_w+cut_w
                        single_coordinate['stack_end_w'] = y*gap_w+img_w-cut_w
                        single_coordinate['patch_start_w'] = cut_w
                        single_coordinate['patch_end_w'] = img_w-cut_w

                if (num_h-1)==0:
                    single_coordinate['stack_start_h'] = 0
                    single_coordinate['stack_end_h'] = whole_h
                    single_coordinate['patch_start_h'] = 0
                    single_coordinate['patch_end_h'] = img_h
                elif (num_h-1)>0:
                    if x == 0:
                        single_coordinate['stack_start_h'] = x*gap_h
                        single_coordinate['stack_end_h'] = x*gap_h+img_h-cut_h
                        single_coordinate['patch_start_h'] = 0
                        single_coordinate['patch_end_h'] = img_h-cut_h
                    elif x == num_h-1:
                        single_coordinate['stack_start_h'] = whole_h-img_h+cut_h
                        single_coordinate['stack_end_h'] = whole_h
                        single_coordinate['patch_start_h'] = cut_h
                        single_coordinate['patch_end_h'] = img_h
                    else:
                        single_coordinate['stack_start_h'] = x*gap_h+cut_h
                        single_coordinate['stack_end_h'] = x*gap_h+img_h-cut_h
                        single_coordinate['patch_start_h'] = cut_h
                        single_coordinate['patch_end_h'] = img_h-cut_h

                if (num_s-1)==0:
                    single_coordinate['stack_start_s'] = 0
                    single_coordinate['stack_end_s'] = whole_s
                    single_coordinate['patch_start_s'] = 0
                    single_coordinate['patch_end_s'] = img_s2
                elif (num_s-1)>0:
                    if z == 0:
                        single_coordinate['stack_start_s'] = z*gap_s2
                        single_coordinate['stack_end_s'] = z*gap_s2+img_s2-cut_s
                        single_coordinate['patch_start_s'] = 0
                        single_coordinate['patch_end_s'] = img_s2-cut_s
                    elif z == num_s-1:
                        single_coordinate['stack_start_s'] = whole_s-img_s2+cut_s
                        single_coordinate['stack_end_s'] = whole_s
                        single_coordinate['patch_start_s'] = cut_s
                        single_coordinate['patch_end_s'] = img_s2
                    else:
                        single_coordinate['stack_start_s'] = z*gap_s2+cut_s
                        single_coordinate['stack_end_s'] = z*gap_s2+img_s2-cut_s
                        single_coordinate['patch_start_s'] = cut_s
                        single_coordinate['patch_end_s'] = img_s2-cut_s

                print('RMBG stack_start_w -----> ',single_coordinate['stack_start_w'], ' patch_start_w -----> ', single_coordinate['patch_start_w'],  \
                'stack_end_w -----> ',single_coordinate['stack_end_w'], ' patch_end_w -----> ',single_coordinate['patch_end_w'])
                print('RMBG stack_start_h -----> ',single_coordinate['stack_start_h'], ' patch_start_h -----> ', single_coordinate['patch_start_h'],  \
                'stack_end_h -----> ',single_coordinate['stack_end_h'], ' patch_end_h -----> ',single_coordinate['patch_end_h'])
                # noise_patch1 = noise_im[init_s:end_s,init_h:end_h,init_w:end_w]
                patch_name = args.datasets_folder+'_x'+str(x)+'_y'+str(y)+'_z'+str(z)
                # train_raw.append(noise_patch1.transpose(1,2,0))
                name_list.append(patch_name)
                # print(' single_coordinate -----> ',single_coordinate)
                coordinate_list[patch_name] = single_coordinate

    return  name_list, noise_im, coordinate_list



def test_preprocess_lessMemoryNoTail_SubImg (args, sub_img):
    img_h = args['img_h']
    img_w = args['img_w']
    img_s2 = args['img_s']
    gap_h = args['gap_h']
    gap_w = args['gap_w']
    gap_s2 = args['gap_s']
    cut_w = (img_w - gap_w)/2
    cut_h = (img_h - gap_h)/2
    cut_s = (img_s2 - gap_s2)/2

    noise_im = sub_img
    noise_im = (noise_im).astype(np.float32)/args['normalize_factor']

    noise_im_ave_single = np.mean(noise_im, axis=0)
    noise_im_ave = np.zeros(noise_im.shape)
    for i in range(0, noise_im.shape[0]):
        noise_im_ave[i,:,:] = noise_im_ave_single
    noise_im = noise_im-noise_im_ave

    whole_w = noise_im.shape[2]
    whole_h = noise_im.shape[1]
    whole_s = noise_im.shape[0]

    num_w = math.ceil((whole_w-img_w+gap_w)/gap_w)
    num_h = math.ceil((whole_h-img_h+gap_h)/gap_h)
    num_s = math.ceil((whole_s-img_s2+gap_s2)/gap_s2)
    # print('int((whole_h-img_h+gap_h)/gap_h) -----> ',int((whole_h-img_h+gap_h)/gap_h))
    # print('int((whole_w-img_w+gap_w)/gap_w) -----> ',int((whole_w-img_w+gap_w)/gap_w))
    # print('int((whole_s-img_s2+gap_s2)/gap_s2) -----> ',int((whole_s-img_s2+gap_s2)/gap_s2))
    name_list = []
    coordinate_list = {}
    for x in range(0,num_h):
        for y in range(0,num_w):
            for z in range(0,num_s):
                single_coordinate={'init_h':0, 'end_h':0, 'init_w':0, 'end_w':0, 'init_s':0, 'end_s':0}
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
                single_coordinate['init_h'] = init_h
                single_coordinate['end_h'] = end_h
                single_coordinate['init_w'] = init_w
                single_coordinate['end_w'] = end_w
                single_coordinate['init_s'] = init_s
                single_coordinate['end_s'] = end_s

                if (num_w-1)==0:
                    single_coordinate['stack_start_w'] = 0
                    single_coordinate['stack_end_w'] = whole_w
                    single_coordinate['patch_start_w'] = 0
                    single_coordinate['patch_end_w'] = img_w
                elif (num_w-1)>0:
                    if y == 0:
                        single_coordinate['stack_start_w'] = y*gap_w
                        single_coordinate['stack_end_w'] = y*gap_w+img_w-cut_w
                        single_coordinate['patch_start_w'] = 0
                        single_coordinate['patch_end_w'] = img_w-cut_w
                    elif y == num_w-1:
                        single_coordinate['stack_start_w'] = whole_w-img_w+cut_w
                        single_coordinate['stack_end_w'] = whole_w
                        single_coordinate['patch_start_w'] = cut_w
                        single_coordinate['patch_end_w'] = img_w
                    else:
                        single_coordinate['stack_start_w'] = y*gap_w+cut_w
                        single_coordinate['stack_end_w'] = y*gap_w+img_w-cut_w
                        single_coordinate['patch_start_w'] = cut_w
                        single_coordinate['patch_end_w'] = img_w-cut_w

                if (num_h-1)==0:
                    single_coordinate['stack_start_h'] = 0
                    single_coordinate['stack_end_h'] = whole_h
                    single_coordinate['patch_start_h'] = 0
                    single_coordinate['patch_end_h'] = img_h
                elif (num_h-1)>0:
                    if x == 0:
                        single_coordinate['stack_start_h'] = x*gap_h
                        single_coordinate['stack_end_h'] = x*gap_h+img_h-cut_h
                        single_coordinate['patch_start_h'] = 0
                        single_coordinate['patch_end_h'] = img_h-cut_h
                    elif x == num_h-1:
                        single_coordinate['stack_start_h'] = whole_h-img_h+cut_h
                        single_coordinate['stack_end_h'] = whole_h
                        single_coordinate['patch_start_h'] = cut_h
                        single_coordinate['patch_end_h'] = img_h
                    else:
                        single_coordinate['stack_start_h'] = x*gap_h+cut_h
                        single_coordinate['stack_end_h'] = x*gap_h+img_h-cut_h
                        single_coordinate['patch_start_h'] = cut_h
                        single_coordinate['patch_end_h'] = img_h-cut_h

                if (num_s-1)==0:
                    single_coordinate['stack_start_s'] = 0
                    single_coordinate['stack_end_s'] = whole_s
                    single_coordinate['patch_start_s'] = 0
                    single_coordinate['patch_end_s'] = img_s2
                elif (num_s-1)>0:
                    if z == 0:
                        single_coordinate['stack_start_s'] = z*gap_s2
                        single_coordinate['stack_end_s'] = z*gap_s2+img_s2-cut_s
                        single_coordinate['patch_start_s'] = 0
                        single_coordinate['patch_end_s'] = img_s2-cut_s
                    elif z == num_s-1:
                        single_coordinate['stack_start_s'] = whole_s-img_s2+cut_s
                        single_coordinate['stack_end_s'] = whole_s
                        single_coordinate['patch_start_s'] = cut_s
                        single_coordinate['patch_end_s'] = img_s2
                    else:
                        single_coordinate['stack_start_s'] = z*gap_s2+cut_s
                        single_coordinate['stack_end_s'] = z*gap_s2+img_s2-cut_s
                        single_coordinate['patch_start_s'] = cut_s
                        single_coordinate['patch_end_s'] = img_s2-cut_s
                '''
                print('RMBG stack_start_w -----> ',single_coordinate['stack_start_w'], ' patch_start_w -----> ', single_coordinate['patch_start_w'],  \
                'stack_end_w -----> ',single_coordinate['stack_end_w'], ' patch_end_w -----> ',single_coordinate['patch_end_w'])
                print('RMBG stack_start_h -----> ',single_coordinate['stack_start_h'], ' patch_start_h -----> ', single_coordinate['patch_start_h'],  \
                'stack_end_h -----> ',single_coordinate['stack_end_h'], ' patch_end_h -----> ',single_coordinate['patch_end_h'])
                '''
                # noise_patch1 = noise_im[init_s:end_s,init_h:end_h,init_w:end_w]
                patch_name = 'x'+str(x)+'_y'+str(y)+'_z'+str(z)
                # train_raw.append(noise_patch1.transpose(1,2,0))
                name_list.append(patch_name)
                # print(' single_coordinate -----> ',single_coordinate)
                coordinate_list[patch_name] = single_coordinate

    return  name_list, noise_im, coordinate_list



def train_preprocess_lessMemory_segBG(args):
    input_pretype = 'mean'
    GT_pretype = 'min'

    name_list = []
    GT_list={}
    input_list={}
    GT_folder = args.datasets_path+'//'+args.datasets_folder+'//'+args.GT_folder
    input_folder = args.datasets_path+'//'+args.datasets_folder+'//'+args.input_folder
    print('input_folder -----> ',input_folder)

    img_num = len(list(os.walk(input_folder, topdown=False))[-1][-1])
    for i in range(0, img_num):
        input_name = list(os.walk(input_folder, topdown=False))[-1][-1][i]
        GT_name = list(os.walk(GT_folder, topdown=False))[-1][-1][i]
        
        input_dir = input_folder+'//'+input_name
        GT_dir = GT_folder+'//'+GT_name
        print('read im_name -----> ',input_dir)
        print('read im_name -----> ',GT_dir)
        input_img = tiff.imread(input_dir)
        GT_img = tiff.imread(GT_dir)
        # im = im.transpose(2,0,1)
        print(input_img.shape)

        # GT_img = np.expand_dims(GT_img, axis=0)
        # input_img = input_img.transpose(2,0,1)
        input_img = input_img.astype(np.float32)/args.normalize_factor

        if input_pretype == 'mean':
            print('input_pretype == mean')
            input_img_ave_single = np.mean(input_img, axis=0)
            input_img_ave = np.zeros(input_img.shape)
            for i in range(0, input_img.shape[0]):
                input_img_ave[i,:,:] = input_img_ave_single
            input_img = input_img-input_img_ave


        if GT_pretype == 'min':
            print('GT_pretype == min')
            GT_img_ave_single = np.min(GT_img, axis=0)
            GT_img_ave = np.zeros(GT_img.shape,dtype='float32')
            for i in range(0, GT_img.shape[0]):
                GT_img_ave[i,:,:] = GT_img_ave_single
            GT_img = GT_img-GT_img_ave

        GT_list[GT_name.replace('.tif','')] = GT_img
        input_list[input_name.replace('.tif','')] = input_img
        name_list.append(GT_name.replace('.tif',''))
    
    num_per_img = math.ceil(args.train_datasets_size/img_num)
    coor_list = []
    for i in range(0, img_num):
        for ii in range(0, num_per_img):
            per_coor = {}
            init_w = np.random.randint(0, GT_img.shape[-2]-args.img_w-1)
            init_h = np.random.randint(0, GT_img.shape[-1]-args.img_h-1)
            init_s = np.random.randint(0, GT_img.shape[0]-args.img_s-1)

            per_coor['name'] = name_list[i]
            per_coor['init_w'] = init_w
            per_coor['init_h'] = init_h
            per_coor['init_s'] = init_s
            coor_list.append(per_coor)
    return  coor_list, GT_list, input_list


def shuffle_datasets_lessMemory(name_list):
    index_list = list(range(0, len(name_list)))
    # print('index_list -----> ',index_list)
    random.shuffle(index_list)
    random_index_list = index_list
    # print('index_list -----> ',index_list)
    new_name_list = list(range(0, len(name_list)))
    for i in range(0,len(random_index_list)):
        new_name_list[i] = name_list[random_index_list[i]]
    return new_name_list