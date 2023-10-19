import argparse
import scipy.io as scio
import os
import time
import datetime
from skimage import io

from DWonder.MN.utils import list2contours, Joint_Mask_List_Simple_Merge_Mask, Joint_Mask_List_Mul, Joint_Mask_List_Simple
from DWonder.MN.utils import listAddcontours_Laplacian_pytorch, list2contours, listAddcontours_Laplacian_pytorch_Merge_Mask
############################################################################################################################################# 
parser = argparse.ArgumentParser()
parser.add_argument('--datasets_path', type=str, default='datasets', help="")
parser.add_argument('--datasets_folder1', type=str, default='test', help="")
parser.add_argument('--datasets_folder2', type=str, default='test', help="")
parser.add_argument('--output_dir', type=str, default='./results_mat', help="output directory")
opt = parser.parse_args()
############################################################################################################################################# 
if not os.path.exists(opt.output_dir): 
    os.mkdir(opt.output_dir)

current_time = datetime.datetime.now().strftime("%Y%m%d-%H%M")
SEG_output_path = opt.output_dir + '//' + 'R_'+opt.datasets_folder1+'_'+opt.datasets_folder2+'_'+current_time
if not os.path.exists(SEG_output_path): 
    os.mkdir(SEG_output_path)
SEG_f_con_output_path = SEG_output_path+ '//' + 'f_con'
if not os.path.exists(SEG_f_con_output_path): 
    os.mkdir(SEG_f_con_output_path)
SEG_mat_con_output_path = SEG_output_path+ '//' + 'mat'
if not os.path.exists(SEG_mat_con_output_path): 
    os.mkdir(SEG_mat_con_output_path)


im_folder1 = opt.datasets_path+'//'+opt.datasets_folder1
im_folder2 = opt.datasets_path+'//'+opt.datasets_folder2

for im_name1 in list(os.walk(im_folder1, topdown=False))[-1][-1]:
    data1 = scio.loadmat(im_folder1+'//'+im_name1)
    final_mask1 = data1['final_mask_list']

for im_name2 in list(os.walk(im_folder2, topdown=False))[-1][-1]:
    data2 = scio.loadmat(im_folder2+'//'+im_name2)
    final_mask2 = data2['final_mask_list']

# print(len(final_mask_list1[0]))
# print(len(final_mask_list2[0]))

final_mask_list1 = [] #1): #
for i in range(0,len(final_mask1[0])):
    neuron = final_mask1[0][i]
    neuron_centroid = neuron['centroid'][0][0][0]
    neuron_position = neuron['position'][0][0]
    neuron_trace = neuron['trace'][0][0][0]

    # print('neuron_position.shape[0]',neuron_position.shape[0])
    neuron_position_list =[] #1): #
    for p_i in range(0, neuron_position.shape[0]):
        now_position = list(neuron_position[p_i,:])
        neuron_position_list.append(now_position)
        # print(now_position[0], now_position[1])

    new_neuron = {}
    new_neuron['centroid'] = neuron_centroid
    new_neuron['position'] = neuron_position_list
    new_neuron['trace'] = neuron_trace

    # print(neuron_centroid)
    # print(neuron_trace)
    # print(neuron_position)
    '''
    for p_i in range(0,1): # len(position)):
        now_position = neuron_position[p_i]
        print(now_position[0], now_position[1])
    '''

    # print(neuron_position)
    # print(neuron_centroid[0])
    # print(neuron_centroid[0][0])
    # print(neuron_centroid[0][0][0])
    final_mask_list1.append(new_neuron)


final_mask_list2 = []
for i in range(0,len(final_mask2[0])):
    neuron = final_mask2[0][i]
    neuron_centroid = neuron['centroid'][0][0][0]
    neuron_position = neuron['position'][0][0]
    neuron_trace = neuron['trace'][0][0][0]

    # print('neuron_position.shape[0]',neuron_position.shape[0])
    neuron_position_list =[] #1): #
    for p_i in range(0, neuron_position.shape[0]):
        now_position = list(neuron_position[p_i,:])
        neuron_position_list.append(now_position)
        # print(now_position)

    new_neuron = {}
    new_neuron['centroid'] = neuron_centroid
    new_neuron['position'] = neuron_position_list
    new_neuron['trace'] = neuron_trace
    # print(neuron['centroid'])
    final_mask_list2.append(new_neuron)

print(len(final_mask_list1))
print(len(final_mask_list2))
# final_mask_list1[0]+final_mask_list2[0]

# final_mask_list = Joint_Mask_List_Simple_Merge_Mask(final_mask_list, [], corr_mark=0.7, area_mark=0.9, active_rate=0, if_coor=False, if_area=True, if_merge=False)
final_mask_list = Joint_Mask_List_Simple(final_mask_list1, final_mask_list2, corr_mark=0.7, area_mark=0.7, active_rate=0, if_coor=False, if_area=True, if_merge=True)

print(len(final_mask_list))

final_mask_list = listAddcontours_Laplacian_pytorch(final_mask_list, 2160, 2560)
final_contours,whole_contours = list2contours(final_mask_list, 2160, 2560)

img_f_contours_name = SEG_f_con_output_path+'//'+opt.datasets_folder1+'_'+opt.datasets_folder2+'_f_con.tif'
final_contours = final_contours.astype('uint16')
io.imsave(img_f_contours_name, final_contours)


mat_save_name = SEG_mat_con_output_path+'//'+opt.datasets_folder1+'_'+opt.datasets_folder2+'.mat'
data = {'a':final_mask_list, 'final_contours':final_contours}
scio.savemat(mat_save_name, {'final_mask_list':data['a'], 'final_contours':data['final_contours']})
'''
'''