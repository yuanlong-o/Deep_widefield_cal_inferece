import os
import sys
import time
'''
flag = sys.argv[1]

###################################################################################################################################################################
# Only train, using 1 GPU and batch_size=1
if flag == 'train':
    os.system('python train.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 0 --batch_size 1 \
                               --fmap 16 \
                               --img_h 200 --img_w 200 --img_s 128 \
                               --train_datasets_size 1000')  

if flag == 'test':
    os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109111033 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 0 --batch_size 1 \
                              --test_datasize 2000')


if flag == 'test_eval':
    os.system('python test_and_eval.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109111033 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --pixel_size 0.8 \
                              --GPU 0 --batch_size 1 \
                              --test_datasize 2000')
   
os.system('python train.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 0 --batch_size 1 \
                               --fmap 16 \
                               --img_h 128 --img_w 128 --img_s 64 \
                               --train_datasets_size 1000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109152243 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 0 --batch_size 1 \
                              --img_h 128 --img_w 128 --img_s 64 \
                              --gap_h 96 --gap_w 96 --gap_s 32 \
                              --test_datasize 200')

os.system('python train.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 0 --batch_size 1 \
                               --fmap 16 \
                               --img_h 128 --img_w 128 --img_s 64 \
                               --train_datasets_size 1000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109162324 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 0 --batch_size 1 \
                              --img_h 128 --img_w 128 --img_s 64 \
                              --gap_h 96 --gap_w 96 --gap_s 32 \
                              --test_datasize 200')

os.system('python train.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 0 --batch_size 1 \
                               --fmap 16 \
                               --img_h 128 --img_w 128 --img_s 64 \
                               --train_datasets_size 1000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109162324 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 1 --batch_size 1 \
                              --img_h 256 --img_w 256 --img_s 64 \
                              --gap_h 192 --gap_w 192 --gap_s 32 \
                              --test_datasize 200')

os.system('python train.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 0 --batch_size 1 \
                               --fmap 16 \
                               --img_h 256 --img_w 256 --img_s 64 \
                               --train_datasets_size 1000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109181041 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 1 --batch_size 1 \
                              --img_h 256 --img_w 256 --img_s 64 \
                              --gap_h 192 --gap_w 192 --gap_s 32 \
                              --test_datasize 2000')

os.system('python train.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 2 --batch_size 1 \
                               --fmap 16 \
                               --img_h 128 --img_w 128 --img_s 192 \
                               --train_datasets_size 1000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109211340 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 2 --batch_size 1 \
                              --img_h 128 --img_w 128 --img_s 192 \
                              --gap_h 96 --gap_w 96 --gap_s 160 \
                              --test_datasize 300')

os.system('python train_GAN.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 2 --batch_size 1 \
                               --fmap 16 \
                               --img_h 256 --img_w 256 --img_s 64 \
                               --train_datasets_size 1000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109221050 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 2 --batch_size 1 \
                              --img_h 256 --img_w 256 --img_s 64 \
                              --gap_h 192 --gap_w 192 --gap_s 32 \
                              --test_datasize 2000')vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_rmbg_cut

os.system('python train_GAN.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                               --n_epochs 50 --GPU 2 --batch_size 1 \
                               --fmap 16 \
                               --img_h 256 --img_w 256 --img_s 64 \
                               --train_datasets_size 1000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109231359 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 2 --batch_size 1 \
                              --img_h 256 --img_w 256 --img_s 64 \
                              --gap_h 192 --gap_w 192 --gap_s 32 \
                              --test_datasize 2000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_rmbg_cut_202109241139 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2_rmbg_cut \
                              --GPU 2 --batch_size 1 \
                              --img_h 256 --img_w 256 --img_s 64 \
                              --gap_h 192 --gap_w 192 --gap_s 32 \
                              --test_datasize 2000')

os.system('python train_GAN_FFD.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --n_epochs 30 --GPU 2 --batch_size 1 \
                               --fmap 16 --lr 0.0002 \
                               --img_h 256 --img_w 256 --img_s 128 \
                               --train_datasets_size 1000')

os.system('python test_FFD.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109271317 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 2 --batch_size 1 \
                              --img_h 256 --img_w 256 --img_s 128 \
                              --gap_h 224 --gap_w 224 --gap_s 96 \
                              --test_datasize 2000')

os.system('python test.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_rmbg_cut_202109241139 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2_rmbg_cut \
                              --GPU 2 --batch_size 1 \
                              --img_h 256 --img_w 256 --img_s 64 \
                              --gap_h 224 --gap_w 224 --gap_s 32 \
                              --test_datasize 2000')

os.system('python train_GAN_FFD.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --n_epochs 30 --GPU 2 --batch_size 1 \
                               --fmap 32 --lr 0.0002 \
                               --img_h 256 --img_w 256 --img_s 128 \
                               --train_datasets_size 1000')

os.system('python test_FFD.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1_202109271949 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_2 \
                              --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h 256 --img_w 256 --img_s 128 \
                              --gap_h 224 --gap_w 224 --gap_s 96 \
                              --test_datasize 2000')

os.system('python train_GAN_FFDMul.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr \
                                --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_2_id_1 \
                               --n_epochs 30 --GPU 2 --batch_size 1 \
                               --fmap 32 --lr 0.0002 \
                               --img_h 256 --img_w 256 --img_s 128 \
                               --train_datasets_size 2000')

os.system('python test_FFD.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_202109282316 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                              --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h 256 --img_w 256 --img_s 128 \
                              --gap_h 224 --gap_w 224 --gap_s 96 \
                              --test_datasize 2000')

os.system('python train_GAN_FFDMul.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                               --n_epochs 30 --GPU 1 --batch_size 1 \
                               --fmap 32 --lr 0.0002 \
                               --img_h 256 --img_w 256 --img_s 128 \
                               --train_datasets_size 2000')

os.system('python test_FFD.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_202109301706 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h 256 --img_w 256 --img_s 128 \
                              --gap_h 224 --gap_w 224 --gap_s 96 \
                              --test_datasize 2000')

os.system('python train_GAN_FFDMul.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                               --n_epochs 30 --GPU 1 --batch_size 1 \
                               --fmap 32 --lr 0.0002 \
                               --img_h 256 --img_w 256 --img_s 128 \
                               --train_datasets_size 4000')
os.system('python test_FFD.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_202110011541 \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --GPU 1 --batch_size 1 --fmap 32 \
                              --img_h 256 --img_w 256 --img_s 128 \
                              --gap_h 224 --gap_w 224 --gap_s 96 \
                              --test_datasize 4000')
os.system('python train_GAN_FFDMul.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                               --n_epochs 100 --GPU 1 --batch_size 1 \
                               --fmap 32 --lr 0.0001 \
                               --img_h 256 --img_w 256 --img_s 256 \
                               --train_datasets_size 4000') 
# vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_202111031743_nomean
# vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_202110011541

os.system('python test_FFD_feature.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_202111031743_nomean \
                              --datasets_folder 23_1_inter \
                              --GPU 1 --batch_size 1 --fmap 32 \
                              --img_h 256 --img_w 256 --img_s 256 \
                              --gap_h 224 --gap_w 224 --gap_s 224 \ vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr
                              --test_datasize 300') vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr

# stime.sleep(3600)  vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr

os.system('python train_GAN_FFDMul.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 16 \
                                 --img_h 32 --img_w 32 --img_s 256 \
                                 --train_datasets_size 2000') 

os.system('python train_GAN_FFDMul.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 6 \
                                 --img_h 64 --img_w 64 --img_s 256 \
                                 --train_datasets_size 2000') 

os.system('python train_GAN_FFDMul.py --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 10 \
                                 --img_h 48 --img_w 48 --img_s 256 \
                                 --train_datasets_size 2000')

i=1
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*128)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=2
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*128)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=4
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*96)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=6
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*80)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=8
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*48)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=10
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*48)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=12
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*32)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=14
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*32)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

i=16
os.system('python test_FFD1.py --denoise_model vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr'+str(i)+' \
                              --datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                              --pth_path pth_step --GPU 2 --batch_size 1 --fmap 32 \
                              --img_h '+str(i*32)+' --img_s 128 \
                              --test_datasize 2000 --realign_rate '+str(i))

os.system('python train_GAN_FFDMul.py --datasets_folder 3_wfto2p \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 2 \
                                 --img_h 128 --img_w 128 --img_s 256 \
                                 --train_datasets_size 2000')

for i in range(1,5):
    os.system('python test_FFD1.py --denoise_model 3_wfto2p_dn4_fm32_rr2_202204250945 \
                                --datasets_folder wd_virus_binning3_'+str(i)+' \
                                --pth_path pth --GPU 2 --batch_size 1 --fmap 32 \
                                --img_h 256 --img_s 128 \
                                --test_datasize 2000 --realign_rate 2')

os.system('python train_GAN_FFDMul.py --datasets_folder vol_1000_400_200_NA_0.30_res_0.40_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 2 \
                                 --img_h 128 --img_w 128 --img_s 256 \
                                 --train_datasets_size 2000')

os.system('python test_FFD1.py --denoise_model vol_3000_800_200_NA_0.30_res_1.60_Hz_10_ndif_0.03_exp_3_dens_10k_arr_dn4_fm32_rr2_202205031605 \
                                --datasets_folder vol_3000_1200_200_NA_0.30_res_1.60_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                --pth_path pth --GPU 2 --batch_size 1 --fmap 32 \
                                --img_h 256 --img_s 128 \
                                --test_datasize 2000 --realign_rate 2')

os.system('python train_GAN_FFDMul.py --datasets_folder vol_1000_400_200_NA_0.30_res_0.40_Hz_10_ndif_0.03_exp_3_dens_10k_arr \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 2 \
                                 --img_h 128 --img_w 128 --img_s 256 \
                                 --train_datasets_size 2000')

os.system('python test_FFD1.py --denoise_model vol_1000_400_200_NA_0.30_res_0.40_nopretrain \
                                --datasets_folder vol_1000_400_200_NA_0.30_res_0.40_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
                                --pth_path pth --GPU 2 --batch_size 1 --fmap 32 \
                                --img_h 256 --img_s 128 \
                                --test_datasize 2000 --realign_rate 2')

os.system('python train_GAN_FFDMul.py --datasets_folder NA_0.60_res_0.60_ndif_0.02_exp_3 \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 100 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 2 \
                                 --img_h 128 --img_w 128 --img_s 256 \
                                 --train_datasets_size 2000')

os.system('python test_FFD1.py --denoise_model NA_0.60_res_0.60_ndif_0.02_exp_3_dn4_fm32_rr2_202205091043 \
                                --datasets_folder NA_0.60_res_0.60_ndif_0.02_exp_3_test \
                                --pth_path pth --GPU 2 --batch_size 1 --fmap 32 \
                                --img_h 256 --img_s 128 \
                                --test_datasize 2000 --realign_rate 2')


os.system('python train_GAN_FFDMul_conti.py --datasets_folder vessel_vol_600_200_NA_0.30_Hz_10_ndif_0.03_exp_5_dens_10k_arr \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --denoise_model vol_600_200_rr2_202205062254 \
                                 --pth_index 76 \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 2 \
                                 --img_h 128 --img_w 128 --img_s 256 \
                                 --train_datasets_size 2000')

os.system('python test_FFD1.py --denoise_model vessel_vol_600_200_NA_0.30_Hz_10_ndif_0.03_exp_5_dens_10k_arr_dn4_fm32_rr2_202205131642 \
                                --datasets_folder vessel_vol_600_200_NA_0.30_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                --pth_path pth --GPU 2 --batch_size 1 --fmap 32 \
                                --img_h 256 --img_s 128 \
                                --test_datasize 2000 --realign_rate 2')

os.system('python train_GAN_FFDMul_big.py --datasets_folder vessel_vol_600_200_NA_0.30_Hz_10_ndif_0.03_exp_5_dens_10k_arr \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --n_epochs 30 --GPU 2 --batch_size 1 \
                                 --fmap 32 --lr 0.0001 --realign_rate 1 \
                                 --img_h 256 --img_w 256 --img_s 128 \
                                 --train_datasets_size 2000')

os.system('python train_GAN_FFDMul_conti.py --datasets_folder NA0.6_train \
                                 --datasets_folder_eval vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_5_dens_10k_arr_test \
                                 --denoise_model vol_600_200_rr2_202205062254 \
                                 --pth_index 76 \
                                 --n_epochs 30 --GPU 1 --batch_size 1 \
                                 --fmap 32 --lr 0.0002 --realign_rate 2 \
                                 --img_h 128 --img_w 128 --img_s 256 \
                                 --train_datasets_size 2000')
'''
os.system('python test_FFD1.py --denoise_model NA_0.60_res_0.60_ndif_0.02_exp_3_dn4_fm32_rr2_nopretrain \
                                --datasets_folder NA0.6_test \
                                --pth_path pth --GPU 1 --batch_size 1 --fmap 32 \
                                --img_h 256 --img_s 128 \
                                --test_datasize 2000 --realign_rate 2')