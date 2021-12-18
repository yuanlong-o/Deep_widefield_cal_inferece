import os
import sys

os.system('python main_RSM.py \
--RMBG_model_folder arr_202110011541 \
--RMBG_model_name E_30_Iter_4009 \
--SEG_model_folder TS3DUnetFFD_datasets_test1_mov_wo_bg_100_64_498_10_2_6000_ic1_oc1_lr0.0001_fm16_20211129-1355 \
--SEG_model_name seg_30 \
--test_datasize 20000 \
--datasets_path datasets \
--datasets_folder vol_600_200_NA_0.30_res_0.80_Hz_10_ndif_0.03_exp_3_dens_10k_arr_test \
--GPU 0')