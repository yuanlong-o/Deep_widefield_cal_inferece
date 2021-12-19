import os
import sys

# train SEG model
os.system('python train_SEG_model.py --img_h 256 --img_w 256 --GPU 0 --f_maps 16 \
--normalize_factor 1000 --input_nc 1 --output_nc 1 --n_epochs 100 --lr 0.0001 --frame_num 64 \
--train_datasets_size 500 --datasets_folder test1_mov_wo_bg_100_64_498_10_2_6000')