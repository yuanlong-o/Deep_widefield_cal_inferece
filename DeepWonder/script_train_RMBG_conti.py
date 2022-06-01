import os
import sys

# train RMBG model
os.system('python train_RMBG_model_conti.py --datasets_folder train_RMBG1 \
                               --n_epochs 2 --GPU 0 --batch_size 1 \
                               --fmap 32 --lr 0.0001 \
                               --denoise_model vol_600_200_rr2 --pth_name E_76_Iter_2004 \
                               --img_h 128 --img_w 128 --img_s 128 \
                               --train_datasets_size 4000')
