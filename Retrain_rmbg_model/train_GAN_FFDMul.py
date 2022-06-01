import os
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.utils.data import DataLoader
import argparse
import time
import datetime
from network import Network_3D_Unet
from data_process import train_preprocess_lessMemory_segBG, shuffle_datasets_lessMemory, trainsetMul,  test_preprocess_lessMemoryNoTail_chooseOne, testset, multibatch_test_save, singlebatch_test_save, train_preprocess_lessMemoryMulStacks_eval
from utils import save_yaml

from torch.utils.tensorboard import SummaryWriter
from Discriminator import NLayerDiscriminator3D
import functools
import numpy as np
from utils import FFDrealign4, inv_FFDrealign4, FFDrealign, inv_FFDrealign

# update: 9/10/2021. 
#   +evaluation
#   +updated tensorboard


from skimage import io
#############################################################################################################################################
if __name__ == '__main__':
    # tensorboard session
    writer = SummaryWriter()

    parser = argparse.ArgumentParser()
    parser.add_argument("--n_epochs", type=int, default=40, help="number of training epochs")
    parser.add_argument('--GPU', type=str, default='0,1', help="the index of GPU you will use for computation")

    parser.add_argument('--batch_size', type=int, default=2, help="batch size")
    parser.add_argument('--img_w', type=int, default=150, help="the width of image sequence")
    parser.add_argument('--img_h', type=int, default=150, help="the height of image sequence")
    parser.add_argument('--img_s', type=int, default=150, help="the slices of image sequence")
    parser.add_argument('--realign_rate', type=int, default=150, help="the slices of image sequence")

    parser.add_argument('--lr', type=float, default=0.00005, help='initial learning rate')
    parser.add_argument("--b1", type=float, default=0.5, help="Adam: bata1")
    parser.add_argument("--b2", type=float, default=0.999, help="Adam: bata2")
    parser.add_argument('--normalize_factor', type=int, default=1, help='normalize factor')
    parser.add_argument('--fmap', type=int, default=16, help='number of feature maps')
    parser.add_argument('--down_num', type=int, default=4, help='number of feature maps')

    parser.add_argument('--output_dir', type=str, default='./results', help="output directory")
    parser.add_argument('--datasets_folder', type=str, default='train', help="A folder containing files for training")
    parser.add_argument('--GT_folder', type=str, default='mov_wo_bg', help="video or image")
    parser.add_argument('--input_folder', type=str, default='mov_w_bg', help="")
    
    parser.add_argument('--datasets_path', type=str, default='datasets', help="dataset root path")
    parser.add_argument('--pth_path', type=str, default='pth', help="pth file root path")
    parser.add_argument('--select_img_num', type=int, default=100000, help='select the number of images used for training')
    parser.add_argument('--train_datasets_size', type=int, default=4000, help='datasets size for training')

    # evaluate
    parser.add_argument('--datasets_folder_eval', type=str, default='train', help="A folder containing files for training")

    opt = parser.parse_args()
    ########################################################################################################################
    opt.img_w = opt.img_w*opt.realign_rate
    opt.img_h = opt.img_h*opt.realign_rate
    # default image gap is 0.75*image_dim
    opt.gap_s=int(opt.img_s*0.5)
    opt.gap_w=int(opt.img_w*0.5)
    opt.gap_h=int(opt.img_h*0.5)
    opt.ngpu=str(opt.GPU).count(',')+1
    print('\033[1;31mTraining parameters -----> \033[0m')
    print(opt)

    ########################################################################################################################
    if not os.path.exists(opt.output_dir): 
        os.mkdir(opt.output_dir)
    current_time = opt.datasets_folder+'_dn'+str(opt.down_num)+'_fm'+str(opt.fmap)+'_rr'+str(opt.realign_rate)+'_'+datetime.datetime.now().strftime("%Y%m%d%H%M")
    output_path = opt.output_dir + '/' + current_time
    pth_path = 'pth//'+ current_time
    if not os.path.exists(pth_path): 
        os.mkdir(pth_path)
    f = open(pth_path+'//loss.txt','w') 
    ########################################################################################################################
    # evaluation datasets, this is optional
    if opt.datasets_folder_eval is not None:
        eval_flag = False # we will do the evaluations
        name_list_eval, noise_img_eval, noise_img2_eval, coordinate_list_eval = train_preprocess_lessMemoryMulStacks_eval(opt)


    ########################################################################################################################
    yaml_name = pth_path+'//para.yaml'
    save_yaml(opt, yaml_name)

    os.environ["CUDA_VISIBLE_DEVICES"] = str(opt.GPU)
    batch_size = opt.batch_size
    lr = opt.lr

    # get training dataset
    coor_list, GT_list, input_list = train_preprocess_lessMemory_segBG(opt)
    print('coor_list -----> ',len(coor_list))
    # print('input_list -----> ',input_list)

    ########################################################################################################################
    # build the network
    L1_pixelwise = torch.nn.L1Loss() # build l1 loss
    L2_pixelwise = torch.nn.MSELoss() # build l2 loss
    Dloss = torch.nn.BCEWithLogitsLoss()

    denoise_generator = Network_3D_Unet(in_channels = opt.realign_rate*opt.realign_rate,
                                        out_channels = opt.realign_rate*opt.realign_rate,
                                        f_maps=opt.fmap,
                                        final_sigmoid = True,
                                        down_num = opt.down_num)
    norm_layer = functools.partial(nn.BatchNorm3d, affine=True, track_running_stats=True)
    Dnet = NLayerDiscriminator3D(input_nc=opt.realign_rate*opt.realign_rate, ndf=16, n_layers=3, norm_layer=norm_layer)

    if torch.cuda.is_available():
        denoise_generator = denoise_generator.cuda()
        Dnet.cuda()
        denoise_generator = nn.DataParallel(denoise_generator, device_ids=range(opt.ngpu))
        print('\033[1;31mUsing {} GPU for training -----> \033[0m'.format(torch.cuda.device_count()))
        L2_pixelwise.cuda()
        L1_pixelwise.cuda()
        Dloss.cuda()
    ########################################################################################################################
    # optimizer
    optimizer_G = torch.optim.Adam(denoise_generator.parameters(),
                                    lr=opt.lr, betas=(opt.b1, opt.b2))
    optimizer_D = torch.optim.Adam(Dnet.parameters(),
                                    lr=opt.lr, betas=(opt.b1, opt.b2))

    ########################################################################################################################
    cuda = True if torch.cuda.is_available() else False
    Tensor = torch.cuda.FloatTensor if cuda else torch.FloatTensor
    prev_time = time.time()

    ########################################################################################################################
    # main training loop

    time_start=time.time()

    # start training
    for epoch in range(0, opt.n_epochs):
        coor_list = shuffle_datasets_lessMemory(coor_list)
        for index in range(len(coor_list)):
            iteration = index
            per_coor = coor_list[index]
            train_im_name = per_coor['name']
            init_w = per_coor['init_w']
            init_h = per_coor['init_h']
            init_s = per_coor['init_s']
            GT_im = GT_list[train_im_name]
            input_im = input_list[train_im_name]
            
            GT_patch =       GT_im[init_s:init_s+opt.img_s, init_w:init_w+opt.img_w, init_h:init_h+opt.img_h]
            input_patch = input_im[init_s:init_s+opt.img_s, init_w:init_w+opt.img_w, init_h:init_h+opt.img_h].copy()

            rand_bg = np.random.randint(0, 5000)
            rand_gama = np.random.randint(1000, 5000)/1000
            input_patch = (input_patch+rand_bg)/rand_gama

            input = torch.from_numpy(np.expand_dims(np.expand_dims(input_patch, 3),0)).cuda()
            input = input.permute([0,4,1,2,3])
            target = torch.from_numpy(np.expand_dims(np.expand_dims(GT_patch, 3),0)).cuda()
            target = target.permute([0,4,1,2,3])
            # print('input ',input.shape)
            # input=torch.from_numpy(np.expand_dims(input_patch, 0))
            # target=torch.from_numpy(np.expand_dims(GT_patch, 0))
            '''
            train_data = trainsetMul(coor_list, GT_list, input_list,opt)
            trainloader = DataLoader(trainsetMul, batch_size=batch_size, shuffle=True, num_workers=3)
    
            for iteration, (input, target) in enumerate(trainloader):
            '''
            input=input.cuda().float()
            target = target.cuda()

            re_type='none'

            real_A = FFDrealign(input, opt.realign_rate)
            real_B = FFDrealign(target, opt.realign_rate)

            real_A = Variable(real_A)

            fake_B = denoise_generator(real_A)

            # print('real_A shape -----> ', real_A.shape)
            # print('real_B shape -----> ',real_B.shape) , ef, df
            # print('fake_B shape -----> ',fake_B.shape)

            L1_loss = L1_pixelwise(fake_B, real_B)
            L2_loss = L2_pixelwise(fake_B, real_B)

            # print('fake_B -----> ',fake_B.shape)
            D_fake_B = Dnet(fake_B)
            valid = Variable(Tensor(D_fake_B.shape).fill_(1.0), requires_grad=False)
            fake = Variable(Tensor(D_fake_B.shape).fill_(0.0), requires_grad=False)
            # print('fake -----> ',fake.shape)
            Dloss1 = Dloss(D_fake_B, fake)
            ################################################################################################################
            optimizer_G.zero_grad() # initialize the gradient
            # Total loss
            # Total_loss =  0.5*L1_loss + 0.5*L2_loss + Dloss1 # 1:1 L1 and L2 loss
            Total_loss = L1_loss + L2_loss + Dloss1 # 1:1 L1 and L2 loss
            Total_loss.backward() # back propagation
            optimizer_G.step() # one step
            ################################################################################################################
            D_fake_B = Dnet(fake_B.detach())
            D_real_A = Dnet(real_A)

            optimizer_D.zero_grad()
            real_loss = Dloss(D_real_A, valid)
            fake_loss = Dloss(D_fake_B, fake)
            d_loss = (real_loss + fake_loss) / 2

            d_loss.backward()
            optimizer_D.step()
            ################################################################################################################
            batches_done = epoch * len(coor_list) + iteration
            batches_left = opt.n_epochs * len(coor_list) - batches_done
            time_left = datetime.timedelta(seconds=int(batches_left * (time.time() - prev_time)))
            prev_time = time.time()

            if iteration%10 == 0:
                time_end=time.time()
                print('\r[Epoch %d/%d] [Batch %d/%d] [Total loss: %.2f, L1 Loss: %.2f, L2 Loss: %.2f] [ETA: %s] [Time cost: %.2d s]        ' 
                % (
                    epoch+1,
                    opt.n_epochs,
                    iteration+1,
                    len(coor_list),
                    Total_loss.item(),
                    L1_loss.item(),
                    L2_loss.item(),
                    time_left,
                    time_end-time_start
                ), end=' ')
                f.write('Dloss1 '+str(Dloss1.item())+' L1_loss '+str(L1_loss.item())+' L2_loss '+str(L1_loss.item())+' \n')
                print('Dloss1 '+str(Dloss1.item())+' L1_loss '+str(L1_loss.item())+' L2_loss '+str(L1_loss.item())+' \n')

                print(real_B.cpu().detach().numpy().max(),' ---> ', real_B.cpu().detach().numpy().min(),' ---> ',\
                fake_B.cpu().detach().numpy().max(),' ---> ',  fake_B.cpu().detach().numpy().min())

            if (iteration+1)%len(coor_list) == 0:
                print('\n', end=' ')

            ################################################################################################################
            # save model
            if (iteration + 1) % (len(coor_list)) == 0:
                model_save_name = pth_path + '//E_' + str(epoch+1).zfill(2) + '_Iter_' + str(iteration+1).zfill(4) + '.pth'
                if isinstance(denoise_generator, nn.DataParallel): 
                    torch.save(denoise_generator.module.state_dict(), model_save_name)  # parallel
                else:
                    torch.save(denoise_generator.state_dict(), model_save_name)         # not parallel

        # tensorboard
        writer.add_scalar('Loss/total', Total_loss.item(),epoch)
        writer.add_scalar('Loss/L1', L1_loss.item(), epoch)
        writer.add_scalar('Loss/L2', L2_loss.item(), epoch)


        if (epoch+1)%1 == 0:
            # output image
            # print('fake_B ---> ',fake_B.shape)
            fake_B_realign = inv_FFDrealign(fake_B, opt.realign_rate)
            real_B_realign = inv_FFDrealign(real_B, opt.realign_rate)
            real_A_realign = inv_FFDrealign(real_A, opt.realign_rate)

            output_img = fake_B_realign.cpu().detach().numpy()
            train_GT = real_B_realign.cpu().detach().numpy()
            train_input = real_A_realign.cpu().detach().numpy()

            train_input = train_input.squeeze().astype(np.float32)*opt.normalize_factor
            train_GT = train_GT.squeeze().astype(np.float32)*opt.normalize_factor
            output_img = output_img.squeeze().astype(np.float32)*opt.normalize_factor
            train_input = np.clip(train_input, 0, 65535)
            train_GT = np.clip(train_GT, 0, 65535).astype('uint16')
            output_img = np.clip(output_img, 0, 65535)
            result_name = pth_path+ '/' + str(epoch)  + '_output.tif'
            noise_img1_name = pth_path + '/' + str(epoch) + '_noise1.tif'
            noise_img2_name = pth_path + '/' + str(epoch)  + '_noise2.tif'
            io.imsave(result_name, output_img)
            io.imsave(noise_img1_name, train_input)
            io.imsave(noise_img2_name, train_GT)

    # close writer
    f.close()
    writer.close()

