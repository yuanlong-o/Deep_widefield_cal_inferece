from DWonder.SEG.model_3DUnet import UNet3D, T_conv_net, TS_UNet3D
import torch.nn as nn
import torch.nn.functional as F
import torch
from DWonder.SEG.model_Unet import Unet, Unet3, Unet4



class SEG_Network_3D_Unet(nn.Module):
    def __init__(self, UNet_type = '3DUNet', in_channels=1, out_channels=1, frame_num=64, final_sigmoid = True,f_maps=64):
        super(SEG_Network_3D_Unet, self).__init__()

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.final_sigmoid = final_sigmoid

        if UNet_type == '3DUNet':
            self.Generator = UNet3D( in_channels = in_channels,
                                     out_channels = out_channels, 
                                     frame_num = frame_num,
                                     final_sigmoid = final_sigmoid,
                                     f_maps = f_maps)
        if UNet_type == 'TS_UNet3D':
            self.Generator = TS_UNet3D( in_channels = in_channels,
                                     out_channels = out_channels, 
                                     frame_num = frame_num,
                                     final_sigmoid = final_sigmoid,
                                     f_maps = f_maps)
    def forward(self, x):
        fake_x = self.Generator(x)
        return fake_x


class Network_TSnet(nn.Module):
    def __init__(self,  in_channels=1, out_channels=1, frame_num=64, f_maps=64):
        super(Network_TSnet, self).__init__()

        self.in_channels = in_channels
        self.out_channels = out_channels

        self.T_conv_net = T_conv_net( in_channels = in_channels, 
                                    frame_num = frame_num,
                                    tc_f_maps = f_maps)
        self.Unet = Unet3( in_ch = f_maps, 
                                out_ch = out_channels,
                                    f_num = f_maps)

        # T_conv_net = T_conv_net(in_channels = 1, frame_num = opt.frame_num, tc_f_maps=opt.tc_f_num)
        # Unet_1 = Unet3(in_ch = opt.tc_f_num, out_ch = 1, f_num = opt.unet_f_num)

    def forward(self, x):
        # print('Initialize TS_net --->')
        pred_patch = self.T_conv_net(x)
        input_u = torch.squeeze(pred_patch, dim=2)
        output_u = self.Unet(input_u)
        # fake_x = self.Generator(x)
        return output_u
