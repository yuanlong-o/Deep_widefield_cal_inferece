import torch
########################################################################################################################
def create_feature_maps(init_channel_number, number_of_fmaps):
    return [init_channel_number * 2 ** k for k in range(number_of_fmaps)]

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


def FFDrealign8(input):
    # batch channel time height width
    realign_input = torch.cuda.FloatTensor(input.shape[0], input.shape[1]*8, int(input.shape[2]/2), int(input.shape[3]/2), int(input.shape[4]/2))
    # print('realign_input -----> ',realign_input.shape)
    # print('input -----> ',input.shape)
    realign_input[:,0,:,:,:] = input[:,:,0::2,0::2, 0::2]
    realign_input[:,1,:,:,:] = input[:,:,0::2,0::2, 1::2]
    realign_input[:,2,:,:,:] = input[:,:,0::2,1::2, 0::2]
    realign_input[:,3,:,:,:] = input[:,:,0::2,1::2, 1::2]

    realign_input[:,4,:,:,:] = input[:,:,1::2,0::2, 0::2]
    realign_input[:,5,:,:,:] = input[:,:,1::2,0::2, 1::2]
    realign_input[:,6,:,:,:] = input[:,:,1::2,1::2, 0::2]
    realign_input[:,7,:,:,:] = input[:,:,1::2,1::2, 1::2]
    return realign_input

def inv_FFDrealign8(input):
    # batch channel time height width
    realign_input = torch.cuda.FloatTensor(input.shape[0], int(input.shape[1]/8), int(input.shape[2]*2), int(input.shape[3]*2), int(input.shape[4]*2))

    realign_input[:,:,0::2,0::2, 0::2] = input[:,0,:,:,:]
    realign_input[:,:,0::2,0::2, 1::2] = input[:,1,:,:,:]
    realign_input[:,:,0::2,1::2, 0::2] = input[:,2,:,:,:]
    realign_input[:,:,0::2,1::2, 1::2] = input[:,3,:,:,:]

    realign_input[:,:,1::2,0::2, 0::2] = input[:,4,:,:,:]
    realign_input[:,:,1::2,0::2, 1::2] = input[:,5,:,:,:]
    realign_input[:,:,1::2,1::2, 0::2] = input[:,6,:,:,:]
    realign_input[:,:,1::2,1::2, 1::2] = input[:,7,:,:,:]
    return realign_input