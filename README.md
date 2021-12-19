# Deep widefield neuron finder (DeepWonder)
<img src="img/deepwonder_logo.png" width="400" align="center">
Implementation for deep widefield neuron finder (DeepWonder)
![PyPI version](https://badge.fury.io/py/segmentation-models-pytorch.svg)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1cluMDiY0G0NR4j62OkhOd8Cs316EbuDc?usp=sharing)
## Overview
Widefield imaging provides optical access to multi-millimeter fields of view in video rate, but cell-resolved calcium imaging has mostly contaminated by tissue scattering and calcium inference suffers multiple times slower processing speed than two-photon microscope in the same scale. We present a deep widefield neuron finder (DeepWonder), which is fueled by numerical recordings but effectively works on experimental data with an order of magnitude of speed improvement and improved precision and accuracy. With DeepWonder, widefield calcium recordings can be quickly converted into background-free movies, and neuronal segments and activites then be demixed out. For more details and results, please see our companion paper at TODO.

<img src="img/Workflow.png" width="600" align="center">

## Run Demo data with Google Colab
We have stored our data in Google Colab, which is a free Jupyter notebook environment that requires no setup and runs entirely in the cloud. A demo script with full processing of DeepWonder on several demo datasets (including NAOMi1p virtual datasets, cropped RUSH datasets, and two-photon validation datasets, all mounted to the Google Colab using Google Drive) is available through Colab via https://colab.research.google.com/drive/1cluMDiY0G0NR4j62OkhOd8Cs316EbuDc?usp=sharing#scrollTo=hE0tSwj8LTg0. 


## Run demo data locally without re-training
This tutorial will show how DeepWonder removes backgrounds in  widefield captures and further segments neurons.
### Environment 
* Ubuntu 16.04 
* Python 3.6
* Pytorch >= 1.3.1
* NVIDIA GPU (24 GB Memory) + CUDA

### 🛠 Environment configuration
* Create a virtual environment and install Pytorch. In the 4th step, please select the correct Pytorch version that matches your CUDA version from https://pytorch.org/get-started/previous-versions/
```
$ conda create -n deepwonder python=3.6
$ source activate deepwonder
$ pip install torch==1.3.1
$ conda install pytorch torchvision cudatoolkit -c pytorch 
```

* Install other dependencies

```
$ conda install -c anaconda matplotlib opencv scikit-learn scikit-image
$ conda install -c conda-forge h5py pyyaml tensorboardx tifffile
```
### Download the source code and data
* The source code can be downloaded from 
* 
```
$ git clone git://github.com/yuanlong-o/Deep_widefield_cal_inferece
$ cd DeepWonder/DeepWonder_pytorch/
```

* The demo data can be downloaded by
* 
### Remove backgrounds
```
TODO
```

### Neuronal segmentations
```
TODO
```


## Work for your own data
### Realistic widefield capture generation
DeepWonder relies on a highly realistic simulation of widefield capture for training a network that removes widefield background. We refered to the Neural Anatomy and Optical Microscopy (NAOMi) package (https://www.sciencedirect.com/science/article/pii/S0165027021001084) to populat the brain tissue with multiple blood vessels, somata, axons, and dendrites. We further modified the NAOMi pipeline such that it could faithfully simulate data acquisition of one-photon excitations with one-photon excitation model and related noise model, which we termed as NAOMi1p. The virtual widefield data by NAOMi1p can be affected by optical parameters, illumination powers, and other parameters like protein concentration. The full key parameters are listed as
1. Optical parameters of microscope: NA, FOV, FN, immersion medium
2. Indicator parameters: expression level, indicator types (e.g. GCaMP6 or GCaMP7)
3. Imaging parameters: session length, frame rate, illumination powers, imaging depth 
As an example, run NAOMi_1p_single.m to generate a 750 x 750 x 1000 frame virtual widefield recording. You can adjust the above parameters based on your own system. To generate multiple training pairs with different tissue distributions, run NAOMi_1p_loop.m

* Software requirement for NAOMi1p:
  * Matlab-compatible version of Windows or Linux (see https://www.mathworks.com/support/requirements/matlab-system-requirements.html)
  * Matlab R2020b or newer
  * Toolboxes: signal processing, image processing, statistics and machine learning, curve fitting, parallel computing
  * Matlab-compatible CUDA driver (>= 10.1 for Matlab R2020a)

### Background removal network training
We train a 3D Unet to remove the background of the network. To achieve that, run the script.py to train a network that removes backgrounds of your system.
```
$ source activate deepcad
$ python script.py train
```
This will generate xxxx.

### Neuronal segmentations
We further use a convolutional neural network (CNN) for neruonal segmentations, and local non-negative matrix factorization (NMF) to readout neuronal signals. Take background-removed files in BMBG, run xxxx to achieves the final results.
(how it trains?)



