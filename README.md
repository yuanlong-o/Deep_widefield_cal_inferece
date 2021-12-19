# Deep widefield neuron finder (DeepWonder)
<img src="img/deepwonder_logo.png" width="800" align="center">
Implementation for deep widefield neuron finder (DeepWonder)

![Imaging modality](https://img.shields.io/badge/Imaging%20modality-Wide--filed%20SIngle--photon-brightgreen)![Purpose](https://img.shields.io/badge/Purpose-Neuron%20analysis-orange)![Datasets](https://img.shields.io/badge/Datasets-Synthetic%20biological%20images-blue)

## 📋 Table of content
 1. [Overview](#Overview)
 2. [Quick start DeepWonder](#Start)
    1. [Environment](#Environment)
    2. [Install dependencies](#Dependencies)
    3. [Download the demo code and data](#Download)
    4. [Run the trained model](#Run)
    5. [Work for your own data](#Own data)
    6. [Run Demo data with Google Colab](#Colab)
 3. [Train DeepWonder](#Train)
    1. [Realistic widefield capture generation](#NAOMI)
    2. [Background removal network training](#Train BR)
    3. [Neuron segmentation network training](#Train NS)
 4. [Other information](#Information)
    1. [Citation](#Citation)
    2. [Email](#Email)

## **📚** Overview <a name="Overview"></a>
Widefield imaging provides optical access to multi-millimeter fields of view in video rate, but cell-resolved calcium imaging has mostly contaminated by tissue scattering and calcium inference suffers multiple times slower processing speed than two-photon microscope in the same scale. We present a deep widefield neuron finder (DeepWonder), which is fueled by numerical recordings but effectively works on experimental data with an order of magnitude of speed improvement and improved precision and accuracy. With DeepWonder, widefield calcium recordings can be quickly converted into background-free movies, and neuronal segments and activites then be demixed out. For more details and results, please see our companion paper at TODO.

<img src="img/Workflow.png" width="600" align="center">

## **⏳** Quick start DeepWonder <a name="Start"></a>
This tutorial will show how DeepWonder removes backgrounds in  widefield captures and further segments neurons.
### **💡** Environment <a name="Environment"></a>
* Ubuntu 16.04 
* Python 3.6
* Pytorch >= 1.3.1
* NVIDIA GPU (24 GB Memory) + CUDA

### **💡** Install dependencies <a name="Dependencies"></a>
* Create a virtual environment and install some dependencies.
```
$ conda create -n deepwonder_env python=3.6
$ source activate deepwonder_env
$ pip install -q torch==1.10.0
$ pip install -q torchvision==0.8.2
$ pip install deepwonder==0.3
$ pip install -q opencv-python==4.1.2.30
$ pip install -q tifffile  
$ pip install -q scikit-image==0.17.2
$ pip install -q scikit-learn==0.24.1
$ pip install -q scipy==1.5.2
$ pip install -q sklearn==0.0
```
### **💡** Download the demo code and data <a name="Download"></a>
```
$ git clone git://github.com/cabooster/Deep_widefield_cal_inferece
$ cd DeepCAD/DeepWonder/
```
We upload three demo data on Google drive: [synthetic wide field data by naomi1p code ](https://drive.google.com/drive/folders/1WiTrL5gRuMUssMYt2uDRDO-5pmmrdNSc?usp=sharing), [Rush data](https://drive.google.com/drive/folders/1CP6CuAmOkAx_hoAhT4h-Pd1o_FTcva9M?usp=sharing) and [ordinary wide field data](https://drive.google.com/drive/folders/1QSqbNWmZTlbctYt0Vh0I529gt-kYNX4w?usp=sharing). You can download them and put them in the *DeepWonder/datasets* folder.
### **💡** Run the trained model <a name="Run"></a>
Run the script.py to analyze the demo data. 
```
$ python script.py 
```
You can get the result in the *DeepWonder/results* folder. The removing background result is in the *DeepWonder/results/XXXX/RMBG* folder. The neuron footprints id in the *DeepWonder/results/XXXX/f_con* folder as a *.tif* file and the *DeepWonder/results/XXXX/mat folder* as a *.mat* file.
### **💡** Work for your own data <a name="Own data"></a>
Because we trained our network on the synthetic data whose pixel size is *0.8 um*, we should first resize your data to *0.8 um* one pixel. Then put your data in the *DeepWonder/datasets* folder and use DeepWonder as described above.
### **💡** Run Demo data with Google Colab <a name="Colab"></a>
We have stored our data in Google Colab, which is a free Jupyter notebook environment that requires no setup and runs entirely in the cloud. A demo script with full processing of DeepWonder on several demo datasets (including NAOMi1p virtual datasets, cropped RUSH datasets, and two-photon validation datasets, all mounted to the Google Colab using Google Drive) is available through Colab via 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1cluMDiY0G0NR4j62OkhOd8Cs316EbuDc?usp=sharing).
You can process your data with the on-line DeepWonder guide by the tutorial.
## **🔁** Train DeepWonder <a name="Train"></a>
### **💡** Realistic widefield capture generation <a name="NAOMI"></a>
DeepWonder relies on a highly realistic simulation of widefield capture for training a network that removes widefield background. We refered to the Neural Anatomy and Optical Microscopy (NAOMi) package (https://www.sciencedirect.com/science/article/pii/S0165027021001084) to populat the brain tissue with multiple blood vessels, somata, axons, and dendrites. We further modified the NAOMi pipeline such that it could faithfully simulate data acquisition of one-photon excitations with one-photon excitation model and related noise model, which we termed as NAOMi1p. The virtual widefield data by NAOMi1p can be affected by optical parameters, illumination powers, and other parameters like protein concentration. The full key parameters are listed as
1. Optical parameters of microscope: NA, FOV, FN, immersion medium
2. Indicator parameters: expression level, indicator types (e.g. GCaMP6 or GCaMP7)
3. Imaging parameters: session length, frame rate, illumination powers, imaging d
  As an example, run NAOMi1p_widefield.m to generate a 750 x 750 x 1000 frame virtual widefield recording. You can adjust the above parameters based on your own system.
### **💡** Background removal network training <a name="Train BR"></a>
Put the GT images generated by NAOMI1p code in the *DeepWonder/datasets/XXXX/GT* and the raw images generated by NAOMI1p code in the *DeepWonder/datasets/XXXX/Input*. To achieve that, run the script_train_RMBG.py to train the background removal network.
```
$ source activate deepwonder_env
$ python script_train_RMBG.py train
```
You will get trained removing background model in *DeepWonder/RMBG_pth* folder.
### **💡** Neuron segmentation network training <a name="Train NS"></a>
You can establish a training datasets for segmentation: input data is the no-background data generated by NAOMI1p code and GT data is generated by biological modeling. Put the GT data in the *DeepWonder/datasets/XXXX/mask* and the raw data in the *DeepWonder/datasets/XXXX/image*. Run the script_train_SEG.py to train the neuron segmentation network.

```
$ source activate deepwonder_env
$ python script_train_SEG.py train
```
You will get trained Neuron segmentation model in *DeepWonder/SEG_pth* folder.
## 🤝 Other information <a name="Information"></a>
### **📝** Citation <a name="Citation"></a>



### **📝** Email <a name="Email"></a>
If you have any questions, you can send email to me zhanggx19@mails.tsinghua.edu.cn.
