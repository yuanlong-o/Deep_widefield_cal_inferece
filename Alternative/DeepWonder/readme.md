
### Calcium Extraction with DeepWonder

DeepWonder offers a streamlined way to extract calcium signatures from recordings. Follow the instructions below to set up your environment and process your data.

#### Setting up the Environment

1. **Create a Virtual Environment**:
   - Using Anaconda, create a new Python 3.9 virtual environment named `DWonder`.
   - Activate the `DWonder` environment by executing:
     ```
     conda activate DWonder
     ```

2. **Install GPU-Version PyTorch (Stable Build)**:
   - Install the GPU-version of PyTorch by referring to the official [PyTorch installation guide](https://pytorch.org/get-started/locally/).
   - Execute the following commands to install PyTorch and its dependencies:
     ```
     pip3 install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu116
     ```
   
3. **Verify PyTorch Installation**:
   - Ensure that PyTorch has been installed successfully by testing its GPU capabilities.
     ```python
     python
     >>> import torch
     >>> torch.cuda.is_available()  # This should return True if successful.
     ```

#### Processing the Data

1. **Prepare Your Data**:
   - Place the video recordings you intend to process in the directory `[Your_experiment_path]`.

2. **Execute DeepWonder Scripts**:
   - Navigate to the `DeepWonder` directory.
   - Run the following commands:
     ```
     python ca_rmbg.py --path [Your_experiment_path]
     python ca_seg.py  --path [Your_experiment_path]
     ```

3. **One-Click Script for Windows Users**:
   - We provide a convenient script for Windows users to execute DeepWonder.
   - Update the `$Paths` variable in `ca_extract.ps1` to your experiment path:
     ```
     $Paths = 'Your_Experiment_Path'
     ```
   - Execute the script using Windows PowerShell.