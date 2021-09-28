# Deep widefield calcium inference (DWCI)

Implementation for deep widefield calcium inference (DWCI)

## Current contents
### A modified NAOMi for widefield imaging, with various adjustable nobs for fitting your system, including
1. Optical parameters of microscope: NA, FOV, FN, immersion medium
2. Indicator parameters: expression level, indicator types (e.g. GCaMP6 or GCaMP7)
3. Imaging parameters: session length, frame rate, illumination powers, imaging d


## Usage
### Single batch generation
  1. Confirm your imaging parameters. Some pre-configed parameters are available. 
  2. run NAOMi_1p_single.m with your parameters.

### Multiple batches generation
  1. Confirm your imaging parameters. Some pre-configed parameters are available. 
  2. run NAOMi_1p_loop.m with your parameters. The program will automatically store multiple copies with different parameters there.


## Notifications
  1. The system NA should not exists 0.5.
  2. Recommanded imaging power is 1mW/mm^2. 
