function mov = applyNoiseModel(clean_mov,noise_params,dynodeResponse)

% function mov = applyNoiseModel(clean_mov,noise_params)
%
% Function to apply the electronics noise model to a two-photon microscopy
% simulated video created with zero noise. Inputs are:
%   - clean_mov    - MxNxT array where each 2-D array is a movie frame with
%                    no noise
%   - noise_params - Struct contaning the parameters for the noise model
%          .mu     - Mean measurement increase per photon (default = 150)
%          .mu0    - Electronics offset (default = 0)
%          .sigma  - Variance increase per photon (default = 5000)
%          .sigma0 - Electronics base noise variance (default = 3.5)
% 
% Output is
%   - mov - Movie with simulated noise 
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run noise model on each frame

if(nargin<3)
  dynodeResponse = [];
end

mov = zeros(size(clean_mov),'single');                                              % Set up movie array

if(~isfield(noise_params,'type')||isempty(noise_params.type)||strcmp(noise_params.type,'poissongauss'))
  for kk = 1:size(clean_mov,3)                                               % Iterate over time steps
    mov(:,:,kk) = PoissonGaussNoiseModel(clean_mov(:,:,kk), noise_params); % Run data through the noise model    
    mov(:,:,kk) = pixel_bleed(mov(:,:,kk), noise_params.bleedp, noise_params.bleedw); % Add bleed-through from previous pixel
  end
elseif(strcmp(noise_params.type,'dynode'))
  if(isempty(dynodeResponse))
    ndynodes = 6;
    dynodeResponse = dynode_chain(6000,ndynodes,1*[0.25 4 3.69 3.69 2.19 2.19 2.19]);
    dynodeResponse = dynodeResponse(2:end);  
  end
  for kk = 1:size(clean_mov,3)
    countprint(kk);
    mov(:,:,kk) = DynodeNoiseModel(clean_mov(:,:,kk), noise_params, dynodeResponse); % Run data through the noise model    
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%