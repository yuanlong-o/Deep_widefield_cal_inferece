function abb = generateZernike(psf_params,offset)

% abb = generateZernike(psf_params,offset)
%  
% Generates the Zernike aberration weights given a input distribution and
% offset. The inputs are
% 
% - psf_params        - Struct contaning the parameters for the PSF
%        .lambda      - Two-photon excitation wavelength (um)
%        .zernikeWt   - Microscope aberration weights (Zernike basis)
%        .zernikeDst  - Microscope aberration weights (Zernike basis) as
%                       a function of position
% - offset            - offset for calculating aberrations at a distance
%                       away from the center of field
%
% The output is
% - abb               - Zernike abberation weights
% 
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<2)
  offset = [0 0];
end

if((~isfield(psf_params,'zernikeWt'))||isempty(psf_params.zernikeWt))
  abb = 0;
else
  if((~isfield(psf_params,'zernikeDst'))||isempty(psf_params.zernikeDst))
    abb = psf_params.zernikeWt;
  else
    abb = psf_params.zernikeDst(offset).*psf_params.zernikeWt;    
  end
end
abb = abb*psf_params.lambda*1e-6;