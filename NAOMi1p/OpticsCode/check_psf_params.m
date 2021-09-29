function psf_params = check_psf_params(psf_params)

% psf_params = check_psf_params(psf_params)
%  
% This function checks the elements of the struct psf_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - psf_params        - Struct contaning the parameters for the PSF
%          .NA          - Numerical aperture of Gaussian beam
%          .n           - Refractive index of propagation material
%          .n_diff      - Shift in refractive index from vessels to tissue
%          .lambda      - Two-photon excitation wavelength (um)
%          .obj_fl      - Objective focal length (mm)
%          .ss          - Subsampling factor for fresnel propagation
%          .sampling    - Spatial sampling for tissue occlusion mask
%          .psf_sz      - Default two-photon PSF size simulated (um)
%          .prop_sz     - Fresnel propagation length outside of volume (um)
%          .blur        - PSF lateral blurring (um)
%          .scatter_sz  - Scattering object sizes (um), column vector
%          .scatter_wt  - Scattering object weights, column vector
%          .zernikeWt   - Microscope aberration weights (Zernike basis)
%
% 2017 - Adam Charles and Alex Song


% modified by YZ. 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(psf_params)
    clear psf_params
    psf_params = struct;
end
dParams.objNA      = 0.3;                                                  % Default objective numerical aperture
dParams.NA         = 0.3;                                        % Default excitation numerical aperture                                                
dParams.n          = 1.35;                                                 % Default index of refraction in tissue
dParams.n_diff     = 0.03;                                                 % Default shift in index of refraction from vessels to tissue. We set 0.03 for better results
dParams.lambda     = 0.488;                                                 % Default one-photon excitation wavelength (microns)
dParams.obj_fl     = 18;                                                  % Default objective focal length (mm), for a 10x objective

% sampling related
dParams.ss         = 2;                                                    % Default subsampling factor for fresnel propagation (from volume voxel size)
dParams.sampling   = 50;                                                   % Default spatial sampling for tissue occlusion mask
dParams.psf_sz     = [36 36 100];                                           % Default one-photon PSF size simulated (microns), this is typically for 200 um range of a 10x objective
dParams.prop_sz    = 10;                                                   % Default fresnel propagation length outside of volume (microns)
dParams.blur       = 3;                                                    % Default psf lateral blurring (microns)
dParams.scatter_sz = [0.51 1.56 4.52 14.78]';                              % Default scattering object sizes (microns), column vector
dParams.scatter_wt = [0.57 0.29 0.19 0.15]';                               % Default scattering object weights, column vector
dParams.zernikeWt  = [0 0 0 0 0 0 0 0 0 0 0];                         % Default microscope aberration weights (Zernike basis). In units of wavelength, a small amount of spherical aberration and astigmatism added as uncorrected "system" aberrations
dParams.taillength = 50;                                                   % Distance from edge of PSF_sz to estimate tailweight (um)    
dParams.type       = 'gaussian';                                           % Default PSF type ('gaussian','vtwins','bessel')          
dParams.scaling    = 'widefield';                                         % Default PSF scaling type ('widefield' only)
dParams.hemoabs    = 0.00674*log(10);                                      % Hemoglobin absorbance scaling factor. Default assumes 150mg/ml Hb, 64500 g/mol Hb, 2.9 (abs/um)/(mol/L) in units of abs/um. Absorbance calculated from Scott Prahl's Hb curve and eGFP emission spectrum
dParams.propcrop   = true;                                                 % Flag to crop scanned beam during optical propagation (default true) 
dParams.fastmask   = true;

psf_params = setParams(dParams, psf_params);
if psf_params.fastmask
    psf_params.FM.sampling = 10;
    psf_params.FM.fineSamp = 2;
    psf_params.FM.ss       = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

