function Uout2 = generateBA(vol_params,psf_params)

% Uout2 = generateBA(vol_params,psf_params)
%
% This function generates a Gaussian back aperture intensity profile with 
% the specified aberrations applied. The inputs are
% 
% - vol_params          - Struct with parameters for the volume generation
%          .vol_sz      - 3-element vector with the size (in um) of the 
%                         volume to generate (default = 100x100x30um)
%          .vres        - resolution to simulate volume at (default = 2
%                         samples/um)
% - psf_params          - Struct contaning the parameters for the PSF
%          .NA          - Numerical aperture of Gaussian beam
%          .objNA       - Numerical aperture of objective lens
%          .n           - Refractive index of propagation material
%          .lambda      - Two-photon excitation wavelength (um)
%          .obj_fl      - Objective focal length (mm)
%          .ss          - Subsampling factor for fresnel propagation
%          .sampling    - Spatial sampling for tissue occlusion mask
%          .zernikeDst  - Microscope aberration weights (Zernike basis) as
%                         a function of position
% The output is
% - Uout2               - Output scalar field
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isfield(vol_params,'vasc_sz'))||isempty(vol_params.vasc_sz)
  vol_params.vasc_sz = gaussianBeamSize(psf_params,vol_params.vol_depth+...
  vol_params.vol_sz(3)/2)+vol_params.vol_sz+[0 0 1]*vol_params.vol_depth;% Get beam size
end

fl = single(psf_params.obj_fl/1000);                                       % focal length [m]
D2 = single(1e-6*(1/vol_params.vres)/psf_params.ss);                       % observation grid spacing [m]
N = single(1e-6*(vol_params.vasc_sz(1:2)-vol_params.vol_sz(1:2))/D2);      % gridsize
D1 = single(max(gaussianBeamSize(psf_params,fl*1e6)/1e6)/min(N));          % source grid spacing [m]
nre = single(psf_params.n);                                                % immersion numerical index
rad = single(tan(asin(psf_params.NA/nre))*fl);                             % source radius [m]
objrad = single(tan(asin(psf_params.objNA/nre))*fl);                       % source radius [m]
k = 2*nre*pi/single(psf_params.lambda*1e-6);                               % optical wavenumber [rad/m]
[X,Y] = meshgrid((-N(1)/2:N(1)/2-1)*D1,(-N(2)/2:N(2)/2-1)*D1);             
Uout = generateGaussianProfile(X,Y,rad,objrad,k,fl);                       % generate gaussian wavefront with apodization 1

imax = round(vol_params.vol_sz(1)/psf_params.sampling)+1;
jmax = round(vol_params.vol_sz(2)/psf_params.sampling)+1;

if(imax*jmax==1||(~isfield(psf_params,'zernikeDst'))||isempty(psf_params.zernikeDst))
  abb = generateZernike(psf_params);
  Uout2 = applyZernike(Uout,X/objrad,Y/objrad,k,abb);
else
  Uout2 = cell(imax,jmax);
  for i = 1:imax
    for j = 1:jmax
      abb = generateZernike(psf_params);
      Uout2{i,j} = applyZernike(Uout,X,Y,k,abb);
    end
  end
end