function [mask,psfs3,psfs,psf2p2] = genCorticalLightPath_1p(vol_params,psf_params,vol,Uin)
 
% [mask,psfs3,psfs] = genCorticalLightPath(vol_params,psf_params,vol_out,Uin)
%
% This function generates a point-by-point map of obscuration for TPM
% imaging setting across a neural volume. The function generates a 3D mask
% that modulates the intensity along the light paths. The inputs are
% 
%   - vol_params      - Struct with parameters for the volume generation
%       .vol_sz       - 3-element vector with the size (in um) of the 
%                       volume to generate (default = 100x100x30um)
%       .vres         - resolution to simulate volume at (default = 2
%                       samples/um)
%       .vol_depth    - Depth of the volume under the brain surface
%   - psf_params      - Struct contaning the parameters for the PSF
%       .n_diff       - Shift in refractive index from vessels to tissue
%       .lambda       - Two-photon excitation wavelength (um)
%       .obj_fl       - Objective focal length (mm)
%       .ss           - Subsampling factor for fresnel propagation
%       .sampling     - Spatial sampling for tissue occlusion mask
%       .psf_sz       - Default two-photon PSF size simulated (um)
%       .prop_sz      - Fresnel propagation length outside of volume (um)
%   - vol             - Simulated volume impacting propagation (vessels)
%   - Uin             - Input scalar field
%
% The outputs are
%   - mask            - 2D mask giving relative two-photon excitation at
%                       each position laterally
%   - psfs3           - Average aberrated PSF across the simulated field
%   - psfs            - All PSFs at each simulated position
%
% 2016 - Alex Song and Adam Charles

% modified by YZ. last update: 12/17/2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for light path obscurations

if(nargout>3)
  psf2p2 = [];
end
vres      = vol_params.vres;
vol_sz    = vol_params.vol_sz;
vol_depth = vol_params.vol_depth;
vasc_sz   = vol_params.vasc_sz;
verbose   = vol_params.verbose;

fl  = single(psf_params.obj_fl/1000); % focal length [m]
ss  = psf_params.ss;
D2  = single(1e-6*(1/vres)/ss); % observation grid spacing [m]
N   = single(1e-6*(vasc_sz(1:2)-vol_sz(1:2))/D2);
D1  = single(max(gaussianBeamSize(psf_params,fl*1e6)/1e6)/min(N)); % source grid spacing [m]
N   = N(1);
nre = single(psf_params.n);
z = single(fl-(vol_depth+vol_sz(3)/2)*1e-6); % propagation distance [m]
wvl = single(psf_params.lambda*1e-6); % optical wavelength [m]
psf_samp = min(psf_params.sampling,1e10);
k = 2*pi/wvl; % optical wavenumber [rad/m]
psfpx = psf_params.psf_sz*vres;
proppx = psf_params.prop_sz*vres;
ndiff = psf_params.n_diff;

z = [0 z];
delta = [D1 D2];
[x1, y1] = meshgrid((-N/2 : N/2-1) * D1);
sg = exp(-(x1/(0.47*N*D1)).^16) .* exp(-(y1/(0.47*N*D1)).^16);
t = repmat(sg, [1 1 2]);
if(~iscell(Uin))
  Uout = fresnel_propagation_multi(Uin, wvl, delta, z, t, nre);
%   Uout = (N/ss)^2*Uout/(sum(abs(Uout(:)).^2));
  Uout = Uout/sqrt(sum(abs(Uout(:)).^2));
end
if verbose == 1
  fprintf('Calculating mask layer...');
elseif verbose >1
  fprintf('Calculating mask layer...\n');
end
imax = round(vol_sz(1)/psf_samp)+1;
jmax = round(vol_sz(2)/psf_samp)+1;
x2 = [];
zA = vres*(vol_depth+vol_sz(3)/2)-psfpx(3)/2;
zB = vres*(vol_depth+vol_sz(3)/2)+psfpx(3)/2;
if(rem(zA,proppx))
  neur_ves_A = imresize(ndiff*k*1e-6/vres*cat(3,groupzproject(single(vol(:,:,1:zA-rem(zA,proppx))),proppx,'sum'), ...
    sum(single(vol(:,:,zA-rem(zA,proppx)+1:zA)),3)),ss);
  zApos = proppx/vres*[(0:size(neur_ves_A,3)-1), zA/(proppx)]*1e-6;
else
  neur_ves_A = imresize(ndiff*k*1e-6/vres*groupzproject(single(vol(:,:,1:zA)),proppx,'sum'),ss);
  zApos = proppx/vres*(0:size(neur_ves_A,3))*1e-6;
end
neur_ves_B = imresize(ndiff*k*1e-6/vres*single(vol(:,:,zA+1:zB)),ss);

if(verbose>1)
  fprintf('Propagating through %d locations:\n',imax*jmax)
end
psfs = cell(imax,jmax);
for i = 1:imax
  for j = 1:jmax
    if(verbose>1)
      tloop = tic;
    end
    
    phzA = neur_ves_A((1:N)+(psf_samp*vres*ss*(i-1)),(1:N)+(psf_samp*vres*ss*(j-1)),:);
    phzA = cat(3,zeros(N,N),phzA);
    phzA = exp(1i*phzA);
    phzA = bsxfun(@times,sg,phzA);
%     tA = repmat(sg,[1 1 size(phzA,3)+1]).*cat(3,ones(N,N),exp(1i*phzA));
%     UoutA = fresnel_propagation_multi(Uout, wvl, D2*ones(size(phzA,3)+1,1), zApos, tA, nre);

    UoutA = fresnel_propagation_multi(Uout, wvl, D2*ones(size(phzA,3),1), zApos, phzA, nre);
    if((~isfield(psf_params,'propcrop'))||psf_params.propcrop)
      N2 = min(max(gaussianBeamSize(psf_params,psfpx(3)/vres/2,3)/1e6)/(1e-6/vres)*2*ss,N);    
    else
      N2 = N;
    end
    phzB = neur_ves_B((N/2-N2/2+1:N/2+N2/2)+(psf_samp*vres*ss*(i-1)),(N/2-N2/2+1:N/2+N2/2)+(psf_samp*vres*ss*(j-1)),:);    
    if(isempty(x2))
      [x2, y2] = meshgrid((-N2/2 : N2/2-1) * D1);
      sg2 = exp(-(x2/(0.47*N2*D1)).^16) .* exp(-(y2/(0.47*N2*D1)).^16);
    end
    
    phzB = cat(3,zeros(N2,N2),phzB);
    phzB = exp(1i*phzB);
    phzB = bsxfun(@times,sg2,phzB);
%     tB = repmat(sg2,[1 1 size(phzB,3)+1]).*cat(3,ones(N2,N2),exp(1i*phzB));
%     [~, UoutAll] = fresnel_propagation_multi(UoutA, wvl, D2*ones(size(phzB,3)+1,1), (0:size(phzB,3))*1e-6/vres, tB, nre);
    UoutA = UoutA(N/2-N2/2+1:N/2+N2/2,N/2-N2/2+1:N/2+N2/2);
    [~, UoutAll] = fresnel_propagation_multi(UoutA, wvl, D2*ones(size(phzB,3)+1,1), (0:size(phzB,3))*1e-6/vres, phzB, nre);

    psf2p = UoutAll(N2/2-psfpx(1)*ss/2+1:N2/2+psfpx(1)*ss/2,N2/2-psfpx(2)*ss/2+1:N2/2+psfpx(2)*ss/2,1:end-1);
    if((~isfield(psf_params,'scaling'))||strcmp(psf_params.scaling,'widefield'))
      psf2p = ss^2*abs(psf2p).^2;
    else
      warning('Needs to be a specified scaling, defaulting to ''widefield''')
      psf2p = ss^2*abs(UoutAll(:,:,1:end-1)).^2;
    end
    psfs{i,j} = ss^2*imresize(imtranslate(psf2p,ss/2-[0.5 0.5]),1/ss)*(vres*(1e6*wvl)^1.5)/(pi*nre);
    if(verbose>1)
      fprintf('Propagation %d finished (%f s)\n',(i-1)*imax+j,toc(tloop));
    end
    if(nargout>3)
      if(isempty(psf2p2))
        psf2p2 = psf2p;
      else
        psf2p2 = psf2p2+psf2p;
      end
    end
  end
end

psfs3 = zeros(size(psfs{1,1}));
for i = 1:imax
  for j = 1:jmax
    psfs3 = psfs3+(abs(psfs{i,j}));
  end
end
psfs3 = psfs3/(imax*jmax);
psfmag = zeros([imax jmax]);
for i = 1:imax
  for j = 1:jmax
    psfmag(i,j) = sum(psfs{i,j}(:));
  end
end

[X,Y] = meshgrid(double(1:vol_sz(1)*vres)-0.5,double(1:vol_sz(1)*vres)-0.5);
[x,y] = meshgrid(double(psf_samp*vres*(0:imax-1)),double(psf_samp*vres*(0:jmax-1)));
mask  = single(griddata(x,y,psfmag,X,Y,'v4'));

fprintf('done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%