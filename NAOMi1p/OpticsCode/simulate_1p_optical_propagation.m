function opt_out = simulate_1p_optical_propagation(vol_params,psf_params,vol_out)

% opt_out = simulate_optical_propagation(vol_params,psf_params,vol_out)
%
% Function to simulate optical propgation through the simulated volume.
% Generates the PSF and appropriate blurring representing aberrations
% within the tissue. The inputs to this function are:
%
%   - vol_params        - Struct with parameters for the volume generation
%       .vol_sz         - 3-element vector with the size (in um) of the 
%                         volume to generate (default = 100x100x30um)
%       .vres           - resolution to simulate volume at (default = 2
%                         samples/um)
%       .vol_depth      - Depth of the volume under the brain surface
%                         (default = 100um)
%   - psf_params        - Struct contaning the parameters for the PSF
%       .scatter_sz     - Scattering object sizes (um), column vector
%       .scatter_wt     - Scattering object weights, column vector
%   - vol_out           - Struct containing simulated volume
%       .neur_ves_all   - Full simulated blood vessels (including
%                         surrounding areas)
%       .neur_ves       - Simulated blood vessels (only within volume)
%
% The output of this function is
%   - opt_out           - Struct containing the PSF-related variables
%       .psf            - 3D array containing the PSF illumination 
%                         profile
%       .mask           - Mask indicating the illumination blockage 
%                         across the field-of-view
%
% 2017- Alex Song and Adam Charles 

%  modified by YZ. last update: 9/27/2021. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

vol_params   = check_vol_params(vol_params);                               % Check volume parameters
psf_params   = check_psf_params(psf_params);                               % Check point spread function parameters
acc_flag     = 0;                                                          % Accuracy flag - set to 0 for light simulation (faster, less memory intensive)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create optical variables

vres   = vol_params.vres;                                                  % Get volume resolution (pixels per micron)
vol_sz = vol_params.vol_sz;                                                % Get volume size (in microns) 

vol_params.vasc_sz = gaussianBeamSize(psf_params,vol_params.vol_depth+ ...
    vol_params.vol_sz(3)/2)+vol_params.vol_sz+[0 0 1]*vol_params.vol_depth;% Get beam size
vasc_sz = vol_params.vasc_sz;

% Setup a temporary vasculature volume if neur_ves_all is the wrong size
TMPvasc = [];
if (isfield(vol_out,'neur_ves_all'))&&~isempty(vol_out.neur_ves_all)
  nvpx = size(vol_out.neur_ves_all);
  vcpx = vol_params.vasc_sz*vres;
  if(~all(nvpx==vcpx))
    TMPvasc = false(vcpx);
    if(all(nvpx<=vcpx));
      dpx = (vcpx-nvpx)/2;
      TMPvasc(1+floor(dpx(1)):end-ceil(dpx(1)),1+floor(dpx(2)):end-ceil(dpx(2)), ... 
        1:end-dpx(3)*2) = vol_out.neur_ves_all;
    elseif(all(nvpx>=vcpx))
      dpx = (nvpx-vcpx)/2;
      TMPvasc = vol_out.neur_ves_all(1+floor(dpx(1)):end-ceil(dpx(1)), ...
        1+floor(dpx(2)):end-ceil(dpx(2)), 1:end-2*dpx(3));
    else
      nvpx = size(vol_out.neur_ves)/vres;
      dpx = (vcpx-nvpx)/2;
      TMPvasc(1+floor(dpx(1)):end-ceil(dpx(1)),1+floor(dpx(2)):end-ceil(dpx(2)), ... 
        1:end-2*dpx(3)) = vol_out.neur_ves;
    end
  end
end



%% Setup scanned volume
if(vol_params.verbose>1)
    tsetup = tic;
end

if(acc_flag)
  % Setup optical wavefront
  if (~isfield(psf_params,'type'))||isempty(psf_params.type)||strcmp(psf_params.type,'gaussian')
    Uin = generateBA(vol_params,psf_params);
  else
    warning('Needs to be a supported type, defaulting to ''gaussian''')
    Uin = generateBA(vol_params,psf_params);
  end

  if (~isfield(vol_out,'neur_ves_all'))||~isequal(size(vol_out.neur_ves_all),vres*vasc_sz)
    vol = false(vres*vasc_sz);
    vol(floor(vres*(vasc_sz(1)-vol_sz(1))/2)+(1:vol_sz(1)*vres),floor(vres*(vasc_sz(2)-vol_sz(1))/2)+(1:vol_sz(2)*vres),:) = vol_out.neur_ves;
  else
    vol = vol_out.neur_ves_all;
  end
  
  if (isfield(psf_params,'scatter_sz'))&&(~isempty(psf_params.scatter_sz))
    if (isfield(psf_params,'scatter_wt'))&&(~isempty(psf_params.scatter_wt))
      if length(psf_params.scatter_wt)==1
        scatter_wt = psf_params.scatter_wt*ones(size(psf_params.scatter_sz,1),1);
      else
        scatter_wt = psf_params.scatter_wt;
      end
    end
    if(size(psf_params.scatter_sz,2)==3)
      scatter_sz = vres*psf_params.scatter_sz;
    elseif(size(psf_params.scatter_sz,2)==1)
      scatter_sz = vres*repmat(psf_params.scatter_sz,1,3);
    else
      error('Scattering sizes must be 1D or 3D (columns)');
    end
    vol = (~vol).*masked_3DGP_v2(vasc_sz*vres,scatter_sz,1,0,[],[],scatter_wt)+vol;
  end
  if(vol_params.verbose>1)
    fprintf('Volume setup computed (%f s)\n',toc(tsetup));
  end
  % Optical propagation
  if ~isstruct(Uin)
      
      % modified 1p here
    [mask,psf] = genCorticalLightPath_1p(vol_params,psf_params,vol,Uin);        % Generate diffraction-limited tissue obscuration mask
    if(psf_params.sampling>max(size(mask)))
      mask       = ones(size(mask),'like',mask);
    end
    mask              = mask/max(mask(:));
    opt_out.mask      = mask;
    opt_out.psf       = psf/mean(mask(:));
  else
    names = fieldnames(Uin);
    opt_out = struct('mask',struct,'psf',struct);
    for i = 1:length(names)
      [mask,psf]   = genCorticalLightPath_1p(vol_params,psf_params,vol,...
        getfield(Uin,names{i}));        % Generate diffraction-limited tissue obscuration mask
      mask         = mask/max(mask(:));
      psf          = psf/mean(mask(:));
      if(psf_params.sampling>max(size(mask)))
        mask       = ones(size(mask),'like',mask);
      end
      opt_out.mask = setfield(opt_out.mask,names{i},mask);
      opt_out.psf  = setfield(opt_out.psf,names{i},psf);
    end
  end
end

if(~acc_flag)
vol_depth = vol_params.vol_depth;
volpx = vol_sz*vres;
vascpx = vasc_sz*vres;
psfsz = psf_params.psf_sz;
psfpx = psfsz*vres;
proppx = psf_params.prop_sz*vres;
  
zA = vres*(vol_depth+vol_sz(3)/2)-psfpx(3)/2;
zB = vres*(vol_depth+vol_sz(3)/2)+psfpx(3)/2;
zAsz = [vasc_sz(1:2)*vres zA/proppx];
if((~isfield(psf_params,'propcrop'))||psf_params.propcrop)
  zBsz = vres*(2*gaussianBeamSize(psf_params,psfsz(3)/2,3)+[vol_params.vol_sz(1:2) psfsz(3)]);
  if(zBsz(1)>zAsz(1))
    zBsz(1:2) = zAsz(1:2);
  end
else
  zBsz = vres*(2*gaussianBeamSize(psf_params,psfsz(3)/2,3)+[vol_params.vol_sz(1:2) psfsz(3)]);
  zBsz(1:2) = zAsz(1:2);
end


if (isfield(psf_params,'taillength'))&&(~isempty(psf_params.taillength))
  zCsz = [vasc_sz(1:2)*vres psf_params.taillength*vres/proppx];
else
  zCsz = [vasc_sz(1:2)*vres (vasc_sz(3)*vres-zB)/proppx];
end

if (isfield(psf_params,'scatter_sz'))&&(~isempty(psf_params.scatter_sz))
    if (isfield(psf_params,'scatter_wt'))&&(~isempty(psf_params.scatter_wt))
        if length(psf_params.scatter_wt)==1
            scatter_wt = psf_params.scatter_wt*ones(size(psf_params.scatter_sz,1),1);
        else
            scatter_wt = psf_params.scatter_wt;
        end
    end
    if (isfield(psf_params,'n_diff_scatter'))&&(~isempty(psf_params.n_diff_scatter))
        scatter_wt = scatter_wt*psf_params.n_diff_scatter/psf_params.n_diff;
    end

    if(size(psf_params.scatter_sz,2)==3)
        scatter_sz = vres*psf_params.scatter_sz;
    elseif(size(psf_params.scatter_sz,2)==1)
        scatter_sz = vres*repmat(psf_params.scatter_sz,1,3);
    else
        error('Scattering sizes must be 1D or 3D (columns)');
    end
    neur_ves_A = masked_3DGP_test(ceil(zAsz),bsxfun(@times,scatter_sz,[1 1 1/proppx]),proppx,0,[],[],scatter_wt);
    if(rem(zAsz(3),1)~=0)
      neur_ves_A(:,:,1) = neur_ves_A(:,:,1)*rem(zAsz(3),1);
%       neur_ves_A(:,:,end) = neur_ves_A(:,:,end)*rem(zAsz(3),1);
    end
    neur_ves_B = masked_3DGP_test(zBsz,scatter_sz,1,0,[],[],scatter_wt);
    neur_ves_C = masked_3DGP_test(ceil(zCsz),bsxfun(@times,scatter_sz,[1 1 1/proppx]),proppx,0,[],[],scatter_wt);
    if(rem(zCsz(3),1)~=0)
      neur_ves_C(:,:,end) = neur_ves_C(:,:,end)*rem(zCsz(3),1);
    end
else
    neur_ves_A = zeros(ceil(zAsz),'single');
    neur_ves_B = zeros(ceil(zBsz),'single');
    neur_ves_C = zeros(ceil(zCsz),'single');
end

% Add in vasculature to propagation volume
if (~isfield(vol_out,'neur_ves_all'))%||~isequal(size(vol_out.neur_ves_all),vres*vasc_sz)
  TMPves = groupzproject(single(vol_out.neur_ves(:,:,1:zA)),proppx,'sum');
  neur_ves_A(floor((zAsz(1)-volpx(1))/2)+(1:volpx(1)), ...
    floor((zAsz(2)-volpx(1))/2)+(1:volpx(2)),:) = TMPves+ ...
    neur_ves_A(floor((zAsz(1)-volpx(1))/2)+(1:volpx(1)), ...
    floor((zAsz(2)-volpx(1))/2)+(1:volpx(2)),:).* ... 
    bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);

  TMPves = single(vol_out.neur_ves(:,:,zA+1:zB));
  neur_ves_B(floor((zBsz(1)-volpx(1))/2)+(1:volpx(1)), ...
    floor((zBsz(2)-volpx(1))/2)+(1:volpx(2)),:) = TMPves+ ...
    neur_ves_B(floor((zBsz(1)-volpx(1))/2)+(1:volpx(1)), ...
    floor((zBsz(2)-volpx(1))/2)+(1:volpx(2)),:).* ... 
    bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves);

  TMPves = groupzproject(single(vol_out.neur_ves(:,:,zB+1:end)),proppx,'sum');
  neur_ves_C(floor((zCsz(1)-volpx(1))/2)+(1:volpx(1)), ...
    floor((zCsz(2)-volpx(1))/2)+(1:volpx(2)),1:size(TMPves,3)) = TMPves+ ...
    neur_ves_C(floor((zCsz(1)-volpx(1))/2)+(1:volpx(1)), ...
    floor((zCsz(2)-volpx(1))/2)+(1:volpx(2)),1:size(TMPves,3)).* ... 
    bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);
elseif(isempty(TMPvasc))
  TMPves = groupzproject(single(vol_out.neur_ves_all(:,:,1:zA)),proppx,'sum');
  neur_ves_A = TMPves + neur_ves_A.* ... 
    bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);

  TMPves = single(vol_out.neur_ves_all(floor((vascpx(1)-zBsz(1))/2)+(1:zBsz(1)), ...
    floor((vascpx(2)-zBsz(2))/2)+(1:zBsz(2)),zA+1:zB));
  neur_ves_B = TMPves + neur_ves_B.* ... 
    bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves);

  TMPves = groupzproject(single(vol_out.neur_ves_all(:,:,zB+1:end)),proppx,'sum');
  if(size(TMPves,3)>size(neur_ves_C,3))
    neur_ves_C = TMPves(:,:,1:size(neur_ves_C,3)) + neur_ves_C.* ... 
      bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves(:,:,1:size(neur_ves_C,3))/proppx);
  elseif(size(TMPves,3)==size(neur_ves_C,3))
    neur_ves_C = TMPves + neur_ves_C.* ... 
      bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);
  else
    neur_ves_C(:,:,1:size(TMPves,3)) = TMPves + neur_ves_C(:,:,1:size(TMPves,3)).* ... 
      bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);    
  end
else
  TMPves = groupzproject(single(TMPvasc(:,:,1:zA)),proppx,'sum');
  neur_ves_A = TMPves + neur_ves_A.* ... 
    bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);

  TMPves = single(TMPvasc(floor((vascpx(1)-zBsz(1))/2)+(1:zBsz(1)), ...
    floor((vascpx(2)-zBsz(2))/2)+(1:zBsz(2)),zA+1:zB));
  neur_ves_B = TMPves + neur_ves_B.* ... 
    bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves);

  TMPves = groupzproject(single(TMPvasc(:,:,zB+1:end)),proppx,'sum');
  if(size(TMPves,3)>size(neur_ves_C,3))
    neur_ves_C = TMPves(:,:,1:size(neur_ves_C,3)) + neur_ves_C.* ... 
      bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves(:,:,1:size(neur_ves_C,3))/proppx);
  elseif(size(TMPves,3)==size(neur_ves_C,3))
    neur_ves_C = TMPves + neur_ves_C.* ... 
      bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);
  else
    neur_ves_C(:,:,1:size(TMPves,3)) = TMPves + neur_ves_C(:,:,1:size(TMPves,3)).* ... 
      bsxfun(@minus,ones(size(TMPves,1),size(TMPves,2),'single'),TMPves/proppx);    
  end
end
 
if(vol_params.verbose>1)
    fprintf('Volume setup computed (%f s)\n',toc(tsetup));
end

wvl = single(psf_params.lambda*1e-6); % optical wavelength [m]
k = 2*pi/wvl; % optical wavenumber [rad/m]
ndiff = psf_params.n_diff;
ss  = psf_params.ss;

%Setup optical wavefront
if (~isfield(psf_params,'type'))||isempty(psf_params.type)||strcmp(psf_params.type,'gaussian')
  Uin = generateBA(vol_params,psf_params);
else
  warning('Needs to be a supported type, defaulting to ''gaussian''')
  Uin = generateBA(vol_params,psf_params);
end


if ~isstruct(Uin)
  if(psf_params.fastmask)
    psf_paramsFM = psf_params;
    FMnames = fieldnames(psf_params.FM);
    for i = 1:length(FMnames)
      psf_paramsFM = setfield(psf_paramsFM,FMnames{i},eval(['psf_paramsFM.FM.' FMnames{i}]));
    end
    if (~isfield(psf_params,'type'))||isempty(psf_params.type)||strcmp(psf_params.type,'gaussian')
      UinFM = generateBA(vol_params,psf_paramsFM);
    else
      warning('Needs to be a supported type, defaulting to ''gaussian''')
      UinFM = generateBA(vol_params,psf_paramsFM);
    end
    ves_A_FM = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_A,psf_paramsFM.ss));
    ves_B_FM = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_B,psf_paramsFM.ss));
    ves_C_FM = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_C,psf_paramsFM.ss));    
    
    [mask,~,~,psfTM,psfBM] = genCorticalLightPathLite_1p(vol_params,...
                          psf_paramsFM, ves_A_FM,ves_B_FM,ves_C_FM,UinFM); % Generate diffraction-limited tissue obscuration mask
    clear ves_A_FM ves_B_FM ves_C_FM
    neur_ves_A = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_A,ss));
    neur_ves_B = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_B,ss));
    neur_ves_C = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_C,ss));

    [~,psf,~,psfT,psfB] = genCorticalLightPathLite_1p(vol_params,psf_params, ...
                                        neur_ves_A,neur_ves_B,neur_ves_C,Uin);        % Generate diffraction-limited tissue obscuration mask
    psfT.mask = psfTM.mask/mean(psfTM.mask(:));
    psfB.mask = psfBM.mask/mean(psfBM.mask(:));
  else
   
    neur_ves_A = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_A,ss));
    neur_ves_B = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_B,ss));
    neur_ves_C = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_C,ss));

    [mask,psf,~,psfT,psfB] = genCorticalLightPathLite_1p(vol_params,psf_params, ...
                                        neur_ves_A,neur_ves_B,neur_ves_C,Uin);        % Generate diffraction-limited tissue obscuration mask
    if(psf_params.sampling>max(size(mask)))
      mask       = ones(size(mask),'like',mask);
    end
  end
  mask              = mask/mean(mask(:));
  opt_out.mask      = mask;
  opt_out.psf       = psf/mean(mask(:));
  opt_out.psfT      = psfT;
  opt_out.psfB      = psfB;
else
  names = fieldnames(Uin);
  opt_out = struct('mask',struct,'psf',struct,'psfT',struct,'psfB',struct);
  if(psf_params.fastmask)
    psf_paramsFM = psf_params;
    FMnames = fieldnames(psf_params.FM);
    for j = 1:length(FMnames)
      psf_paramsFM = setfield(psf_paramsFM,FMnames{j},eval(['psf_paramsFM.FM.' FMnames{j}]));
    end
    ves_A_FM = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_A,psf_paramsFM.ss));
    ves_B_FM = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_B,psf_paramsFM.ss));
    ves_C_FM = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_C,psf_paramsFM.ss));
    neur_ves_A = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_A,ss));
    neur_ves_B = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_B,ss));
    neur_ves_C = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_C,ss));
    if (~isfield(psf_params,'type'))||isempty(psf_params.type)||strcmp(psf_params.type,'gaussian')
      UinFM = generateBA(vol_params,psf_paramsFM);
    else
      warning('Needs to be a supported type, defaulting to ''gaussian''')
      UinFM = generateBA(vol_params,psf_paramsFM);
    end
    for i = 1:length(names)
      [mask,~,~,psfTM,psfBM] = genCorticalLightPathLite_1p(vol_params,psf_paramsFM, ...
        ves_A_FM,ves_B_FM,ves_C_FM,getfield(UinFM,names{i}));        % Generate diffraction-limited tissue obscuration mask
      
      [~,psf,~,psfT,psfB] = genCorticalLightPathLite_1p(vol_params,psf_params, ...
        neur_ves_A,neur_ves_B,neur_ves_C,getfield(Uin,names{i}));        % Generate diffraction-limited tissue obscuration mask
      psfT.mask = psfTM.mask/mean(psfTM.mask(:));
      psfB.mask = psfBM.mask/mean(psfBM.mask(:));
      opt_out.mask = setfield(opt_out.mask,names{i},mask);
      opt_out.psf  = setfield(opt_out.psf,names{i},psf);
      opt_out.psfT = setfield(opt_out.psfT,names{i},psfT);
      opt_out.psfB = setfield(opt_out.psfB,names{i},psfB);
    end
  else
    neur_ves_A = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_A,ss));
    neur_ves_B = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_B,ss));
    neur_ves_C = exp(1i*imresize(ndiff*k*1e-6/vres*neur_ves_C,ss));
    for i = 1:length(names)
      [mask,psf,~,psfT,psfB] = genCorticalLightPathLite_1p(vol_params,psf_params, ...
        neur_ves_A,neur_ves_B,neur_ves_C,getfield(Uin,names{i}));        % Generate diffraction-limited tissue obscuration mask
      mask         = mask/mean(mask(:));
      psf          = psf/mean(mask(:));
      if(psf_params.sampling>max(size(mask)))
        mask       = ones(size(mask),'like',mask);
      end
      opt_out.mask = setfield(opt_out.mask,names{i},mask);
      opt_out.psf  = setfield(opt_out.psf,names{i},psf);
      opt_out.psfT = setfield(opt_out.psfT,names{i},psfT);
      opt_out.psfB = setfield(opt_out.psfB,names{i},psfB);
    end
  end
end
end

%% Calculate the collection mask (absorption on collection side)
colpx = vascpx;
colpx(3) = ceil(colpx(3)/proppx);

if (~isfield(vol_out,'neur_ves_all'))
  TMPvasc = zeros(colpx,'single');
  TMPvasc(floor((colpx(1)-volpx(1))/2)+(1:volpx(1)), floor((colpx(2)-volpx(1))/2)+(1:volpx(2)),:) = groupzproject(single(vol_out.neur_ves(:,:,1:size(vascpx,3))),proppx,'sum');  
elseif(isempty(TMPvasc))
  TMPvasc = groupzproject(single(vol_out.neur_ves_all),proppx,'sum');  
else
  TMPvasc = groupzproject(single(TMPvasc),proppx,'sum');
end

coldist = vres*tan(asin(psf_params.objNA/psf_params.n))*((vol_depth+vol_sz(3)/2)-(proppx/vres*((1:size(TMPvasc,3))-0.5)));
colmask = zeros(colpx(1),colpx(2),'single');
[X,Y] = meshgrid(-colpx(1)/2:colpx(1)/2-1,-colpx(2)/2:colpx(2)/2-1);
rho = sqrt(X.^2+Y.^2);
for i = 1:length(coldist)
  if(coldist(i)>0)
    se = autocrop(single(rho<=coldist(i)));
    se = se/sum(se(:));
    colmask = colmask+conv2(TMPvasc(:,:,i),se,'same');
  end
end
colmask = colmask(floor((colpx(1)-volpx(1))/2)+(1:volpx(1)), floor((colpx(2)-volpx(2))/2)+(1:volpx(2)));

opt_out.colmask = 10.^(-colmask/vres*psf_params.hemoabs);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
