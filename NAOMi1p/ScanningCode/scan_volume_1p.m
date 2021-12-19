function varargout = scan_volume_1p(neur_vol, PSF_struct, neur_act, ...
                    vol_params, scan_params, noise_params, spike_opts,  wdm_params, pixel_size, exp_level, output_dir)

% mov = scan_volume(neur_vol, PSF, neur_act, scan_params, varargin)
%
% Scan a volume and create a simulated two-photon calcium imaging movie.
% The inputs to this function are: 
%   neur_vol     - A struct that contains the output of the volume
%                  generating code. 
%   PSF          - A struct representing the point-spread function of the
%                  microscopy setup
%   neur_act     - Struct containing the simulated fluorescence activity 
%                  for each of the K neurons for all nt time-steps
%    .soma       - Kx(nt) matrix where each row is the activity of each 
%                  cell at the soma 
%    .dend       - Kx(nt) matrix where each row is the activity of each 
%                  cell at the dendrites [can be optional]
%    .bg         - N_bgx(nt) matrix where each row is the activity of each 
%                  background/neuropil component [can be optional]
%   scan_params  - Struct contaning the scanning parameters. Includes
%    .scan_avg   - Sampling rate of the scanning in terms of how many
%                  granular pixels to scan into one pixel (default = 2)
%    .motion     - True/false option for whether or not to simulate
%                  motion while performing scanning (default = true)
%    .scan_buff  - Number of granular pixels to keep as a buffer from the
%                  edge of the volume (default = 10)
%    .verbose    - Level of verbosity in the output during the volume
%                  generation. Can be 0,1,2. 0 = no text updates, 1 = 
%                  some text outputs. 2 = detailed text outputs (default
%                  = 1)
%   noise_params - Struct containing parameters for the noise model
%    .mu         - Mean measurement increase per photon (default = 0.5)
%    .mu0        - Electronics offset (default = 0)
%    .sigma      - Variance increase per photon (default = 0.3)
%    .sigma0     - Electronics base noise variance (default = 0.1)
%
% The outputs of this function is: 
%   mov         - A 3D array where each 2D slice is one frame in the video
%                 recroding of the simulated microscopy data
%
% 2016 - Adam Charles and Alex Song

%  modified by YZ. last update: 9/28/2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs

if ~isfield(PSF_struct, 'psf')
    error('Must provide PSF to scan!')
else
    PSF    = PSF_struct.psf;                                               % Extract point spread function from PSF struct
end
if ~isfield(PSF_struct, 'mask')
    t_mask = [];                                                           % No masking if mask is not supplied
else
    t_mask = PSF_struct.mask;                                              % Extract mask from PSF struct
    t_thresh = 1e-5;
    t_mask(t_mask<t_thresh) = t_thresh;
end
if isfield(PSF_struct, 'colmask')
    t_mask = t_mask.*PSF_struct.colmask;                                   % No masking if mask is not supplied
end

if ~isfield(PSF_struct, 'g_blur')
    g_blur = [];                                                            % No additional blurring if transversal blur function not supplied
else
    g_blur = PSF_struct.blur;                                              % Extract point transversal blur function from PSF struct
end

start_Nt = 100;    

scan_params  = check_imaging_params(scan_params);                             % Check the scanning parameter struct for missing elements
noise_params = check_noise_params(noise_params);                           % Check the noise parameter struct for missing elements
spike_opts   = check_spike_opts(spike_opts);                               % Check spike/fluorescence simulation parameters
wdm_params   = check_wdm_params(wdm_params);

if(floor(scan_params.sfrac)~=scan_params.sfrac)
    warning('Binning parameter is not an integer - this causes local interpolation of pixel positions');
end

vol_sz_px             = size(neur_vol.neur_vol);                           % Get image size (equivalent to vol_params.vol_sz*vol_params.vres

%%% CURRENT SIGSCALE NEEDS TO BE MODIFIED TO BETTER ACCOUNT FOR DIFFERENT
%%% VOLUME SIZES   
% calculate signal scale
noise_params.sigscale = wdmSignalscale(wdm_params,[], vol_params)*(spike_opts.dt)* ...
                          (scan_params.sfrac^2);
                      
if ~isstruct(neur_act)                                                     % Make sure that neur_act is a struct
    if size(neur_act,3)==1
        neur_act = struct('soma',neur_act);                                % Can be a Kx(nt)x1 array if only soma activity is provided
    elseif size(neur_act,3)==2
        neur_act = struct('soma',neur_act(:,:,1),'dend',neur_act(:,:,2));  % Can be a Kx(nt)x2 array if soma and dendrite activity is provided
    end
end
if (~isfield(neur_act,'dend'))||(isempty(neur_act.dend))
    neur_act.dend = neur_act.soma;                                         % Default to dendrites having the same activity as somas, if no dendrite activity is provided
end
if (~isfield(neur_act,'bg'))||(isempty(neur_act.bg))
    neur_act.bg = zeros(1,size(neur_act,2),'single');                      % Default to no background if background is not provided
end
if (isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1)
    neur_act.nuc  = neur_act.soma;
    neur_act.soma = 0*neur_act.soma;
    neur_act.dend = 0*neur_act.dend;
    neur_act.bg   = 0*neur_act.bg;
end

if (~isfield(scan_params,'fsimPath'))
  scan_params.fsimPath = [];
end

if (~isfield(scan_params,'fsimCleanPath'))
  scan_params.fsimCleanPath = [];
end

if (~isfield(scan_params,'saveBlocksize'))||(isempty(scan_params.saveBlocksize))
  scan_params.saveBlocksize = 1000;
end

if (~isfield(scan_params,'movout'))||(isempty(scan_params.movout))
  scan_params.movout = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize data

scan_buff = scan_params.scan_buff;                                         % Get scanning side buffer number (how much to stay away from the edges of the volume)
mot_opt   = scan_params.motion;                                            % Get scanning motion option (true or false) to decide whether to simulate tissue-FOV motion
scan_avg  = scan_params.scan_avg;                                          % Get scanning stepping amount (how many sub-resolution steps between each sample point in the FOV)
sfrac     = scan_params.sfrac;                                             % Subsampling factor

if (~isfield(scan_params,'vol_sz'))||(isempty(scan_params.vol_sz))
    N1 = size(neur_vol.neur_vol,1);                                        % Get length dimension
    N2 = size(neur_vol.neur_vol,2);                                        % Get width dimension
    N3 = size(neur_vol.neur_vol,3);                                        % Get depth dimension
else
    N1 = scan_params.vol_sz(1);                                            % Get length dimension
    N2 = scan_params.vol_sz(2);                                            % Get width dimension
    N3 = scan_params.vol_sz(3);                                            % Get depth dimension
end
Nt = size(neur_act.soma,2);                                                % Get number of time steps
end_Nt =  Nt  - start_Nt;
if ~isstruct(PSF)                                                          % Get size of the PSF
    [Np1,Np2,Np3] = size(PSF);
elseif isfield(PSF,'left')
    [Np1,Np2,Np3] = size(PSF.left);                                        % Reserved for vTwINS testing
else
    error('Unknown input configuration!')
end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform some error checking

if (N1 < Np1)||(N2 < Np2)
    error('PSF extent is bigger than the volume!')                         % Check that the PSF fits inside the volume (transversally)
end
if (N3 < Np3)
    error('PSF depth is larger than the volume depth!')                    % Check that the PSF fits inside the volume (axially)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize movie and temporary volume

scale_ratio = vol_params.vres / pixel_size;

mov_w_bg = zeros(floor([(N1)/sfrac * scale_ratio, (N2)/sfrac * scale_ratio, Nt]), 'single'); % Initialize movie array
mov_wo_bg = zeros(floor([(N1)/sfrac * scale_ratio, (N2)/sfrac * scale_ratio, Nt]), 'single'); % Initialize movie array

% initial volume
WMPvol = zeros(N1,N2,N3,'single');                                         % Initialize temporary volume
if (~isfield(scan_params,'zoffset'))||(isempty(scan_params.zoffset))
  z_base = floor(0.5*(N3-Np3));
else
  z_base = floor(0.5*(N3-Np3))+scan_params.zoffset;
end
z_loc  = z_base;                                                           % Set scan depth
x_loc  = floor(scan_buff+1);                                               % Set initial scan location in x
y_loc  = floor(scan_buff+1);                                               % Set initial scan location in y

if mot_opt                                                                 % If motion is requested, initialize the motion simulation variables
    zmaxdiff = 2;%5;
    d_stps   = [-1,1,zeros(1,5)];                                          % Set the vector of steps to sample from when choosing how much to move the volume relative to the FOV at each frame
    d_stpsZ  = [-1,1,zeros(1,100)];                                        % Set the vector of steps to sample from when choosing how much to move the volume (axially) relative to the FOV at each frame
    d_stps2  = -3:3;                                                       % Set the list of admissible steps (in granular resolution steps) that the volume can move relative to the FOV due to larger motion jumps
    p_jump   = 0.05;                                                       % Set the probability of a bigger jump in the FOV
    maxshear = 1/200;                                                      % fraction of max shearing that can be from fast y motion
    if nargout > 2
        mot_hist = zeros(3,Nt);                                            % If needed, initialize an array that will eventually store all the motion information over time. 
    end
else                                                                       % ----- 
    zmaxdiff = 0;
    d_stps   = [0,0,0];                                                    %   | 
    d_stpsZ  = [0,0,0];                                                    % If no motion is requested, set all the motion simulation variables to zero (this will stop all motion)
    d_stps2  = [0,0,0];                                                    %   |
    p_jump   = 0;                                                          %   |
    maxshear = 0;                                                          % -----
end

if (isfield(PSF_struct,'psfT'))&&(~isempty(PSF_struct.psfT))
    psfT          = PSF_struct.psfT;
    psfB          = PSF_struct.psfB;
    psfT.mask     = psfT.mask/mean(psfT.mask(:));
    psfB.mask     = psfB.mask/mean(psfB.mask(:));
    psfT.convmask = psfT.convmask/sum(psfT.convmask(:));
    psfB.convmask = psfB.convmask/sum(psfB.convmask(:));
    psfT.freq_psf = psf_fft(size(WMPvol), psfT.convmask);  
    psfB.freq_psf = psf_fft(size(WMPvol), psfB.convmask);  
else
    psfT = [];
    psfB = [];
end

freq_psf = psf_fft(size(WMPvol), PSF, scan_avg);                           % Pre-calculate the FFT of the PSF for faster scanning
sigscale = noise_params.sigscale;                                          % scaling factor for signal
somaVol  = cell(size(neur_vol.gp_vals,1),2);                               % --- 
dendVol  = cell(size(neur_vol.gp_vals,1),2);                               %  |
for i = 1:size(neur_vol.gp_vals,1)                                         %  |
    somaVol{i,1} = neur_vol.gp_vals{i,1}(neur_vol.gp_vals{i,3});           % Pre-alocate separate soma/dendrite indexing for faster activity modulation in the volume
    dendVol{i,1} = neur_vol.gp_vals{i,1}(~neur_vol.gp_vals{i,3});          %  |
    if(isempty(t_mask))
        somaVol{i,2} = neur_vol.gp_vals{i,2}(neur_vol.gp_vals{i,3});       %  |
        dendVol{i,2} = neur_vol.gp_vals{i,2}(~neur_vol.gp_vals{i,3});      %  |
    else
        somaVol{i,2} = neur_vol.gp_vals{i,2}(neur_vol.gp_vals{i,3}).* ...
                                    t_mask(max(1,mod(somaVol{i,1},N1*N2)));%  |
        dendVol{i,2} = neur_vol.gp_vals{i,2}(~neur_vol.gp_vals{i,3}).* ...
                                    t_mask(max(1,mod(dendVol{i,1},N1*N2)));%  |      
    end
end                                                                        % ---
if(~isempty(neur_vol.bg_proc))                                             % t_mask multiplied by gp_vals{i,2}
    if(isempty(t_mask))
        axonVol = neur_vol.bg_proc;
    else
        axonVol = cell(size(neur_vol.bg_proc));
        for i = 1:size(neur_vol.bg_proc,1)
            axonVol{i,1} = neur_vol.bg_proc{i,1};
            axonVol{i,2} = neur_vol.bg_proc{i,2}.* ...
                 t_mask(max(1,mod(axonVol{i,1},N1*N2)));
        end
    end
else
    axonVol = [];
end
if(isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1)
    if(isempty(t_mask))
        nucVol = neur_vol.gp_nuc;
    else
        nucVol = cell(size(neur_vol.gp_nuc));
        for i = 1:size(neur_vol.gp_nuc,1)
            nucVol{i,1} = neur_vol.gp_nuc{i,1};
            nucVol{i,2} = t_mask(max(1,mod(nucVol{i,1},N1*N2)));
        end
    end
end

%%
%% Setup the temporal signals for the scanning volume
%  list these three since the calcium signals are generated for them

% activities of soma
soma_act = single(neur_act.soma);

% activities of dendrites
dend_act = single(neur_act.dend);

bg_act = single(neur_act.bg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalize soma activities
somaVol_act = [];
window_1 = 20;
window_2 = 200; 
for ll = 1:size(soma_act,1) % for neuron number
    buf_act = soma_act(ll,1 : end);
    if max(buf_act(:)) == 0
        somaVol_act(ll) = 0;
    else
        somaVol_act(ll) = calcu_std_and_smooth(buf_act, window_1, window_2, 1);
    end
end

median_val = median(somaVol_act) ;
somaVol_act(somaVol_act > median_val) = median_val;
somaVol_act(somaVol_act < median_val * 0.1) = 0;

% find amplification such that all neurons are uniformly firing
for ll = 1:size(soma_act,1) % for neuron number
    buf_act = soma_act(ll,1 : end);
    buf2 = calcu_std_and_smooth(buf_act, window_1, window_2, 1);
    if buf2 == 0
        soma_act(ll,:) = 0;
    else
        soma_act(ll,:) = soma_act(ll,:) / buf2 *  somaVol_act(ll);
    end
end

% disable apical dendrites

%% reduce the variance of background fluctuations
% min activity
% 
% cut off the too-high baseline

bg_act_max = max(bg_act,[],2);
soma_max = max(soma_act,[],2);
dend_max = max(dend_act,[],2);
% median_min_bg995 = prctile(bg_act_min, 99.5);
median_soma = median(soma_max );
median_bg= median(bg_act_max);


% background dendirites
% for ll = 1:size(bg_act,1) % for neuron number
%     buf = bg_act(ll,:) +median_bg * 0.2;
%     if max(buf) > 0
%         bg_act(ll,:) = buf / max(buf(:)) * bg_act_max(ll);
%     end
% end


% % dendrites
for ll = 1:size(dend_act,1) % for neuron number
    buf = dend_act(ll,:) + max(dend_max);
    if max(buf) > 0
        dend_act(ll,:) = buf / max(buf(:)) * median_soma * 1e-2;
    end
end

% disable apical dendrites
dend_act(vol_params.N_neur + 1 : vol_params.N_neur + vol_params.N_den, :) = 0;
% statistics
max(soma_act(:))
max(bg_act(:))
max(dend_act(:))

max(std(soma_act, 0, 2))
max(std(dend_act, 0, 2))
max(std(bg_act, 0, 2))

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preset a base volume that is motified with the activity of each frame.. 
% This saves a LOT of computational time 

fprintf('Initializing the base volume...')
soma_min = min(soma_act,[],2);
soma_act = bsxfun(@minus,soma_act,soma_min);

dend_min = min(dend_act,[],2);
dend_act = bsxfun(@minus,dend_act,dend_min);

bg_min = min(bg_act,[],2);
bg_act = bsxfun(@minus,bg_act,bg_min);

f0vol = zeros(N1,N2,N3,'single');
f0vol_wo_bg = zeros(N1,N2,N3,'single');
for ll = 1:size(soma_min,1)
    array_SubSubTest(f0vol,somaVol{ll,1},somaVol{ll,2},soma_min(ll));      % Iteratively add in each neuron's soma activity
end
for ll = 1:size(soma_min,1)
    array_SubSubTest(f0vol_wo_bg,somaVol{ll,1},somaVol{ll,2},soma_min(ll));      % Iteratively add in each neuron's soma activity
end

for ll = 1:size(neur_vol.gp_nuc,1)
    if (numel(neur_vol.gp_nuc{ll,2}) == 1)&&(neur_vol.gp_nuc{ll,2}(1) == 0)% Nothing to do
    elseif numel(neur_vol.gp_nuc{ll,2}) == numel(neur_vol.gp_nuc{ll,1})
        array_SubSubTest(f0vol,(neur_vol.gp_nuc{ll,1}), ...
            single(neur_vol.gp_nuc{ll,2}),single(soma_min(ll)));           % Iteratively add in each neuron's nucleus activity
    elseif (numel(neur_vol.gp_nuc{ll,2})==1)&&(neur_vol.gp_nuc{ll,2}(1)~=0)
        f0vol(neur_vol.gp_nuc{ll,1}) = neur_vol.gp_nuc{ll,2}*soma_min(ll); % Iteratively add in each neuron's nucleus activity
    end
end
for ll = 1:size(dend_min,1)                                                % For large volumes, mex functions decrease runtimes a lot
    array_SubSubTest(f0vol,dendVol{ll,1},dendVol{ll,2},dend_min(ll));      % Iteratively add in each neuron's soma activity
end
for ll = 1:size(bg_min,1)
    if(~isempty(axonVol{ll,1}))
        array_SubModTest(f0vol, (axonVol{ll,1}),axonVol{ll,2},bg_min(ll)); % Iteratively add in each background component's activity
    end
end
if(isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1)
    nuc_act = single(neur_act.nuc);
    nuc_min = min(nuc_act,[],2);
    nuc_act = bsxfun(@minus,nuc_act,bg_min);
    nuc_act(nuc_act<cutoff) = 0;  
    for ll = 1:size(nucVol,1)
      array_SubSubTest(f0vol,nucVol{ll,1},nucVol{ll,2},nuc_min(ll));       % Iteratively add in each neuron's soma activity
    end
end
fprintf('base initialized.\n')

img_w0 = single_scan(f0vol(:,:,z_loc:(z_loc+Np3-1)),...
                                     size(PSF), freq_psf, scan_avg, true); 
img_w0 = imresize(img_w0, [floor((N1)/sfrac * scale_ratio), floor((N2)/sfrac * scale_ratio)]);
        
saveastiff(im2uint16(img_w0 / max(img_w0(:))), sprintf('%s\\f0img.tiff', output_dir))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteratively scan the volume

if (~isempty(scan_params.fsimPath))
  [tifLinkFsim,tagFsim] = tifinitialize(scan_params.fsimPath, ...
              floor([(N1)/sfrac, (N2)/sfrac,]) - 2*scan_buff/sfrac);
end
if (~isempty(scan_params.fsimCleanPath))
  [tifLinkFsimClean,tagFsimClean] = tifinitialize(scan_params.fsimCleanPath, ...
              floor([(N1)/sfrac, (N2)/sfrac,]) - 2*scan_buff/sfrac);
end

if scan_params.verbose == 1
    fprintf('Scanning...')
elseif scan_params.verbose > 1
    fprintf('Scanning...\n')
end

%% main scan loop
if ~isstruct(PSF)
    for kk = start_Nt + 1 : Nt  % Iterate over time steps
        if scan_params.verbose >1
            tic
        end
        if rand(1) > p_jump                                                % Randomly choose whether to jump
            x_loc = min(max(1,x_loc+randsample(d_stps2,1)), 2*scan_buff+1);% Update random jump for x-axis (slow) scan start 
            y_loc = min(max(1,y_loc+randsample(d_stps2,1)), 2*scan_buff+1);% Update random jump for y-axis (fast) scan start
        end
        %% initilaize motions
        x_pos = min(max(1,x_loc+randsample(d_stps, 1)), 2*scan_buff+1);    % Update random drift for x-axis (slow) scan start
        y_pos = min(max(1,y_loc+randsample(d_stps, 1)), 2*scan_buff+1);    % Update random drift for y-axis (fast) scan start
        z_loc = min(max(z_base-zmaxdiff,z_loc+...
                       randsample(d_stpsZ, 1)), z_base+zmaxdiff);          % Update random drift for depth scan start
        z_loc = min(max(1,z_loc+randsample(d_stps, 1)), N3-Np3+1);         % Update random drift for depth scan start
        if nargout > 2
            mot_hist(:,kk) = [x_pos;y_pos;z_loc];                          % If needed store the current scan location to the output array
        end
        y_shr = [zeros(randsample(floor(2*N1/5),1),1); ...
                linspace(0,1,round(rand*3*N1/5))'*(2*(rand-0.5))*maxshear*N1];    % Create shearing vector
        y_shr = cat(1,y_shr, y_shr(end)*ones(max(0,N1-numel(y_shr)),1));   % Make sure vector is the full size (N1)
        y_off = vec(min(max(1,y_pos+y_shr + ...
                        randsample(d_stps, N1, true)'),2*scan_buff+1));    % Randomly select offset for each row (+/-0.5um), in addition to the shearing vector
        clear TMPvol; 
        %% add volumes
        WMPvol_w_bg = zeros(N1,N2,N3,'single');                                 % Reset the temporary volume - It's an order of magnitde faster to scrap and start over                                                
        WMPvol_w_bg_pure = zeros(N1,N2,N3,'single');
        
        for ll = 1:size(soma_act,1)                                        % For large volumes, mex functions decrease runtimes a lot
            if(soma_act(ll,kk)>0)&&(~isempty(somaVol{ll,2}))
                array_SubSubTest(WMPvol_w_bg,somaVol{ll,1},...
                                        somaVol{ll,2}, exp_level * soma_act(ll,kk));    % Iteratively add in each neuron's soma activity
            end
        end 
        clean_img_wo_bg = single_scan(WMPvol_w_bg(:,:,...
                                        z_loc:(z_loc+Np3-1)) + f0vol_wo_bg(:,:,z_loc:(z_loc+Np3-1)), ...
                                        size(PSF), freq_psf, scan_avg, true); 
        
        
        % nuc
        if(isfield(scan_params,'nuc_label'))&&(scan_params.nuc_label>=1)
            for ll = 1:size(neur_vol.gp_nuc,1)
                if(nuc_act(ll,kk)>0)&&(~isempty(nucVol{ll,2}))
                    array_SubSubTest(WMPvol_w_bg,nucVol{ll,1},...
                                        nucVol{ll,2},nuc_act(ll,kk));    % Iteratively add in each neuron's soma activity
                    array_SubSubTest(WMPvol_w_bg_pure,nucVol{ll,1},...
                                        nucVol{ll,2},nuc_act(ll,kk));    % Iteratively add in each neuron's soma activity                
                end            
            end
        end
        
%         % dend
        for ll = 1:size(dend_act,1)                                        % For large volumes, mex functions decrease runtimes a lot
            if (dend_act(ll,kk)>0)&&(~isempty(dendVol{ll,1}))
                array_SubSubTest(WMPvol_w_bg,dendVol{ll,1},...
                                           dendVol{ll,2},dend_act(ll,kk)); % Iteratively add in each neuron's soma activity
                array_SubSubTest(WMPvol_w_bg_pure,dendVol{ll,1},...
                                           dendVol{ll,2},dend_act(ll,kk)); % Iteratively add in each neuron's soma activity
            end
        end
        
        % bg
        for ll = 1:size(bg_act,1)
            if(~isempty(axonVol{ll,1}) && bg_act(ll,kk)>0)
                array_SubModTest(WMPvol_w_bg_pure, axonVol{ll,1}, ...
                                           axonVol{ll,2},bg_act(ll,kk )); % Iteratively add in each background component's activity
                array_SubSubTest(WMPvol_w_bg_pure,dendVol{ll,1},...
                                           dendVol{ll,2},dend_act(ll,kk)); % Iteratively add in each neuron's soma activity
            end
        end
        
        %% get clean image
        clean_img_w_bg = (sigscale/(2*sfrac^2))*single_scan(WMPvol_w_bg(:,:,...
                  z_loc:(z_loc+Np3-1))+f0vol(:,:,z_loc:(z_loc+Np3-1)),...
                                     size(PSF), freq_psf, scan_avg, true); % Scan a single frame, normalization for pixel size and scanning code
        
        %% contaminations of clean image
        % further blur ness
        if(~isempty(g_blur))                                               % Keep clean_img with no additive blur
            clean_img_w_bg = clean_img_w_bg + conv2(clean_img_w_bg,g_blur,'same');        % Add in additional blurring for dimmer areas
        end
        
        % top and bottom masks
        if(~isempty(psfT))
            if(~isempty(t_mask))
                top_mask = 1./t_mask;
                bot_mask = 1./t_mask;
            else
                top_mask = ones(size(clean_img_w_bg,1),size(clean_img_w_bg,2),'single');
                bot_mask = ones(size(clean_img_w_bg,1),size(clean_img_w_bg,2),'single');            
            end
            if(isfield(psfT,'mask'))
                top_mask = top_mask.*psfT.mask;
                bot_mask = bot_mask.*psfB.mask;
            end
            % top image
            if(isempty(1:z_loc-1))
                % so the top img is still based on the whole volume
                top_img = blurredBackComp2(WMPvol_w_bg_pure,1:size(WMPvol_w_bg_pure,3),...
                           psfT.freq_psf,psfT.weight, top_mask, 1, [], f0vol);          
            else
                top_img = blurredBackComp2(WMPvol_w_bg_pure,1:z_loc-1,psfT.freq_psf,...
                                         psfT.weight, top_mask, 1, psfT.psfZ, f0vol);
            end
            if(isempty(z_loc+Np3:size(WMPvol_w_bg_pure,3)))
                bot_img = blurredBackComp2(WMPvol_w_bg_pure,1:size(WMPvol_w_bg_pure,3),...
                           psfB.freq_psf,psfB.weight, bot_mask, 1, [], f0vol);                      
            else
                bot_img = blurredBackComp2(WMPvol_w_bg_pure,z_loc+Np3:size(WMPvol_w_bg_pure,3),...
                           psfB.freq_psf,psfB.weight, bot_mask, 1, psfB.psfZ, f0vol);            
            end
            top_img = top_img*(sigscale/(sfrac^2));
            bot_img = bot_img*(sigscale/(sfrac^2));
               
            clean_img_w_bg = clean_img_w_bg + bot_img * 0.1+ top_img * 0.1;

        end
        
        %% motions 
        if scan_params.motion % typically without motion
            clean_img_wo_bg = imgSubRowShift(clean_img_wo_bg, scan_buff, x_pos, ...
                                                            round(y_off)); % Sub-select a portion of an image with row shifting
        end
        clean_img_wo_bg = imresize(clean_img_wo_bg, [floor((N1)/sfrac * scale_ratio), floor((N2)/sfrac * scale_ratio)]);
        
        if scan_params.motion
            clean_img_w_bg = imgSubRowShift(clean_img_w_bg, scan_buff, x_pos, ...
                                                            round(y_off)); % Sub-select a portion of an image with row shifting
        end
        clean_img_w_bg = imresize(clean_img_w_bg , [floor((N1)/sfrac * scale_ratio), floor((N2)/sfrac * scale_ratio)]);
        
        mov_wo_bg(:,:,kk) = clean_img_wo_bg;
        %% add noise
        
        samp_img  = PoissonGaussNoiseModel_camera(clean_img_w_bg, noise_params);        % Run data through the noise model        
        mov_w_bg(:,:,kk) = samp_img;  
        
        % debug mode
        saveastiff(im2uint16(samp_img / 65535), sprintf('%s\\debug\\%d.tiff', output_dir, kk - start_Nt))
        saveastiff(im2uint16(clean_img_wo_bg / 65535), sprintf('%s\\debug_gt\\%d.tiff', output_dir, kk - start_Nt))
        %% saving tiff, optional for large files
        if (~isempty(scan_params.fsimPath))
            tifLinkFsim = tifappend(tifLinkFsim,samp_img,tagFsim, ...
              scan_params.saveBlocksize,kk,scan_params.fsimPath);
        end
        if (~isempty(scan_params.fsimCleanPath))
            tifLinkFsimClean = tifappend(tifLinkFsimClean,clean_img_w_bg, ...
              tagFsimClean,scan_params.saveBlocksize,kk,scan_params.fsimCleanPath);
        end        
        if scan_params.verbose == 1
            fprintf('.');
        elseif scan_params.verbose >1
            Ttmp = toc;
            fprintf('Scanned frame %d (%f s)\n',kk - start_Nt,Ttmp)
        end
    end
elseif isfield(PSF,'left')                                                 % Check if separate beam v_twins configuration is being used
    error('Supported PSF configuration!')    
else
    error('Unknown input configuration!')
end

if (~isempty(scan_params.fsimPath))
  tifLinkFsim.close();
end
if (~isempty(scan_params.fsimCleanPath))
  tifLinkFsimClean.close();
end


if scan_params.verbose >= 1
    fprintf('done.\n')
end

% save
mov_w_bg = mov_w_bg(:, :, start_Nt + 1 : end);
mov_wo_bg = mov_wo_bg(:, :, start_Nt + 1: end);
saveastiff(im2uint16(mov_w_bg / max(mov_w_bg(:))), sprintf('%s\\mov_w_bg.tiff', output_dir))
saveastiff(im2uint16(mov_wo_bg / max(mov_wo_bg(:))), sprintf('%s\\mov_wo_bg.tiff', output_dir))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Output parsing

if nargout > 0
    varargout{1} = mov_w_bg;                                               % Output full movie
end
if nargout > 1
    varargout{2} = mot_hist;                                               % Output history of motion (x/y/z) (if requested)
end




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% utility function
function activity = calcu_std_and_smooth(input, window_1, window_2, poly_index)
    input = input(:);
    smoothed_value = smoothdata(input, 1,'movmean', window_1);

    detrend_data = input - smoothed_value;
    smoothed_value2 = smoothdata(detrend_data, 1,'movmean', window_2);
    tmp = smoothed_value2.^poly_index;
    activity =std(tmp);
end