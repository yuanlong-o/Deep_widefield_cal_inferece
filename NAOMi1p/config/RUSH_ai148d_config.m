%% RUSH configuration

vol_params.vol_sz    = [300,300,150];   % Volume size to sample (in microns)
vol_params.vol_depth = 50;  % Set the depth of imaging (depth at the middle of the simulated volume)
vol_params.neur_density = 3e4;
pixel_size = 0.8; % your system pixel size
vol_params.vres         = ceil(pixel_size);  % pixel size
neur_params.avg_rad     = 1*5.9;  % for slightly larger neuron


%% widefield system parameters
FN = 18; % FN of objective, Olympus is 18, Nikon is 20, and Zeiss is 16.
M = 10; % system magnification

obj_immersion = 'air'; % or 'water'
psf_params.obj_fl = FN / M;
psf_params.objNA     = 0.3;  % emission NA
psf_params.NA        = 0.3;  % excitation NA
psf_params.lambda = 0.488; % excitation wavelength

if strcmp(obj_immersion, 'water')
    psf_params.zernikeWt  = [0 0 0 0 0 0 0 0 0 0 0]; % system aberrations
    wdm_params.nidx = 1.33;
elseif strcmp(obj_immersion, 'air')
    % add some SA
    psf_params.zernikeWt  = [0 0 0 0.1 0 0 0 0 0 0 0]; % 4th SA. 
    wdm_params.nidx = 1; % collection medium
else
    error('Require water or air as immersion medium!')
end
% z_range = min(vol_params.vol_sz(3), 150); % we maximize the size of 1P for better background simulations
% x_range = round(z_range / 4) - mod(round(z_range / 4), 4);
psf_params.psf_sz = [36, 36, 100];


%% widefield system parameters
spike_opts.prot = 'GCaMP6';
mode = 'w_dend'; % choose the simulation to be with dendrites or noe
if strcmp(mode, 'w_dend')
    spike_opts.dendflag = 1; % if the data is soma confined, set it to be 0;
elseif strcmp(mode, 'wo_dend')
    spike_opts.dendflag = 0; 
end

spike_opts.nt = 1000;                                              % Set number of time step
spike_opts.nt = spike_opts.nt  + 200;
frate = 10;
spike_opts.dt = 1/frate;   % frame rate
spike_opts.rate  = 1e-3; % 0.25 for hawk, and 1e-3 for others
spike_opts.smod_flag = 'hawke'; 

%% widefield system parameters
wdm_params.lambda = 0.532; % emission wavelength
wdm_params.pavg = 0.5;                                                 % power in units of mW, for whole FOV
wdm_params.qe  = 0.7; % sensor QE
exp_level = 1; % control the brightness of neurons 

scan_params.motion = false; %  motion simulation flat=g
scan_params.verbose = 2; % details



%% check those parameters
vol_params   = check_vol_params(vol_params);                               % Check volume parameters
vasc_params  = check_vasc_params([]);                                      % Make default set of vasculature parameters
neur_params  = check_neur_params(neur_params);                                      % Make default set of neuron parameters
dend_params  = check_dend_params([]);                                      % Make default set of dendrite parameters
axon_params  = check_axon_params([]);                                      % Make default set of axon parameters
bg_params    = check_bg_params([]);                                        % Make default set of background parameters
spike_opts   = check_spike_opts(spike_opts);                               % Check spike/fluorescence simulation parameters
noise_params = check_noise_params([]);                                     % Make default noise parameter struct for missing elements
psf_params   = check_psf_params(psf_params);                               % Check point spread function parameters
scan_params  = check_imaging_params(scan_params);                                      % Check the scanning parameter struct
wdm_params   = check_wdm_params(wdm_params);                               % Check the auxiliary two-photon imaging parameter struct 
