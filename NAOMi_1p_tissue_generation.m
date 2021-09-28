clc, clear
close all

%% this file is used to generate NAOMi1p data
%  last update: 9/27/2021. YZ

%% add path
installNAOMi

%% volume  parameters
vol_params.vol_sz    = [300,300,150];   % Volume size to sample (in microns)
vol_params.vol_depth = 50;  % Set the depth of imaging (depth at the middle of the simulated volume)
vol_params.neur_density = 1e4;
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
exp_level = 3; % control the brightness of neurons 

scan_params.motion = false; %  motion simulation flat=g
scan_params.verbose = 2; % details
% scan_params.fsimPath

% noise_params

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


% output path

output_dir = sprintf('Z:\\YZ_personal_storage\\deep_widefield_calcium_inference\\data2\\%s\\res_%.2f\\vol_%d_%d_NA_%.2f_Hz_%d_exp_%d_d_%dk_pw_%.2f', ...
                                            mode, pixel_size, vol_params.vol_sz(1), vol_params.vol_sz(3), ...
                                            psf_params.objNA,frate, exp_level,vol_params.neur_density / 1e3, ...
                                            wdm_params.pavg);

mkdir(output_dir)
%% generate volume
tic
[vol_out,vol_params,neur_params,vasc_params,dend_params,bg_params, ...
    axon_params] = simulate_neural_volume(vol_params, neur_params, ...
            vasc_params, dend_params, bg_params, axon_params, psf_params); % Draw a random volume - this takes the longest amound of time
        
% assign all targets with fluroescences
spike_opts.K = size(vol_out.gp_vals,1); % Read off the number of neurons              
        
fprintf('Simulated neural volume in %f seconds.\n', toc);
save(sprintf('%s\\vol_out.mat', output_dir), 'vol_out', '-v7.3')
save(sprintf('%s\\prams.mat', output_dir), 'vol_params', 'neur_params', 'vasc_params', 'dend_params', 'bg_params', 'axon_params');
%% save necessary volumes
neur_vol = vol_out.neur_vol;
saveastiff(im2uint16(neur_vol / max(neur_vol(:))), sprintf('%s\\neur_vol.tiff', output_dir));

neur_ves = vol_out.neur_ves_all;
saveastiff(im2uint16(neur_ves / max(neur_ves(:))), sprintf('%s\\neur_ves_all.tiff', output_dir));
%% generate PSFs
tic
PSF_struct = simulate_1p_optical_propagation(vol_params,psf_params,vol_out);  % Create the point-spread function and mask for scanning
fprintf('Simulated optical propagation in %f seconds.\n', toc); 
save(sprintf('%s\\PSF_struct.mat', output_dir), 'PSF_struct')
%% save necessary PSF files
saveastiff(im2uint16(PSF_struct.psf / max(PSF_struct.psf, [], 'all')), sprintf('%s\\psf.tiff', output_dir))
saveastiff(im2uint16(PSF_struct.mask / max(PSF_struct.mask, [], 'all')), sprintf('%s\\psf_mask.tiff', output_dir))
saveastiff(im2uint16(PSF_struct.colmask / max(PSF_struct.colmask, [], 'all')), sprintf('%s\\psf_colmask.tiff', output_dir))
saveastiff(im2uint16(PSF_struct.psfB.mask / max(PSF_struct.psfB.mask, [], 'all')), sprintf('%s\\psfB_mask.tiff', output_dir))
saveastiff(im2uint16(PSF_struct.psfT.mask/ max(PSF_struct.psfT.mask, [], 'all')), sprintf('%s\\psfT_colmask.tiff', output_dir))


%% generate neurons
tic                                           
[neur_act,spikes] = generateTimeTraces(spike_opts,[],vol_out.locs);        % Generate time traces using AR-2 process
fprintf('Simulated temporal activity in %f seconds.\n', toc); 

%% plot traces
figure('position', [100, 100, 400, 800]), imagesc(neur_act.soma(:, : )),  title('soma'), colormap(othercolor('BuGn7'))
saveas(gcf, sprintf('%s\\soma_heat.jpg', output_dir)), close
figure('position', [100, 100, 400, 800]), imagesc(neur_act.bg(:, : )),  title('bg'), colormap(othercolor('BuGn7'))
saveas(gcf, sprintf('%s\\bg_heat.jpg', output_dir)), close
figure('position', [100, 100, 400, 800]), imagesc(neur_act.dend(:, : )), title('dendrites'), colormap(othercolor('BuGn7'))
saveas(gcf, sprintf('%s\\dend_heat.jpg', output_dir)), close

figure('position', [100, 100, 400, 800]), temporal_trace_render(neur_act.soma, 'k'), title('soma')
saveas(gcf, sprintf('%s\\soma.jpg', output_dir)),close
figure('position', [100, 100, 400, 800]), temporal_trace_render(neur_act.bg, 'k'), title('bg')
saveas(gcf, sprintf('%s\\bg.jpg', output_dir)),close
figure('position', [100, 100, 400, 800]), temporal_trace_render(neur_act.dend, 'k'), title('dendrites')
saveas(gcf, sprintf('%s\\dend.jpg', output_dir)),close

%% peform imaging    
clc
scan_volume_1p(vol_out, PSF_struct, neur_act, ...
                       vol_params, scan_params, noise_params, spike_opts, wdm_params, pixel_size, 1, output_dir); % Perform the scanning simulation

%%


