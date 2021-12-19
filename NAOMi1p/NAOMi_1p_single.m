clc, clear
close all

%% this file is used to generate NAOMi1p data
%  last update: 12/17/2021. YZ

%% add path
installNAOMi1p


%% load pre-defined data
% pre-defined configuration file for RUSH
RUSH_ai148d_config

%% RUSH 148d config
FOV_sz = 600;% FOV, um
nt = 1000;  % frames
fn = 10; % frame rate
pavg = 0.5; %mW per mm^2


%% parser
vol_params.vol_sz(1) = FOV_sz;
vol_params.vol_sz(2) = FOV_sz;
spike_opts.nt = nt + 200;
spike_opts.dt = 1 / fn;
wdm_params.pavg = pavg;


%% output path
output_dir = sprintf('Z:\\YZ_personal_storage\\deep_widefield_calcium_inference\\data2\\%s\\res_%.2f\\vol_%d_%d_NA_%.2f_Hz_%d_exp_%d_d_%dk_pw_%.2f', ...
                                            mode, pixel_size, vol_params.vol_sz(1), vol_params.vol_sz(3), ...
                                            psf_params.objNA,frate, exp_level,vol_params.neur_density / 1e3, ...
                                            wdm_params.pavg);
% make sub folders
buf = true;
id = 1;
while buf
    if exist(sprintf('%s\\%d', output_dir, id), 'dir') == 7
        id = id +1;
    else
        buf = false;
    end 
end
output_dir = sprintf('%s\\%d', output_dir, id);
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
% [neur_act,spikes] = generateTimeTraces(spike_opts,[],vol_out.locs);        % Generate time traces using AR-2 process
fun_time_trace_generation(vol_out, nt, fn, output_dir)
fprintf('Simulated temporal activity in %f seconds.\n', toc); 
%% plot traces
figure('position', [100, 100, 400, 800]), imagesc(neur_act.soma(:, : )),  title('soma'), colormap(othercolor('BuGn7'))
saveas(gcf, sprintf('%s\\soma_heat.jpg', output_dir)), close
figure('position', [100, 100, 400, 800]), imagesc(neur_act.bg(:, : )),  title('bg'), colormap(othercolor('BuGn7'))
saveas(gcf, sprintf('%s\\bg_heat.jpg', output_dir)), close
figure('position', [100, 100, 400, 800]), imagesc(neur_act.dend(:, : )), title('dendrites'), colormap(othercolor('BuGn7'))
saveas(gcf, sprintf('%s\\dend_heat.jpg', output_dir)), close

%% peform imaging    
clc
% 
% vol_out = importdata(sprintf('%s\\vol_out.mat', output_dir));
% PSF_struct = importdata(sprintf('%s\\PSF_struct.mat',output_dir));     
spike_opts = importdata(sprintf('%s\\firing_rate_0.001_smod_flag_other\\spikes_opts.mat', output_dir));
neur_act = importdata(sprintf('%s\\firing_rate_0.001_smod_flag_other\\S.mat', output_dir));

exp_level = 5;
scan_volume_1p(vol_out, PSF_struct, neur_act, ...
                       vol_params, scan_params, noise_params, spike_opts, wdm_params, pixel_size, exp_level, output_dir); % Perform the scanning simulation



