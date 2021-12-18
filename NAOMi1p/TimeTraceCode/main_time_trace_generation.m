clc, clear
close all
%% this file is to allocate the temporal signals to the generated sample.
%  note there is no noise contaminations here.
%  all traces looks perfectly great.

% note the trace is highly related to the neuron volume, so I will send it
% back to the volume dir.
% last update: 8/20/2020. YZ
%% spike parameters
% choose which volume you would like to use
input_dir = '..\\Volume_simulation_all\\volume_size_600_600_200_depth_0_res_1_neuron_num_720_1';
input_file = 'volume_output_size_600_600_200_depth_0_res_1_neuron_num_720.mat';
all_output = importdata(sprintf('%s\\%s', input_dir, input_file));



n_locs = all_output.vol_out.locs;
spike_opts.K      = size(all_output.vol_out.gp_vals,1);                               % Read off the number of neurons
spike_opts.rate   = 1e-3;   
% The average rate of firing for all the components can be modulated using
% this parameter. default is 1e-3
spike_opts.mu = 0;  % Default the mean of the normal r.v. used in the log-normal to zero
spike_opts.sig = 1; % Default the standard-deviation of the normal r.v. used in the log-normal to one
spike_opts.dyn_type = 'Ca_DE'; % Default to a single-compartment Ca2+ model with double-exponential smoothing, double exponential

spike_opts.rate_dist = 'gamma';  % Default to a gamma distribution of firing rates
spike_opts.dt = 1/6;   % Default sampling rate is 30 Hz
spike_opts.nt = 1000; % Default number of time-steps to simulate is 1000
% spike_opts.rate = 1e-3; % Default inverse average of 1s between bursts (1000 1/100 s units) WAS: spike_opts.rate = 0.16;
spike_opts.N_bg = 0; % Default to only one background component
spike_opts.prot = 'GCaMP6'; 
spike_opts.alpha = 1; % Default the Gamma distribution to an Exponential distribution (alpha=1 for exponential distribuion)
spike_opts.burst_mean = 10; % Default the mean of the Poisson distribution for number of spikes per burst to be 4

% spike_opts.smod_flag = 'hawkes'; % Default simulation model is the Hawke's model
spike_opts.smod_flag = 'other';
% spike_opts.p_off = 0.2;
spike_opts.p_off = -1;
spike_opts.selfact = 1.2; 
spike_opts.min_mod = [0.4 2.53];
spike_opts.spikeflag = 1;  
spike_opts.dendflag = 1; % Flag to simulate dendrite traces per component
spike_opts.axonflag = 1;  % Flag to simulate axon traces per component
                          
outdir = sprintf('%s\\firing_rate_%g_smod_flag_%s', input_dir, spike_opts.rate, spike_opts.smod_flag);
mkdir(outdir)                        
% Check time-trace generation options
if (isfield(spike_opts,'spikeflag')) && spike_opts.spikeflag               % Check whether or not to save the spiking process itself
    spikeflag = 1;                                                         % Default to saving it
else
    spikeflag = 0;                                                         % Alternative is not to save it 
end

if (isfield(spike_opts,'dendflag')) && ~spike_opts.dendflag                % Check the dendritic simulation flag (might not be set for nuclear imaging)
    dendflag = 0;                                                          % If not set, don't do extra work
else
    dendflag = 1;                                                          % If set, also simulate dendrites
end

if (isfield(spike_opts,'axonflag')) && ~spike_opts.axonflag                % Check the axonal simulation flag (might not be set for nuclear imaging)
    axonflag = 0;                                                          % If not set, don't do extra work
else
    axonflag = 1;                                                          % If set, also simulate axons
end

if (spike_opts.N_bg>0)&&(axonflag)                                         % Check if both kind of background is set (axonal neuropil models or the older Gaussian Process neuropil models)
    error('background must be either axons or GPs!')                       % Error out if both are set as it is not clear what to do
end

spike_opts.nt = ceil(spike_opts.nt*100*spike_opts.dt);                     % Spiking is simulated at 100Hz, so make the time traces longer so that later downsampling gives the desired # outputs
prot          = spike_opts.prot;                                           % Extract the type of protein to simulate
l_buff        = 500;                                                       % Create a buffer to allow differential equations to hit steady-state



if (spike_opts.N_bg>0)&&(axonflag)                                         % Check if both kind of background is set (axonal neuropil models or the older Gaussian Process neuropil models)
    error('background must be either axons or GPs!')                       % Error out if both are set as it is not clear what to do
end

%% not sure what is the S_times
S_times = [];
if(iscell(S_times))
    S_times = single(spikesCellToVec(S_times))*7.6e-6;
end

if (~isempty(S_times))&&((size(S_times,1) ~= spike_opts.K)||...
                                   (size(S_times,2) ~= spike_opts.nt)) % Check if there is a correctly sized pre-determined set of spiking activity for the cells
    warning('Incompatible size: Make sure number of components and time traces length is OK')
end
if (~isempty(S_times))
    if(size(S_times,1) ~= spike_opts.K)
        spike_opts.K = size(S_times,1);
    end
    if(size(S_times,2) ~= spike_opts.nt)
        spike_opts.nt = size(S_times,2);
    end
    n_desired = floor(spike_opts.nt/(100*spike_opts.dt));              % Store the official value for later use
else
    n_desired = spike_opts.nt;                                         % Store the official value for later use
end

%% other parameters
disc_opt = true; % ????
mod_vals = []; % ????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate spike times
if isempty(S_times)
    switch spike_opts.smod_flag
        case 'hawkes'
            spike_opts2    = spike_opts;                                   % Copy over spike opts
            spike_opts2.dt = 1/100;                                        % reset the dt to be at 1ms sampling rate
            S_cell    = genCorrelatedSpikeTrains2(spike_opts2,[],n_locs,...
                                                          disc_opt, true); % Samples spikes from a Hawkes process (samples somas and background together)
            S_times   = single(S_cell.soma);                               % Returns a KxN matrix of spike times
            S_cell    = S_cell.bg;
            S_cell.bg = S_cell;
        otherwise
            S_times = gen_burst_spike_times(spike_opts);                   % Returns a KxN matrix of spike times
    end
    if strcmp(spike_opts.dyn_type,'AR1')||...
                                         strcmp(spike_opts.dyn_type,'AR2') % If using the calcium dynamics differential equations...
        S_times(S_times==1) = (1+rand)*exp(spike_opts.mu +...
                        spike_opts.sig*randn(size(S_times(S_times==1))));  % Modulates the height of the KxN matrix of spike times
    elseif strcmp(spike_opts.dyn_type,'single')||...
                                strcmp(spike_opts.dyn_type,'Ca_DE')||...
                                      strcmp(spike_opts.dyn_type,'double') % If using the calcium dynamics differential equations...
        S_times       = [zeros(spike_opts.K,l_buff), S_times];             % ...buffer time traces to let the dynamics reach steady-state before spiking starts
        S_times       = (7.6e-6)*(S_times);                                % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
    end
else 
    S_times = [zeros(spike_opts.K,l_buff), S_times];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate underlying calcium fluorescence levels from the spikes
switch spike_opts.dyn_type
    case 'AR1'
        % Generate impulse response
        h_ca    = make_calcium_impulse(0.9,1);                             % Returns an AR-1 impulese response to convolve with the spike times
        h_ca    = 0.5*h_ca/max(h_ca);                                      % Appropriately normalize the impulse response
        b_cell  = abs(1+0.1*randn(spike_opts.K,1));                        % Create baseline values
        S_somas = bsxfun(@plus, bsxfun(@times, 2.5*convn(S_times, h_ca,...
                           'same'), b_cell*(0.5+0.5*rand(1))), b_cell(:)); % Convolve and add baselines
    case 'AR2'
        % Generate impulse response
        h_ca    = make_calcium_impulse(0.9,[1,1]);                         % Returns an AR-2 impulese response to convolve with the spike times
        h_ca    = 0.5*h_ca/max(h_ca);                                      % Appropriately normalize the impulse response
        b_cell  = abs(1+0.1*randn(spike_opts.K,1));                        % Create baseline values
        S_somas = bsxfun(@plus, bsxfun(@times, 2.5*convn(S_times, h_ca,...
                           'same'), b_cell*(0.5+0.5*rand(1))), b_cell(:)); % Convolve and add baselines
    case 'single'
        cal_params.sat_type = 'single';                                    % Set the calcium dynamics to 'single' mode (simpler dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 800;
        S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));        % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
        [~,~,S_somas] = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time               
        S_somas       = S_somas(:,l_buff+1:end);                           % Remove buffering period
    case 'Ca_DE'
        cal_params.sat_type = 'Ca_DE';                                     % Set the calcium dynamics to 'single' mode (simpler dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 292.3;                                       % Somas needs an even higher extrusion rate
        cal_params.t_on     = 0.8535;
        cal_params.t_off    = 98.6173;
        cal_params.ca_amp   = 76.1251;                                     % Somas needs an even higher extrusion rate
        [~,~,S_somas]  = calcium_dynamics(S_times, cal_params, prot);      % Simulate the dynamics and recover the fluoresence over time               
        S_somas        = S_somas(:,l_buff+1:end);                          % Remove buffering period
    case 'double'
        warning('The "double" method is not fully vetted: use at your own risk.')
        cal_params.sat_type = 'double';                                    % Set the calcium dynamics to 'double' mode (complex dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 800;
        S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));        % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
        
        [~,~,S_somas] = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time
        S_somas       = S_somas(:,l_buff+1:end);                           % Remove buffering period
    otherwise
        error('Unknown dynamics type to simulate!')
end
S.soma = single(S_somas);                                                  % Store the soma activity
clear S_somas                                                              % Free some memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate a dendrite-specific time-traces for each cell
if(dendflag)
switch spike_opts.dyn_type
    case 'AR1'
        % Generate impulse response
        h_ca   = make_calcium_impulse(0.8,1);                              % Returns an AR-1 impulese response to convolve with the spike times
        h_ca   = 0.5*h_ca/max(h_ca);                                       % Appropriately normalize the impulse response
        b_cell = abs(1+0.1*randn(spike_opts.K,1));                         % Create baseline values
        S_dend = bsxfun(@plus, bsxfun(@times, 2.5*convn(S_times, ...
                     h_ca, 'same'), b_cell*(0.5+0.5*rand(1))), b_cell(:)); % Convolve and add baselines
    case 'AR2'
        % Generate impulse response
        h_ca   = make_calcium_impulse(0.8,[1,1]);                          % Returns an AR-2 impulese response to convolve with the spike times
        h_ca   = 0.5*h_ca/max(h_ca);                                       % Appropriately normalize the impulse response
        b_cell = abs(1+0.1*randn(spike_opts.K,1));                         % Create baseline values
        S_dend = bsxfun(@plus, bsxfun(@times, 2.5*convn(S_times, ...
                     h_ca, 'same'), b_cell*(0.5+0.5*rand(1))), b_cell(:)); % Convolve and add baselines
    case 'single'
        cal_params.sat_type = 'single';                                    % Set the calcium dynamics to 'single' mode (simpler dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 400;
        S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));        % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
        [~,~,S_dend]  = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time               
        S_dend        = S_dend(:,l_buff+1:end);                            % Remove buffering period\
    case 'Ca_DE'
        cal_params.sat_type = 'Ca_DE';                                     % Set the calcium dynamics to 'single' mode (simpler dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 292.3/2;                                     %
        cal_params.t_on     = 0.8535;
        cal_params.t_off    = 98.6173;
        cal_params.ca_amp   = 76.1251;     
        [~,~,S_dend]  = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time               
        S_dend        = S_dend(:,l_buff+1:end);                            % Remove buffering period
    case 'double'
        warning('The "double" method is not fully vetted: use at your own risk.')
        cal_params.sat_type = 'double';                                    % Set the calcium dynamics to 'double' mode (complex dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 1000;
        S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));        % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
        [~,~,S_dend]  = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time
        S_dend        = S_dend(:,l_buff+1:end);                            % Remove buffering period
    otherwise
        error('Unknown dynamics type to simulate!')
end
else
  S_dend = [];
end
S.dend = single(S_dend);                                                   % Store the dendrite activity
clear S_dend                                                               % Free some memory
if(spikeflag)
    spikes.somas = S_times(:,l_buff+1:end);                                % Save the spike times as a potential output
else
    spikes = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate an axon-specific time-traces for each cell

if(axonflag)
switch spike_opts.dyn_type
    case 'AR1'
        % Generate impulse response
        h_ca   = make_calcium_impulse(0.8,1);                              % Returns an AR-1 impulese response to convolve with the spike times
        h_ca   = 0.5*h_ca/max(h_ca);                                       % Appropriately normalize the impulse response
        b_cell = abs(1+0.1*randn(spike_opts.K,1));                         % Create baseline values
        S_axon = bsxfun(@plus, bsxfun(@times, 2.5*convn(S_times, ...
                     h_ca, 'same'), b_cell*(0.5+0.5*rand(1))), b_cell(:)); % Convolve and add baselines
    case 'AR2'
        % Generate impulse response
        h_ca   = make_calcium_impulse(0.8,[1,1]);                          % Returns an AR-2 impulese response to convolve with the spike times
        h_ca   = 0.5*h_ca/max(h_ca);                                       % Appropriately normalize the impulse response
        b_cell = abs(1+0.1*randn(spike_opts.K,1));                         % Create baseline values
        S_axon = bsxfun(@plus, bsxfun(@times, 2.5*convn(S_times, ...
                     h_ca, 'same'), b_cell*(0.5+0.5*rand(1))), b_cell(:)); % Convolve and add baselines
    case 'single'
        cal_params.sat_type = 'single';                                    % Set the calcium dynamics to 'single' mode (simpler dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 400;
        S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));        % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
        [~,~,S_axon]  = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time               
        S_axon        = S_axon(:,l_buff+1:end);                            % Remove buffering period\
    case 'Ca_DE'
        cal_params.sat_type = 'Ca_DE';                                     % Set the calcium dynamics to 'single' mode (simpler dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 292.3/4;                                     % 
        cal_params.t_on     = 0.8535;
        cal_params.t_off    = 98.6173;
        cal_params.ca_amp   = 76.1251;     
        [~,~,S_axon]  = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time               
        S_axon        = S_axon(:,l_buff+1:end);                            % Remove buffering period
    case 'double'
        warning('The "double" method is not fully vetted: use at your own risk.')
        cal_params.sat_type = 'double';                                    % Set the calcium dynamics to 'double' mode (complex dynamics)
        cal_params.dt       = 1/100;                                       % Set the calcium dynamics framerate
        cal_params.ext_rate = 1000;
        S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));        % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
        [~,~,S_axon]  = calcium_dynamics(S_times, cal_params, prot);       % Simulate the dynamics and recover the fluoresence over time
        S_axon        = S_axon(:,l_buff+1:end);                            % Remove buffering period
    otherwise
        error('Unknown dynamics type to simulate!')
end
else
  S_axon = [];
end
S.bg = single(S_axon);                                                   % Store the dendrite activity
clear S_dend                                                               % Free some memory
if(spikeflag)
    spikes.bg = S_times(:,l_buff+1:end);                                % Save the spike times as a potential output
else
    spikes = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate a background time-trace (if requested)

if spike_opts.N_bg>0                                                       % If background activity is requested, calculate it
    bgscale      = 0.5;                                                    % Scale the background to ~ 1/2 power
    opts_bg      = spike_opts;                                             % Transfer all options down to the background time-trace generator
    if opts_bg.N_bg>0
        opts_bg.K = opts_bg.N_bg;                                          % Set number of background components 
    end
    opts_bg.rate = 0.25*ones(opts_bg.N_bg,1);                              % Set "spiking" rate for background/neuropil
    opts_bg.sig  = 0.2;                                                    % Set variance in background/neuropil generation to a static 0.2 (can be made adjustable in later iteration)
    switch spike_opts.smod_flag
        case 'hawkes'
            S_times = S_cell.bg;                                           % Returns a N_bgxN matrix of spike times
        otherwise
            S_times = gen_burst_spike_times(opts_bg);                      % Returns a N_bgxN matrix of spike times 
    end
    if strcmp(opts_bg.dyn_type,'AR1')||strcmp(opts_bg.dyn_type,'AR2')      % If using the calcium dynamics differential equations...
        S_times(S_times==1) = (1+rand)*exp(opts_bg.mu +...
                        opts_bg.sig*randn(size(S_times(S_times==1))));     % Modulates the height of the KxN matrix of spike times
    elseif strcmp(opts_bg.dyn_type,'single')||...
                                   strcmp(opts_bg.dyn_type,'Ca_DE')||...
                                         strcmp(opts_bg.dyn_type,'double') % If using the calcium dynamics differential equations...
        S_times = [zeros(opts_bg.K,l_buff), S_times];                      % ...buffer time traces to let the dynamics reach steady-state before spiking starts
        S_times = (7.6e-6)*(S_times);                                      % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
    end
    switch spike_opts.dyn_type
        case 'AR1'
            % Generate impulse response
            h_ca = make_calcium_impulse(0.8,1);                            % Returns an AR-1 impulese response to convolve with the spike times
            h_ca = 0.5*h_ca/max(h_ca);                                     % Appropriately normalize the impulse response
            b_cell       = abs(1+0.25*randn(spike_opts.N_bg,opts_bg.nt));  % Initialize the array of background/neuropil activity
            S_bg = bgscale*convn(S_times, h_ca,'same')+b_cell;             % Convolve and add baselines
            S_bg = S_bg/mean(S_bg(:));                                     % Normalize the background activity to the mean
        case 'AR2'
            % Generate impulse response
            h_ca = make_calcium_impulse(0.8,[1,1]);                        % Returns an AR-2 impulese response to convolve with the spike times
            h_ca = 0.5*h_ca/max(h_ca);                                     % Appropriately normalize the impulse response
            b_cell       = abs(1+0.25*randn(spike_opts.N_bg,opts_bg.nt));  % Initialize the array of background/neuropil activity
            S_bg = bgscale*convn(S_times, h_ca,'same')+b_cell;             % Convolve and add baselines
            S_bg = S_bg/mean(S_bg(:));                                     % Normalize the background activity to the mean
        case 'single'
            cal_params.sat_type = 'single';                                % Set the calcium dynamics to 'single' mode (simpler dynamics)
            cal_params.dt       = 1/100;                                   % Set the calcium dynamics framerate
            cal_params.ext_rate = 1600;
            S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));    % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
            S_times       = [zeros(spike_opts.N_bg,l_buff), S_times];      % Buffer time traces to let the dynamics reach steady-state before spiking starts
            [~,~,S_bg]    = calcium_dynamics(S_times, cal_params, prot);   % Simulate the dynamics and recover the fluoresence over time               
            S_bg          = S_bg(:,l_buff+1:end);                          % Remove buffering period
        case 'Ca_DE'
            cal_params.sat_type = 'Ca_DE';                                 % Set the calcium dynamics to 'single' mode (simpler dynamics)
            cal_params.dt       = 1/100;                                   % Set the calcium dynamics framerate
            cal_params.ext_rate = 2800;                                    % Background needs an even higher extrusion rate
            cal_params.t_on     = 3.1295;
            cal_params.t_off    = 0.020137 ;
            cal_params.ca_amp   = 130.917;
            [~,~,S_bg]  = calcium_dynamics(S_times, cal_params, prot);     % Simulate the dynamics and recover the fluoresence over time               
            S_bg        = S_bg(:,l_buff+1:end);                            % Remove buffering period
        case 'double'
            warning('The "double" method is not fully vetted: use at your own risk.')
            cal_params.sat_type = 'double';                                % Set the calcium dynamics to 'double' mode (complex dynamics)
            cal_params.dt       = 1/100;                                   % Set the calcium dynamics framerate
            cal_params.ext_rate = 200;
            S_times       = (7.6e-6)*(S_times/min(S_times(S_times>0)));    % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
            S_times       = [zeros(spike_opts.N_bg,l_buff), S_times];      % Buffer time traces to let the dynamics reach steady-state before spiking starts
            [~,~,S_bg]    = calcium_dynamics(S_times, cal_params, prot);   % Simulate the dynamics and recover the fluoresence over time
            S_bg          = S_bg(:,l_buff+1:end);                          % Remove buffering period
        otherwise
            error('Unknown dynamics type to simulate!')
    end
    S.bg = single(S_bg);                                                   % Store the dendrite activity
    clear S_bg                                                             % Free some memory
    if(spikeflag)
        spikes.bg = S_times(:,l_buff+1:end);                               % Save the spike times as a potential output
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do some post processing

if mod(floor(1./spike_opts.dt), 10) == 0                                   % For efficiency, get lowest common denomenator factors
    r1 = 10;                                                               % Set up/down sampling rates for resampling
    r2 = floor(1./spike_opts.dt)/10;
else
    r1 = 100;                                                              % Set up/down sampling rates for resampling
    r2 = floor(1./spike_opts.dt);
end

if spike_opts.dt ~= 1/100                                                  % Only do this if the sampling rate is not 100Hz
    buff = 100;
    buff2  = floor(1./spike_opts.dt);
    
    soma_min = min(S.soma,[],2);    
    S.soma = [S.soma(:,1)*ones(1,buff),S.soma,S.soma(:,end)*ones(1,buff)]; % Buffer the data to avoid ringing from the resampling (for somas)
    S.soma = resample(double(S.soma.'), r2, r1).';                         % Resample time-traces to the correct sampling rate (for somas)
    S.soma = S.soma(:, (buff2+1):(end-buff2));                             % Extract resampled time-traces (for somas)
    S.soma = bsxfun(@max, S.soma, soma_min(:));

    if(~isempty(S.dend))
      dend_min = min(S.dend,[],2);
      S.dend = [S.dend(:,1)*ones(1,buff),S.dend,S.dend(:,end)*ones(1,buff)];% Buffer the data to avoid ringing from the resampling (for dendrites)
      S.dend = resample(double(S.dend.'), r2, r1).';                       % Resample time-traces to the correct sampling rate (for dendrites)
      S.dend = S.dend(:, (buff2+1):(end-buff2));                           % Extract resampled time-traces (for dendrites)
      S.dend = bsxfun(@max, S.dend, dend_min(:));
    end
    
    if (spike_opts.N_bg > 0)||(axonflag)
      bg_min   = min(S.bg,[],2);
      S.bg   = [S.bg(:,1)*ones(1,buff),S.bg,S.bg(:,end)*ones(1,buff)];     % Buffer the data to avoid ringing from the resampling (for background)
      S.bg   = resample(double(S.bg.'), r2, r1).';                         % Resample time-traces to the correct sampling rate (for background)
      S.bg   = S.bg(:, (buff2+1):(end-buff2));                             % Extract resampled time-traces (for background)
      S.bg   = bsxfun(@max, S.bg, bg_min(:));
    end
end

if size(S.soma,2) > n_desired
    S.soma = S.soma(:,1:n_desired);                                        % If necessary, cut the sizes of the time courses down (for somas)
end

if(~isempty(S.dend))
    if size(S.dend,2) > n_desired
        S.dend = S.dend(:,1:n_desired);                                        % If necessary, cut the sizes of the time courses down (for dendrites)
    end
end

if spike_opts.N_bg > 0
    if size(S.bg,2) > n_desired
        S.bg = 0.4*S.bg(:,1:n_desired);                                    % If necessary, cut the sizes of the time courses down (for background)  
    end
end

if(isempty(mod_vals))
    mod_vals = expression_variation(size(S.soma,1), spike_opts.p_off, spike_opts.min_mod);    % Create a set of modulatory factors for each cell
end
S.soma = bsxfun(@times, S.soma, mod_vals);                                 % Modulate the activity of each cell's soma

if(~isempty(S.dend))
    S.dend = bsxfun(@times, S.dend, mod_vals);                             % Modulate the activity of each cell's dendrites
end
if spike_opts.N_bg > 0
    S.bg = expression_variation(S.bg, 0, 0.8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing 

% if nargout == 1
%     varargout{1} = S;                                                      % If only one output: just return the fluorescence traces
% elseif nargout == 2
spikes.bg = sparse(spikes.bg);
spikes.somas = sparse(spikes.somas);

save(sprintf('%s\\S.mat', outdir), 'S')     % return the fluorescence traces
save(sprintf('%s\\spikes.mat', outdir), 'spikes')  % ... and the spiking activity
save(sprintf('%s\\spikes_opts.mat', outdir), 'spike_opts')                                          

figure, imagesc(S.soma(:, : )), saveas(gca, sprintf('%s\\figure of soma', outdir)), title('soma')
figure, imagesc(S.bg(:, : )), saveas(gca, sprintf('%s\\figure of bg', outdir)), title('bg')
figure, imagesc(S.dend(:, : )), saveas(gca, sprintf('%s\\figure of dend', outdir)), title('dendrites')

%% detailed plot
ts = zscore(S.soma(1 : 200, : ), 0, 2);
% ts = zscore(global_T_deconv_filtered( 10500 : 10550, 1 : end) , 0, 2);
for i = 1 : size(ts, 1)
    ts(i, :) = smooth(ts(i, :).', 3);
end
y_shift = 2;
clip = false;
% rng(10021)
sel = 1:size(ts,1);

nixs = 1:size(ts,1);
sel_nixs = nixs(sel);

figure('Position', [100, 100, 600, 800]);
title(['Temporal activity' ' - timeseries, z-scored'], 'Interpreter', 'none');
hold on
for n_ix = 1:numel(sel_nixs)
    ax = gca();
    ax.ColorOrderIndex = 1;
    loop_ts = ts(sel_nixs(n_ix),:);
    if clip
        loop_ts(loop_ts > 3*y_shift) = y_shift;
        loop_ts(loop_ts < -3*y_shift) = -y_shift;
    end
    t = (0:size(ts,2)-1);
    if max(loop_ts) - mean(loop_ts) > 1 * y_shift
        
        plot1 = plot(t, squeeze(loop_ts) + y_shift*(n_ix-1), 'color', 'k');
        plot1.Color(4) = 0.5;
    else
        % with alpha
        plot1 = plot(t, squeeze(loop_ts) + y_shift*(n_ix-1), 'b');
        plot1.Color(4) = 0.3;
    end
    
%     text(30, y_shift*(n_ix-1), num2str(sel_nixs(n_ix)));
end
xlabel('Frame');
xlim([min(t) max(t)]);
axis off
hold off;
axis tight;
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gca, sprintf('%s\\stack_trace_1_100_soma.png', outdir))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
