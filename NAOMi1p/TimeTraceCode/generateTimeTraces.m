function varargout = generateTimeTraces(spike_opts, varargin)

% S = generateTimeTraces(spike_opts, S_times) 
% 
% Generate time traces for a number of simulated neurons. Each neuron fires
% with as am independent Poisson process with a given rate. The firing
% rates are Gamma-distributed, and the firing strength is log-normal
% distributed. The neural spikes are then either convolved with an AR-1 or
% AR-2 impulse response to model the calcium fluorescence dynamics, or
% directly used in a wither a single- or quad-occupancy model of Ca2+ and
% protein interactions to create more biophysiologically accurate
% fluorescence traces.
% The inputs to this function is:
%   -spike_opts   - Struct containing the main options used for generating
%                   time traces. Includes
%     .nt         - Number of time-steps to simulate (default = 1000)
%     .K          - Number of neurons to generate time traces for (default =
%                   30)
%     .rate       - Spiking rate for each neuron (default = 0.001 t/spike)
%     .dt         - Time difference between time-points (inverse of
%                   frame-rate to generate samples at) (default = 1/30 s) 
%     .mu         - (default = 0)
%     .sig        - (default = 1)
%     .dyn_type   - Type of dynamics to simulate. Can be 'AR1'
%                   (autoregressive model with 1 time lag), 'AR2' (AR model
%                   with 2 time lags), 'single' (single-compartment Ca2+ 
%                   <-> protein dynamics), 'Ca_DE' (single-compartment Ca2+
%                    <-> double-expo protein dynamics) (default = 'Ca_DE')
%     .N_bg       - Number of background/neuropil time-traces to generate
%                   (default = 0)
%     .prot       - Protein selection for Ca2+ dynamics and fluorescence
%                   simulation. Options are 'GCaMP3' and 'GCaMP6' (default
%                   = 'GCaMP6')
%     .alpha      - Gamma distribution parameter that modulates the rate
%                   distribution. Defaults to an to an Exponential
%                   distribution (alpha=1) 
%     .burst_mean - Mean of the Poisson distribution for number of spikes
%                   per burst to be (default = 4)
%   -spike_predef - OPTIONAL second argument that lets a user put in a
%                   Kx(nt) pre-defined spike-train for each neuron. In this
%                   case the neuropil/background components are still 
%                   generated independently.
% 
% The ouputs of this function is
%  - S      - Struct containing the simulated fluorescence activity for
%             each of the K neurons for all nt time-steps
%     .soma - Kx(nt) matrix where each row is the activity of each cell at
%             the soma 
%     .dend - Kx(nt) matrix where each row is the activity of each cell at
%             the dendrites [can be optional]
%     .bg   - N_bgx(nt) matrix where each row is the activity of each 
%             background/neuropil component [can be optional]
% 
% 2016 - Adam Charles & Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

spike_opts    = check_spike_opts(spike_opts);                              % Check time-trace generation options
spike_opts.nt = ceil(spike_opts.nt*100*spike_opts.dt);                     % Spiking is simulated at 100Hz, so make the time traces longer so that later downsampling gives the desired # outputs
prot          = spike_opts.prot;                                           % Extract the type of protein to simulate
l_buff        = 500;                                                       % Create a buffer to allow differential equations to hit steady-state

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

if nargin > 1
    S_times = varargin{1};
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
        n_desired = floor(spike_opts.nt/(100*spike_opts.dt));              % Store the official value for later use
    end
else
    n_desired = spike_opts.nt;                                             % Store the official value for later use
    S_times   = [];
end

if nargin > 2
    n_locs = varargin{2};
else
    n_locs = [];
end

if nargin > 3
    disc_opt = varargin{3};
    if isempty(disc_opt); disc_opt = true; end
else
    disc_opt = true;
end
if isempty(disc_opt); disc_opt = true; end

if nargin > 4
    mod_vals = varargin{4};
else
    mod_vals = [];
end

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
%         cal_params.ext_rate = 292.3;                                       % Somas needs an even higher extrusion rate
        [~,~,S_somas]  = calcium_dynamics(S_times, cal_params, prot, 1, 1);% Simulate the dynamics and recover the fluoresence over time               
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
%         cal_params.ext_rate = 292.3/2;                                     %
        [~,~,S_dend]  = calcium_dynamics(S_times, cal_params, prot, 1,1/2);% Simulate the dynamics and recover the fluoresence over time               
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
%         cal_params.ext_rate = 292.3/4;                                     % 
        [~,~,S_axon]  = calcium_dynamics(S_times, cal_params, prot, 1,1/4);% Simulate the dynamics and recover the fluoresence over time               
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

if size(S.soma,2) > n_desired;     S.soma = S.soma(:,1:n_desired);   end   % If necessary, cut the sizes of the time courses down (for somas)

if(~isempty(S.dend))
    if size(S.dend,2) > n_desired; S.dend = S.dend(:,1:n_desired);   end   % If necessary, cut the sizes of the time courses down (for dendrites)
end

if spike_opts.N_bg > 0
    if size(S.bg,2) > n_desired;   S.bg = 0.4*S.bg(:,1:n_desired);   end   % If necessary, cut the sizes of the time courses down (for background)  
else
    if size(S.bg,2) > n_desired;   S.bg = S.bg(:,1:n_desired);       end   % If necessary, cut the sizes of the time courses down (for background)  
end

if(isempty(mod_vals))
    mod_vals = expression_variation(size(S.soma,1), spike_opts.p_off, ...
                                                      spike_opts.min_mod); % Create a set of modulatory factors for each cell
end
S.soma = bsxfun(@times, S.soma, mod_vals);                                 % Modulate the activity of each cell's soma
if(~isempty(S.dend)); S.dend = bsxfun(@times, S.dend, mod_vals);  end      % Modulate the activity of each cell's dendrites

if ~isempty(S.bg);    S.bg = bsxfun(@times, S.bg, mod_vals);               % Modulate the activity of each cell's dendrites
else;                 S.bg = S.dend;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing 

if nargout == 1
    varargout{1} = S;                                                      % If only one output: just return the fluorescence traces
elseif nargout == 2
    varargout{1} = S;                                                      % If two outputs: return the fluorescence traces...
    varargout{2} = spikes;                                                 % ... and the spiking activity
else
    varargout{1} = S;                                                      % If more than two outputs: return the fluorescence traces...
    varargout{2} = spikes;                                                 % ... and the spiking activity...
    varargout{3} = mod_vals;                                               % ... and the modulation values
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
