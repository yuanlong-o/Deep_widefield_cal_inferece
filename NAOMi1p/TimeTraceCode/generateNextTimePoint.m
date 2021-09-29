function [resampStruct,spikes,fluor] = generateNextTimePoint(spike_opts, cal_params, n_locs, mod_vals, resampStruct)

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


if(~strcmp(spike_opts.smod_flag,'hawkes'))
  warning('Using hawkes process....')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate spike times

spike_opts.dt = 1/100;

if(isempty(resampStruct))
  resampStruct = struct;
end

[resampStruct,spikes]  = genNextSpikeTimepoint(spike_opts, n_locs, resampStruct);        % Samples spikes from a Hawkes process (samples somas and background together)
spikes       = single((7.6e-6)*spikes);                                % Normalize spike amplitudes to be reasonable calcium concentrations (in M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate underlying calcium fluorescence levels from the spikes

if(~strcmp(cal_params.sat_type,'Ca_DE'))
  warning('Using Ca_DE....')
end

cal_params.dt        = 1/100;                                       % Set the calcium dynamics framerate

[resampStruct,fluor]  = genNextCalciumDynamics(spikes, cal_params, resampStruct);      % Simulate the dynamics and recover the fluoresence over time
fluor = fluor.*mod_vals;
