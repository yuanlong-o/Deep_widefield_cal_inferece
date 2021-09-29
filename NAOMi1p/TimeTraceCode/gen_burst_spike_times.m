function S = gen_burst_spike_times(spike_opts)

% S = gen_burst_spike_times(spike_opts)
%
% This function creates a matrix of spike times for K neurons over nt
% time-steps. The input to this function is a struct with the following
% properties: 
% 
%   spike_opts - Struct containing the main options used for generating
%                time traces. Includes
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
%
% The rates for each neuron are calculated by drawing a random Gamma sample
% from a distribution with the given spike_opts.alpha parameter and a beta 
% parameter given by spike_opts.rate. The Gamma distribution is
% 
%      P(x) = (x^(a-1)e^(-x/b))/(Gamma(a)*b^a)
%
% Each neuron then has spikes chosen by running a Poisson process with the
% inverse rate selected. The out put is
% 
% S - spike_opts.K by spike_opts.nt matrix with ones at spike locations. 
%
% 2016 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs

spike_opts = check_spike_opts(spike_opts);                                 % Check the options struct

if numel(spike_opts.rate) == 1
    switch spike_opts.rate_dist
        case 'uniform'        
            spike_opts.rate = spike_opts.rate.*ones(spike_opts.K,1);       % Sample each neuron's firing rate from a Gamma distribution
        case 'gamma'
            rate_tmp = gamrnd(spike_opts.alpha,spike_opts.rate,...
                                                        [spike_opts.K,1]); % Sample each neuron's firing rate from a Gamma distribution
            spike_opts.rate = max(min(rate_tmp, spike_opts.rate*10), ...
                                                      spike_opts.rate/10); % Restrict the range of possible rates
        otherwise
            spike_opts.rate = spike_opts.rate.*ones(spike_opts.K,1);       % Sample each neuron's firing rate from a Gamma distribution
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create spiking times

spike_opts.rate = 1./spike_opts.rate;                                      % Transform rate to average inter-spike time
ref_time = 5;                                                              % Number of samples for a refractory rate. 5ms (base sampling rate is at 1ms)

S = zeros(spike_opts.K, spike_opts.nt);                                    % Initialize spike matrix

for kk = 1:spike_opts.K
    T_tot = 0;                                                             % Initialize the time counter to zero
    while T_tot<spike_opts.nt
        t_arr = -1*spike_opts.rate(kk)*log(rand(1));                       % Take an independent exponential step (def of Poisson process)
        S(kk,min(ceil(T_tot+t_arr),spike_opts.nt+1)) = 1;                  % Set a spike down at that location - If the spike is outside the [1,nt] range, set it to be at nt+1 so that a huge matrix is not generated
        T_tot = T_tot + t_arr;                                             % Keep counting time (increment time counter)
        if spike_opts.burst_mean > 0                                       % Check if bursting behavior is warrented given the mean number of spikes per burst (spike_opts.burst_mean)
            num_in_burst = 1 + poissrnd(spike_opts.burst_mean);            % Get number of spikes in the burst
            if num_in_burst > 1
                for ll = 2:num_in_burst
                    t_arr = ref_time + 2*rand(1);                          % Jitter the actual times of each spike in the burst
                    S(kk,min(ceil(T_tot+t_arr),spike_opts.nt+1)) = 1;      % Set a spike down at that location - If the spike is outside the [1,nt] range, set it to be at nt+1 so that a huge matrix is not generated
                    T_tot = T_tot + t_arr;                                 % Keep counting time (increment time counter)
                end
            end
        end
    end
end
S = S(:,1:spike_opts.nt);                                                  % Take only the specified number of time-points (just in case there is overflow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
