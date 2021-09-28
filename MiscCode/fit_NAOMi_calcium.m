function cal_params = fit_NAOMi_calcium(spike_time, t_frame, fmean_norm, t_ephys, varargin)

% function out_val = fit_NAOMi_calcium(spike_time, t_frame, fmean_norm, t_ephys, params, {OptType})
% 
% Function to set up the comparison used for learning the parameters for
% the calcium imaging time-trace generation. Specifically this function
% fits the parameters: ext_rate, t_on, t_off, and ca_amp.
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 4; optType = varargin{1};
else;          optType = 'fmincon';                                        % 'fmincon' or 'BayesOpt': default is 'fmincon'
end

if nargin > 5; nRestart = varargin{2};
else;          nRestart = 1000;                                            % How many times to restart the parameter optimization
end

if nargin > 6; protType = varargin{3};
else;          protType = 'GCaMP6';                                        % What indicator to train for
end

if isempty(nRestart); nRestart = 1000; end

if strcmp(optType, 'BayesOpt')
    ext_rate = optimizableVariable('ext_rate', [100,1500]);
    t_on     = optimizableVariable('t_on',     [1e-3,10]);
    t_off    = optimizableVariable('t_off',    [1e-3,10]);
    ca_amp   = optimizableVariable('ca_amp',   [0,100]);


    cal_fun = @(z) calciumFun(spike_time, t_frame, fmean_norm, ...
                           t_ephys, [z.ext_rate,z.t_on,z.t_off,z.ca_amp]); % Set up function for BayesOpt optimization
    zout    = bayesopt(cal_fun, [ext_rate, t_on, t_off,ca_amp], ...
                   'AcquisitionFunctionName','expected-improvement-plus'); % Run Bayes opt. I've also used " 'UseParallel',true " here to try and speed things up
               
    cal_params.ext_rate = zout.XAtMinObjective.ext_rate;                   % Extract the "ext_rate" parameter from the optimal point
    cal_params.t_on     = zout.XAtMinObjective.t_on;                       % Extract the "t_on" parameter from the optimal point
    cal_params.t_off    = zout.XAtMinObjective.t_off;                      % Extract the "t_off" parameter from the optimal point
    cal_params.ca_amp   = zout.XAtMinObjective.ca_amp;                     % Extract the "ca_amp" parameter from the optimal point
else
    
    cal_fun2 = @(z) calciumFun(spike_time, ...
                         t_frame(1:length(fmean_norm)), fmean_norm,...
                         max(t_frame(1:length(fmean_norm))), z, protType); % Set up function for fmincon optimization
    options = optimoptions('fmincon','Display','off', ...
                                              'OptimalityTolerance',1e-8); % Create options structure for fmincon optimization
    l_bound = [200,0,1e-3,0];                                              % Set appropriate lower bounds for the parameters
    u_bound = [300,1,100,100];                                            % Set appropriate upper bounds for the parameters
    fout    = Inf;                                                         % Initialize cost best function at infinity
    zout   = nan(4,1);                                                     % Initialize best parameters to nans 
    fprintf('.')
    for kk = 1:nRestart      
        [ztmp, ftmp] = fmincon(cal_fun2, ...
                       l_bound+rand(1,4).*(u_bound-l_bound),[],[],[],...
                                           [],l_bound,u_bound,[],options); % Run fmincon to find the best parameters with a random initialization (uniform between lower and upper bound)
        if ftmp < fout
            fout  = ftmp;
            zout = ztmp;
        end
        if mod(kk,50)==0; fprintf('.'); end
    end
    fprintf('\n')
    cal_params.ext_rate = zout(1);                                         % Extract the "ext_rate" parameter from the optimal point
    cal_params.t_on     = zout(2);                                         % Extract the "t_on" parameter from the optimal point
    cal_params.t_off    = zout(3);                                         % Extract the "t_off" parameter from the optimal point
    cal_params.ca_amp   = zout(4);                                         % Extract the "ca_amp" parameter from the optimal point
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions

function out_val = calciumFun(spike_time, t_frame, fmean_norm, t_ephys, params, varargin)

% function out_val = fit_NAOMi_calcium(spike_time, t_frame, fmean_norm, t_ephys, params)
% 
% Function to set up the comparison used for learning the parameters for
% the calcium imaging time-trace generation.
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up parameters

if nargin > 5
    protType = varargin{1};
else
    protType = 'GCaMP6';
end


stimes  = ceil(spike_time*100);                                            % 100Hz estimate of spike times
spikes2 = zeros(1,round(max(t_ephys)));
spikes2(stimes) = 7.6e-6;

spike_opts          = check_spike_opts([]);                                % Set default spiking options
cal_params.sat_type = 'Ca_DE';                                             % Set the calcium dynamics to 'single' mode (simpler dynamics)
cal_params.dt       = 1/100;                                               % Set the sampling frequency

cal_params.ext_rate = params(1);                                           % --------- 
cal_params.t_on     = params(2);                                           % Pass through the parameters to the calcium 
cal_params.t_off    = params(3);                                           %     parameter structure.
cal_params.ca_amp   = params(4);                                           % --------- 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the fluorescence

[~,~,simFluor]  = calcium_dynamics(spikes2, cal_params, protType);         % Simulate the fluorescence time-traces
simFluor        = simFluor/min(simFluor);                                  % Normalize the minimum value to one

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare to the real data

fmean_norm = vec(double(fmean_norm/median(fmean_norm)));                   % Ensure that the comparison data is appropriately normalized and a vector
fmean_norm = [fmean_norm(1);fmean_norm];                                   % Append an anchor point at zero
t_frame    = [0;vec(t_frame)];                                             % Append an anchor point at zero

[cal_out, t_rs] = resample(fmean_norm, t_frame, 100);                      % Resample data to a grid the same sampling frequency as the simulated time-trace
t_max           = min(length(simFluor), numel(t_rs));                      % Find the number of samples to compare (as many samples as are in both vectors)
out_val         = norm(cal_out(1:t_max) - vec(simFluor(1:t_max)));         % Calcualte l2 norm of the error

out_val = double(out_val);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
