function cal_params = check_cal_params(cal_params, prot_type)

% vol_params = check_cal_params(vol_params)
%  
% This function checks the elements of the struct vol_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - cal_params  - Struct with parameters for the calcium simulation
%     .ext_rate   - Extrusion rate for the (default = 265.73)
%     .ca_bind    - Calcium binding constant (default = 110)
%     .ca_rest    - Resting-state calcium concentration (default = 50e-9)
%     .ind_con    - Indicator concentration (default = 200e-6)
%     .ca_dis     - Calcium disassociation constant (default = 290e-9)
%     .ca_sat     - Optional calcium saturation parameter(default = 1)
%     .sat_type   - Type of dynamics to simulate (default = 'double')                                                                 
%     .dt         - Time-trace sampling rate - should be at least the video
%                   sampling rate (default = 1/30) 
%     .ca_amp     - Calcium transient amplitude (default = 130.917 for 
%                   GCaMP6; default = 0.05 for GCaMP3) 
%     .t_on       - Rising time-constant for calcium transient (default =
%                   3.1295 for GCaMP6; default = 1 for GCaMP3) 
%     .t_off      - Falling time-constant for calcium transient (default =
%                   0.020137 for GCaMP6; default = 1 for GCaMP3) 
%     .a_bind     - Binding rate for more detailed simulation (default =
%                   3.5) 
%     .a_ubind    - Unbinding rate for more detailed simulation (default =
%                   7) 
% 
% 2017 - Adam Charles and Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(cal_params)                                                     % Make sure that vol_params is a struct
    clear cal_params
    cal_params = struct;
end

% if (isfield(cal_params,'sat_type'))&&strcmp(cal_params.sat_type,'Ca_DE')
%   cal_params.ext_rate  = 265.73;                                      % Somas needs an even higher extrusion rate
%   cal_params.t_on      = 3.1295;
%   cal_params.t_off     = 0.020137;
%   cal_params.ca_amp    = 130.917;                                       % Somas needs an even higher extrusion rate
% end

if (~isfield(cal_params,'ext_rate'))||isempty(cal_params.ext_rate)         % Default extrusion rate
    cal_params.ext_rate = 265.73;                                          % Default for mouse layer 2/3: from Helmchen & Tank 2011 
end

if (~isfield(cal_params,'ca_bind'))||isempty(cal_params.ca_bind)           % Default calcium binding ratio
    cal_params.ca_bind = 110;                                              % Default for mouse layer 2/3: from Helmchen & Tank 2011
end

if (~isfield(cal_params,'ca_rest'))||isempty(cal_params.ca_rest)           % Default calcium resting level
    cal_params.ca_rest = 50e-9;                                            % ~30-100 n mol: Helmchen et al 1996, Maravall et al. 2000 via Helmchen & Tank 2011
end

if (~isfield(cal_params,'ind_con'))||isempty(cal_params.ind_con)           % Default indicator concentration
    cal_params.ind_con =  200e-6;
end

if (~isfield(cal_params,'ca_dis'))||isempty(cal_params.ca_dis)             % Default calcium dissociation constant
    cal_params.ca_dis = 290e-9;                                            % Default for GCaMP6f: from Badura et al. 2014
end

if (~isfield(cal_params,'ca_sat'))||isempty(cal_params.ca_sat)             % Default to no saturation
    cal_params.ca_sat = 1;
end

if (~isfield(cal_params,'sat_type'))||isempty(cal_params.sat_type)         % Default to a single cacium binding/unbinding equation
    cal_params.sat_type = 'double';
end

if (~isfield(cal_params,'dt'))||isempty(cal_params.dt)                     % Default to a a sampling rate of 30 Hz
    cal_params.dt = 1/30;                                                  % Default is 30Hz
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bind-unbind constants for full simulation

if (~isfield(cal_params,'a_bind'))||isempty(cal_params.a_bind)             % Default binding rate of 3.5
    cal_params.a_bind = 3.5;
end
if (~isfield(cal_params,'a_ubind'))||isempty(cal_params.a_ubind)           % Default unbinding rate of 7
    cal_params.a_ubind = 7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Double exponential parameters for ca + double exp simulation

if (~isfield(cal_params,'ca_amp'))||isempty(cal_params.ca_amp)             % Default Amplitude of 1
    switch prot_type                                                       % First calculate the dF/F curve (as in Badura et al.)
        case {'GCaMP6','gcamp6'}
            cal_params.ca_amp = 130.917;
        case {'GCaMP7','gcamp7'}
            cal_params.ca_amp = 230.917;
        case {'GCaMP3','gcamp3'}
            cal_params.ca_amp = 0.05;
        otherwise
    end
end
if (~isfield(cal_params,'t_on'))||isempty(cal_params.t_on)                 % Default on rate of 0.1
    switch prot_type                                                       % First calculate the dF/F curve (as in Badura et al.)
        case {'GCaMP6','gcamp6'}
            cal_params.t_on = 3.1295;
        case {'GCaMP7','gcamp7'}
            cal_params.t_on = 3.1295;
        case {'GCaMP3','gcamp3'}
            cal_params.t_on = 1;
        otherwise
    end
end
if (~isfield(cal_params,'t_off'))||isempty(cal_params.t_off)               % Default off rate of 1
    switch prot_type                                                       % First calculate the dF/F curve (as in Badura et al.)
        case {'GCaMP6','gcamp6'}
            cal_params.t_off = 0.020137;
        case {'GCaMP7','gcamp7'}
            cal_params.t_off = 0.020137;
        case {'GCaMP3','gcamp3'}
            cal_params.t_off = 1;
        otherwise
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%