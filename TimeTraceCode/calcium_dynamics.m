function varargout = calcium_dynamics(S, varargin)

% [CB,C,F] = calcium_dynamics(S, cal_params)
%  
% Function to simulate the underlying calcium dynamics with either a
% single-occupancy or quad-occupancy model of protein <=> Ca2+
% interactions. The inputs to this model are:
%  - S           - A Kx(nt) array of the spiking activity at each time-step
%                  and for each neuron 
%  - cal_params  - Struct with parameters for the calcium simulation
%    .ext_rate   - Extrusion rate for the (default = 1800)
%    .ca_bind    - Calcium binding constant (default = 110)
%    .ca_rest    - Resting-state calcium concentration (default = 50e-9)
%    .ind_con    - Indicator concentration (default = 200e-6)
%    .ca_dis     - Calcium disassociation constant (default = 290e-9)
%    .ca_sat     - Optional calcium saturation parameter(default = 1)
%    .sat_type   - Type of dynamics to simulate (default = 'double')                                                                 
%    .dt         - Time-trace sampling rate - should be at least the video
%                  sampling rate (default = 1/30) 
%    .ca_amp     - Calcium transient amplitude (default = 0.09 for GCaMP6;
%                  default = 0.05 for GCaMP3) 
%    .t_on       - Rising time-constant for calcium transient (default =
%                  0.1 for GCaMP6; default = 1 for GCaMP3) 
%    .t_off      - Falling time-constant for calcium transient (default =
%                  1.5 for GCaMP6; default = 1 for GCaMP3) 
%    .a_bind     - Binding rate for more detailed simulation (default =
%                  3.5) 
%    .a_ubind    - Unbinding rate for more detailed simulation (default =
%                  7) 
% 
% The ouputs from this model are
%  - CB - A Kx(nt) array of the bound calcium concentration at each
%         time-step and for each neuron
%  - C  - A Kx(nt) array of the total calcium concentration at each
%         time-step and for each neuron
%  - F  - A Kx(nt) array of the overall fluorescence at each time-step and
%         for each neuron 
% 
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

if nargin < 2;    cal_params = struct;
else;             cal_params = varargin{1};
end

if nargin < 3;    prot_type = 'GCaMP6f';
else;             prot_type = varargin{2};
end

if nargin < 4;     over_samp = 1;
else;              over_samp = varargin{3};
end

if isempty(over_samp); over_samp = 1; end

if nargin < 5;     ext_mult = 1;
else;              ext_mult = varargin{4};
end


cal_params = check_cal_params(cal_params, prot_type);                      % Check param struct and fill in with defaults
% Extract necessary params
ext_rate = ext_mult*cal_params.ext_rate;
ca_bind  = cal_params.ca_bind;
ca_rest  = cal_params.ca_rest;
ind_con  = cal_params.ind_con;
ca_dis   = cal_params.ca_dis;
ca_sat   = cal_params.ca_sat;
sat_type = cal_params.sat_type;
dt       = cal_params.dt;
a_bind   = cal_params.a_bind;
a_ubind  = cal_params.a_ubind;
ca_amp   = cal_params.ca_amp;
t_on     = cal_params.t_on;
t_off    = cal_params.t_off;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the calcium iteratively
if(over_samp>1)
    S = reshape([S;zeros((over_samp-1)*size(S,1), size(S,2))], ...
                                                          size(S,1), []);  % Over-sample the spike trains
end
C      = zeros(size(S,1),size(S,2),'single');                              % Initialize the calcium array to be the same size as the event-array
C(:,1) = max(ca_rest,S(:,1));                                              % Initialize the first column of the calcium array

if strcmp(sat_type,'single')
    oversampFlag = 1;
    a  = a_bind(1); b = a_ubind(1); 
    r  = 0;                                                                % Try to calculate an initial guess for the steady-state rest values
    CB = 0*C;                                                              % Initialize the vlaues of the bound calcium concentrations
    CB(:,1) = min(r(r>=0));                                                % Set the values at zero to the calculated initial value
    
    for kk = 2:size(S,2)
        C(:,kk) = C(:,kk-1) + dt*b*CB(:,kk-1) + (-dt*ext_rate*(C(:,kk-1)...
                     - CB(:,kk-1) - ca_rest) + S(:,kk))./(1 + ca_bind + ...
                                (ind_con*ca_dis)./(C(:,kk-1) + ca_dis).^2);% Update the calcium based on the model dynamics
        if (ca_sat < 1)&&(ca_sat >= 0)
            C(:,kk) = min(C(:,kk), ca_dis*ca_sat./(1-ca_sat));             % Saturate the concentration, if necessary
        end
        CB(:,kk) = CB(:,kk-1) + dt*(-b*CB(:,kk-1) + ...
                         a*(C(:,kk-1)-CB(:,kk-1)).*(ind_con-CB(:,kk-1)));  % Equation for calcium interaction with the protiens
    end
elseif strcmp(sat_type,'Ca_DE')
    oversampFlag = 0;
    a  = a_bind(1)*100*dt; b = a_ubind(1)*100*dt;                          % Decay rates are optimized of 100Hz sampling. the "100*dt' multiplicative factor modifies the decay rates to the actua, desired, sampling 
    for kk = 2:size(S,2)
        C(:,kk) = C(:,kk-1) + (-dt*ext_rate*(C(:,kk-1) ...
                       - ca_rest) + S(:,kk))./(1 + ca_bind + ...
                       (ind_con*ca_dis)./(C(:,kk-1) + ca_dis).^2);         % Update the calcium based on the model dynamics
        if (ca_sat < 1)&&(ca_sat >= 0)
            C(:,kk) = min(C(:,kk), ca_dis*ca_sat./(1-ca_sat));             % Saturate the concentration, if necessary
        end
    end
    clear S
    h_ca = single(mk_doub_exp_ker(t_on, t_off, ca_amp, dt));
    TMP  = convn(C(1,:)-ca_rest, h_ca(:).', 'full')+ca_rest;
    TMP  = TMP(:,1:over_samp:end);
    CB   = zeros(size(C,1),size(TMP,2),'like',C);
    for kk = 1:size(C,1)
      TMP      = convn(C(kk,:)-ca_rest, h_ca(:).', 'full') + ca_rest;      % Apply the impulse response and convolve
      CB(kk,:) = TMP(1:over_samp:end);
    end
    C  = C(:,1:over_samp:end);
    CB = CB(:,1:size(C,2));
elseif strcmp(sat_type,'double')
    oversampFlag = 1;
    a = a_bind; b = a_ubind;                                             
    if numel(a) == 1; a = [a,a]; else; a = a(1:2); end
    if numel(b) == 1; b = [b,b]; else; b = b(1:2); end
    CB1 = 0*C;                                                             % Initialize the vlaues of the bound calcium concentrations: type 1 (fast) binding
    CB2 = 0*C;                                                             % Initialize the vlaues of the bound calcium concentrations: type 2 (slow) binding
    for kk = 2:size(S,2)
        C(:,kk) = C(:,kk-1) + dt*(b(1)*CB1(:,kk-1) + b(2)*CB2(:,kk-1)) ...
            + (-dt*ext_rate*(C(:,kk-1) - CB1(:,kk-1) - CB2(:,kk-1) ...
                       - ca_rest) + S(:,kk))./(1 + ca_bind + ...
                       (ind_con*ca_dis)./(C(:,kk-1) + ca_dis).^2);         % Update the calcium based on the model dynamics
        if (ca_sat < 1)&&(ca_sat >= 0)
            C(:,kk) = min(C(:,kk), ca_dis*ca_sat./(1-ca_sat));             % Saturate the concentration, if necessary
        end
        CB1(:,kk) = CB1(:,kk-1) + dt*(-b(1)*CB1(:,kk-1)+a(1)*(C(:,kk-1)...
                         - CB1(:,kk-1)-CB2(:,kk-1)).*(ind_con...
                                               -CB1(:,kk-1)-CB2(:,kk-1))); % Equation for calcium interaction with the protiens
        CB2(:,kk) = CB2(:,kk-1) + dt*(-b(2)*CB2(:,kk-1)+a(2)*(C(:,kk-1)...
                         - CB1(:,kk-1)-CB2(:,kk-1)).*(ind_con...
                                               -CB1(:,kk-1)-CB2(:,kk-1))); % Equation for calcium interaction with the protiens
    end
    CB = CB1 + CB2;                                                        % Total fluorescing indicator is the sum of both types of bound indicator
else
    error('Unknown model!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if(oversampFlag)
    C  = C(:,1:over_samp:end);
    CB = CB(:,1:over_samp:end);
end

if nargout == 1
    varargout{1} = CB;
elseif nargout == 2
    varargout{1} = CB;
    varargout{2} = C;
elseif nargout == 3
    varargout{1} = CB;
    varargout{2} = C;
    if strcmp(sat_type,'single')
        CB = CB + ca_rest + (b/a)*CB./(ind_con-CB);
    end
    F = sat_nonlin(CB, prot_type);                                         % If requested, also output the fluoresence levels
    varargout{3} = F;
end
end

function F = sat_nonlin(CB, prot_type)
    switch lower(prot_type)                                                % First calculate the dF/F curve (as in Badura et al.)
        case {'gcamp6','gcamp6f'}
            F0 = 1;
            F = 25.2*(1./(1 + (290e-9./CB).^2.7));                         % Hill equation values taken from Badura et al. 2014
        case 'gcamp6s'
            F0 = 1;
%            F = 53.8*(1./(1 + (147e-9./CB).^2.45));                        % Hill equation values taken from Dana et al. 2019
             F = 27.2*(1./(1 + (147e-9./CB).^2.45));                        % Hill equation values taken from Dana et al. 2019
        case 'gcamp3'
            F0 = 2;
            F = 12*(1./(1 + (287e-9./CB).^2.52));                          % Hill equation values taken from Badura et al. 2014
        case {'ogb1','ogb-1'}
            F0 = 1;
            F = 14*(1./(1 + 250e-9./CB));                                  % Hill equation values taken from Badura et al. 2014
        case {'gcamp6-rs09','gcamp6rs09'}
            F0 = 1.4;
            F = 25*(1./(1 + (520e-9./CB).^3.2));                           % Hill equation values taken from Badura et al. 2014
        case {'gcamp6-rs06','gcamp6rs06'}
            F0 = 1.2;
            F = 15*(1./(1 + (320e-9./CB).^3));                             % Hill equation values taken from Badura et al. 2014
        case 'jgcamp7f'
            F0 = 1;
            F = 30.2*(1./(1 + (174e-9./CB).^2.3));                         % Hill equation values taken from Dana et al. 2019
        case 'jgcamp7s'
            F0 = 1;
            F = 40.4*(1./(1 + (68e-9./CB).^2.49));                         % Hill equation values taken from Dana et al. 2019
        case 'jgcamp7b'
            F0 = 1;
            F = 22.1*(1./(1 + (82e-9./CB).^3.06));                         % Hill equation values taken from Dana et al. 2019
        case 'jgcamp7c'
            F0 = 1;
            F = 145.6*(1./(1 + (298e-9./CB).^2.44));                       % Hill equation values taken from Dana et al. 2019
        otherwise 
            warning('Unknown protien type. Defaultin to GCaMP6f...\n')
            F0 = 1;
            F = 25.2*(1./(1 + (290e-9./CB).^2.7));                         % Hill equation values taken from Badura et al. 2014
    end
    F = F0 + F0*F;                                                         % Translate dF/F into fluorescence values by solving (F-F0)/F0 = f([CaB])                                                     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%