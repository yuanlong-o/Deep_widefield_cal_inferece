function [sp,F] = genNextCalciumDynamics(S, cal_params, sp)

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
% 2019 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate the calcium iteratively

% Harded-coded changes reflect changes in fitted parameters - adjust the
% h_ca response to be appropriate for the actual changes in dt
dt        = 3*cal_params.dt;                                              % Adjust the delta-time appropriately
ca_amp    = 1.185*3*cal_params.ca_amp;

if (~isfield(sp,'h_ca'))
  sp.h_ca = single(mk_doub_exp_ker(cal_params.t_on, cal_params.t_off, ca_amp, dt))';
end

if (~isfield(sp,'C'))
  sp.C = zeros(size(S,1),length(sp.h_ca),'single');
  sp.C(:,1) = max(cal_params.ca_rest,S);                                           % Initialize the first column of the calcium array
end

C = sp.C(:,1);
C = C + (-dt*cal_params.ext_rate*(C - cal_params.ca_rest) + S)./(1 + cal_params.ca_bind + ...
  (cal_params.ind_con*cal_params.ca_dis)./(C + cal_params.ca_dis).^2);         % Update the calcium based on the model dynamics

if (cal_params.ca_sat < 1)&&(cal_params.ca_sat >= 0)
  C = min(C, cal_params.ca_dis*cal_params.ca_sat./(1-cal_params.ca_sat));             % Saturate the concentration, if necessary
end

sp.C = [C sp.C(:,1:end-1)];
sp.CB = (sp.C-cal_params.ca_rest)*sp.h_ca+cal_params.ca_rest;                 % Apply the impulse response and convolve
F = sat_nonlin(sp.CB, cal_params.prot_type);                                         % If requested, also output the fluoresence levels

end

function F = sat_nonlin(CB, prot_type)
switch prot_type                                                       % First calculate the dF/F curve (as in Badura et al.)
  case {'GCaMP6','gcamp6'}
    F0 = 1;
    F = 25.2*(1./(1 + (290e-9./CB).^2.7));                         % Hill equation values taken from Badura et al. 2014
  case {'GCaMP3','gcamp3'}
    F0 = 2;
    F = 12*(1./(1 + (287e-9./CB).^2.52));                          % Hill equation values taken from Badura et al. 2014
  case {'ogb1','ogb-1','OGB-1','OGB1'}
    F0 = 1;
    F = 14*(1./(1 + 250e-9./CB));                                  % Hill equation values taken from Badura et al. 2014
  case {'GCaMP6-RS09','GCaMP6RS09','gcamp6-rs09','gcamp6rs09'}
    F0 = 1.4;
    F = 25*(1./(1 + (520e-9./CB).^3.2));                           % Hill equation values taken from Badura et al. 2014
  case {'GCaMP6-RS06','GCaMP6RS06','gcamp6-rs06','gcamp6rs06'}
    F0 = 1.2;
    F = 15*(1./(1 + (320e-9./CB).^3));                             % Hill equation values taken from Badura et al. 2014
  case 'GCaMP7'
    F0 = 1;
    F = 25.2*(1./(1 + (200e-9./CB).^1.9));                         % Hill equation values made up
  otherwise
    warning('Unknown protien type. Defaultin to GCaMP6...\n')
    F0 = 1;
    F = 25.2*(1./(1 + (290e-9./CB).^2.7));                         % Hill equation values taken from Badura et al. 2014
end
F = F0 + F0*F;                                                         % Translate dF/F into fluorescence values by solving (F-F0)/F0 = f([CaB])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%