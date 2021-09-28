function h_out = make_calcium_impulse(ca_scale,varargin)

% h_out = make_calcium_impulse(ca_scale, dt)
% 
% This function generates a temporal calcium impulse response. It takes the
% time constants provided in ca_scale and creates the characteristic
% polynomial to be used with MATLAB's impulse function, which simulates the
% impulse function that is returned in h_out. The sampling rate dt can be
% input as an optional second input. 
% 
% 2015: Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if max(ca_scale) <= 0
    error('Need positive time-scales for on/off decay functions.')
end

if nargin > 1
    dt = varargin{1};
else
    dt = 1/30;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the impulse response

a     = poly(exp(-ca_scale));                                              % Use the time-scale values to create a list of AR coefficients
sys_a = arima('AR',num2cell(-a(2:end)),'Constant',0);                      % Make an ARIMA system to use with MATLAB's impulse response function
h_out = impulse(sys_a,10/dt);                                              % Get the impulse response for the AR system

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%