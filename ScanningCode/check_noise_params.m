function noise_params = check_noise_params(noise_params)

% noise_params = check_noise_params(noise_params)
%  
% This function checks the elements of the struct noise_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - noise_params - Struct contaning the parameters for the noise model
%          .mu       - Mean measurement increase per photon (default = 100)
%          .mu0      - Electronics offset (default = 0)
%          .sigma    - Variance increase per photon (default = 2300)
%          .sigma0   - Electronics base noise variance (default = 2.7)
%          .sigscale - (default = 2.0000e-07)
%          .bleedp   - (default = 0.3000)
%          .bleedw   - (default = 0.4000)
%    
% 2017 - Adam Charles and Alex Song

%  modified by YZ for widefield imaging.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(noise_params)
    clear noise_params
    noise_params = struct;
end

if (~isfield(noise_params,'mu'))||isempty(noise_params.mu)                 % Default mean measurement increase per photon is 100
    noise_params.mu = 100;
end
if (~isfield(noise_params,'mu0'))||isempty(noise_params.mu0)               % Default electronics offset is 0
    noise_params.mu0 = 0;
end
if (~isfield(noise_params,'sigma'))||isempty(noise_params.sigma)           % Default variance increase per photon is 2300
    noise_params.sigma = 2300;
end
if (~isfield(noise_params,'sigma0'))||isempty(noise_params.sigma0)         % Default electronics base noise variance is 2.7
    noise_params.sigma0 = 2.7;
end
if (~isfield(noise_params,'darkcount'))||isempty(noise_params.darkcount)   % Default CMOS
    noise_params.darkcount = 0.05;
end
if (~isfield(noise_params,'sigscale'))||isempty(noise_params.sigscale)     % Default signal magnitude scale is 2e-7
    noise_params.sigscale = 2e-7;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
