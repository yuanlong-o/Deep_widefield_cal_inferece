function noise_params = check_noise_params(noise_params)

%  modified by YZ for widefield imaging.
%  last update: 12/17/2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(noise_params)
    clear noise_params
    noise_params = struct;
end

if (~isfield(noise_params,'offset'))||isempty(noise_params.offset)                
    noise_params.offset = 100; % Default mean measurement increase per photon is 100
end
if (~isfield(noise_params,'gamma'))||isempty(noise_params.gamma)              
    noise_params.gamma = 2.2; % Default electronics offset is 0
end
if (~isfield(noise_params,'sigma'))||isempty(noise_params.sigma)           
    noise_params.sigma = 2300; % Default variance increase per photon is 2300
end
if (~isfield(noise_params,'sigma_p'))||isempty(noise_params.sigma_p)         
    noise_params.sigma_p = 200; 
end

if (~isfield(noise_params,'sigscale'))||isempty(noise_params.sigscale)    
    noise_params.sigscale = 0.5;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
