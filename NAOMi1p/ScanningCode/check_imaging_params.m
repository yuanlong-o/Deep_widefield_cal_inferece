function scan_params = check_imaging_params(scan_params)

% scan_params = check_scan_params(scan_params)
%  
% This function checks the elements of the struct scan_params to ensure
% that all fields are set. Fields not set are set to a list of default
% parameters. The struct checked is:
% 
%   - scan_params  - Struct contaning the parameters for the scanning
%       .scan_avg  - Sampling rate of the scanning in terms of how many
%                    granular pixels to scan into one pixel (default = 2) 
%       .motion    - True/false option for whether or not to simulate
%                    motion while performing scanning (default = true) 
%       .scan_buff - Number of granular pixels to keep as a buffer from the
%                    edge of the volume (default = 10) 
%       .verbose   - Level of verbosity in the output during the volume
%                    generation. Can be 0,1,2. 0 = no text updates, 1 = 
%                    some text outputs. 2 = detailed text outputs (default
%                    = 1)
%
% 2017 - Adam Charles and Alex Song

%  slightly modified by YZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the checks

if isempty(scan_params)                                                    % Make sure that scan_params is a struct
    clear scan_params
    scan_params = struct;
end

if (~isfield(scan_params,'scan_buff'))||isempty(scan_params.scan_buff)     % Default motion buffer
    scan_params.scan_buff = 10;
end
if (~isfield(scan_params,'motion'))||isempty(scan_params.motion)           % Default to having motion simulated
    scan_params.motion = false;
end
if (~isfield(scan_params,'scan_avg'))||isempty(scan_params.scan_avg)       % Default scan type is diffraction limited
    scan_params.scan_avg = 2;
end
if (~isfield(scan_params,'sfrac'))||isempty(scan_params.sfrac)             % Scan subsampling factor
    scan_params.sfrac = 1;
end
if (~isfield(scan_params,'verbose'))||isempty(scan_params.verbose)         % Default verbose level is 1
    scan_params.verbose = 1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%