function tiff_writer(outputFileName,mov,varargin)

% function Y = tiff_writer(name,mov)
%
% Code to write tiff stack.
%
% 2015 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs and set variables

if nargin < 2
    error('Need to give a name and a file to write!')
end

if nargin > 2;    saveOpt = varargin{1};
else;             saveOpt = 'standard';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write images

switch lower(saveOpt)
    case 'standard'
        for kk=1:length(mov(1, 1, :))
           imwrite(double(mov(:, :, kk)), outputFileName, 'WriteMode', 'append','compression','none');
        end
    case '16bit'
        for kk=1:length(mov(1, 1, :))
           imwrite(double(mov(:, :, kk)), outputFileName, 'WriteMode', 'append','compression','none');
        end
    otherwise
        error('Bad save option given!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
