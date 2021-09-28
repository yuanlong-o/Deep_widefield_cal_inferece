function idList = isolateVisibleSomas(vol,psf,vol_params, neur_params, varargin)

% idList = isolateVisibleSomas(vol,psf,vol_params, neur_params, {thresh})
% Compare all locations of neural centers in vol to the main imaging 
% plane (center of the volume). The distance threshold should be related
% to the falloff of the PSF
%
% 2020 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 4
    thresh = varargin{1};
else
    thresh = 0;
end

vMid         = vol_params.vol_sz(3)/2;                                   % Get midpoint of volume (plane of imaging)
nNeur        = size(vol.gp_nuc,1);                                       % Get number of neural somas in the volume
zLocs        = vol.locs(1:nNeur,end);                                    % Get the z (axial) locations of the neurons
psfHalfWidth = widthestimate3D(psf.psf);                                 % Calculate the half-width of the point-spread function
thresh       = thresh + psfHalfWidth(end)/2 + neur_params.avg_rad(1);    % The threshold is calculated based on the psf half-width + neural radius. 
idList       = find(abs(zLocs - vMid)<thresh);                           % Find all the somas no more than 'thresh' away from the scan plane

end


