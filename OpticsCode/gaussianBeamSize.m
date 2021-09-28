function gauss_sz = gaussianBeamSize(psf_params,dist,apod)

% gauss_sz = gaussianBeamSize(psf_params,dist,apod)
%
% Function to calculate an overestimate of the total Gaussian beam waist at
% a large distance from the focal point. This overestimate is approximately
% equivalent to an apodization of 5, and is returned as a 3-vector
% corresponding to the 3 dimension of the beam (axial given last).
% The inputs are
%
%   - psf_params - Struct contaning the parameters for the PSF
%          .NA     - Numerical aperture of Gaussian beam
%          .n      - Refractive index of propagation material
%
%   - dist         - Distance away from focal point of Gaussian beam
%   - apod         - Scaling factor (linear)
%
% The output is
%   - gauss_sz     - 3-vector (X,Y,Z) corresponding to maximum interaction 
%                    size in X and Y, 0 in Z
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin<3)
  apod = 2;
end

gauss_sz = ...
        ceil(tan(asin(psf_params.objNA/psf_params.n))*dist*1.5)*apod*[1 1 0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%