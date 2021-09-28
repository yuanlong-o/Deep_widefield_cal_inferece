function Uout = applyZernike(Uin,X,Y,k,abb)

% Uout = applyZernike(Uin,X,Y,k,abb)
%  
% Applies the Zernike aberration given an input field, position and
% aberrations. The inputs are:
% 
%   - Uin            - Input scalar field
%   - X              - X field positions
%   - Y              - Y field positions
%   - k              - Wavenumber [1/m]
%   - abb            - Zernike aberrations
%
% The outpus is:
%   - Uout           - Output scalar field
%
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phase = zeros(size(X),'single');
idxs = find(abb~=0);
if(~isempty(idxs))
  for i = idxs
    phase = phase+abb(i)*zernike(i,X,Y);
  end
end
phase = phase-mean(phase(:));
Uout = Uin.*exp(1i*k*phase);