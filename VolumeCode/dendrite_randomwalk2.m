function path_out = dendrite_randomwalk2(M,root,ends,distsc,maxlength,fillweight,maxel,minlength)

% path_out = dendrite_randomwalk2(M,root,ends,distsc,maxlength,fillweight,
%            maxel,minlength)
% 
% Wrapper for dendrite_randomwalk_cpp.cpp, a fast function to produce a
% random walk through the neural volume. The inputs for this function are
%
%   - M          - matrix that indicates the difficulty to occupy a
%                  particular location
%   - root       - 1x3 start location for the random walk path
%   - ends       - Nx3 directed ending locations for the random walk
%   - distsc     - weighting parameter greater than 0 setting relative
%                  weight of directed walk
%   - maxlength  - maximum length for the random path
%   - fillweight - weighting value for a single step against occupancy rate
%   - maxel      - maximum number of components within a single voxel
%   - minlength  - minimum length for the random path
%
% The output of this function is
%   - path_out   - Output paths of the random walks
% 
% 2017 - Alex Song

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing 

if(~isa(M,'single'))
  warning('M must be a single, casting as a single')
  M = single(M);
end

if(~isa(root,'int32'))
  root = int32(root);
end

if(~isa(ends,'int32'))
  ends = int32(ends);
end

if(~isa(distsc,'single'))
  distsc = single(distsc);
end

if(~isa(maxlength,'int32'))
  maxlength = int32(maxlength);
end

if(~isa(fillweight,'single'))
  fillweight = single(fillweight);
end

if(~isa(maxel,'int32'))
  maxel = int32(maxel);
end

if(~isa(minlength,'int32'))
  minlength = int32(minlength);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,path_out] = dendrite_randomwalk_cpp(M,root,ends,distsc,maxlength,...
                               fillweight,maxel,minlength,int32(size(M))); % Run fast random-walk code

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%