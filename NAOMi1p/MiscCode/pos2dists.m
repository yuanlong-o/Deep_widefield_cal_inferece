function distances = pos2dists(positions)

% distances = pos2dists(positions)
%
% Function to convert a series of points into a matrix of distances between
% all the points. Input is the NxD list of N D-dimensional points. Output
% is the NxN matrix of distances between any two points (position i,j in
% the output is the distance between points i and j).
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dims      = size(positions,2);
numpoints = size(positions,1);
distances = zeros(numpoints);

for i = 1:numpoints
    distVec = zeros(numpoints,1);
    for j = 1:dims
        distVec = distVec+(positions(:,j)-positions(i,j)).^2;
    end
    distances(:,i) = sqrt(distVec);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%