function varargout = corrTimeAndSpace(t1,t2,s1,s2,varargin)

% varargout = corrTimeAndSpace(t1,t2,s1,s2,varargin)
%
% Correlate two decompositions using both the timecourses (t1 and t2) and 
% spatial profiles (s1 and s2).
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 4
    thresh_t = varargin{1};
else
    thresh_t = 0.1;
end

if nargin > 5
    thresh_s = varargin{2};
else
    thresh_s = 0.5;
end

if nargin > 6
    thresh_pix = varargin{3};
else
    thresh_pix = 0.1;
end

s1 = reshape(s1,[],size(s1,ndims(s1)));                                    % Make sure s1 is in matrix form
s2 = reshape(s2,[],size(s2,ndims(s2)));                                    % Make sure s2 is in matrix form

if (size(s1,2)~=size(t1,2))||(size(s2,2)~=size(t2,2))
    error('Must have matching number of spatial and temporal components.')
end

if (size(t1,1)~=size(t1,1))||(size(s1,1)~=size(s2,1))
    error('Must have matching dimensions for spatial/temporal components.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate temporal correlations and 

s1  = bsxfun(@gt,s1,thresh_pix*max(s1,[],1));                              % Threshold small pixels
s2  = bsxfun(@gt,s2,thresh_pix*max(s2,[],1));                              % Threshold small pixels
tp1 = sum(s1,1);
tp2 = sum(s2,1);

olap = (s1.')*single(s2);                                                            % Calculate number of overlaping pixels between spatial maps
C = corrcoef([t1,t2]);                                                     % Calculate (temporal) correlation matrix
C = C(1:size(t1,2),(size(t1,2)+1:end));                                    % Extract only the correlations needed

C(C<thresh_t) = NaN;                                                       % Threshold all correlations that are too small
C((bsxfun(@rdivide, olap, tp2)<thresh_s)&...
                        ((bsxfun(@rdivide, olap, tp1.')<thresh_s))) = NaN; % Threshold correlations for non-sufficiently-overlapping profiles

if nargout > 2
    varargout{3} = C;
end

corr_vec  = [];                                                            % Initialize the correlation storage vector
all_pairs = [];                                                            % Initialize the set of pairings
while any(~isnan(C(:)))
    cm        = max(C(~isnan(C)));                                         % Find the current max correlation 
    [I,J]     = find(C==cm,1,'first');                                     % Find the location of the current max correlation 
    corr_vec  = cat(1,corr_vec,cm);                                        % Record the max correlation 
    all_pairs = cat(1,all_pairs,[I,J]);                                    % Add in the new pair to the pair list
    C(:,J)    = NaN;                                                       % No longer pair this column
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output parsing

if nargout == 1
    tmp.corrs    = corr_vec;
    tmp.pairs    = all_pairs;
    varargout{1} = tmp;
elseif nargout > 1
    varargout{1} = corr_vec;
    varargout{2} = all_pairs;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% C(olap<thresh_s) = NaN;                                                    % Threshold correlations for non-sufficiently-overlapping profiles
% olap = bsxfun(@rdivide, olap, sum(s2,1));                                  % Normalize to the number of total pixels in s2
% C    = t1.'*t2;                                                            % Calculate (temporal) correlation matrix

