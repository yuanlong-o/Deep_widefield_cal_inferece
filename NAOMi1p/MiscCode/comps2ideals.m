function ideal = comps2ideals(comps, baseim, varargin)

% ideal = comps2ideals(comps,varargin)
% 
% Calculate the SNR-adjusted ideal components frm the full components.
% 
% 2018 - Alex Song & Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if nargin > 2
    k = varargin{1};
else
    k = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

minNumEl = 5;                                                              % Set the minimum number of 
compsz   = size(comps);                                                    % Get size of components
if(length(compsz)==2)
  compsz = [compsz 1];
end
ideal    = reshape(bsxfun(@rdivide,comps,baseim),[],compsz(3));

for i = 1:size(comps,3)
    
%     cutoff = 1/(1+2/max(ideal(:,i)));
%     cutoff = 1/(k+2/prctile(ideal(:,i),95));
    sort_i = sort(ideal(:,i),'descend');
    sort_i = sort_i(~isnan(sort_i));
    cutoff = 1/(k+2/mean(sort_i(1:minNumEl)));
    
    if(sum(ideal(:,i)>cutoff))
        rp = regionprops(reshape(ideal(:,i)>cutoff,compsz(1),[]),'PixelIdxList');
        if(~isempty(rp))
            [~,val] = max(cellfun(@length,{rp.PixelIdxList}));
            temp    = zeros(size(ideal,1),1);
            if(length(rp(val).PixelIdxList)>=minNumEl)
                temp(rp(val).PixelIdxList) = ideal(rp(val).PixelIdxList,i);
            end
            ideal(:,i) = temp;
        end
    else
      ideal(:,i) = zeros(size(ideal,1),1);
    end
end

ideal = reshape(ideal,compsz);                                             % Reshape ideal components back to the sizes of the components
ideal = (ideal>0).*comps;                                                  % Extract only the portions of the components that exceed the SNR threshold

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%