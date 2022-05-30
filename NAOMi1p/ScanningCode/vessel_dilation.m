
function vessel_img = vessel_dilation(valid_vessel, N_cut, seed, max_disk, dilate_trace, t)
%% this file is used to simulate the vessel dilations
%  N_cut: FOV divisions
%  seed: rand seed for repeating experiments
%  max_disk: dilation strength
%  dilate_trace: dilation amplitude

%  last update: 4/24/2022.

% load vessel file

% random seletect parts.
[size_h, size_w, size_z] = size(valid_vessel);

%  random cut point
rng(seed)
cut_point = rand(3, N_cut);


% 
cut_point_array = [];
for i = 1 : 3
    for j = 1 : N_cut
        tmp = j * 1 / (N_cut + 1) + (cut_point(i, j) - 0.5) / 2 / N_cut;
        cut_point_array(i, j) = min(max(round(tmp * size(valid_vessel, i)), 1), size(valid_vessel, i)); 
    end    
end


rng(seed + 1)

vessel_img = zeros(size_h, size_w, size_z, 'single');
for i = 1 : N_cut + 1
    dim1_array = [1, cut_point_array(1, :), size_h];
    for j = 1 : N_cut + 1
        dim2_array = [1, cut_point_array(2, :), size_w];
        for k = 1 : N_cut + 1
            dim3_array = [1, cut_point_array(3, :), size_z];

            curr_vessel = valid_vessel(dim1_array(i): dim1_array(i + 1),...
                                       dim2_array(j): dim2_array(j + 1), ...
                                       dim3_array(k): dim3_array(k + 1));

            % do dilate
            curr_ind = sub2ind([N_cut+1, N_cut+1, N_cut+1], i, j, k);
            se = strel('sphere', round(max(max_disk * dilate_trace(curr_ind, t), 1)));
%                 tic
            curr_vessel = imdilate(curr_vessel>0, se );
            curr_vessel = curr_vessel > 0;
%                 toc
            vessel_img(dim1_array(i): dim1_array(i + 1),...
                       dim2_array(j): dim2_array(j + 1), ...
                       dim3_array(k): dim3_array(k + 1)) = curr_vessel;
        end
    end
end
vessel_img = 1 - vessel_img;

end
