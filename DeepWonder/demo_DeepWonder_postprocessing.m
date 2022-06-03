clc, clear
close all

%% a demo script to show discard too small/too dim components after network segmentation
%  you can supply any other criteria (like throw away boundary components, apply vessel masks, etc)
%  all threshold can be adjusted.
%  last update: 5/30/2020. YZ

%%
network_result = importdata(sprintf('%s\\network\\mat\\seg_29_wd.mat', target_dir));
% generate A
N_comp = length(network_result.final_mask_list);
network_A = zeros(size(network_result.final_contours, 1), size(network_result.final_contours, 2), N_comp);
network_C = [];
network_A_center = [];
for i = 1 : N_comp
    buf = zeros(size(network_A, 1), size(network_A, 2));
    valid_ind = sub2ind(size(buf), network_result.final_mask_list{i}.position(:, 1) + 1, network_result.final_mask_list{i}.position(:, 2) + 1);
    buf(valid_ind) = 1;
    network_A(:, :, i) = buf;
    network_C = [network_C; network_result.final_mask_list{i}.trace];
end
% throw away components that are too small
area_threshold = 50;
invalid_ind = [];
for i = 1 : N_comp
    buf = network_A(:, :, i);
    % calculate area
    curr_area = sum(buf, 'all');
     
    if curr_area < area_threshold
        invalid_ind(i) = 1;
    else
        invalid_ind(i) = 0;
    end
end
network_A(:, :, find(invalid_ind)) = [];
network_C(find(invalid_ind), :) = [];


%%
network_A_readout = [];
buf = dir(sprintf('%s\\network\\RMBG\\*.tif', target_dir));
network_raw = loadtiff(sprintf('%s\\network\\RMBG\\%s', target_dir, buf(1).name));
network_raw = single(network_raw);

for i = 1 : size(network_A ,3)
    i
    buf = network_A(:, :, i);
    curr_net_sig = squeeze(mean(bsxfun(@times, network_raw, buf), [1, 2]));

    network_A_readout(i, :) = curr_net_sig;
end
%% throw away those too dim neurons
threshold = 0.05;
std_network_A_readout = std(network_A_readout, 0, 2);
max_val = max(std_network_A_readout);
valid_ind = std_network_A_readout > max_val * threshold;
valid_ind = find(valid_ind);

%% output filtered components
network_A_center = com(reshape(network_A, [], size(network_A, 3)), size(network_A, 1), size(network_A, 2));
network_A_center_filt = network_A_center(valid_ind, :);
network_A_filt= network_A(:, :, valid_ind);
network_C_filt = network_C(valid_ind, :);



function cm = com(A,d1,d2,d3)
% cm short for center of mass
% center of mass calculation
% inputs:
% A: d X nr matrix, each column in the spatial footprint of a neuron
% d1, d2, d3: the dimensions of the 2-d (or 3-d) field of view

% output:
% cm: nr x 2 (or 3) matrix, with the center of mass coordinates

    if nargin < 4
        d3 = 1;
    end
    if d3 == 1
        ndim = 2;
    else
        ndim = 3;
    end

    nr = size(A,2);
    Coor.x = kron(ones(d2*d3,1),double(1:d1)');
    Coor.y = kron(ones(d3,1),kron(double(1:d2)',ones(d1,1)));
    Coor.z = kron(double(1:d3)',ones(d2*d1,1));
    cm = [Coor.x, Coor.y, Coor.z]'*A/spdiags(sum(A)',0,nr,nr);
    cm = cm(1:ndim,:)';
end