function noisy_out = PoissonGaussNoiseModel_camera(clean_in, noise_params)

% noise model for camera

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

% noise_params = check_noise_params(noise_params);                           % Check the noise parameter struct for missing elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw noisy outputs
noisy_out = noise_params.gamma * poissrnd(clean_in);
noisy_out = noisy_out  + noise_params.offset;
noisy_out = noisy_out+ normrnd(0, noise_params.sigma_p, size(noisy_out, 1), size(noisy_out, 2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
