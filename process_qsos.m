% process_qsos: run DLA detection algorithm on specified objects

% load C4 catalog
Full_catalog = ...
    load(sprintf('%s/catalog', processed_directory(training_release)));

if (ischar(prior_ind))
  prior_ind = eval(prior_ind);
end

prior.z_qsos  = Full_catalog.all_zqso(prior_ind);
prior.c4_ind = prior_ind;
% My prior_ind here is already those OK sight of lines that have CIV


% filter out CIVs from prior catalog corresponding to region of spectrum below
% Ly-alpha QSO rest

prior.z_c4 = Full_catalog.all_z_c4(prior_ind);

for i = size(prior.z_c4)
  if (observed_wavelengths(civ_1548_wavelength , prior.z_c4(i)) < ...
      observed_wavelengths(lya_wavelength, prior.z_qsos(i)))
    prior.c4_ind(i) = false;
  end
end

prior = rmfield(prior, 'z_c4');

% load QSO model from training release
variables_to_load = {'rest_wavelengths', 'mu', 'M'};
load('data/dr7/processed/learned_model_Cooskey_all_qso_catalog_norm_1310-1325.mat',variables_to_load{:});

% load DLA samples from training release
variables_to_load = {'offset_samples', 'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples', processed_directory(training_release)), ...
     variables_to_load{:});

% load redshifts from catalog to process
catalog = load(sprintf('%s/catalog', processed_directory(release)));

% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask'};
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_load{:});

% enable processing specific QSOs via setting to_test_ind
if (ischar(test_ind))
  test_ind = eval(test_ind);
end

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);

z_qsos = catalog.all_zqso(test_ind);

num_quasars = numel(z_qsos);

% preprocess model interpolants
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');
M_interpolator = ...
    griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');
% log_omega_interpolator = ...
%     griddedInterpolant(rest_wavelengths,        log_omega, 'linear');

% initialize results
min_z_c4s                   = nan(num_quasars, 1);
max_z_c4s                   = nan(num_quasars, 1);
log_priors_no_c4           = nan(num_quasars, 1);
log_priors_c4              = nan(num_quasars, 1);
log_likelihoods_no_c4      = nan(num_quasars, 1);
sample_log_likelihoods_c4  = nan(num_quasars, num_C4_samples);
log_likelihoods_c4         = nan(num_quasars, 1);
log_posteriors_no_c4       = nan(num_quasars, 1);
log_posteriors_c4          = nan(num_quasars, 1);

% c_0   = exp(log_c_0);
% tau_0 = exp(log_tau_0);
% beta  = exp(log_beta);

for quasar_ind = 1:num_quasars
  tic;

  z_qso = z_qsos(quasar_ind);

  fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...', ...
         quasar_ind, num_quasars, z_qso);

  this_wavelengths    =    all_wavelengths{quasar_ind};
  this_flux           =           all_flux{quasar_ind};
  this_noise_variance = all_noise_variance{quasar_ind};
  this_pixel_mask     =     all_pixel_mask{quasar_ind};
  
  this_wavelengths    =    this_wavelengths';  
  this_flux           =           this_flux';         
  this_noise_variance = this_noise_variance';
  this_pixel_mask     =     this_pixel_mask';   

  % convert to QSO rest frame
  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);

  unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
                 (this_rest_wavelengths <= max_lambda);
  % keep complete copy of equally spaced wavelengths for absorption
  % computation
  this_unmasked_wavelengths = this_wavelengths(unmasked_ind);

  % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
  % in read_spec_dr7.m
  ind = unmasked_ind & (~this_pixel_mask);

  this_wavelengths      =      this_wavelengths(ind);
  this_rest_wavelengths = this_rest_wavelengths(ind);
  this_flux             =             this_flux(ind);
  this_noise_variance   =   this_noise_variance(ind);

%   this_lya_zs = ...
%       (this_wavelengths - lya_wavelength) / ...
%       lya_wavelength;

  % c4 existence prior
  less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));

  this_num_c4    = nnz(prior.c4_ind(less_ind));
  this_num_quasars = nnz(less_ind);
  this_p_c4 = this_num_c4 / this_num_quasars;

  log_priors_c4(quasar_ind) = ...
      log(                   this_num_c4) - log(this_num_quasars);
  log_priors_no_c4(quasar_ind) = ...
      log(this_num_quasars - this_num_c4) - log(this_num_quasars);

  fprintf('\n');
  fprintf(' ...     p(   CIV | z_QSO)  mvn      : %0.3f\n',     this_p_c4);
  fprintf(' ...     p(no CIV | z_QSO)        : %0.3f\n', 1 - this_p_c4);
  
  % interpolate model onto given wavelengths
  this_mu = mu_interpolator( this_rest_wavelengths);
  this_M  =  M_interpolator({this_rest_wavelengths, 1:k});

%   this_log_omega = log_omega_interpolator(this_rest_wavelengths);
%   this_omega2 = exp(2 * this_log_omega);

%   this_scaling_factor = 1 - exp(-tau_0 .* (1 + this_lya_zs).^beta) + c_0;

%   this_omega2 = this_omega2 .* thiscaling_factor.^2;

  % baseline: probability of no DLA model
%   disp(size(this_M))
%   disp(size(this_flux))
%   disp(size(this_noise_variance))
  log_likelihoods_no_c4(quasar_ind) = ...
      log_mvnpdf_low_rank(this_flux, this_mu, this_M, ...
          this_noise_variance);

  log_posteriors_no_c4(quasar_ind) = ...
      log_priors_no_c4(quasar_ind) + log_likelihoods_no_c4(quasar_ind);

  fprintf(' ... log p(D | z_QSO, no CIV)     : %0.2f\n', ...
                log_likelihoods_no_c4(quasar_ind));

  min_z_c4s(quasar_ind) = min_z_c4(this_wavelengths, z_qso);
  max_z_c4s(quasar_ind) = max_z_c4(this_wavelengths, z_qso);

  sample_z_c4 = ...
       min_z_c4s(quasar_ind) +  ...
      (max_z_c4s(quasar_ind) - min_z_c4s(quasar_ind)) * offset_samples;

  % ensure enough pixels are on either side for convolving with
  % instrument profile
  padded_wavelengths = ...
      [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
                log10(min(this_unmasked_wavelengths)) - pixel_spacing,         ...
                width)';                                                       ...
       this_unmasked_wavelengths;                                              ...
       logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,         ...
                log10(max(this_unmasked_wavelengths)) + width * pixel_spacing, ...
                width)'                                                        ...
      ];

  % [mask_ind] to retain only unmasked pixels from computed absorption profile
  % this has to be done by using the unmasked_ind which has not yet
  % been applied this_pixel_mask.
  ind = (~this_pixel_mask(unmasked_ind));

  % compute probabilities under DLA model for each of the sampled
  % (normalized offset, log(N HI)) pairs
  parfor i = 1:num_C4_samples
    % absorption corresponding to this sample
    absorption = voigt(padded_wavelengths, sample_z_c4(i), ...
                       nciv_samples(i), num_lines);

    absorption = absorption(ind);

    c4_mu     = this_mu     .* absorption;
    c4_M      = this_M      .* absorption;
%     dla_omega2 = this_omega2 .* absorption.^2;

    sample_log_likelihoods_c4(quasar_ind, i) = ...
        log_mvnpdf_low_rank(this_flux, c4_mu, c4_M, ...
             this_noise_variance);
  end

  % compute sample probabilities and log likelihood of DLA model in
  % numerically safe manner
  max_log_likelihood = max(sample_log_likelihoods_c4(quasar_ind, :));
  sample_probabilities = ...
      exp(sample_log_likelihoods_c4(quasar_ind, :) - ...
          max_log_likelihood);
  log_likelihoods_c4(quasar_ind) = ...
      max_log_likelihood + log(mean(sample_probabilities));

  log_posteriors_c4(quasar_ind) = ...
      log_priors_c4(quasar_ind) + log_likelihoods_c4(quasar_ind);

  fprintf(' ... log p(D | z_QSO,    CIV)     : %0.2f\n', ...
                log_likelihoods_c4(quasar_ind));
  fprintf(' ... log p(CIV | D, z_QSO)        : %0.2f\n', ...
                log_posteriors_c4(quasar_ind));
  % fprintf(' ... Num_CIV                      : %d\n ', ... 
  %               Full_catalog.all_Num_c4_sys(quasar_ind))
  % % fprintf('... FilterFlag                    : %d\n ', filter)
  fprintf(' took %0.3fs.\n', toc);
end

% compute model posteriors in numerically safe manner
max_log_posteriors = ...
    max([log_posteriors_no_c4, log_posteriors_c4], [], 2);

model_posteriors = ...
    exp([log_posteriors_no_c4, log_posteriors_c4] - max_log_posteriors);

model_posteriors = model_posteriors ./ sum(model_posteriors, 2);

p_no_c4 = model_posteriors(:, 1);
p_c4    = 1 - p_no_c4;

% save results
variables_to_save = {'training_release', 'training_set_name', ...
                     'c4_catalog_name', 'prior_ind', 'release', ...
                     'test_set_name', 'test_ind', 'prior_z_qso_increase', ...
                     'max_z_cut', 'num_lines', 'min_z_c4s', 'max_z_c4s', ...
                     'log_priors_no_c4', 'log_priors_c4', ...
                     'log_likelihoods_no_c4', 'sample_log_likelihoods_c4', ...
                     'log_likelihoods_c4', 'log_posteriors_no_c4', ...
                     'log_posteriors_c4', 'model_posteriors', 'p_no_c4', ...
                     'p_c4'};
cc
filename = sprintf('%s/processed_qsos_%s', ...
                   processed_directory(release), ...
                   test_set_name);

save(filename, variables_to_save{:}, '-v7.3');
