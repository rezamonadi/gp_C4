fprintf('Setting paramters ...\n')
set_parameters

train_ratio =0.5;


fprintf('Building catalogs ...\n')
build_catalog
fprintf('Preloading QSOs ...\n')

preload_qsos
%load('data/dr7/processed/preloaded_qsos.mat')
fprintf('preparing voigt.c ...\n')

cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
mex voigt.c -lcerf
fprintf('preparing testing, training, and prior indeces ...\n')

f = load(sprintf('%s/filter_flags', processed_directory(training_release)));
filter_flags= f.filter_flags;

half_ID = randsample(all_QSO_ID, int32(train_ratio*numel(all_QSO_ID)));
test_ind = ((~ismember(all_QSO_ID, half_ID)) & (filter_flags==0));

prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
    (ismember(all_QSO_ID, half_ID)));
train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
   ismember(all_QSO_ID, half_ID) );

fprintf('Learning model ...\n')

learn_qso_model
%variables_to_load = {'training_release', 'train_ind', 'max_noise_variance', ...
%                    'minFunc_options', 'rest_wavelengths', 'mu', ...
%                     'initial_M', 'M',  'log_likelihood', ...
%                     };

%load(sprintf('%s/learned_model_%s_norm_%d-%d',             ...
%             processed_directory(training_release), ...
%             training_set_name, ...
%	     normalization_min_lambda, normalization_max_lambda), ...
%     variables_to_load{:});

fprintf('Generating samples for integrating out parameters in the model...\n')
generate_c4_samples


%variables_to_load = {'uniform_min_log_nciv', 'uniform_max_log_nciv', ...
%                     'fit_min_log_nciv', 'fit_max_log_nciv', 'alpha', ...
%                     'extrapolate_min_log_nciv', ...
%                     'offset_z_samples', 'offset_sigma_samples', 'log_nciv_samples', 'nciv_samples'};
%load(sprintf('%s/civ_samples-%s', processed_directory(training_release), training_set_name), ...
%     variables_to_load{:});
c4_catalog_name = 'cooksey';
fprintf('Processing ...\n')
process_qsosL2
% % % variables_to_load = {'training_release', 'training_set_name', ...
% % %     'c4_catalog_name', 'prior_ind', 'release', ...
% % %     'test_set_name', 'test_ind', 'prior_z_qso_increase', ...
% % %     'max_z_cut', 'num_lines', 'min_z_c4s', 'max_z_c4s', ...
% % %     'log_priors_no_c4', 'log_priors_c4', ...
% % %     'log_likelihoods_no_c4', 'sample_log_likelihoods_c4', ...
% % %     'log_likelihoods_c4', 'log_posteriors_no_c4', ...
% % %     'log_posteriors_c4', 'model_posteriors', 'p_no_c4', ...
% % %     'p_c4', 'map_z_c4', 'map_N_c4'};

% % % filename = sprintf('%s/b-sample-processed_qsos_%s', ...
% % %     processed_directory(release), ...
% % %     test_set_name);

% % % load(filename, variables_to_load{:});


for quasar_ind=1604:1700
   fprintf('This ind:%d\n',quasar_ind)
   plot_QSO
end

