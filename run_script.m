c4_set_parameters
c4_build_catalog
c4_preload_qsos
% cd minFunc_2012
% addpath(genpath(pwd))
% mexAll
% cd ..
% training_release = 'dr7';
% filter_flags = load(sprintf('%s/filter_flags', processed_directory(training_release)));
% tif = (filter_flags==0);
% prior_ind = randsample(tif, int32(0.5*numel(tif)));
% test_ind = (~prior_ind);
% train_ind = (~ismember(all_QSO_ID(prior_ind), c4_QSO_ID));
% % learn_qso_model
% % generate_c4_samples
