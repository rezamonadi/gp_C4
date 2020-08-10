c4_set_parameters
c4_build_catalog
% c4_preload_qsos
% cd minFunc_2012
% addpath(genpath(pwd))
% mexAll
% cd ..
training_release = 'dr7';
f = load(sprintf('%s/filter_flags', processed_directory(training_release)));
filter_flags= f.filter_flags;

half_ID = randsample(all_QSO_ID, int32(0.5*numel(all_QSO_ID)));
test_ind = ((~ismember(all_QSO_ID, half_ID)) & (filter_flags==0));

prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
               (ismember(all_QSO_ID, half_ID)));  
train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
                ismember(all_QSO_ID, half_ID) );

c4_inds = ismember(all_QSO_ID, c4_QSO_ID);
% % learn_qso_model
% % generate_c4_samples
