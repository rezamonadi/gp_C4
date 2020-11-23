c4_set_parameters
min_lambda = 1315;
max_lambda = 1550;
dlambda = 0.05;
train_ratio =0.5;



c4_build_catalog
c4_preload_qsos
% load('data/dr7/processed/preloaded_qsos.mat')
cd minFunc_2012
addpath(genpath(pwd))
mexAll
cd ..
f = load(sprintf('%s/filter_flags', processed_directory(training_release)));
filter_flags= f.filter_flags;

half_ID = randsample(all_QSO_ID, int32(train_ratio*numel(all_QSO_ID)));
test_ind = ((~ismember(all_QSO_ID, half_ID)) & (filter_flags==0));

prior_ind = ((ismember(all_QSO_ID, c4_QSO_ID)) & (filter_flags==0) & ...
   (ismember(all_QSO_ID, half_ID)));
train_ind = (~ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & ...
   ismember(all_QSO_ID, half_ID) );

learn_qso_model
generate_c4_samples
test_set_name= 'cooksey_half-2';
c4_catalog_name = 'cooksey-half-1';
process_qsos

y_score = p_c4;
y_test = ind_c4(test_ind);
fname = sprintf('ROC-l%d-%d-dl%f-tr-%f.mat',min_lambda,max_lambda,...
    dlambda, train_ratio);
save(fname, 'y_score', 'y_test', 'mu', 'M', 'rest_wavelengths');
save('indeces.mat', 'test_ind', 'train_ind', 'prior_ind');


