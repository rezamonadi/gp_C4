c4_set_parameters
c4_build_catalog
c4_preload_qsos
cd minFunc_2012
addpath(genpath(pwd))
mexAll
train_release = 'dr7'
catalog = load(sprintf('%s/catalog', processed_directory(training_release)));
train_ind_full= [' ~ismember(all_QSO_ID, c4_QSO_ID) & ' '(catalog.filter_flags==0)'];
tif = eval(train_ind_all);
train_ind = randsample(tif, int32(0.5*numel(tif)));
learn_qso_model
