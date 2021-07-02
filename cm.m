% clc
clear
set_parameters;
build_catalog;
training_set_name = 'b-cond-sample-10000-L12-Nciv-1352-1650-dz-0-zCut-5000.mat'
% training_set_name ='sample-20000-L12-Nciv-14';
filename = sprintf('%s/processed_qsos_%s', ...
    processed_directory(release), ...
    training_set_name);
load(filename, 'test_ind', 'p_L1','p_no_c4', 'map_z_c4L2', 'training_set_name');
training_set_name
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_zqso                = cooksey_catalog{4};
ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0) & (all_RATING==3);

% cc=0;
% for quasar_ind=1:num_quasars
%     if ((all_zqso(quasar_ind) - 5000/speed_of_light) < all_z_c4(quasar_ind))
%         test_ind(quasar_ind) = 0;
%         cc=cc+1;

%     end
% end
a = 1e-6;
p_c4 = 1- p_no_c4 -a*p_L1;
% z_c4_condition = (all_zqso(test_ind) - 5e6/speed_of_light) > all_z_c4(test_ind);
z_c4_prior = (all_zqso(test_ind) - 5e6/speed_of_light) > map_z_c4L2;
% p_c4(~z_c4_prior) = 0; 
true_c4s = ind_has_c4(test_ind);% &  z_c4_condition;
%  all_zqso(test_ind)*(-5e6/speed_of_light)
for tr=[0.5];%, 0.8, 0.9, 0.95, 0.99]
    
    TP = nnz(true_c4s & p_c4>tr);
    TN = nnz(~true_c4s & p_c4<tr );
    FN = nnz(true_c4s & p_c4<tr);
    FP = nnz(~true_c4s & p_c4>tr)
    P = nnz(true_c4s );
    N = nnz(~true_c4s );
    confusion_matrix=[TP/P, FN/P; FP/N, TN/N];
    Accuracy = (TP+TN)/(P+N);
    ErrorRate = (FP+FN)/(P+N);
    % Sensitivity = TP/P;
    % Specificity = TN/N;
    fprintf('p:%.2f\nCM:[%.4f, %.4f\n    %.4f, %.4f]\nAccuracy:%.4f\nError Rate:%.4f\n'...
    ,tr, TP/P, FN/P, FP/N,TN/N,Accuracy,ErrorRate);
    fprintf('--------------------\n\n')
end


% tr=0.9;
% figure();
% % histogram(p_no_c4(ind_has_c4(test_ind) & p_c4<tr))
% % hold on
% histogram(p_L1(ind_has_c4(test_ind) & p_c4<0.1))
% % hold on
% % histogram(p_c4(ind_has_c4(test_ind) & p_c4<tr))
% % legend('pNull','pL1', 'pC4')

% xlabel('p(single line)')
% title('C13 detected a doublet but p(doublet)<0.1')
% saveas(gcf,'p_single.png')
