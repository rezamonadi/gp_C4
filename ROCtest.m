% clear
% clc 
% set_parameters
% build_catalog
% filename = sprintf('%s/processed_qsos_%s.mat', ...
%     processed_directory(release), ...
%     training_set_name);
% load(filename);
% mkdir(sprintf('ROC-%s', training_set_name));

for a =[1, 0, 0.000001, 0.00001, 0.0001, 0.001]
    fig=figure();
    p_c4New = 1- p_no_c4 - a*p_L1;
    y_score = p_c4New;
    y_true = all_ind_c4(test_ind);
    [X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
    plot(X,Y)
    legend(sprintf('a=%.2e, AUC=%.5f',a, AUC))
    set(get(gca, 'YLabel'), 'String', 'TPR');
    set(get(gca, 'XLabel'), 'String', 'FPR');
    exportgraphics(fig, sprintf('ROC-%s/a-%.10f.pdf', training_set_name, a)...
                                    ,'ContentType','vector')
end
