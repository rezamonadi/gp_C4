% variables_to_load = {'p_c4', 'test_ind'};
% release='dr7';

% c4_catalog_name = 'cooksey-50';
% filename = sprintf('%s/b-processed_qsos_%s', ...
set_parameters;
build_catalog;
training_set_name = 'sample-25000-L12-Nciv-1352-1550-dz-125-zCut-5000'

filename = sprintf('%s/processed_qsos_%s.mat', ...
    processed_directory(release), ...
    training_set_name);
load(filename);

load(sprintf('%s/catalog', processed_directory(release)));
% variables_to_load = {'offset_z_samples', 'log_nciv_samples'};
% load(sprintf('%s/civ_samples-%s', processed_directory(training_release), training_set_name), ...
%         variables_to_load{:});
% catalog = load(sprintf('%s/catalog', processed_directory(release)));
% % load preprocessed QSOs

variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask'};
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_load{:});


% enable processing specific QSOs via setting to_test_ind
if (ischar(test_ind))
  test_ind = eval(test_ind);
end
test_ind = test_ind & filter_flags==0;
all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);


num_quasars = sum(test_ind);
ID = all_QSO_ID(test_ind);
Z_Cooksey_compare=[0];
NCIV_Cooksey_compare = [0];
jj=0;
for quasar_ind=1:num_quasars

    % this_wavelengths    =    all_wavelengths{quasar_ind};
    % this_wavelengths    =              this_wavelengths';
    % this_pixel_mask     =     all_pixel_mask{quasar_ind};
    % this_pixel_mask     =               this_pixel_mask';
    this_ID = ID{quasar_ind};
    this_systems = ismember(c4_QSO_ID, this_ID);
    % convert to QSO rest frame
    z_qsos = all_zqso(test_ind);
    % this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qsos(quasar_ind));
    % unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
    %  (this_rest_wavelengths <= max_lambda);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    % this_unmasked_wavelengths = this_wavelengths(unmasked_ind);

    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    % ind = unmasked_ind & (~this_pixel_mask);

    % this_wavelengths      =      this_wavelengths(ind);
    % this_rest_wavelengths = this_rest_wavelengths(ind);


    % sample_z_c4 = ...
    % min_z_c4s(quasar_ind) +  ...
    % (max_z_c4s(quasar_ind) - min_z_c4s(quasar_ind))*offset_z_samples;

    if(sum(this_systems)>0)
        this_c4s = NCIV(this_systems);
        if (this_c4s>0)
            jj=jj+1;

            this_Zs  = Z_c4(this_systems);
            % [~, maxind] = nanmax(sample_log_likelihoods_c4(quasar_ind, :));
            % map_log_c4  = log_nciv_samples(maxind);
            % map_z_c4    = sample_z_c4(maxind);        
            matched_Z_c4= this_Zs(abs(this_Zs -map_z_c4L2(quasar_ind))==min(abs(this_Zs-map_z_c4L2(quasar_ind))));
            matched_N_c4= this_c4s(abs(this_Zs -map_z_c4L2(quasar_ind))==min(abs(this_Zs-map_z_c4L2(quasar_ind))));
            NCIV_Cooksey_compare(jj) = matched_N_c4 - map_N_c4L2(quasar_ind);
            % if(NCIV_Cooksey_compare(jj)<-2)
            %     jj
            %     matched_Z_c4
            %     map_z_c4

            %     matched_N_c4
            %     map_log_c4
            % end
            Z_Cooksey_compare(jj) = matched_Z_c4 - map_z_c4L2(quasar_ind);
        end
    end
end
fig=figure();
histogram(NCIV_Cooksey_compare, 50)
set(get(gca, 'XLabel'), 'String', 'N(C13) - N(MAP)');
exportgraphics(fig, sprintf('%s-N(13)-N(MAP).pdf',...
                            training_set_name),...
                            'ContentType','vector')
fig=figure();
histogram(Z_Cooksey_compare, 50)
set(get(gca, 'XLabel'), 'String', 'z(C13) -z(MAP)');
exportgraphics(fig, sprintf('%s-z(13)-z(MAP).pdf',...
                            training_set_name),'ContentType','vector')


a=1e-10;
p_c4 = 1 - p_no_c4  - a*p_L1;                           
fig=figure();
y_score = p_c4;
y_true = all_ind_c4(test_ind);
[X,Y,T,AUC] =perfcurve(y_true, y_score, 'true');
plot(X,Y)
legend(sprintf('a=%.2e, AUC=%.5f',a, AUC))
set(get(gca, 'YLabel'), 'String', 'TPR');
set(get(gca, 'XLabel'), 'String', 'FPR');
exportgraphics(fig, sprintf('%s-ROC.pdf', training_set_name)...
                                    ,'ContentType','vector')

