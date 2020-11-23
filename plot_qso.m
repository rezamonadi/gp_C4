
% plot_a_processed_qso.m : plot sample of spectrum with model dla
% and the sample likelihoods
%
% Usage:
% ----
% % first load the catalogues and spectra
% addpath dr16q/
% load_processed_catalogs
% % then process on the chosen spectra
% selected_thing_ids = [43880646]; % 23097883 43880646 355787041 352241122
% process_a_qso_multiple_dlas_meanflux
% % plot the selected thingID
% plot_a_processed_qso

% load QSO model from training release
variables_to_load = {'rest_wavelengths', 'mu', 'M'};
load('data/dr7/processed/learned_model_Cooskey_all_qso_catalog_norm_1310-1325.mat',variables_to_load{:});

% load DLA samples from training release
variables_to_load = {'offset_samples', 'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples', processed_directory(training_release)), ...
     variables_to_load{:});
     
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

% create the absorber redshifts
min_z_c4s(quasar_ind) = min_z_c4(this_wavelengths, z_qso);
max_z_c4s(quasar_ind) = max_z_c4(this_wavelengths, z_qso);

sample_z_c4 = ...
   min_z_c4s(quasar_ind) +  ...
  (max_z_c4s(quasar_ind) - min_z_c4s(quasar_ind)) * offset_samples;

% preprocess model interpolants
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');


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

% build the model
this_mu = mu_interpolator( this_rest_wavelengths);


unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
               (this_rest_wavelengths <= max_lambda);
% 
% ind = (~this_pixel_mask(unmasked_ind));
% available_ind = find(~ind & ~this_pixel_mask);
% ind(min(available_ind(available_ind > find(ind, 1, 'last' )))) = true;
% ind(max(available_ind(available_ind < find(ind, 1, 'first')))) = true;

fig = figure('visible', 'off');
clf();
hold on
%% plotting the sample_log_likelihoods_dla in z_dla and log_nhi axes
%     subplot('position', [0.05 0.49 0.90 0.45]);
%     hold on
       
    [~, maxind] = nanmax(sample_log_likelihoods_c4(quasar_ind, :));
    map_z_c4    = sample_z_c4(maxind);        
    map_log_c4  = nciv_samples(maxind);

    % construct dla_mu_map
    absorption = voigt(padded_wavelengths, map_z_c4, ...
        nciv_samples(maxind), num_lines);
    s= size(absorption)*[1;0];
    absorption = absorption(ind);

    c4_mu = this_mu(1:s).* absorption;
    

%     subplot('position', [0.05 0.05 0.90 0.38]);

%     hold on
        this_z_c4 = (this_wavelengths / 1550) - 1;
        plot(this_z_c4, this_flux, 'Color', [0.01, 0.2, 0.8]);
        hold on
        xlim([min(sample_z_c4) max(sample_z_c4)]);
        ylim([-1     5]);
        xlabel('(observed wavelengths $\lambda$ (\AA) / 1550 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
        ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
        if( p_c4(quasar_ind)>0.7)
            plot(this_z_c4(1:s), c4_mu, 'Color', [1,0.2, .5]);
            
            hold on
        end
    
        plot( (1 + z_qso) * (rest_wavelengths / 1550) - 1, mu, 'Color', [0.01, 0.98, 0.05]);
    hold on

ID = all_QSO_ID{test_ind};
test_ind_c4 = ind_c4(test_ind);
tit = sprintf('ID: %s, C4:%d, P(CIV):%.2f', ID, test_ind_c4(quasar_ind), p_c4(quasar_ind));
title(tit);

fid = num2str(quasar_ind);
% saveas(fig, fid, 'jpg');
fidpdf = sprintf('FP/%s.pdf',fid);
exportgraphics(fig, fidpdf,'ContentType','vector')