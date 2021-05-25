
% plot_a_processed_qso.m : plot sample of spectrum with model dla
% and the sample likelihoods
%
% Usage:
% ----
% % first load the catalogues and spectra
% addpath dr16q/
% load_processed_catalogsset_parameters;
set_parameters;
build_catalog;


% variables_to_load = {'training_release', 'training_set_name', ...
%     'c4_catalog_name', 'prior_ind', 'release', ...
%     'test_ind', 'prior_z_qso_increase', ...
%     'max_z_cut', 'num_lines', 'min_z_c4s', 'max_z_c4s', ...
%     'log_priors_no_c4', 'log_priors_c4', ...
%     'log_likelihoods_no_c4', 'sample_log_likelihoods_c4', ...
%     'log_likelihoods_c4', 'log_posteriors_no_c4', ...
%     'log_posteriors_c4', 'model_posteriors', 'p_no_c4', ...
%     'p_c4', 'map_z_c4', 'map_N_c4', 'map_sigma_c4'};

filename = sprintf('%s/processed_qsos_%s', ...
    processed_directory(release), ...
    training_set_name);
processed=matfile(filename);


load('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');
% % then process on the chosen spectra
% selected_thing_ids = [43880646]; % 23097883 43880646 355787041 352241122
% process_a_qso_multiple_dlas_meanflux
% % plot the selected thingID
% plot_a_processed_qso

% load QSO model from training release
variables_to_load = {'rest_wavelengths', 'mu', 'M'};
load(sprintf('%s/learned_model_%s_norm_%d-%d',...
processed_directory(processed.training_release),processed.training_set_name,...
normalization_min_lambda, normalization_max_lambda),variables_to_load{:});

% load DLA samples from training release
variables_to_load = {'offset_sigma_samples', 'offset_z_samples', 'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples-%s', processed_directory(processed.training_release), processed.training_set_name), ...
     variables_to_load{:});
     
% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask'};
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_load{:});

sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;


test_ind = processed.test_ind;

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);

this_wavelengths    =    all_wavelengths{quasar_ind};
this_wavelengths    =              this_wavelengths';
this_flux           =           all_flux{quasar_ind}; 
this_flux           =                     this_flux';
this_noise_variance = all_noise_variance{quasar_ind};
this_noise_variance =           this_noise_variance';
this_pixel_mask     =     all_pixel_mask{quasar_ind};
this_pixel_mask     =               this_pixel_mask';

% convert to QSO rest frame
catalog = load(sprintf('%s/catalog', processed_directory(release)));
z_qsos = catalog.all_zqso(test_ind);
this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qsos(quasar_ind));
    
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

% interpolate model onto given wavelengths
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');
this_mu = mu_interpolator( this_rest_wavelengths);


min_z_c4s = min_z_c4(this_wavelengths, z_qsos(quasar_ind));
max_z_c4s = max_z_c4(this_wavelengths, z_qsos(quasar_ind));

sample_z_c4 = ...
    min_z_c4s +  ...
    (max_z_c4s - min_z_c4s) * offset_z_samples;


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
ID = all_QSO_ID(test_ind);
this_ID = ID{quasar_ind};
this_systems = ismember(c4_QSO_ID, this_ID);
load(sprintf('%s/filter_flags', processed_directory(release)), ...
     'filter_flags');
p_c4 = processed.p_c4;
sample_log_likelihoods_c4 = processed.sample_log_likelihoods_c4;
if(sum(this_systems)==0 & p_c4(quasar_ind)>0.9 & filter_flags(quasar_ind)==0)
    fprintf('This FP:%d\n',quasar_ind);
    fig = figure('visible', 'off');
    clf();
    hold on
    % plotting the sample_log_likelihoods_dla in z_dla and log_nhi axes
    % subplot('position', [0.05 0.49 0.90 0.45]);
    % hold on
        
    [~, maxind] = nanmax(sample_log_likelihoods_c4(quasar_ind, :));
    map_z_c4    = sample_z_c4(maxind);        
    map_log_c4  = nciv_samples(maxind);
    map_sigma_c4 = sample_sigma_c4(maxind);

    % construct dla_mu_map
    absorption = voigt(padded_wavelengths, map_z_c4, ...
        nciv_samples(maxind), num_lines, map_sigma_c4);
    % Testing Voigt profile 
    
    all_c4_NCIV_test=all_c4_NCIV(test_ind);
    all_c4_Z_test=all_z_c4(test_ind);
    this_c4s = NCIV(this_systems);
    this_Zs  = Z_c4(this_systems);
    matched_Z_c4= this_Zs(abs(this_Zs -map_z_c4)==min(abs(this_Zs-map_z_c4)));
    matched_N_c4= this_c4s(abs(this_Zs -map_z_c4)==min(abs(this_Zs-map_z_c4)));
    % absorption = voigt(padded_wavelengths, this_Zs(1), ...
        % 10^this_c4s(1), num_lines, 1661700.9357883865);
        
    absorption = absorption(ind);
    c4_mu = this_mu.* absorption;


    this_z_c4 = (this_wavelengths / 1550) - 1;

    % xlim([]);
    % ylim([-0.1   5]);
    xlabel('(observed wavelengths $\lambda$ (\AA) / 1549 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
    ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
    % if( p_c4(quasar_ind)>0.7)
        plot(this_z_c4, c4_mu, 'Color', 'b');
    
    hold on
    plot(this_z_c4, this_flux, 'Color', 'r');
    hold on

    plot( (1 + z_qsos(quasar_ind)) * (rest_wavelengths / 1549) - 1, mu, 'Color', 'g');



    test_ind_c4 = all_z_c4>0;
    test_ind_c4= test_ind_c4(test_ind);
    tit = sprintf('ID:%s\nMAP(N):%.3f,  N:%.2f, zqso=%.2f\nMAP(z):%.3f, Z:%.2f, sigma:%.2f, P:%.2f',...
                ID{quasar_ind}, log10(map_log_c4),  matched_N_c4, z_qsos(quasar_ind), ...
                map_z_c4, matched_Z_c4, map_sigma_c4/1e5, p_c4(quasar_ind));
    % tit = sprintf('ID:%s\nN:%.2f, zqso=%.2f, Z:%.2f, sigma:%.2f, P:%.2f',...
    %             ID{quasar_ind}, this_c4s(1), z_qsos(quasar_ind), ...
    %             this_Zs(1), map_sigma_c4 p_c4(quasar_ind));
    title(tit);
    % % text(1.7,2.5,tit)
    ylim([-1 3])
    % fid = num2str(quasar_ind);
    dir = sprintf('FP-%s',training_set_name);
    mkdir(dir);
    fid = sprintf('FP-%s/FP-ID%s.pdf', training_set_name,  ID{quasar_ind});
    % saveas(fig, fid, 'jpg');
    % fidpdf = sprintf('FP/%s.pdf',fid);
    exportgraphics(fig, fid,'ContentType','vector')
end