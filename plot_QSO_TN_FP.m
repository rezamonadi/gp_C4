
% plot_a_processed_qso.m : plot sample of spectrum with model dla
% and the sample likelihoods
%
% Usage:
% ----
% % first load the catalogues and spectra
% addpath dr16q/
% load_processed_catalogsset_parameters;
clc
clear
set_parameters;
build_catalog;

variables_to_load = {'training_release', 'training_set_name', ...
    'c4_catalog_name', 'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4', 'sample_log_likelihoods_c4L1', ...
    'sample_log_likelihoods_c4L2','log_likelihoods_c4L1', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L1',...
    'model_posteriors', 'p_no_c4', ...
    'p_L1', 'map_z_c4L1', 'map_N_c4L1', 'map_sigma_c4L1', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2'};

filename = sprintf('%s/processed_qsos_%s', ...
    processed_directory(release), ...
    training_set_name);
processed=matfile(filename);
load('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

% load QSO model from training release
variables_to_load = {'rest_wavelengths', 'mu', 'M'};
load(sprintf('%s/learned_model_%s_norm_%d-%d',...
processed_directory(processed.training_release),processed.training_set_name,...
normalization_min_lambda, normalization_max_lambda),variables_to_load{:});

% load C4 samples from training release
variables_to_load = {'offset_sigma_samples', 'offset_z_samples', 'log_nciv_samples', 'nciv_samples'};
load(sprintf('%s/civ_samples-%s', processed_directory(processed.training_release), processed.training_set_name), ...
     variables_to_load{:});
     
% load preprocessed QSOs
variables_to_load = {'all_wavelengths', 'all_flux', 'all_noise_variance', ...
                     'all_pixel_mask'};
load(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_load{:});

% building samples-> Z and sigma
sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;




test_ind = processed.test_ind;
all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
catalog = load(sprintf('%s/catalog', processed_directory(release)));
z_qsos = catalog.all_zqso(test_ind);
load(sprintf('%s/filter_flags', processed_directory(release)), ...
        'filter_flags');


p_L1 = processed.p_L1;
p_c4s = 1- processed.p_no_c4 -(1e-5)*p_L1;
sample_log_likelihoods_c4L2s = processed.sample_log_likelihoods_c4L2;
map_z_c4L2s = processed.map_z_c4L2;
map_N_c4L2s = processed.map_N_c4L2;
map_sigma_c4L2s = processed.map_sigma_c4L2;

map_z_c4L1s = processed.map_z_c4L1;
map_N_c4L1s = processed.map_N_c4L1;
map_sigma_c4L1s = processed.map_sigma_c4L1;

        % Testing Voigt profile 
% dir = sprintf('TN-Posterior-%s',training_set_name);
dir = sprintf('FP-Posterior-N-sigma-%s',training_set_name);
mkdir(dir);
count =0;
for quasar_ind=1:11000
    if count>100
       break
    end
    % fprintf('quasar_ind:%d\n',quasar_ind);
    this_wavelengths    =    all_wavelengths{quasar_ind};
    this_wavelengths    =              this_wavelengths';
    this_flux           =           all_flux{quasar_ind}; 
    this_flux           =                     this_flux';
    this_noise_variance = all_noise_variance{quasar_ind};
    this_noise_variance =           this_noise_variance';
    this_pixel_mask     =     all_pixel_mask{quasar_ind};
    this_pixel_mask     =               this_pixel_mask';
    % convert to QSO rest frame
  
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
    sample_z_c4 = min_z_c4s + (max_z_c4s - min_z_c4s) * offset_z_samples;

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
    
    % sample_log_likelihoods_c4 = sample_log_likelihoods_c4L2s(quasar_ind);
    num_systems = sum(this_systems);
    if(num_systems==0 & p_c4s(quasar_ind)==1) % FP
        % if(num_systems==0 & p_c4s(quasar_ind)<0.2 & filter_flags(quasar_ind)==0) % TN
        count=count+1;
        fprintf('FP:%d\n',count);
        fig = figure('visible', 'off', 'Position', [0,0,2024,1800]);
        clf();
        % subplot('position', [0.05 0.49 0.90 0.45]);
% construct dla_mu_map
        num_lines=2;
        absorptionL2 = voigt(padded_wavelengths, map_z_c4L2s(quasar_ind), ...
            10^map_N_c4L2s(quasar_ind), num_lines, map_sigma_c4L2s(quasar_ind));
        
        num_lines = 1;
        absorptionL1 = voigt(padded_wavelengths, map_z_c4L1s(quasar_ind), ...
            10^map_N_c4L1s(quasar_ind), num_lines, map_sigma_c4L1s(quasar_ind));
        % Testing Voigt profile 
        
        all_c4_NCIV_test=all_c4_NCIV(test_ind);
        all_c4_Z_test=all_z_c4(test_ind);
        this_c4s = NCIV(this_systems); % cooksey's found C4s
        this_Zs  = Z_c4(this_systems); % Cooksey's found Zs
        % s2 = scatter(this_Zs, this_c4s, 50, 'r', 'filled', 'DisplayName', 'C13');
        % s2.MarkerFaceAlpha=1;
        % s2.Marker= 'X';
        % s2.MarkerEdgeColor='k';
        % legend('GP', 'C13')
        % hold on
        matched_Z_c4L1= this_Zs(abs(this_Zs -map_z_c4L1s(quasar_ind))==min(abs(this_Zs-map_z_c4L1s(quasar_ind))));
        matched_N_c4L1= this_c4s(abs(this_Zs -map_z_c4L1s(quasar_ind))==min(abs(this_Zs-map_z_c4L1s(quasar_ind))));

        matched_Z_c4L2= this_Zs(abs(this_Zs -map_z_c4L2s(quasar_ind))==min(abs(this_Zs-map_z_c4L2s(quasar_ind))));
        matched_N_c4L2= this_c4s(abs(this_Zs -map_z_c4L2s(quasar_ind))==min(abs(this_Zs-map_z_c4L2s(quasar_ind))));
        % absorption = voigt(padded_wavelengths, this_Zs(1), ...
            % 10^this_c4s(1), num_lines, 1661700.9357883865);
            
        absorptionL1 = absorptionL1(ind);
        absorptionL2 = absorptionL2(ind);
        c4_muL1 = this_mu.* absorptionL1;
        c4_muL2 = this_mu.* absorptionL2;
               
        this_z_c4 = (this_wavelengths / 1550) - 1;

        % xlim([]);
        % ylim([-0.1   5]);
        subplot(3,1,1);
        xlabel('(observed wavelengths $\lambda$ (\AA) / 1549 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
        ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
        
        plot(this_z_c4, c4_muL1, 'Color', 'b');
        hold on
        plot(this_z_c4, c4_muL2, 'Color', 'g');
        hold on
        plot(this_z_c4, this_flux, 'Color', 'r');
        hold on
        legend( 'L1', 'L2', 'Raw')
        % plot( (1 + z_qsos(quasar_ind)) * (rest_wavelengths / 1549) - 1, mu, 'Color', 'g');
        % alpha(0.7)


        test_ind_c4 = all_z_c4>0;
        test_ind_c4= test_ind_c4(test_ind);
        % tit = sprintf('ID:%s\nMAP(N):%.3f,  N:%.2f, zqso=%.2f\nMAP(z):%.3f, Z:%.2f, sigma:%.2f, P:%.2f',...
        %             ID{quasar_ind}, map_N_c4L2s(quasar_ind),  matched_N_c4, z_qsos(quasar_ind), ...
        %             map_z_c4s(quasar_ind), matched_Z_c4, map_sigma_c4s(quasar_ind)/1e5, p_c4s(quasar_ind));
        sub_title = sprintf('MAP(NL1):%.2f  MAP(NL2):%.2f\nMAP(ZL1):%.2f,    MAP(ZL2):%.2f\nPL2:%.2f, PL1:%.2f\n',...
           map_N_c4L1s(quasar_ind),...
           map_N_c4L2s(quasar_ind),...
           map_z_c4L1s(quasar_ind),...
           map_z_c4L2s(quasar_ind),...
           p_c4s(quasar_ind), p_L1(quasar_ind));
        % tit = sprintf('ID:%s\nN:%.2f, zqso=%.2f, Z:%.2f, sigma:%.2f, P:%.2f',...
        %             ID{quasar_ind}, this_c4s(1), z_qsos(quasar_ind), ...
        %             this_Zs(1), map_sigma_c4 p_c4(quasar_ind));
           
        [t,s] = title(sprintf('ID:%s Zqso:%.2f', ID{quasar_ind}, z_qsos(quasar_ind)), sub_title); 
        t.FontSize=12;
        s.FontSize=10;
        
        % text(1.7,2.5,tit)
        % ylim([-1 3]);
        % fid = num2str(quasar_ind);
        xlim([min(this_z_c4), max(this_z_c4)])

        xlabel('Z')
        ylabel('Normalized Flux')
        hold on 



        subplot(3,1,2);
        hold on
        norm_sample_log_likelihoods = processed.sample_log_likelihoods_c4L2(quasar_ind, :);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - max(norm_sample_log_likelihoods);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

        s=scatter(sample_sigma_c4, log_nciv_samples,20,...
                     norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        s.MarkerFaceAlpha = 0.4;
        hcb = colorbar('southoutside');
        hcb.Label.String= 'Likelihood L2';
        % title(sprintf('thingID = %d, zQSO = %.2f', selected_thing_ids, z_qso), 'FontSize', 20, 'Interpreter','latex');
        % xlim([min(sample_z_dlas) max(sample_z_dlas)]);
        % xlabel('$z_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        xlabel('$\sigma$', 'FontSize', 20, 'Interpreter','latex');
        ylabel('$\log N_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        % xlim([min(this_z_c4), max(this_z_c4)])
        % caxis([-100 0]);
        hold on
        
        subplot(3,1,3);
        hold on
        norm_sample_log_likelihoods = processed.sample_log_likelihoods_c4L1(quasar_ind, :);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - max(norm_sample_log_likelihoods);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

        s=scatter(sample_sigma_c4, log_nciv_samples,20,...
                     norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        s.MarkerFaceAlpha = 0.4;
        hcb = colorbar('southoutside');
        hcb.Label.String= 'Likelihood L1';
        % title(sprintf('thingID = %d, zQSO = %.2f', selected_thing_ids, z_qso), 'FontSize', 20, 'Interpreter','latex');
        % xlim([min(sample_z_dlas) max(sample_z_dlas)]);
        % xlabel('$z_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        xlabel('$\sigma$', 'FontSize', 20, 'Interpreter','latex');
        ylabel('$\log N_{CIV}$', 'FontSize', 20, 'Interpreter','latex');
        % xlim([min(this_z_c4), max(this_z_c4)])

        % caxis([-100 0]);
     
        
        fid = sprintf('FP-Posterior-N-sigma-%s/FP-ID%s.pdf', training_set_name,  ID{quasar_ind});
        % fid = sprintf('TN-Posterior-%s/TN-ID%s.pdf', training_set_name,  ID{quasar_ind});
        % saveas(fig, fid, 'jpg');
        % fidpdf = sprintf('FP/%s.pdf',fid);
        exportgraphics(fig, fid,'ContentType','vector')
    end
end
