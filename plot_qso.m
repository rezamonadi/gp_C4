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

figure(1);
clf();
hold on
    %% plotting the sample_log_likelihoods_dla in z_dla and log_nhi axes
    subplot('position', [0.05 0.49 0.90 0.45]);
    hold on
        norm_sample_log_likelihoods = this_sample_log_likelihoods_dla(:, 1);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - max(norm_sample_log_likelihoods);
        norm_sample_log_likelihoods = norm_sample_log_likelihoods - log(sum(exp(norm_sample_log_likelihoods)));

        s = scatter(sample_z_dlas, log_nhi_samples, 40, norm_sample_log_likelihoods, 'filled', 'DisplayName', 'sample likelihoods');
        s.MarkerFaceAlpha = 0.75;
        colorbar('southoutside');
        title(sprintf('quasar ind = %d, z qso = %.2f', quasar_ind, z_qso)); 
        xlim([min(sample_z_dlas) max(sample_z_dlas)]);
        xlabel('z_{dla1}');
        ylabel('log NHI1');
        caxis([-100 0]);
        num_dlas = find( model_posteriors(quasar_ind, :) >= 0.5) - 1 - 1; % -1 for null model -1 for sub_dlas

        if num_dlas >= 1
            roman_z_dlas   = MAP_z_dlas(quasar_ind,   num_dlas, 1:num_dlas);
            roman_log_nhis = MAP_log_nhis(quasar_ind, num_dlas, 1:num_dlas);

            for i = 1:numel(roman_z_dlas)
                s_roman = scatter(roman_z_dlas(i), roman_log_nhis(i), 80, 'b', 'filled', 'd', 'DisplayName', 'MAP predictions');
            end

            legend( [s s_roman], {'$p(\mathbf{y} \mid z_{DLA}, \log N_{HI}, \mathcal{M}_{DLA})$', 'This work'},...
                'Interpreter','latex', 'FontSize', 14);
        end
    hold off
    
    if num_dlas >= 1
        map_inds =  MAP_inds(quasar_ind, num_dlas, 1:num_dlas);

        % construct dla_mu_map
        i = map_inds(1);
        absorption = voigt(padded_wavelengths, sample_z_dlas(i), ...
            nhi_samples(i), num_lines);
        for j = 2:num_dlas
            i = map_inds(j);
            absorption = absorption .* ...
                voigt(padded_wavelengths, sample_z_dlas(i), ...
                    nhi_samples(i), num_lines);
        end

        absorption = absorption(mask_ind);
        % [beyond lya] set Lya absorption to 1 if beyond lya
        indicator  = this_rest_wavelengths > lya_wavelength;
        absorption(indicator) = 1;

        dla_mu = this_mu .* absorption;
    
    else
        dla_mu = this_mu;
    end

    subplot('position', [0.05 0.05 0.90 0.38]);
    
    hold on
        this_z_dlas = (this_wavelengths / lya_wavelength) - 1;
        p_flux = plot(this_z_dlas, this_flux);
        xlim([min(sample_z_dlas) max(sample_z_dlas)]);
        ylim([-1     5]);
        xlabel('(observed wavelengths $\lambda$ (\AA) / 1216 (\AA)) - 1', 'FontSize', 14, 'Interpreter','latex');
        ylabel('normalized flux $\mathbf{y}$',                            'FontSize', 14, 'Interpreter','latex');
        p_dla = plot(this_z_dlas, dla_mu);

        plot( (1 + z_qso) * (rest_wavelengths / lya_wavelength) - 1, mu);
    hold off
hold off

% plot the full spectrum
% predict the continuum
% select y1, y2, mu1, mu2, M1, M2, d1, d2
% 1: Hydrogen absorption region
% 2: Metal-line region
ind_1    = this_rest_wavelengths <= lya_wavelength;
y2       = this_flux(~ind_1);
x1       = this_rest_wavelengths(ind_1);
this_mu1 = this_mu( ind_1);
this_mu2 = this_mu(~ind_1);
this_M1  = this_M( ind_1, :);
this_M2  = this_M(~ind_1, :);
d1       = this_noise_variance( ind_1) + this_omega2( ind_1);
d2       = this_noise_variance(~ind_1) + this_omega2(~ind_1);

addpath continuum/
[mu1, Sigma11] = conditional_mvnpdf_low_rank(y2, ...
    this_mu1, this_mu2, this_M1, this_M2, d1, d2);

figure(2);
hold on
    f     = plot(this_rest_wavelengths, this_flux);
    dm    = plot(this_rest_wavelengths, dla_mu);
    tm    = plot(this_rest_wavelengths, this_mu);
    conmu = plot(x1, mu1, 'Color', '#4DBEEE');
    ylim([-1 5])
    xlabel('Rest-wavelengths (\AA)', 'FontSize', 14, 'Interpreter','latex');
    ylabel('Normalized Flux $\mathbf{y}$', 'FontSize', 14, 'Interpreter','latex');

    legend( [f dm tm conmu], {'observed flux', '$\mu$(DLA)', '$\mu$(Null)', '$\mu$(Continuum)'},...
    'Interpreter','latex', 'FontSize', 14);
hold off

figure(3)
hold on

    norm_sample_kim_log_likelihoods = sample_kim_log_likelihoods(1, :);
    norm_sample_kim_log_likelihoods = norm_sample_kim_log_likelihoods - max(norm_sample_kim_log_likelihoods);
    norm_sample_kim_log_likelihoods = norm_sample_kim_log_likelihoods - log(sum(exp(norm_sample_kim_log_likelihoods)));

    % scatter(tau_0_samples, beta_samples, 40, norm_sample_kim_log_likelihoods, 'filled');

    % s = scatter(tau_0_map, beta_map, 80, 'b', 'filled', 'd', 'DisplayName', 'MAP predictions');

    scatter(tau_0_samples, norm_sample_kim_log_likelihoods, 40, 'filled');

    s = scatter(tau_0_map, max(norm_sample_kim_log_likelihoods), 80, 'b', 'filled', 'd', 'DisplayName', 'MAP predictions');

    xlim([min(tau_0_samples) max(tau_0_samples)])
    xlabel('$\tau_o$', 'FontSize', 14, 'Interpreter','latex');
    ylabel('Normalised Likelihoods', 'FontSize', 14, 'Interpreter','latex');

    legend( [s], {'MAP effective optical depth'},...
    'Interpreter','latex', 'FontSize', 14);
hold off