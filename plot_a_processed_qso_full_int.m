% plot_a_processed_qso_full_int.m : plot sample of spectrum with model civ

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
        s.MarkerFaceAlpha = 0.25;
        colorbar('southoutside');
        title(sprintf('thingID = %d, zQSO = %.2f', selected_thing_ids, z_qso), 'FontSize', 20, 'Interpreter','latex');
        xlim([min(sample_z_dlas) max(sample_z_dlas)]);
        xlabel('$z_{DLA1}$', 'FontSize', 20, 'Interpreter','latex');
        ylabel('$\log N_{HI1}$', 'FontSize', 20, 'Interpreter','latex');
        caxis([-100 0]);
    %     num_dlas = find( model_posteriors(quasar_ind, :) >= 0.5) - 1 - 1; % -1 for null model -1 for sub_dlas

    %     if num_dlas >= 1
    %         roman_z_dlas   = MAP_z_dlas(quasar_ind,   num_dlas, 1:num_dlas);
    %         roman_log_nhis = MAP_log_nhis(quasar_ind, num_dlas, 1:num_dlas);

    %         for i = 1:numel(roman_z_dlas)
    %             s_roman = scatter(roman_z_dlas(i), roman_log_nhis(i), 80, 'b', 'filled', 'd', 'DisplayName', 'MAP predictions');
    %         end

    %         legend( [s s_roman], {'$p(\mathbf{y} \mid z_{DLA}, \log N_{HI}, \mathcal{M}_{DLA})$', 'This work'},...
    %             'Interpreter','latex', 'FontSize', 20);
    %     end
    % hold off

    % if num_dlas >= 1
    %     map_inds =  MAP_inds(quasar_ind, num_dlas, 1:num_dlas);

    %     % construct dla_mu_map
    %     i = map_inds(1);
    %     absorption = voigt(padded_wavelengths, sample_z_dlas(i), ...
    %         nhi_samples(i), num_lines);
    %     for j = 2:num_dlas
    %         i = map_inds(j);
    %         absorption = absorption .* ...
    %             voigt(padded_wavelengths, sample_z_dlas(i), ...
    %                 nhi_samples(i), num_lines);
    %     end

    %     absorption = absorption(mask_ind);
    %     % [beyond lya] set Lya absorption to 1 if beyond lya
    %     indicator  = this_rest_wavelengths > lya_wavelength;
    %     absorption(indicator) = 1;

    %     lya_absorption =  all_lya_absorption(:, i);

    %     dla_mu     = this_mu     .* absorption .* lya_absorption;
    %     dla_M      = this_M      .* absorption .* lya_absorption;
    %     dla_omega2 = this_omega2 .* absorption .* lya_absorption.^2;          

    % else
    %     [~, i] = max(sample_kim_log_likelihoods(quasar_ind, :));

    %     lya_absorption =  all_lya_absorption(:, i);

    %     dla_mu     = this_mu     .* lya_absorption;
    %     dla_M      = this_M      .* lya_absorption;
    %     dla_omega2 = this_omega2 .* lya_absorption.^2;          
    % end

    subplot('position', [0.05 0.06 0.90 0.38]);
    
    hold on
        this_z_dlas = (this_wavelengths / lya_wavelength) - 1;
        p_flux = plot(this_z_dlas, this_flux);
        xlim([min(sample_z_dlas) max(sample_z_dlas)]);
        ylim([-1     5]);
        xlabel('(Observed Wavelengths $\lambda$ (\AA) / 1216 (\AA)) - 1', 'FontSize', 20, 'Interpreter','latex');
        ylabel('Normalized Flux $\mathbf{y}$',                            'FontSize', 20, 'Interpreter','latex');
        p_dla = plot(this_z_dlas, dla_mu);

        plot( (1 + z_qso) * (rest_wavelengths / lya_wavelength) - 1, mu);
    hold off
hold off