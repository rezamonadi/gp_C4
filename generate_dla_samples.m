% generate_dla_samples: generates DLA parameter samples from training
% catalog

% load training catalog
catalog = load(sprintf('%s/catalog', processed_directory(training_release)));

% generate quasirandom samples from p(normalized offset, log₁₀(N_CIV))
rng('default');
sequence = scramble(haltonset(2), 'rr2');

% the first dimension can be used directly for the uniform prior over
% offsets
offset_samples  = sequence(1:num_c4_samples, 1)';

% we must transform the second dimension to have the correct marginal
% distribution for our chosen prior over column density, wc4ch is a
% mixture of a uniform distribution on log₁₀ N_c4 and a distribution
% we fit to observed data

% uniform component of column density prior
u = makedist('uniform', ...
             'lower', uniform_min_log_nc4, ...
             'upper', uniform_max_log_nc4);

% extract observed log₁₀ N_c4 samples from catalog
all_log_nc4s = catalog.log_nc4s(dla_catalog_name);
ind = cellfun(@(x) (~isempty(x)), all_log_nc4s);
log_nc4s = cat(1, all_log_nc4s{ind});

% make a quadratic fit to the estimated log p(log₁₀ N_c4) over the
% specified range
x = linspace(fit_min_log_nc4, fit_max_log_nc4, 1e3);
kde_pdf = ksdensity(log_nc4s, x);
f = polyfit(x, log(kde_pdf), 2);

% convert this to a PDF and normalize
unnormalized_pdf = @(nc4) (exp(polyval(f, nc4)));
Z = integral(unnormalized_pdf, fit_min_log_nc4, 25.0);

% create the PDF of the mixture between the unifrom distribution and
% the distribution fit to the data
normalized_pdf = @(nc4) ...
          alpha  * (unnormalized_pdf(nc4) / Z) + ...
     (1 - alpha) * (pdf(u, nc4));

cdf = @(nc4) (integral(normalized_pdf, fit_min_log_nc4, nc4));

% use inverse transform sampling to convert the quasirandom samples on
% [0, 1] to appropriate values
log_nc4_samples = zeros(1, num_dla_samples);
for i = 1:num_dla_samples
  log_nc4_samples(i) = ...
      fzero(@(nc4) (cdf(nc4) - sequence(i, 2)), 20.5);
end

% precompute N_c4 samples for convenience
nc4_samples = 10.^log_nc4_samples;

variables_to_save = {'uniform_min_log_nc4', 'uniform_max_log_nc4', ...
                     'fit_min_log_nc4', 'fit_max_log_nc4', 'alpha', ...
                     'offset_samples', 'log_nc4_samples', 'nc4_samples'};
save(sprintf('%s/dla_samples', processed_directory(training_release)), ...
     variables_to_save{:}, '-v7.3');
