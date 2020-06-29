% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
variables_to_load = {'all_QSO_ID','all_zqso', 'all_bal_flags'};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});

num_quasars = numel(all_zqso);

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);

% building MJD, Plate, and Fiber_ID from all_QSO_ID
all_plates    = zeros(num_quasars,1);
all_mjds      = zeros(num_quasars,1);
all_fiber_ids = zeros(num_quasars,1);
for i=1:num_quasars
    mpf = all_QSO_ID{i};
    all_mjds(i) = str2double(mpf(1:5));
    all_plates(i)      = str2double(mpf(7:10));
    all_fiber_ids(i) = str2double(mpf(12:end));
end
% 
% for i = 1:num_quasars
% %   if (filter_flags(i) > 0)bal
% %     continue;
% %   end
%   if (all_bal_flags(i)>0)
%     continue;
%   end
%   
%   [this_wavelengths, this_flux, this_noise_variance, this_pixel_mask] ...
%       = file_loader(all_mjds(i), all_plates(i),all_fiber_ids(i));
% 
%   % do not normalize flux: this is done in the learning and processing code.
% 
%   ind = (this_wavelengths >= (min_lambda * (z_qso_cut) + 1)) & ...
%         (this_wavelengths <= (max_lambda * (z_qso_training_max_cut) + 1)) & ...
%         (~this_pixel_mask);
% 
% %   % bit 3: not enough pixels available
% %   if (nnz(ind) < min_num_pixels)
% %     filter_flags(i) = bitset(filter_flags(i), 4, true);
% %     continue;
% %   end
% 
%   all_wavelengths{i}    =    this_wavelengths;%no longer (ind)
%   all_flux{i}           =           this_flux;
%   all_noise_variance{i} = this_noise_variance;
%   all_pixel_mask{i}     =     this_pixel_mask;
% 
%   fprintf('loaded quasar %i of %i (%i/%i/%04i)\n', ...
%           i, num_quasars, all_plates(i), all_mjds(i), all_fiber_ids(i));
% end

variables_to_save = {'loading_min_lambda', 'loading_max_lambda', ...
                     'normalization_min_lambda', 'normalization_max_lambda', ...
                     'min_num_pixels', 'all_wavelengths', 'all_flux', ...
                     'all_noise_variance', 'all_pixel_mask', ...
                     'all_normalizers'};
save(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_save{:}, '-v7.3');

% write new filter flags to catalog
save(sprintf('%s/catalog', processed_directory(release)), ...
     'filter_flags', '-append');
