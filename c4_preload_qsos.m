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
    mpf = all_QSO_ID{i}; % mpf is the string id of qso MJD-Plate-FiberID
    all_mjds(i) = str2double(mpf(1:5));
    all_plates(i)      = str2double(mpf(7:10));
    all_fiber_ids(i) = str2double(mpf(12:end));
end

for i = 1:num_quasars


  if (filter_flags(i)~=0)
    continue;
  end
  
  
  [this_wavelengths, this_flux, this_noise_variance, this_pixel_mask] ...
      = file_loader(all_mjds(i), all_plates(i),all_fiber_ids(i));
	% Here file_loader uses dr7 spectrum reader function and given mpf to read
	% spectrum 
      
  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, all_zqso(i));

  % normalize flux
  
  ind = (this_rest_wavelengths >= normalization_min_lambda) & ...
        (this_rest_wavelengths <= normalization_max_lambda) & ...
        (~this_pixel_mask);

  this_median = nanmedian(this_flux(ind));
  
  % bit 2: cannot normalize (all normalizing pixels are masked)
  if (isnan(this_median))
    filter_flags(i) = bitset(filter_flags(i), 3, true);
    continue;
  end
  
  ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda) & ...
        (~this_pixel_mask);
        
  % bit 3: not enough pixels available
  if (nnz(ind) < min_num_pixels)
    filter_flags(i) = bitset(filter_flags(i), 4, true);
    continue;
  end

  all_normalizers(i) = this_median;

  this_flux           = this_flux           / this_median;
  this_noise_variance = this_noise_variance / this_median^2;
 
  all_wavelengths{i}    =    this_wavelengths;%no longer (ind)
  all_flux{i}           =           this_flux;
  all_noise_variance{i} = this_noise_variance;
  all_pixel_mask{i}     =     this_pixel_mask;

  fprintf('loaded quasar %i of %i (%i/%i/%04i)\n', ...
          i, num_quasars, all_plates(i), all_mjds(i), all_fiber_ids(i));
end

variables_to_save = {'loading_min_lambda', 'loading_max_lambda', ...
                     'normalization_min_lambda', 'normalization_max_lambda', ...
                     'min_num_pixels', 'all_wavelengths', 'all_flux', ...
                     'all_noise_variance', 'all_pixel_mask', ...
                     'all_normalizers'};
save(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_save{:});

% write new filter flags to catalog
save(sprintf('%s/filter_flags', processed_directory(release)), ...
     'filter_flags');
