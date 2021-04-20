% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest

% load QSO catalog
tic;
variables_to_load = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7',...
'all_QSO_ID','all_zqso'};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});

num_quasars = numel(all_zqso);

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);

for i = 1:num_quasars


  if (filter_flags(i)~=0)
    continue;
  end
  
  [this_wavelengths, this_flux, this_noise_variance, this_pixel_mask] ...
      = file_loader(all_mjd_dr7(i), all_plate_dr7(i),all_fiber_dr7(i));
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
 
  % add one pixel on either side
  available_ind = find(~ind & ~this_pixel_mask);
  ind(min(available_ind(available_ind > find(ind, 1, 'last' )))) = true;
  ind(max(available_ind(available_ind < find(ind, 1, 'first')))) = true;

  all_wavelengths{i}    =    this_wavelengths(ind);%no longer (ind)
  all_flux{i}           =           this_flux(ind);
  all_noise_variance{i} = this_noise_variance(ind);
  all_pixel_mask{i}     =     this_pixel_mask(ind);

  fprintf('loaded quasar %i of %i (%i/%i/%04i)\n', ...
          i, num_quasars, all_plate_dr7(i), all_mjd_dr7(i), all_fiber_dr7(i));
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
toc;
