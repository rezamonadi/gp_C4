% set_parameters: sets various parameters for the CIV detection 
% pipeline
% Desined for using DR16 spectra  
%flags for changes
extrapolate_subdla = 0; %0 = off, 1 = on
add_proximity_zone = 0;
integrate          = 1;
optTag = [num2str(integrate), num2str(extrapolate_subdla), num2str(add_proximity_zone)];

% physical constants
lya_wavelength = 1215.6701;                   % Lyman alpha transition wavelength  Å
lyb_wavelength = 1025.7223;                   % Lyman beta  transition wavelength  Å
lyman_limit    =  911.7633;                   % Lyman limit wavelength             Å
civ_1548_wavelength = 1548.2049;		 % CIV transition wavelength  Å
civ_1550_wavelength =  1550.77845; 		 % CIV transition wavelength  Å
speed_of_light = 299792458;                   % speed of light                     m s⁻¹

% converts relative velocity in km s^-1 to redshift difference
kms_to_z = @(kms) (kms * 1000) / speed_of_light;

% utility functions for redshifting
emitted_wavelengths = ...
    @(observed_wavelengths, z) (observed_wavelengths / (1 + z));

observed_wavelengths = ...
    @(emitted_wavelengths,  z) ( emitted_wavelengths * (1 + z));

release = 'dr7';
% download Cooksey's dr7 spectra from this page: 
% http://www.guavanator.uhh.hawaii.edu/~kcooksey/SDSS/CIV/index.html 
% go to table: "SDSS spectra of the sightlines surveyed for C IV."
file_loader = @(mjd, plate, fiber_id) ...
  (read_spec_dr7(sprintf('data/dr7/spectro/1d_26/%04i/1d/spSpec-%05i-%04i-%03i.fit',...
  plate, mjd, plate, fiber_id)));
training_release  = 'dr7';
training_set_name = 'Cooskey_all_qso_catalog';

% file loading parameters
loading_min_lambda = 1350;          % range of rest wavelengths to load  Å
loading_max_lambda = 1570;                    
% The maximum allowed is set so that even if the peak is redshifted off the end, the
% quasar still has data in the range

% preprocessing parameters
%z_qso_cut      = 2.15;                   % filter out QSOs with z less than this threshold
z_qso_cut      = 1.7;                      % according to Cooksey z>1.7                      
z_qso_training_max_cut = 4.5;
%z_qso_training_max_cut = 5;                   % roughly 95% of training data occurs before this redshift; assuming for normalization purposes (move to set_parameters when pleased)
min_num_pixels = 400;                         % minimum number of non-masked pixels

% normalization parameters
% I use 1216 is basically because I want integer in my saved filenames%
%normalization_min_lambda = 1216 - 40;              % range of rest wavelengths to use   Å
normalization_min_lambda = 1310; 
%normalization_max_lambda = 1216 + 40;              %   for flux normalization
normalization_max_lambda = 1325; 
% null model parameters
min_lambda         =  1216;                    % range of rest wavelengths to       Å
max_lambda         = 1600;                    %   model
dlambda            = 0.25;                    % separation of wavelength grid      Å
k                  = 20;                      % rank of non-diagonal contribution
max_noise_variance = 4^2;                     % maximum pixel noise allowed during model training

% optimization parameters
minFunc_options =               ...           % optimization options for model fitting
    struct('MaxIter',     10000, ...
           'MaxFunEvals', 10000);

% C4 model parameters: parameter samples (for Quasi-Monte Carlo)
num_C4_samples           = 10000;                  % number of parameter samples
alpha                    = 0.9;                    % weight of KDE component in mixture
uniform_min_log_nciv     = 13.0189;                   % range of column density samples    [cm⁻²]
uniform_max_log_nciv     = 16;                   % from uniform distribution
fit_min_log_nciv         = 13.0189;                   % range of column density samples    [cm⁻²]
fit_max_log_nciv         = 15.8;                   % from fit to log PDF
extrapolate_min_log_nciv = 13.0189;               % normalization range for the extrapolated region
min_sigma                = 5e5;                   % cm/s -> b/sqrt(2) -> min Doppler par from Cooksey
max_sigma                = 40e5;                   % cm/s -> b/sqrt(2) -> max Doppler par from Cooksey

% model prior parameters

prior_z_qso_increase = kms_to_z(30000);       % use QSOs with z < (z_QSO + x) for prior

% instrumental broadening parameters
width = 3;                                    % width of Gaussian broadening (# pixels)
pixel_spacing = 1e-4;                         % wavelength spacing of pixels in dex

% DLA model parameters: absorber range and model
num_lines = 2;                                % number of members of CIV series to use

max_z_cut = kms_to_z(3000);                   % max z_DLA = z_QSO - max_z_cut
max_z_c4 = @(wavelengths, z_qso) ...         % determines maximum z_DLA to search
    (max(wavelengths)/civ_1548_wavelength - 1) - max_z_cut;

min_z_cut = kms_to_z(3000);                   % min z_DLA = z_Ly∞ + min_z_cut
min_z_c4 = @(wavelengths, z_qso) ...         % determines minimum z_DLA to search
    max(min(wavelengths) / civ_1548_wavelength - 1,                          ...
        observed_wavelengths(min_lambda, z_qso) / civ_1548_wavelength - 1 + ...
        min_z_cut);

% base directory for all data
base_directory = 'data';
% utility functions for identifying various directories
distfiles_directory = @(release) ...
   sprintf('%s/%s/distfiles', base_directory, release);

% spectra_directory   = 'data/dr16/spectra';

processed_directory = @(release) ...
   sprintf('%s/%s/processed', base_directory, release);

c4_catalog_directory = @(name) ...
   sprintf('%s/C4_catalogs/%s/processed', base_directory, name);

   
% replace with @(varargin) (fprintf(varargin{:})) to show debug statements
% fprintf_debug = @(varargin) (fprintf(varargin{:}));
% fprintf_debug = @(varargin) ([]);

