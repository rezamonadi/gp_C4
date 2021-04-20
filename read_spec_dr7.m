% in the same units as the spectrum),
% the forth row is the mask array.
%
% The spectra are binned log-linear. Units are 10^(-17) erg/cm^2/s/Ang
% http://classic.sdss.org/dr7/products/spectra/read_spSpec.html
% Notice that the wavelength vector is not contained in the image,
% but must be generated from parameters in the header. SDSS spectra
% are binned in constant Log(Λ) and the wavelength can be obtained
% from the header parameters COEFF0 and COEFF1 (or alternatively
% CRVAL1 and CD1_1) as follows:
%
% Λ = 10**(COEFF0 + COEFF1*i), where i denotes the (zero indexed) pixel number.
%
% Remember that these are vacuum wavelengths.
%

function [wavelengths, flux, noise_variance, pixel_mask] = read_spec_dr7(filename)

	measurements = fitsread(filename);
	%           'binarytable',  1, ...
	%           'tablecolumns', 1:4);

	% the Primary table contains spectra

	% acquire un-continuum subtracted spectrum
	noise = measurements(3,:);
	
	flux = measurements(1,:);

	% noise (standard deviation)

	noise_ratio = noise./flux;
	noise_variance = noise.^2;
	% mask array
	and_mask = measurements(4,:);

	% acquire rest wavelength
	header = fitsinfo(filename);
	s = size(header.PrimaryData.Keywords);
	for i=1:s(1)
	    if (size(header.PrimaryData.Keywords{i,1})==[1,6])
		if (header.PrimaryData.Keywords{i,1}=='COEFF0')
		    coeff0=header.PrimaryData.Keywords{i,2};
		end
		if (header.PrimaryData.Keywords{i,1}=='COEFF1')
		    coeff1=header.PrimaryData.Keywords{i,2};
		    
		end
	    end
	end
	%coef0 = header.PrimaryData.Keywords{208,2} % Center wavelength (log10) of first pi
	%coef1 = header.PrimaryData.Keywords{209,2} % Log10 dispersion per pixel

	length = numel(flux);

	wavelengths = 10.^(linspace(coeff0, coeff0 + coeff1*(length + 1), length));
	
	pixel_mask =  (noise_ratio>=4) | (and_mask==hex2dec('0x800000')) | (noise == 0); 
end
