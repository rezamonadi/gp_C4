'''
read DR7 SPspec

fits description:
http://classic.sdss.org/dr7/dm/flatFiles/spSpec.html

Quoted
----
The spectrum.
The first row is the spectrum,
the second row is the continuum subtracted spectrum,
the third row is the noise in the spectrum (standard deviation,
in the same units as the spectrum),
the forth row is the mask array.

The spectra are binned log-linear. Units are 10^(-17) erg/cm^2/s/Ang
http://classic.sdss.org/dr7/products/spectra/read_spSpec.html
Notice that the wavelength vector is not contained in the image,
but must be generated from parameters in the header. SDSS spectra
are binned in constant Log(Λ) and the wavelength can be obtained
from the header parameters COEFF0 and COEFF1 (or alternatively
CRVAL1 and CD1_1) as follows:

Λ = 10(COEFF0 + COEFF1*i), where i denotes the (zero indexed) pixel number.

Remember that these are vacuum wavelengths.
'''
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

filename = "spSpec-51630-0266-053.fit"

hdu = fits.open(filename)

# the Primary table contains spectra
data = hdu[0].data

# acquire un-continuum subtracted spectrum
flux = data[0, :]

# noise (standard deviation)
noise_variance = data[2, :] ** 2

# mask array
and_mask = data[3, :]

# acquire rest wavelength
coef0 = hdu[0].header['COEFF0'] # Center wavelength (log10) of first pi
coef1 = hdu[0].header['COEFF1'] # Log10 dispersion per pixel

length = len(flux)

wavelengths = 10**(
    np.linspace(coef0, coef0 + coef1 * (length + 1), length) )

# plot it
plt.plot(wavelengths, flux, label="flux")
plt.plot(wavelengths, noise_variance, label="noise")
