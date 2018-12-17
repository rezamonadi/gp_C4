# -*- coding: utf-8 -*-
"""
Module to compute the CDDF (and dN/dX, Omega_DLA) for DLAs from the catalogue Roman Garnett gave me.

Important Functions:

log_p_nhi(lnhi_min, lnhi_max, filehandle, sampfilehandle):
- Computes the likelihood that each spectrum has log(NHI) within the bin given by lnhi_bin, given the DLA model and the data,
    ie: P(NHI && M_DLA | D)
path_length(min_z, max_z, filehandle):
- Computes the integrated comoving path length, dX. This has a marginal dependence on Omega_m,
but the error on Omega_m now is small enough that it won't matter.

We also want dN/dX, which is the total expected number of DLAs in the catalogue, thus:

SUM( Pr(M_DLA | z_QSO, D) ) / dX

and Omega_DLA which is the total expected NHI of DLAs, so should be:

SUM( P(NHI && M_DLA | D) ) / dX

NaN in any of log_priors_dla, log_likelihoods_no_dla and log_likelihoods_dla means the spectrum was cut for some reason.
(NOT that there is zero probability of a DLA)
This is about half the sample of spectra.
"""

import math
#Complex number
import cmath
import operator
import h5py
import numpy as np
import scipy.integrate as integrate
from scipy.stats import poisson
import matplotlib.pyplot as plt

class DLACatalogue(object):
    """Class to contain the DLA catalogue and hold the files containing the data"""
    def __init__(self, processed_file = "processed_qsos_dr7q.mat", sample_file = "dla_samples.mat", raw_file = "preloaded_qsos_dr7.mat", snrs_file = "snrs_qsos_dr7.mat", snr = -2, lowzcut=False):
        #Should we include the second DLA?
        self.second_dla = False
        #Spectra with a DLA probability below this value are assumed to have p = 0, as an optimization.
        #Can be set as high as 0.1 without changing results much.
        #Can be increased, but never decreased
        self.p_thresh_spec = 5e-2
        #This excludes *samples* whose probability is below this value
        self.p_thresh_sample = 1e-4
        #p value to switch from the Poisson approximation to direct summation.
        #0.25 is the value given in Le Cam 1960. In practice 0.5 seems not terrible.
        self.p_switch = 0.25
        #Exclude spectra closer to the DLA than this, which has fewer DLAs than average.
        self.lowzcut = lowzcut
        self.proximity_zone = 0.1
        self.raw_file = raw_file
        self.processed_file = processed_file
        self.tophat_prior = False

        #Load data from the file
        self.filehandle = h5py.File(processed_file,'r')
        #First load small arrays
        self._z_min = self.filehandle["min_z_dlas"][0]
        self._z_max = self.filehandle["max_z_dlas"][0]
        #Probability of at least one DLA in each spectrum
        self.p_dla = self.filehandle["p_dlas"][0]
        #Probability of exactly two DLAs in each spectrum
        if self.second_dla:
            self.p_dla_2 = self.filehandle["model_posteriors"][2]
        #Index of each spectrum in the file containing the flux: raw_file
        self.real_index = np.where(self.filehandle["test_ind"][0] != 0)[0]
        #umber of bins of dNdX or Omega_DLA to plot per unit z interval
        self.bins_per_z = 6
        #Exclude things which have a low SNR. This is tested to be converged on DR7.
        self.filter_noisy_pixels = False
        self.noise_thresh = 0.5**2
        ff = h5py.File(snrs_file,'r')
        self.snrs = np.array(ff["snrs"])
        if self.filter_noisy_pixels:
            self.pixel_noise = np.array(ff["pixel_noise"])
        ff.close()
        self.set_snr(snr)
        self.do_resample = False
        #This allows us to filter by quasar redshift later
        self.condition = np.ones_like(self.z_min, dtype=np.bool)
        #Now load big arrays.
        #First do the DLA1 likelihoods
        #Load normalization constant for the DLA likelihoods
        self.log_norm_like_cache = {}
        dla_ind = self.filter_dla_spectra(second=False)
        if len(np.shape(self.filehandle["sample_log_likelihoods_dla"])) > 2:
            log_norm_like = self.filehandle["sample_log_likelihoods_dla"][0]
        else:
            log_norm_like = self.filehandle["sample_log_likelihoods_dla"]
        #Normalize by the total likelihood of a DLA in each spectrum, so that sum_spectrum ( like) == 1
        #Each DLA in a spectrum is a different column
        log_dla_like = self.filehandle["log_likelihoods_dla"][0]
        #log_norm_like -= (log_dla_like + np.log(np.shape(self.log_norm_like)[0]))
        for spec in dla_ind[0]:
            self.log_norm_like_cache[spec] = np.array(log_norm_like[:,spec] - (log_dla_like[spec] + np.log(np.shape(log_norm_like)[0])))
        del log_norm_like
        del log_dla_like

        #Now build caches for the DLA2 likelihoods and base_sample values
        if self.second_dla:
            dla_ind_2 = self.filter_dla_spectra(second=True)
            #First the log_likelihood of DLA2
            self.log_norm_like_2_cache = {}
            log_norm_like_2 = self.filehandle["sample_log_likelihoods_dla"][1]
            for spec in dla_ind_2[0]:
                self.log_norm_like_2_cache[spec] = self._do_norm_log_norm_like_2(log_norm_like_2[:,spec], spec)
            del log_norm_like_2
            #Build a cache for the base_sample_ind values we will use
            base_sample_inds = np.array(self.filehandle["base_sample_inds"])
            self.base_sample_inds_cache = {}
            for spec in dla_ind_2[0]:
                self.base_sample_inds_cache[spec] = np.array(base_sample_inds[:,spec])-1
            del base_sample_inds

        #Load samples
        samplefilehandle = h5py.File(sample_file,'r')
        #Get the redshift of each sample
        self.z_offsets = samplefilehandle["offset_samples"][:,0]
        #Get the value of NHI at each sample: we do not want to include samples with a column density below the cut.
        self.lnhi_vals = samplefilehandle["log_nhi_samples"][:,0]
        samplefilehandle.close()

    def resample(self, do_it=True, nspec=0):
        """Generate a new sample (with replacement) of the same size as the original."""
        assert not self.second_dla  #not implemented
        assert not self.filter_noisy_pixels #not implemented
        #z_max, z_min, p_dla, snrs and log_norm_like will now be sampled from the new set.
        self.do_resample = do_it
        #Stop if we aren't resampling
        if not do_it:
            return
        #Get the new sample set
        if nspec == 0:
            nspec = np.size(self.p_dla)
        self._resample = np.empty(nspec,dtype=int)
        #Find the redshift above which there are only 5 DLAs,
        #so that we don't have overly small sized bins
        newmax = np.max(self._z_max) - 0.2
        while np.sum(self._z_max > newmax)*nspec/np.size(self.p_dla) < 10:
            newmax -=0.2
        newmin = np.min(self._z_min) + 0.2
        while np.sum(self._z_min > newmin)*nspec/np.size(self.p_dla) < 10:
            newmin +=0.2
        #This extends the last bin over a wider redshift range
        z_bins = np.linspace(newmin, newmax, 10)
        z_bins[-1] = np.max(self._z_max)
        z_bins[0] = np.min(self._z_min)
        #Roughly preserve the redshift distribution of the quasars.
        #Because high redshift quasars are quite rare,
        #if we just resample entirely randomly we could end up with very few of them.
        total = 0
        for zm,zp in zip(z_bins[:-1],z_bins[1:]):
            ii = np.where(np.logical_and(self._z_max > zm,self._z_max <= zp))
            nthisbin = np.min([int(np.floor(np.size(ii)/np.size(self.p_dla)*nspec)),nspec - total])
            assert nthisbin >= 10
            rand = np.random.randint(0,nthisbin, nthisbin)
            self._resample[total:total+nthisbin] = ii[0][rand]
            total += nthisbin
        assert total == nspec

    def get_sample_errors(self, *, z_min=2, z_max=5, nsample=5):
        """Do a number of resamplings to get error bars on omega_dla and dNdX."""
        dndx_sample = []
        om_sample = []
        self.resample(True)
        for _ in range(nsample):
            (_, dNdX, _, _, _) = self.line_density(z_min=z_min, z_max=z_max)
            (_, omega_dla, _, _, _) =  self.omega_dla_cddf(z_min=z_min, z_max=z_max,lnhi_nbins=15.)
            om_sample.append(1000*omega_dla)
            dndx_sample.append(dNdX)
            self.resample(True)
        self.resample(False)
        dndx_sample = np.array(dndx_sample)
        om_sample = np.array(om_sample)
        self.dndx_68_sample = np.array((np.percentile(dndx_sample, 100-32/2,axis=0), np.percentile(dndx_sample, 32/2,axis=0)))
        assert np.shape(self.dndx_68_sample)[1] == np.shape(dNdX)[0]
        self.dndx_95_sample = np.array((np.percentile(dndx_sample, 100-5/2,axis=0), np.percentile(dndx_sample, 5/2,axis=0)))
        self.omega_68_sample = np.array((np.percentile(om_sample, 100-32/2,axis=0), np.percentile(om_sample, 32/2,axis=0)))
        self.omega_95_sample = np.array((np.percentile(om_sample, 100-5/2,axis=0), np.percentile(om_sample, 5/2,axis=0)))
        self.omega_sample = np.median(om_sample, axis=0)
        self.dndx_sample = np.median(dndx_sample, axis=0)

    def plot_dndx_sample_errors(self, *, z_min=2, z_max=5, nsample=5):
        """Plot the sample errors"""
        try:
            self.dndx_68_sample
        except AttributeError:
            self.get_sample_errors(z_min=z_min, z_max=z_max, nsample=nsample)
        (z_cent, dNdX, dndx68, dndx95, xerrs) = self.line_density(z_min=z_min, z_max=z_max)
        plt.fill_between(z_cent_fill(z_cent, xerrs), dndx95[:,0], dndx95[:,1], color="grey", alpha=0.5)
        yerr = (dNdX-dndx68[:,0], dndx68[:,1] - dNdX)
        plt.errorbar(z_cent, dNdX, yerr=yerr, xerr = xerrs, fmt='o', label="Total")
        yerr = (self.dndx_sample-self.dndx_68_sample[0,:], self.dndx_68_sample[1,:] - self.dndx_sample)
        plt.errorbar(z_cent, self.dndx_sample, yerr=yerr, xerr = xerrs, fmt='o', label="Resampled")
        plt.xlabel(r'z')
        plt.ylabel(r'dN/dX')
        plt.xlim(z_min, z_max)

    def plot_omega_sample_errors(self, *, z_min=2, z_max=5, nsample=5):
        """Plot the sample errors"""
        try:
            self.omega_68_sample
        except AttributeError:
            self.get_sample_errors(z_min=z_min, z_max=z_max, nsample=nsample)
        (z_cent, omega_dla, omega68, omega95, xerrs) = self.omega_dla_cddf(z_min=z_min, z_max=z_max)
        plt.fill_between(z_cent_fill(z_cent, xerrs), 1000*omega95[:,0], 1000*omega95[:,1], color="grey", alpha=0.5)
        yerr = (1000*omega_dla-1000*omega68[:,0], 1000*omega68[:,1] - 1000*omega_dla)
        plt.errorbar(z_cent, 1000*omega_dla, yerr=yerr, xerr = xerrs, fmt='o', label="Total")
        yerr = (self.omega_sample-self.omega_68_sample[0,:], self.omega_68_sample[1,:] - self.omega_sample)
        plt.errorbar(z_cent, self.omega_sample, yerr=yerr, xerr = xerrs, fmt='o', label="Resampled")
        plt.xlabel(r'z')
        plt.ylabel(r'$10^3 \times \Omega_\mathrm{DLA}$')
        plt.xlim(z_min, z_max)


    def _base_sample_inds(self, spec):
        """Load the base_sample index to look up NHI for the second DLA, for spectrum spec"""
        try:
            return self.base_sample_inds_cache[spec]
        except KeyError:
            #base_sample_inds starts off zero indexed and needs to be 1-indexed.
            self.base_sample_inds_cache[spec]= np.array(self.filehandle["base_sample_inds"][:,spec])-1
            return self.base_sample_inds_cache[spec]

    def _log_norm_like(self, spec, *, second=False):
        """Get the probability (normalised likelihood) values for the samples in a particular spectrum from the disc"""
        #Loading this from the disc each time is unreasonably slow
        if self.do_resample:
            spec = self._resample[spec]
        if not second:
            try:
                return self.log_norm_like_cache[spec]
            except KeyError:
                if len(np.shape(self.filehandle["sample_log_likelihoods_dla"])) > 2:
                    log_norm_like = self.filehandle["sample_log_likelihoods_dla"][0,:,spec]
                else:
                    log_norm_like = self.filehandle["sample_log_likelihoods_dla"][:,spec]
                #Normalize by the total likelihood of a DLA in each spectrum, so that sum_spectrum ( like) == 1
                #Each DLA in a spectrum is a different column
                log_dla_like = self.filehandle["log_likelihoods_dla"][0,spec]
                log_norm_like -= (log_dla_like + np.log(np.shape(log_norm_like)[0]))
                self.log_norm_like_cache[spec] = log_norm_like
                assert 0.95 < np.sum(np.exp(log_norm_like)) < 1.05
                return log_norm_like
        # Or get for the second DLA:
        # We will want P(DLA @ q = (N,z)) = P(n_DLA >= 1) P(DLA1 @ q) + P(n_DLA == 2) P(DLA2 @ q)
        # and P(DLA2 @ q ) = sum(DLA1 @ q') P(DLA1 @ q' and DLA2 @ q | data ) P(DLA1 @ q' | data)
        #                  = sum(DLA1 @ q') P(data | DLA1 @ q' and DLA2 @ q ) P( data | DLA1 @ q' ) P(DLA2 @ q | DLA1 @ q') P(DLA1 @ q')
        # P( data | DLA1 @ q' )  is sample_log_likelihood_dla[0] P(data | DLA1 @ q' and DLA2 @ q ) is sample_log_likelihood_dla[1].
        # Note that the parameters for DLA1 sample_log_likelihood_dla[1] are the same as sample_log_likelihood_dla[0]
        # So this sum is over all DLA1 samples.
        # then we have P(DLA2 @ q | DLA1 @ q') == P(DLA1 @ q') == 1/Nsample
        # The parameters of DLA1 are sample j are nhi[j], z[j]
        # The parameters of DLA2 are spectrum dependent and given by nhi[base_sample_inds[i,j]], z[base_sample_inds[i, j]]
        # Mask out nan values by making them very low probability: these correspond to samples where the DLAs are too close.
        try:
            return self.log_norm_like_2_cache[spec]
        except KeyError:
            log_nhi_like = self.filehandle["sample_log_likelihoods_dla"][1,:,spec]
            self.log_norm_like_2_cache[spec] = self._do_norm_log_norm_like_2(log_nhi_like, spec)
            return self.log_norm_like_2_cache[spec]

    def _do_norm_log_norm_like_2(self,log_nhi_like, spec):
        """Compute the normalized probabilities for DLA2 samples from the likelihood values for a spectrum."""
        log_nhi_like[np.isnan(log_nhi_like)] = -1e30
        log_norm_like_2 = log_nhi_like + self._log_norm_like(spec, second=False)
        #Normalize so that the sum of these likelihoods is unity.
        #First add something so we don't underflow our floating points.
        #This has the bonus that for peaked distributions, the normalization constant will be basically one already.
        log_norm_like_2 -= np.max(log_norm_like_2)
        norm = np.logaddexp(log_norm_like_2)
        assert np.isfinite(norm)
        log_norm_like_2 -= norm
        return log_norm_like_2

    def filter_dla_spectra(self, *, second=False):
        """
        Find the spectra we are not interested in, because the probability of a DLA is below the desired threshold.
        Or because the SNR is insufficient
        """
        return np.where((self._p_dla(second=second) > self.p_thresh_spec)*self._filter_snr_spectra())

    def _filter_snr_spectra(self):
        """Helper function to get SNR mask."""
        snrs = self.snrs
        if self.do_resample:
            snrs = self.snrs[self._resample]
        return (snrs > self.snr_thresh)*self.condition

    def filter_snr_spectra(self):
        """Remove spectra whose SNR is below snr_thresh"""
        return np.where(self._filter_snr_spectra())

    def set_snr(self, snr_thresh):
        """Set the value of SNR to be used, loading the SNR array if needed"""
        self.snr_thresh = snr_thresh

    def _p_dla(self, *, second=False):
        """Get the probability of a DLA. If second=False, return the probabilities of at least one DLA in each spectrum.
        If second=True, return the probability of exactly two DLAs in each spectrum.
        """
        if not second:
            if self.do_resample:
                return self.p_dla[self._resample]
            return self.p_dla
        else:
            return self.p_dla_2

    def z_max(self, spec=None):
        """Returns the maximum redshift of the quasar spectrum."""
        if spec is None:
            if self.do_resample:
                return self._z_max[self._resample]
            return self._z_max
        else:
            if self.do_resample:
                return self._z_max[self._resample[spec]]
            return self._z_max[spec]

    def z_min(self, spec=None):
        """Returns the minimum redshift of the quasar spectrum."""
        if spec is None:
            if self.do_resample:
                return self._z_min[self._resample]
            return self._z_min
        else:
            if self.do_resample:
                return self._z_min[self._resample[spec]]
            return self._z_min[spec]

    def path_length(self, z_min, z_max):
        """Compute the path length, dX, over which we looked for DLAs.
        Exclude any paths beyond min_z or max_z and any pixels with pixel_noise > thresh^2
        To compute dX use dz from

        max_z_dlas - min_z_dlas

        and
        dX = (1+z)^2 H_0 / H(z) dz
        """
        assert z_min < z_max
        #Make a clean copy
        #Filter spectra that don't make the SNR cut
        ind = self.filter_snr_spectra()
        max_z_dlas = np.array(self.z_max())[ind]
        min_z_dlas = np.array(self.z_min())[ind]
        #Increase the minimum redshift to remove spectra contaminated by the lyman beta forest.
        if self.lowzcut:
            max_z_dlas = np.max([np.min([max_z_dlas, self.proximity(max_z_dlas)],axis=0), min_z_dlas],axis=0)
        assert np.all(max_z_dlas - min_z_dlas >= 0)
        #Filter spectra that aren't in our redshift range
        i2 = np.where(np.logical_and(min_z_dlas < z_max, max_z_dlas > z_min))
        max_z_dlas = max_z_dlas[i2]
        min_z_dlas = min_z_dlas[i2]
        total = 0
        #Shortcut for spectra which cross the whole bin
        whole_bin = np.logical_and(max_z_dlas > z_max, min_z_dlas < z_min)
        #Find spectra where all pixels pass noise cuts
        if self.filter_noisy_pixels:
            pixel_noise = self.pixel_noise[ind][i2]
            no_filters = np.array([np.all(ftrns < self.noise_thresh) for ftrns in pixel_noise])
            whole_bin = np.logical_and(whole_bin, no_filters)
        i3 = np.where(whole_bin)
        (tbin,err) = integrate.quad(path_length_int, z_min, z_max)
        total += np.size(i3) * tbin
        #Integrate only remaining spectra
        i3 = np.where(np.logical_not(whole_bin))
        max_z_dlas = max_z_dlas[i3]
        min_z_dlas = min_z_dlas[i3]
        if not self.filter_noisy_pixels:
            for (zmin, zmax) in zip(min_z_dlas, max_z_dlas):
                assert zmin <= zmax
                #Do the spectra
                pathzmax = np.min([z_max, zmax])
                pathzmin = np.max([z_min, zmin])
                (ans, err) = integrate.quad(path_length_int, pathzmin, pathzmax)
                total += ans
                assert err < 1e-6
        else:
            total += self._do_filtered_path(z_max, z_min, min_z_dlas, max_z_dlas, pixel_noise, no_filters, i3)
        #The total dX for the path length we looked in
        return total

    def _do_filtered_path(self, z_max, z_min, min_z_dlas, max_z_dlas, pixel_noise, no_filters, i3):
        """Compute the path length for spectra where certain pixels have been filtered due to their SNR."""
        total = 0.
        pixel_noise = pixel_noise[i3]
        no_filters = no_filters[i3]
        #Clamp remaining max and min to limits
        #max_z_dlas[np.where(max_z_dlas > z_max)] = z_max
        #min_z_dlas[np.where(min_z_dlas < z_min)] = z_min
        for (zmin, zmax, pn, nf) in zip(min_z_dlas, max_z_dlas, pixel_noise, no_filters):
            assert zmin < zmax
            #Do the spectra that have good noise properties
            pathzmax = np.min([z_max, zmax])
            pathzmin = np.max([z_min, zmin])
            if nf:
                (ans, err) = integrate.quad(path_length_int, pathzmin, pathzmax)
            #Do the others
            else:
                zzs = zmin+(zmax-zmin)*np.arange(np.size(pn))/(np.size(pn)-1)
                #This will contain a list of contiguous regions with good noise properties
                regions = []
                #Find the first pixel within the redshift range which has good noise.
                ii = np.where(np.logical_and(zzs >= pathzmin, pn < self.noise_thresh))
                if np.size(ii) == 0:
                    continue
                ii = ii[0][0]
                #As long as there is more spectrum to look at within our redshift range
                while np.logical_and(ii < np.size(pn)-1, zzs[ii] <= pathzmax):
                    #Find the next pixel which exceeds the noise bound
                    ie = np.where(np.logical_and(pn[ii:] > self.noise_thresh, zzs[ii:] < pathzmax))
                    #If no more pixels exceed the noise bound, exit the loop
                    if np.size(ie) == 0:
                        regions+=[(zzs[ii], pathzmax)]
                        break
                    #If this pixel exists, mark it as the end of the region
                    ie = ie[0][0]+ii
                    regions+=[(zzs[ii], zzs[ie-1])]
                    #Find the start of the next regions with low noise
                    ind = np.where(pn[ie:] < self.noise_thresh)
                    #If it doesn't exist, exit the loop
                    if np.size(ind) == 0:
                        break
                    ii = ind[0][0]+ie
                ans=0
                err=0
                #Do it piecewise: first argument is the start of each bin, second is the end.
                for zrr in regions:
                    (a1, e1) = integrate.quad(path_length_int, zrr[0], zrr[1])
                    ans+=a1
                    err+=e1
            total += ans
            assert err < 1e-6
        return total

    def column_density_function(self, z_min=1., z_max=6., lnhi_nbins=30, lnhi_min=20.,lnhi_max=23.):
        """This computes the column density function, which is the number
            of absorbers per sight line with HI column densities in the interval
            [NHI, NHI+dNHI] at the absorption distance X.

            So we have f(N) = d n_DLA/ dN dX
            and n_DLA(N) = expected number of absorbers per sightline in each column density bin.
            ie, f(N) = n_DLA / ΔN / ΔX
            Note f(N) has dimensions of cm^2, because N has units of cm^-2 and X is dimensionless.
            Returns:
                (NHI, f_N_table) - N_HI (binned in log) and corresponding f(N)
        """
        #Get the NHI bins
        l_nhi = np.linspace(lnhi_min, lnhi_max, num=lnhi_nbins+1)
        #Get the mean and variance of the probability distribution of DLAs.
        (ndlas, l68, l95) = self._get_confidence_intervals(q_bins=l_nhi, lred=z_min, ured=z_max, lnhi_min=lnhi_min, nhi=True)
        dX = self.path_length(z_min, z_max)
        dN = np.array([10**lnhi_x - 10**lnhi_m for (lnhi_m, lnhi_x) in zip(l_nhi[:-1], l_nhi[1:])])
        cddf = np.array(ndlas) / dX / dN
        #Broadcasting failure
        cddf68 = np.array(l68) / dX / np.vstack([dN,dN]).T
        cddf95 = np.array(l95) / dX / np.vstack([dN,dN]).T
        l_Ncent = np.array([(lnhi_x +lnhi_m)/2. for (lnhi_m, lnhi_x) in zip(l_nhi[:-1], l_nhi[1:])])
        xerrs = (10**l_Ncent - 10**l_nhi[:-1],  10**l_nhi[1:] - 10**l_Ncent)
        return (l_Ncent, cddf, cddf68, cddf95, xerrs)

    def plot_cddf(self, zmin=1., zmax=6., label="GP", color=None, moment=False, twosigma=True):
        """Plot the column density function"""
        (l_N, cddf, cddf68, cddf95, xerrs) = self.column_density_function(z_min=zmin, z_max=zmax)
        if moment:
            cddf *= 10**l_N
            for x in (0,1):
                cddf68[:,x] *= 10**l_N
                cddf95[:,x] *= 10**l_N
        #2 sigma contours.
        if twosigma:
            plt.fill_between(10**l_N, cddf95[:,0], cddf95[:,1], color="grey", alpha=0.5)
        yerr = (cddf-cddf68[:,0], cddf68[:,1] - cddf)
        ii = np.where(cddf68[:,0] > 0.)
        if np.size(ii) > 0:
            plt.errorbar(10**l_N[ii], cddf[ii], yerr=(yerr[0][ii],yerr[1][ii]), xerr = (xerrs[0][ii], xerrs[1][ii]), fmt='o', label=label, color=color)
        i2 = np.where(cddf68[:,0] == 0)
        if np.size(i2) > 0:
            plt.errorbar(10**l_N[i2], cddf[i2]+yerr[1][i2], yerr=yerr[1][i2]/2., xerr = (xerrs[0][i2], xerrs[1][i2]), fmt='o', label=None, uplims=True, color=color,lw=2)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel(r'$N_\mathrm{HI}$ (cm$^{-2}$)')
        plt.ylabel(r"$f(N_\mathrm{HI})$")
        return (l_N, cddf, cddf68, cddf95)

    def line_density(self, z_min=2, z_max=4):
        """Compute the line density of DLAs as a function of redshift
        Default bins chosen to match Noterdaeme 2012"""
        #Get the redshifts
        nbins = np.max([int((z_max-z_min)*self.bins_per_z),1])
        z_bins = np.linspace(z_min, z_max, nbins+1)
        #Get the mean and variance of the probability distribution of DLAs.
        (maxlike, l68, l95) = self._get_confidence_intervals(q_bins=z_bins, lred=z_min, ured=z_max, lnhi_min=20.3, nhi=False)
        #Check the outputs are reasonably ordered.
        dX = np.array([self.path_length(z_m, z_x) for (z_m, z_x) in zip(z_bins[:-1], z_bins[1:])])
        ii = np.where(dX > 0)
        dX = dX[ii]
        dNdX = np.array(maxlike)[ii]/dX
        dndx68 = np.array(l68)[ii] / np.vstack([dX,dX]).T
        dndx95 = np.array(l95)[ii] / np.vstack([dX,dX]).T
        z_cent = np.array([(z_x +z_m)/2. for (z_m, z_x) in zip(z_bins[:-1], z_bins[1:])])
        xerrs = (z_cent[ii] - z_bins[:-1][ii],  z_bins[1:][ii] - z_cent[ii])
        return (z_cent[ii], dNdX, dndx68, dndx95, xerrs)

    def plot_line_density(self, zmin=2, zmax=4, label="GP"):
        """Plot the line density as a function of redshift"""
        (z_cent, dNdX, dndx68, dndx95, xerrs) = self.line_density(z_min=zmin, z_max=zmax)
        #2 sigma contours.
        plt.fill_between(z_cent_fill(z_cent,xerrs), dndx95[:,0], dndx95[:,1], color="grey", alpha=0.5)
        yerr = (dNdX-dndx68[:,0], dndx68[:,1] - dNdX)
        plt.errorbar(z_cent, dNdX, yerr=yerr, xerr = xerrs, fmt='o', label=label)
        plt.xlabel(r'z')
        plt.ylabel(r'dN/dX')
        plt.xlim(zmin, zmax)
        return (z_cent, dNdX, dndx68, dndx95)

    def omega_dla_cddf(self, z_min=2, z_max=4, hubble = 0.7, lnhi_nbins=30):
        """
            Compute Omega_dla, the sum of the mass in a given absorber,
            divided by the volume of the spectra, divided by the critical density.
            This is computed by summing the column density function, rather than directly by summing
            columns. Should be the same as omega_dla.

            So we get omega_dla = m_P H_0 / (c rho_c) int dN N f(N)
        """
        nbins = np.max([int((z_max-z_min)*self.bins_per_z),1])
        z_bins = np.linspace(z_min, z_max, nbins+1)
        protonmass=1.67262178e-24
        #H0 in 1/s units
        h100=3.2407789e-18*hubble
        #Speed of light in cm/s
        light = 2.99e10
        omega_dla = np.array([])
        omega_dla_68 = np.array([]).reshape(0,2)
        omega_dla_95 = np.empty_like(omega_dla_68)
        xerrs = np.empty_like(omega_dla_68)
        z_cent = np.array([])
        conversion = protonmass/light*h100/rho_crit(hubble)
        lnhi_bins = np.linspace(20.3, 23, num=lnhi_nbins+1)
        for zz in range(nbins):
            dX = self.path_length(z_bins[zz], z_bins[zz+1])
            if dX == 0.:
                continue
            (nhi_like, nhi_68,nhi_95) = self._get_omega_confidence_intervals(lnhi_bins=lnhi_bins, lred=z_bins[zz], ured=z_bins[zz+1])
            #Check the outputs are reasonably ordered.
            assert nhi_95[0] <= nhi_68[0] <= nhi_like
            assert nhi_95[1] >= nhi_68[1] >= nhi_like
            #The 1+z factor converts lightspeed to comoving
            omega_dla = np.append(omega_dla, conversion*nhi_like / dX)
            omega_dla_68 = np.append(omega_dla_68, conversion*np.array(nhi_68).reshape(1,2) / dX,axis=0)
            omega_dla_95 = np.append(omega_dla_95, conversion*np.array(nhi_95).reshape(1,2) / dX,axis=0)
            z_c = (z_bins[zz]+z_bins[zz+1])/2.
            z_cent = np.append(z_cent, z_c)
            xerrs = np.append(xerrs, np.array([z_c - z_bins[zz],  z_bins[zz+1] - z_c]).reshape(1,2),axis=0)
        assert np.shape(omega_dla_68) == (np.shape(omega_dla)[0],2)
        return (z_cent, omega_dla, omega_dla_68, omega_dla_95, xerrs.T)

    def _get_omega_confidence_intervals(self, lnhi_bins, lred=2., ured=4., tailprob=5e-4):
        """
        Get the confidence interval on the total abundance of HI in DLAs in a given redshift range (this should be called for each bin in Omega_DLA).
        We do this be computing the CDDF in NHI bins and then summing the PDFs for each one.
        Returns: (maximum a posteriori likelihoods, lower 68 % confidence levels, upper 68% confidence levels, lower and upper 95 % confidence levels)
        """
        (probs, poissons) = self._split_distributions(lnhi_bins, lred=lred, ured=ured, lnhi_min=lnhi_bins[0], lnhi_max=lnhi_bins[-1], nhi=True)
        #probs[i] now contains a list of arrays
        #Now we have built a list of probabilities in each z bin of interest and we want to solve for the Poisson binomial coefficients.
        #to get each combined pdf.
        #Empty pdf: P(NHI=0) = 1
        pdf_comb = np.ones(1)
        nhi_comb = np.zeros(1)
        #We could probably get more accuracy by doing some sort of interpolation and then integrating...
        nhi_cent = 10**np.array([(lnhi_x +lnhi_m)/2. for (lnhi_m, lnhi_x) in zip(lnhi_bins[:-1], lnhi_bins[1:])])
        #Loop over bins in the column density function
        for (pp, pmean, nhi_cc) in zip(probs, poissons, nhi_cent):
            pdf = get_poisson_binomial_pdf(pp)
            #Get the pdf for this NHI bin
            (pdf_one_bin, offset_one_bin) = self._get_combined_levels(pdf, pmean)
            #If the last CDDF bin is consistent with zero, stop.
            if self.tophat_prior:
                (lowtest, _) = interval(np.cumsum(pdf_one_bin), 0.68)
                if lowtest < 1:
                    continue
            #Store the PDFs
            (dlow, dhigh) = interval(np.cumsum(pdf_one_bin), 1-1e-4)
            #We want to include dhigh, as long as it is in the array
            maxr = np.min([dhigh+1,np.size(pdf_one_bin)])
            pdf_comb = np.ravel(np.array([[pdf_comb[j]*pdf_one_bin[i] for i in range(dlow,maxr)] for j in range(np.size(pdf_comb))]))
            #Store the NHI values corresponding to each PDF
            nhi_comb = np.ravel(np.array([[nhi_comb[j]+(offset_one_bin+i)*nhi_cc for i in range(dlow,maxr)] for j in range(np.size(nhi_comb))]))
            assert 1.01 > math.fsum(pdf_comb) > 0.99
            #Sort the pdf by increasing NHI
            sort = np.argsort(nhi_comb)
            nhi_comb = nhi_comb[sort]
            pdf_comb = pdf_comb[sort]
            #Now we want to shrink the arrays a little, by combining options within 1% of each other, as well as merging low-probability tails.
            #If we don't do this the array quickly gets out of hand.
            #First do tails.
            cdf = np.cumsum(pdf_comb)
            t1 = np.where(cdf < tailprob)
            t2 = np.where(cdf > 1-tailprob)
            #Replace the last few points in this distribution with the sum of their pdfs
            if np.size(t2) > 0:
                pdf_comb = np.append(pdf_comb[:t2[0][0]], np.sum(pdf_comb[t2]))
                nhi_comb = np.append(nhi_comb[:t2[0][0]], np.min(nhi_comb[t2]))
            if np.size(t1) > 0:
                pdf_comb = np.insert(pdf_comb[t1[0][-1]+1:], 0,np.sum(pdf_comb[t1]))
                nhi_comb = np.insert(nhi_comb[t1[0][-1]+1:], 0,np.max(nhi_comb[t1]))
            assert 1.01 > math.fsum(pdf_comb) > 0.99
            #Now find options which are indistinguishable for all reasonable purposes
            new_pdf = [pdf_comb[0],]
            new_nhi = [nhi_comb[0],]
            #Here we need a 'real' for loop
            low_ind = 1
            while low_ind < np.size(pdf_comb)-1:
                i3 = np.where(np.logical_and(nhi_comb[low_ind:-1]/nhi_comb[low_ind] < 1+1e-3, np.cumsum(pdf_comb[low_ind:-1]) < pdf_comb[low_ind:-1]+0.04) )
                new_pdf.append(math.fsum(pdf_comb[low_ind:-1][i3]))
                new_nhi.append(np.median(nhi_comb[low_ind:-1][i3]))
                low_ind += i3[0][-1]+1
            #Add the last sample
            if np.size(pdf_comb) > 1:
                new_pdf.append(pdf_comb[-1])
                new_nhi.append(nhi_comb[-1])
            assert np.size(new_pdf) == np.size(new_nhi)
            assert np.abs(math.fsum(new_pdf) -  math.fsum(pdf_comb)) < 1e-4
            pdf_comb = np.array(new_pdf)
            nhi_comb = np.array(new_nhi)
        #Unpack maximum likelihoods and 68/95% contours
        (maxlikes, levels68, levels95) = pdf_confidence(pdf_comb, 0)
        #Edge case
        if levels95[1] >= np.size(nhi_comb):
            levels95=(levels95[0], levels95[1]-1)
        return (nhi_comb[maxlikes], (nhi_comb[levels68[0]], nhi_comb[levels68[1]]), (nhi_comb[levels95[0]], nhi_comb[levels95[1]]))

    def omega_dla(self, z_min=2, z_max=4, hubble=0.7, lnhi_max=23., lnhi_min=20.3):
        """
        Compute the matter density of DLAs as a function of redshift, by summing DLAs.
        This gives us:
            Omega_DLA = m_P H_0 / (c rho_c) * sum(NHI) / dX
        """
        #Get the redshifts
        nbins = np.max([int((z_max-z_min)*self.bins_per_z),1])
        z_bins = np.linspace(z_min, z_max, nbins+1)
        #Get the mean and variance of the probability distribution of DLAs.
        (mean, variance) = self._get_z_nhi_hist(q_bins=z_bins, lred=z_min, ured=z_max, lnhi_min=lnhi_min, lnhi_max=lnhi_max, nhi=False, moment=True)
        #This returns the total matter in DLAs at each redshift in atoms/cm^2.
        #Need to turn this into g/cm^2, divide by path length in (comoving) cm, and then divide by rho_crit.
        #proton mass in g
        protonmass=1.67262178e-24
        dX = np.array([self.path_length(z_m, z_x) for (z_m, z_x) in zip(z_bins[:-1], z_bins[1:])])
        #H0 in 1/s units
        h100=3.2407789e-18*hubble
        #Speed of light in cm/s
        light = 2.99e10
        conversion = protonmass * h100 / light / dX / rho_crit()
        omega_DLA = mean * conversion
        err = np.sqrt(variance) * conversion
        z_cent = np.array([(z_x +z_m)/2. for (z_m, z_x) in zip(z_bins[:-1], z_bins[1:])])
        return (z_cent, omega_DLA, err, z_bins)

    def plot_omega_dla_var(self, zmin=2, zmax=4, label="GP", color=None):
        """Plot omega_DLA as a function of redshift, with errors given by (an approximation to) the distribution variance"""
        (z_cent, omega_DLA, err, z_bins) = self.omega_dla(z_min=zmin, z_max=zmax)
        xerrs = (z_cent - z_bins[:-1],  z_bins[1:] - z_cent)
        plt.errorbar(z_cent, 1000*omega_DLA, yerr=1000*err, xerr = xerrs, fmt='s', label=label,color=color)
        plt.xlabel(r'z')
        plt.ylabel(r'$10^3 \times \Omega_\mathrm{DLA}$')

    def plot_omega_dla(self, zmin=2, zmax=4, label="GP", color=None, twosigma=True):
        """Plot omega_DLA as a function of redshift, with full Bayesian errors"""
        (z_cent, omega_dla, omega_dla_68, omega_dla_95, xerrs) =  self.omega_dla_cddf(z_min=zmin, z_max=zmax)
        if twosigma:
            plt.fill_between(z_cent_fill(z_cent,xerrs), 1000*omega_dla_95[:,0], 1000*omega_dla_95[:,1], color="grey", alpha=0.5)
        omega_dla*=1000
        yerr = (omega_dla-1000*omega_dla_68[:,0], 1000*omega_dla_68[:,1] - omega_dla)
        plt.errorbar(z_cent, omega_dla, yerr=yerr, xerr = xerrs, fmt='s', label=label)
        plt.xlabel(r'z')
        plt.ylabel(r'$10^3 \times \Omega_\mathrm{DLA}$')
        plt.xlim(zmin, zmax)
        return (z_cent, omega_dla, omega_dla_68, omega_dla_95)

    def _get_sample_params(self,spec,*,second=False):
        """Get the (n,z) values for each sample in this spectrum. spec is the spectrum number,
        second denotes whether to return the parameters of the second DLA."""
        #Compute redshift of each sample
        redshifts = self.z_min(spec) + (self.z_max(spec) - self.z_min(spec)) * self.z_offsets
        lnhi_vals = self.lnhi_vals
        #Get N,z values for this spectrum
        if second:
            base_sample = self._base_sample_inds(spec)
            lnhi_vals = lnhi_vals[base_sample]
            redshifts = redshifts[base_sample]
        return (lnhi_vals, redshifts)

    def _get_prob_dla_this_bin(self, spec, index, *, second=False):
        """Get the probability of a DLA with the samples specified in index."""
        return np.exp(self._log_norm_like(spec,second=second)[index]) * self._p_dla(second=second)[spec]

    def _split_distributions(self, q_bins, lred=2., ured=4., lnhi_min=20.3, lnhi_max=23., *, nhi=False):
        """Split the distributions for both the first and second DLA, in turn"""
        (probs, poissons) = self._split_distributions_single(q_bins, lred=lred, ured=ured, lnhi_min=lnhi_min, lnhi_max=lnhi_max, nhi=nhi, second=False)
        if self.second_dla:
            (probs2, poissons2) = self._split_distributions_single(q_bins, lred=lred, ured=ured, lnhi_min=lnhi_min, lnhi_max=lnhi_max, nhi=nhi, second=True)
            #List addition is concatenation, but we want a list of lists.
            probs = list(map(operator.add, probs, probs2))
            #Array addition is element-wise addition
            poissons += poissons2
        return (probs, poissons)

    def lymanbeta(self, zqso):
        """Compute the redshift at which the lyman beta forest at the redshift of the quasar will show up."""
        waveratios = 1026.72/1215.67
        zlyb = (1+zqso) * waveratios - 1
        return zlyb

    def proximity(self, zqso):
        """Remove a redshift range close to the quasar"""
        dz = self.proximity_zone
        return zqso - dz

    def _split_distributions_single(self, q_bins, lred=2., ured=4., lnhi_min=20.3, lnhi_max=23., *, nhi=False, second=False):
        """
            Split the sampled probabilities (in the desired bin) into two sets; those with small probabilities, for which we just keep the mean and sum of squares
            and will model with a Poisson distribution, and those with large probabilities, which we keep exactly for further computation.
        """
        #A list of probabilities for each redshift bin
        probs = [list() for _ in q_bins[:-1]]
        poisson_list = [list() for _ in q_bins[:-1]]
        dla_ind = self.filter_dla_spectra(second=second)
        for spec in dla_ind[0]:
            #Compute redshift of each sample
            (lnhi_vals, redshifts) = self._get_sample_params(spec, second=second)
            #The low cutoff redshift.
            upper_z = ured
            if self.lowzcut:
                upper_z = np.min([self.proximity(self.z_max(spec)), ured])
            #Select only samples with a DLA value, within the redshift we want.
            desired_samples = (lnhi_vals > lnhi_min)*(lnhi_vals < lnhi_max)*(redshifts < upper_z)*(redshifts > lred)
            if self.filter_noisy_pixels:
                #Exclude pixels which have too large noise within them
                #These are the indexes of the samples in the pixel noise vector
                pn = self.pixel_noise[spec]
                pind = np.array((redshifts-self.z_min(spec))/(self.z_max(spec)-self.z_min(spec))*np.size(pn),dtype=np.int)
                desired_samples *=(pn[pind] < self.noise_thresh)
            ind = np.where(desired_samples)
            if np.size(ind) == 0:
                continue
            #Find the probability that we have a DLA from this spectrum in each redshift bin
            p_dla_each_bin = self._get_prob_dla_this_bin(spec, ind[0], second=second)
            ind2 = np.where(p_dla_each_bin > self.p_thresh_sample)
            if np.size(ind2) == 0:
                continue
            #If this is computing the CDDF, use lnhi_vals. Otherwise use redshift for dN/dX and omega_DLA
            if nhi:
                quantity = lnhi_vals[ind]
            else:
                quantity = redshifts[ind]
            for iz in range(np.size(q_bins)-1):
                p_dla_this_z = p_dla_each_bin[ind2][np.where((quantity[ind2] > q_bins[iz])*(quantity[ind2] < q_bins[iz+1]))]
                if np.size(p_dla_this_z) == 0:
                    continue
#                 assert np.all(p_dla_this_z > 1e-4)
                #Add small probability events to the Poisson approximation: use a stable sum as this is probably *very* unstable.
                ipois = np.where(p_dla_this_z < self.p_switch)
                if np.size(ipois) > 0:
                    poisson_list[iz].append(math.fsum(p_dla_this_z[ipois]))
                #Add large probability events to the direct compute chain
                idla = np.where(p_dla_this_z >= self.p_switch)
                if np.size(idla) > 0:
                    probs[iz].append(p_dla_this_z[idla])
        poissons= np.array([math.fsum(pl) for pl in poisson_list])
        #Check that the Poisson approximation is a reasonable one; in practice this seems pretty good.
        #poissonsquare= np.array([math.fsum(pl**2) for pl in poisson_list])
        #assert np.all(poissonsquare/poissons < 0.2)
        return probs, poissons

    def _get_combined_levels(self, pdf_pb, pmean):
        """Get the combined pdf of a Poisson binomial process and a Poisson distribution with parameter pmean"""
        cdf_dla = np.cumsum(pdf_pb)
        #Properties of a zero poisson distribution are not defined.
        if pmean == 0.:
            return (pdf_pb, 0)
        weak = poisson(pmean)
        #So now we have the PDF of the likely DLAs (which may not be Poisson). Add in the PDF of the Poisson process describing the others
        #Neglect the tails where either CDF is < 1e-4
        (plow, phigh) = weak.interval(1-1e-4)
        plow=int(plow)
        phigh=int(phigh)
        (dlow, dhigh) = interval(cdf_dla, 1-1e-4)
        #print(pmean, plow, phigh, np.argmax(pdf_pb), dlow, dhigh)
        #Note that in practice a not terrible approximation is just to sum the confidence intervals.
        #But that marginally overestimates the errors!
        pdf_comb = np.array([math.fsum([weak.pmf(N-i)*pdf_pb[i] for i in range(dlow,np.min([dhigh+1,np.size(pdf_pb)]))]) for N in range(plow+dlow,phigh+dhigh+1)])
        assert 1.00 > math.fsum(pdf_comb) > 0.99
        return (pdf_comb, plow+dlow)

    def _get_confidence_intervals(self, q_bins, lred=2., ured=4., lnhi_min=20.3, lnhi_max=23., nhi=False):
        """
        Get the confidence interval on the number of DLAs in a given redshift (and column density) bin.
        The number of DLAs is the sum of n binomial processes, and so given by a likelihood looking like:
        P(N=n) = sum(all subsets of n) prod (1-p_i) * prod p_i where the first product is over all non-DLA spectra and the second over all DLA spectra.
        This function is too complex to be evaluated directly, but can be solved using an FFT for p large and small N.
        For all p < p_switch we approximate the distribution as a Poisson distribution using Le Cam's (1960) theorem; the error from this is bounded by D_2 (sum(p_j^2)/sum(p_j),
        where D_2 < 16 and is probably ~ 1 here.
        Returns: (maximum a posteriori likelihoods, lower 68 % confidence levels, upper 68% confidence levels, lower and upper 95 % confidence levels)
        """
        (probs, poissons) = self._split_distributions(q_bins, lred=lred, ured=ured, lnhi_min=lnhi_min, lnhi_max=lnhi_max, nhi=nhi)
        #probs[i] now contains a list of arrays
        #Now we have built a list of probabilities in each z bin of interest and we want to solve for the Poisson binomial coefficients.
        maxlikes = []
        levels68=[]
        levels95 = []
        for (pp, pmean) in zip(probs, poissons):
            pdf = get_poisson_binomial_pdf(pp)
            (pdf_comb, offset) = self._get_combined_levels(pdf, pmean)
            (maxlike, ll68, ll95) = pdf_confidence(pdf_comb, offset)
            #Check correctly ordered
            assert ll95[0] <= ll68[0] <= maxlike
            assert ll95[1] >= ll68[1] >= maxlike
            #Unpack maximum likelihoods and 68/95% contours
            maxlikes.append(maxlike)
            levels68.append(ll68)
            levels95.append(ll95)
        return (maxlikes, levels68, levels95)

    def _get_z_nhi_hist(self, q_bins, lred=2., ured=4., lnhi_min=20.3, lnhi_max=23., nhi=False, moment=False):
        """
        Estimate the mean and standard deviation on the number of DLAs in a given redshift bin.
        Since each DLA has some probability of being in a given bin, p_dla * p_in_this_bin,
        each DLA is a binomial process, and the sum is a binomial poisson process.
        Thus the mean is sum(p_dla * p_in_this_bin) and the variance sum[p(1-p)]
        Ignore spectra with p_DLA < p_thresh, as an optimization.
        """
        dla_ind = self.filter_dla_spectra()
        means = np.zeros(np.size(q_bins)-1)
        variances = np.zeros(np.size(q_bins)-1)
        for spec in dla_ind[0]:
            #Compute redshift of each sample
            (lnhi_vals, redshifts) = self._get_sample_params(spec)
            #Select only samples with a DLA value, within the redshift we want.
            ind = np.where((lnhi_vals > lnhi_min)*(lnhi_vals < lnhi_max)*(redshifts < ured)*(redshifts > lred))
            if np.size(ind) == 0:
                continue
            #Find the probability that we have a DLA from this spectrum in each redshift bin
            p_dla_each_bin = self._get_prob_dla_this_bin(spec, ind[0], second=False)
            #Multiply by the column density to get total amount of HI instead of the number of DLAs
            if moment:
                weight = 10**lnhi_vals[ind]
            else:
                weight = 1.
            #If this is computing the CDDF, use lnhi_vals. Otherwise use redshift for dN/dX and omega_DLA
            if nhi:
                quantity = lnhi_vals[ind]
            else:
                quantity = redshifts[ind]
            #These are the means
            (t_hist, _) = np.histogram(quantity, bins=q_bins, weights=weight*p_dla_each_bin)
            means += t_hist
            #These are the variances
            (t_var, _) = np.histogram(quantity, bins=q_bins, weights=weight*weight*(1-p_dla_each_bin)*p_dla_each_bin)
            variances += t_var
        #Don't forget Poisson term from sample variance.
        #The variance before this indicates the uncertainty arising from our imperfect knowledge of the properties of the DLAs
        #in our spectra; this term indicates our imperfect *sampling* of the total population
        #If we had one spectrum which we were certain contained a DLA, this would estimate the error.
        variances += means
        return means, variances

    def find_delta_NHI(self, nspec):
        """Find the range of NHI values in nspec with a likelihood 1/2e times the max.
        This is an easily calculable value which is the 2-sigma contour if the likelihood is Gaussian"""
        likes = self._log_norm_like(nspec)
        mlike = np.max(likes)
        nvals = self.lnhi_vals[np.where(likes > mlike-2)]
        return np.max(nvals) - np.min(nvals)

    def find_delta_z(self, nspec):
        """Find the range of redshift values in nspec with a likelihood 1/2e times the max.
        This is an easily calculable value which is the 2-sigma contour if the likelihood is Gaussian"""
        likes = self._log_norm_like(nspec)
        mlike = np.max(likes)
        nvals = (self.z_max(nspec) - self.z_min(nspec)) * self.z_offsets[np.where(likes > mlike-2)] + self.z_min(nspec)
        return np.max(nvals) - np.min(nvals)

    def find_max_like(self, nspec, *, second=False):
        """Find the maximum likelihood values of NHI and redshift"""
        likes = self._log_norm_like(nspec, second=second)
        mlike = np.argmax(likes)
        (lnhi_vals, redshifts) = self._get_sample_params(nspec, second=second)
        return lnhi_vals[mlike], redshifts[mlike]

    def find_real(self, nspec, *, field = "flux"):
        """Find the index of a quasar in the raw datafile"""
        #Load the indices of the quasars we have data for in the raw file
        nspec_real = self.real_index[nspec]
        hh = h5py.File(self.raw_file,'r')
        flux = hh[hh["all_"+field][0][nspec_real]][0]
        nflux = np.size(flux)
        zzs = (self.z_max(nspec) - self.z_min(nspec))*range(nflux)/nflux+self.z_min(nspec)
        hh.close()
        return zzs, flux

def find_snr(nspec, real_index, raw_file, zmin, zmax):
    """Find the signal to noise ratio, according to the definition where it is the flux/s.d. noise."""
    #Get noise variance
    _ = zmin
    nspec_real = real_index[nspec]
    hh = h5py.File(raw_file,'r')
    wavelengths = hh[hh["all_wavelengths"][0][nspec_real]][0]
    #ipix = np.where(np.logical_and(wavelengths > 1215.67*(1+ zmin), wavelengths < 1215.67*(1+zmax)))
    ipix = np.where(wavelengths > 1215.67*(1+zmax))
    flux = np.array(hh[hh["all_flux"][0][nspec_real]][0])[ipix]
    try:
        norm = hh["all_normalizers"][0][nspec_real]
        #This is so that we don't have an unrealistically low noise threshold inside of absorbers.
        flux[np.where(flux/norm < 0.1)] = norm*0.1
    except KeyError:
        flux[np.where(flux < 0.1)] = 0.1
    noise_var = np.array(hh[hh["all_noise_variance"][0][nspec_real]][0])[ipix]
    hh.close()
    return 1/np.median(np.sqrt(noise_var)/np.abs(flux))

def find_pixel_noise(nspec,real_index, raw_file, zmin, zmax):
    """Find pixels where the absolute value of the noise is below thresh a particular value.
    So we want pixels with: all_noise_variance/all_normalizers^2 < thresh^2
    where all_noise_variance is the noise and defined in preloaded_qsos."""
    nspec_real = real_index[nspec]
    hh = h5py.File(raw_file,'r')
    norm = hh["all_normalizers"][0][nspec_real]
    wavelengths = hh[hh["all_wavelengths"][0][nspec_real]][0]
    ipix = np.where(np.logical_and(wavelengths > 1215.67*(1+ zmin), wavelengths < 1215.67*(1+zmax)))
    noise_var = np.array(hh[hh["all_noise_variance"][0][nspec_real]][0])[ipix]
    hh.close()
    return noise_var/norm**2

def find_pixel_snr(nspec,real_index, raw_file, zmin, zmax):
    """Find pixels where the absolute value of the noise is below thresh a particular value.
    So we want pixels with: all_noise_variance/all_normalizers^2 < thresh^2
    where all_noise_variance is the noise and defined in preloaded_qsos."""
    nspec_real = real_index[nspec]
    hh = h5py.File(raw_file,'r')
    wavelengths = hh[hh["all_wavelengths"][0][nspec_real]][0]
    ipix = np.where(np.logical_and(wavelengths > 1215.67*(1+ zmin), wavelengths < 1215.67*(1+zmax)))
    flux = np.array(hh[hh["all_flux"][0][nspec_real]][0])[ipix]
    noise_var = np.array(hh[hh["all_noise_variance"][0][nspec_real]][0])[ipix]
    try:
        norm = hh["all_normalizers"][0][nspec_real]
        #This is so that we don't have an unrealistically low noise threshold inside of absorbers.
        flux[np.where(flux/norm < 0.1)] = norm*0.1
    except KeyError:
        flux[np.where(flux < 0.1)] = 0.1
    hh.close()
    return np.sqrt(noise_var)/np.abs(flux)


def compute_all_snrs(*, raw_file="preloaded_qsos.mat", processed_file="processed_qsos_dr12q_lyb_lya.mat", save_file="snrs_qsos_dr12.mat"):
    """Compute the SNR for all spectra and save to a separate file"""
    ff = h5py.File(processed_file,'r')
    real_index = np.where(ff["test_ind"][0] != 0)[0]
    min_z_dla = np.array(ff["min_z_dlas"][0])
    max_z_dla = np.array(ff["max_z_dlas"][0])
    ff.close()
    snrs = np.array([find_snr(nn, real_index, raw_file, min_z_dla[nn], max_z_dla[nn]) for nn in range(np.size(real_index))])
    f = h5py.File(save_file, 'w')
    f["snrs"] = snrs
#     dt = h5py.special_dtype(vlen=np.dtype('float64'))
#     dset = f.create_dataset('pixel_noise', (np.size(real_index),), dtype=dt)
#     for nn in range(np.size(real_index)):
#         dset[nn] = find_pixel_noise(nn, real_index, raw_file, min_z_dla[nn], max_z_dla[nn])
#     dset = f.create_dataset('pixel_snr', (np.size(real_index),), dtype=dt)
#     for nn in range(np.size(real_index)):
#         dset[nn] = find_pixel_snr(nn, real_index, raw_file, min_z_dla[nn], max_z_dla[nn])
    f.close()

def HubbleByH0(z, Omega_m=0.279):
    """Hubble function divided by H0, H/H0(z).
        H/H0(z)**2 = Omega_m/a^3 + Omega_lambda
        We neglect curvature and radiation, and assume Omega_lambda = 1- Omega_m.
        Omega_m is WMAP 9 by default
    """
    return math.sqrt(Omega_m* (1+z)**3 + (1 - Omega_m))

def interval(cdf, level, offset=0):
    """Return a tuple with the confidence interval at level for the given cdf.
    level should be between 0 and 1. Larger values mean wider intervals"""
    if np.size(cdf) == 1:
        return (offset, offset)
    ii = np.where((cdf <= 0.5+level/2.)*(cdf >= 0.5-level/2.))
    #This can happen when all the cdf is in one bin, and it is on the edge.
    if True or np.size(ii) == 0:
        high=1+offset
        low = offset
        idown = np.where(cdf < 0.5-level/2)
        if np.size(idown) != 0:
            low += idown[0][-1]+1
        iup = np.where(cdf > 0.5+level/2)
        if np.size(iup) != 0:
            high += iup[0][0]
        else:
            high = np.size(cdf)
        return (low, high)
    return (ii[0][0]+offset, 1+ii[0][-1]+offset)

def pdf_confidence(pdf_comb, offset):
    """Get the maximum likelihood value, and 68 and 95 % confidence levels of a probability density function.
    offset is a value to add to all derived quantities
    """
    #Find the cumulative distribution function
    cdf_comb = np.cumsum(pdf_comb)
    #Find the maximum a posteriori likelihood.
    maxlike = interval(cdf_comb, 0., offset=offset)[0]
    ll68 = interval(cdf_comb, 0.68, offset=offset)
    ll95 = interval(cdf_comb, 0.95,offset=offset)
    assert maxlike >= ll68[0] >= ll95[0]
    assert maxlike <= ll68[1] <= ll95[1]
    return maxlike, ll68, ll95

def get_poisson_binomial_pdf(pp):
    """Get the (exact) PDF of a poisson binomial process from an array listing probabilities"""
    #Check input is reasonable
    if np.size(pp) == 0:
        return np.ones(1)
    ppa = np.array(np.concatenate(pp))
    assert ppa.dtype == np.float64
    assert np.size(np.shape(ppa)) == 1
    Nsamp = np.size(ppa)
    #Compute the coefficients of a DFT which we will use to find the poisson-binomial coefficients.
    #See Fernandez, M. and Williams, S. 2010 or wikipedia
    nco = lambda nn: cmath.exp( - 2*math.pi*nn*1j/(Nsamp+1) ) -1
    #Use the symmetry of the fourier transform to only compute the first half of the array: we know that the pdf is real.
    coeffs = np.array([stable_complex_product(1+ppa*nco(nn) ) for nn in range((Nsamp+1)//2 + 1)])
    #Check for roundoff; should be ok as all the coefficients that are multiplied are within the unit circle
    #Almost all coeffs should be complex...
    assert np.any(np.absolute(coeffs) > 0)
    #Do the FFT
    pdf = np.fft.irfft(coeffs, n=Nsamp+1)
    #Make sure we got a reasonable answer
    assert np.all(np.logical_not(np.isinf(pdf)))
    #Check correctly normalized
    assert np.abs(math.fsum(pdf) - 1.) < 1e-7
    return pdf

def stable_complex_product(iterable):
    """
    Compute the product of a large array of complex numbers using logs and python's stable summation routine.
    Some shenanigans are necessary because we want to avoid taking a complex logarithm.
    We say that if z = r e ^i theta , then
    prod(z) = prod(r) exp(sum i theta) = exp( sum log (r) + sum i theta )
    Returns a long double result so we can store very small values!
    """
    rr = np.absolute(iterable)
    theta = np.angle(iterable)
    return np.exp( math.fsum(np.log(rr)) + 1j*math.fsum(theta),dtype=np.complex256)

def path_length_int(z, Omega_m=0.279):
    """Integrand function for the path length integral above.
    dX = (1+z)^2 H_0 / H(z) dz
    returns dX
    """
    return (1+z)**2 / HubbleByH0(z, Omega_m)

def rho_crit(hubble=0.7):
    """Get the critical density at z=0 in units of g cm^-3"""
    #H in units of 1/s
    #h * 100 km/s/Mpc in h/s
    h100=3.2407789e-18*hubble
    gravcgs =  6.674e-8
    rho_c=3*h100**2/(8*math.pi*gravcgs)
    return rho_c

def z_cent_fill(z_cent, xerrs):
    """Small function to find an expanded filling region.
    Used in plotting 2 sigma contours."""
    filler = np.array(z_cent)
    filler[0] -= xerrs[0][0]
    filler[-1] += xerrs[-1][-1]
    return filler

