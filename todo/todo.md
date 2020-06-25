# Redshift Code -> CIV code

- [x] Redshift code has a learning script without absorption.

- [ ] Redshift code has out-data model (the i.i.d. Gaussian models for data outside the range of modelling window), which is not needed for CIV code.

- [ ] Need a MATLAB version of `read_spec.m` for DR7. Python version in [tutorials/read_spec_dr7.py](tutorials/read_spec_dr7.py)

- Redshift code doesn't have Voigt.c, which should be added with CIV profile (assumed to be a doublet voigt).
  - [ ] add `voigt.c` and change to doublet.

- Redshift code doesn't build absorber catalogue, which should be added.
  - [ ] `build_catalogs.m` add CIV hash table.

- Redshift code doesn't process absorber model, which should be added.
  - [ ] `process_qsos.m` add CIV model log likelihood.

- The prior for CIV column density should be added (previously was
generate_dla_samples.m in Garnett's code).
  - [ ] `generate_dla_samples.m` add CIV column density prior.

## Redshift Code descriptions

- [ ] (can be removed) `log_mvnpdf_iid.m`, calculate log likelihood of a multivariate normal pdf without off-diagonal covariance.
- [ ] (can be removed) `read_data.m`, Leah's script for DLA finding + zQSO estimation.
- [ ] (can be removed) `process_spectra.m`, Leah's script for ROC curve of DLA finding + zQSO estimation.
- [ ] (can be removed) `test_logpdf_iid.m`, test file for `log_mvnpdf_iid.m`.
- [ ] (can be removed) `test_outdata_logpdf.m`, test file for `log_mvnpdf_iid.m`.
