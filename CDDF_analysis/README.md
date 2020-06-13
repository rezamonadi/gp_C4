# Basic usage of the QSOLoader class

Also refer to `tutorials/tutorial.py`.

This is how to instantiate this little class:

```python
from CDDF_analysis.qso_loader import QSOLoaderZ

# in python
qsos = QSOLoaderZ(
        preloaded_file="preloaded_zqso_only_qsos.mat",
        catalogue_file="zqso_only_catalog.mat",
        learned_file="learned_model_outdata_dr9q_minus_concordance_norm_1176-1256.mat",
        processed_file="processed_zqso_only_qsos_dr12q-100_uniformprior.mat",
        dla_concordance="dla_catalog",
        los_concordance="los_catalog",
        sample_file="dla_samples.mat",
        occams_razor=False)
```

- `preloaded_file`, `learned_file`, `processed_file` and `catalog_file` are in the data product folder: http://tiny.cc/gp_zestimation_catalogue.
- For concordance catalogues, they are in the path `data/dla_catalogs/dr9q_concordance/processed/dla_catalog` and `data/dla_catalogs/dr9q_concordance/processed/los_catalog` if you have run this before:

```bash
# in shell
cd data/scripts
./download_catalogs.sh
```

The most useful feature is to plot a given spectrum with the GP mean prior:

```python
# the index of the catalogue
nspec = 1

qsos.plot_this_mu(nspec)
```

There are other routines in `make_zqso_plots.py` to reproduce the plots in the paper.

```python
from CDDF_analysis import make_zqso_plots

# plot GP mu and kernel
make_zqso_plots.do_procedure_plots(qsos)

# plot figure 6 & 7 in the paper
make_zqso_plots.do_plot_thing_ids(qsos, selected_thing_ids=[544031279, 27885089])

# reproduce z_map - z_PCA plot
make_zqso_plots.do_velocity_dispersions(qsos, dr12q_fits='data/dr12q/distfiles/DR12Q.fits')
```
