# Azimuthal anisotropic ambient noise tomography

This contains a collection of Python scripts that can be used to reproduce the ambient-noise based Rayleih phase dispersion measurements, used in

**Anisotropic reversible-jump MCMC shear-velocity tomography of the eastern Alpine crust (Kästle and Tilmann, 2024), submitted to G3**

The scripts are not fully documented and some of the paths will have to be adapted. For an overview of the method and the chosen parameters, please see the above cited article and the references therein. Feel free to contact Emanuel Kaestle in case of question or if you need help with any of the scripts.

If you plan to use these scripts in your work, please cite the article above.

## Usage

You need a Python (3.xx) installation and a number of Python packages, above all Obspy. For recommended installation options, please see the Obspy website (https://docs.obspy.org/). You will further need the following packages: numpy, scipy, matplotlib, mpi4py, pandas and sqlite (possibly others).

The scripts should be run according to the numbering. The most important options are explained in each script.

__01_massdownloader.py__

Uses the Obspy massdownloader routine to download the necessary raw data. At the time of the study, most of the data was only accessible using an EIDA token. Most of the data should, however, be publicly available by now. So that the script can be adapted for a use without token. The script can be run repeatedly, to check for updated data, download missing files or additional station data. It should check automatically which files are downloaded already.

__02_preprocess_files.py__

This script takes the previously downloaded data and runs the preprocessing steps (downsampling, filtering, removal of instrument response). This script can be executed while the download is still running. Depending on the user choice, the original, raw files are left untouched, deleted or overwritten.

__03_create_ccs.py__

Creates cross correlations between all available station pairs. For bookkeeping, a database file is created. The processing allows to extend the existing cross correlation files, if more data is added at a later time. Optionally, it also saves monthly cross correlations.

__04_process_spectra_monthly_zz.py__

This uses the cross correlation files to extract pair-wise phase dispersion measurements. 

__noise.py__

This contains several functions that are necessary to run the other scripts.

## Further information

Plase see also the other packages on my github page (https://github.com/ekaestle/). The amb noise tools package contains a similar collection of scripts with a few simple application examples. Scripts that can be used to create tomographic models can be found in the other sections.
