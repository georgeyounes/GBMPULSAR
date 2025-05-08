# gbmpulsar - utilizing Fermi GBM as a pulsar timing tool

## Introduction

This library provides tools to process gbm daily continuous time-tagged event (TTE) data for pulsar analysis. There are 
several steps required for this task: (1) Download the data, (2) barycentering POSHIST, (3) filter the data for 
maximizing signal-to-noise and barycenter correct the events TIME, (4) phase-calculation, and (5) visualization. 
Continuous TTE files are quite large, and we shall be dealing with over 10^8 rows of "good" events after filtering the 
data per day. For this purpose, **Dask** is utilized to perform almost all operations lazily before computing the 
output. Testing this on my 32 GB, 10-core, M1 machine, I manage to run a 1-day analysis in about 2 minutes utilizing 8 
threads, compared to about 18 minutes when not running in parallel. The tools are scalable so the user should be able to
speed things up even further with more cores and available memory. 

## Acknowledgements

I am grateful to the discussions I have had with Paul Ray on data barycentering and pulsar timing. This work is 
inspired by the long-term RXTE, Swift/XRT, and NICER magnetar and pulsar monitoring work.

## Installation

At the moment, gbmpulsar can be installed locally after cloning the directory or simply downloading 
and untarring it. Then from the GBMPULSAR root directory:

```bash
  python -m pip install .
```

You may add the '-e' option to pip to install as an editable package.

The code has been tested on Python 3.9.16, and it requires matplotlib, pandas, scipy, astropy, numpy, and dask.

## Detailed description

Below, I provide a short description of each step in the process of building the final products:

#### (1) Download continuous TTE data

A command line tool **getgbmdailydata** is available to download the required data, e.g., as follows

```bash
  getgbmdailydata 2017/11/01 -t tte
```
which will download all the tte files for all detectors into the current directory. Optional arguments allow for a 
different output directory (which will be created if it does not exist). To download the POSHIST file

```bash
  getgbmdailydata 2017/11/01 -t poshist
```

Note that the TTE files are broken into chuncks of 24 files for each hour of the day, and there are 12 NaI detectors. 
Each seperate file is about 20 MB large (total of several GB per day). Hence, make sure you are on a fast internet 
before you go down this route. Moreover, I highly recommend that all 13+ years worth of data is downloaded (about 30 TBs
or so) to a cluster or any data storage space. You can always point the code to that root directory for analysis. 
Of course, if your analysis can be done with portion of the dates, the above is not a necessity.

#### (2) Barycentring

Unfortunately the GBM team does not provide a barycentering tool for their data, and building one 
from scratch requires considerable thought and, especially, testing. Fermi LAT team, on the other hand, provide a fast 
barycentering tool **gtbary**, which is based on the well tested and documented tools axbary and barycorr of HEASoft. I 
opted to hack the barycentering of GBM data. I note that the barycentering is done on the GBM POSHIST time stamps which 
are 1-second ticks, i.e., about 86400 per each daily POSHIST file. This is greatly faster than barycentering the full 
array of times. The latter is done below via interpolation. 

For barycentering the POSHIST file, we utilize a dummy Fermi LAT event file, and the actual Fermi LAT spacecraft file. 
We swap the TIME column in the former with the TIME column of the POSHIST file, along with the necessary header 
keywords. **gtbary** is then run on that file, and the resulting barycentered TIMES are added to the POSHIST file as a 
new column (named BARY_TIME). This can be performed with the command-line tool **barycenterposhist** as follows:

```bash
  barycenterposhist glg_poshist_all_171101_v01.fit 83.633218 22.014464 lat_evtFile_placeHolderTTE.fits lat_spacecraft_merged.fits
```
where *glg_poshist_all_171101_v01.fit* is the POSHIST file to be barycentered, *83.633218 22.014464* are the RA and DEC
of the pulsar of interest (in this case the Crab), *lat_evtFile_placeHolderTTE.fits* is the dummy placeholder LAT event
file, and *lat_spacecraft_merged.fits* is the LAT spacecraft file. Make sure the latter is up-to-date for your analysis 
(you can in fact have one that runs through the full Fermi mission lifetime which can be found 
[here](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/).

Finally, you of course need to have the [fermitools](https://github.com/fermi-lat/Fermitools-conda/wiki) installed in 
order to have access to **gtbary**.

#### (3) Filtering and event barycentering

In order to maximize signal-to-noise, we perform a series of orbital filtering to minimize the background contribution 
to the source events. These are detector-to-source angle which dictates the detector effective area to the source flux 
(default $\leq50^{\circ}$), detector-to-zenith which ensures that the detector is looking far enough from earth limb
(default $\leq60^{\circ}$), and source-to-zenith which ensures the source is above earth limb for Fermi (i.e., not earth 
blocked, default $\leq75^{\circ}$). Once this filtering is done, we utilize the POSHIST barycentered times to perform
the barycentering of the "clean" events through an interpolation - the barycentered times are fit to a spline and that
is utilized to predict the correction for each event timestamp. This interpolation is accurate to less than 10 
microsecond which is sufficient for timing any pulsar.

#### (4) Phase assignment

Rotational phases are assigned utilizing [PINT](https://nanograv-pint.readthedocs.io/en/latest/) according to a .par 
file. Polycos, which is implemented in PINT, is used to perform the task to speed up the process. Alternatively, the 
user can use [CRIMP](https://github.com/georgeyounes/CRIMP) which is natively fast, yet this is only valid for isolated 
pulsars! All of the above work only on barycentered data. If the pulsar you are interested in have proper motion or 
parallax variation, you will need to take care of that separately, e.g., by performing the daily barycentering with a 
varying position accordingly. 

Steps (2) and (3) can be performed with the primary task of the pipeline as follows:

```bash
  gbmpulsarpipe glg_poshist_all_171101_v01.fit all 171101 83.633218 22.014464 crab_171104.par crab ../171101
```
where *glg_poshist_all_171101_v01.fit* is the barycentered poshist file for the day (can add path to it if desired), 
*all* is a list of detectors to be analyzed (or comma-separated list), 171101 is the date of interest in format YYMMDD, 
*83.633218 22.014464* are the RA and DEC of the pulsar of interest, *crab_171104.par* is the .par timing file, *crab* 
is the pulsar name, *../171101* is the input directory for where the TTE files reside. There are several optional 
arguments that will be discussed with full documentation. The user can run **gbmpulsarpipe -h** for a full list of 
arguments.

Once finished, pulse profiles with 120 bin/cycle resolution and in several PHA bands will be created in .csv format, 
and dropped in an output directory folder (**-od** or **--OUTPUTDIR**, default **args.date** or ./171101 in the example 
case above). The PHA bins are as follows:

```bash
phabins = [4, 7, 22, 33, 51, 73, 86, 104, 128]  # corresponding to
energybins = [8, 10, 30, 50, 100, 200, 300, 500, 1000]  # keV (bins are [start, end))
channels = [0, 1, 2, 3, 4, 5, 6, 7]  # corresponding channels
```

#### (5) visualization

Once the code has been run on a variety of dates, visualization can be done with the command-line tool 
**gbmmergepulseprofiles** which merges the pulse profiles in channel (energy) and time. A waterfall plot is created. The
script can be run as follows:

```bash
gbmmergepulseprofiles 2017-11-01 2017-11-06 -pn crab_171101-171106
```
which merges all the pulse profiles from the above dates and all channels together, but does not perform any further 
binning per cycle or time. Note that the script is expecting the input directory to the .csv files to be the current 
directory and dates of interest (i.e., ./171101, ./171102, etc.; the root directory can be changed with -rt --ROOT_DIR). 
An exmaple plot for the Crab data on the dates 2017-11-01 to 2017-11-06 can be found in the folder data. 
**gbmmergepulseprofiles** also affords multiple optional input parameters for varying binning.

## Expected data structure

Your data directory should be organized as follows:

    root_dir/
        171126/
            pulseprofile_<srcname>_4-7_pha.csv
            pulseprofile_<srcname>_7-22_pha.csv
            ...
        171127/
            pulseprofile_<srcname>_4-7_pha.csv
            ...

This should be automatically the case if you follow with the pipeline. An example of such can be found in the folder 
data.

## To do

Creating pulse profile is only part of the fun. One would want to calibrate those so that any pulsed flux measurement 
can be converted to a pulsed flux. This is not implmented yet, but will be in a future release. Thankfully, the Crab 
is a natural calibration source and happens to be exquisitely bright in GBM. This makes the job a little easier. Note
that the pipeline already properly trakcs the good time interval for each detector per source, hence so the measured 
rates are already accurate. 

Another wish is to make the library compatible with BGO detectors. As it stands currently, it only works (properly) with
NaI detectors.

## Disclaimer

This code is distributed in the hope that it will be useful, but WITHOUT ANY EXPRESSED OR 
IMPLIED WARRANTY. Use this code solely if you understand what it is doing. Feel free to 
reach out with any questions.

## License

[MIT](https://choosealicense.com/licenses/mit/)
