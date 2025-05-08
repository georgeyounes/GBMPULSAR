import multiprocessing as mp

# Use fork start method before any child processes
mp.set_start_method("fork", force=True)

# Suppress PINT extended-precision RuntimeWarnings globally and polycos debug registration messages globally
import logging
import warnings

logging.getLogger("pint.polycos").setLevel(logging.WARNING)
warnings.filterwarnings(
    "ignore",
    message="This platform does not support extended precision floating-point"
)

import argparse
import sys
import time
import numpy as np
import shutil
import os
import pandas as pd

import dask.dataframe as dd
from dask.distributed import Client

from gdt.missions.fermi.time import Time

from gbmpulsar.filterttedataforsource import (
    gtifromposhistpersource,
    mergettefullday,
    readposhistforttebarycenter,
    timefilterttedf,
    barycenter_ttetimes,
    definegbmdetectors
)
from gbmpulsar.gbmpulsar_logging import get_logger
from gbmpulsar.binphases import binphases
from crimp.pulseprofile import plotpulseprofile
from crimp.calcphase import calcphase

import pint.models
from pint.polycos import Polycos

pint.logging.setup(level="WARNING")

# Configure gbmpulsar module-level logger
logger = get_logger(__name__)

# Avoid writing .pyc files
sys.dont_write_bytecode = True


def filter_bary_mjd(df, gti, poshist_df):
    df2 = timefilterttedf(df, gti)
    if df2.empty:
        return df2
    df2.loc[:, 'BARY_TIME'] = barycenter_ttetimes(poshist_df, df2['TIME'])
    df2.loc[:, 'BARY_TIME_MJD'] = Time(df2['BARY_TIME'], format='fermi', scale='tdb').mjd
    assert isinstance(df2, pd.DataFrame), "filter_bary_mjd did not return a DataFrame"
    return df2


def eval_phases_pint(df, ptable):
    if df.empty:
        df.loc[:, 'PHASES'] = np.array([], dtype=float)
        return df
    phases = ptable.eval_phase(df['BARY_TIME_MJD'].to_numpy())
    phases[phases < 0] += 1.0
    df.loc[:, 'PHASES'] = phases
    return df


def eval_phases_crimp(df, parfile):
    if df.empty:
        df.loc[:, 'PHASES'] = np.array([], dtype=float)
        return df
    _, phases = calcphase(df['BARY_TIME_MJD'], parfile)
    df.loc[:, 'PHASES'] = phases
    return df


def buildpulseprofiles(df, exposure, outdirectory, srcname, nbrpulsebins=120):
    # Corresponding to [8-10, 10-30, 30-50, 50-100, 100-200, 200-300, 300-500, 500-1000] keV (bins are [start, end))
    phabins = np.array([4, 7, 22, 33, 51, 73, 86, 104, 128])

    # Bin the PHA column using phabins and group by each bin - right=False implies [start, end)
    df['phabin'] = pd.cut(df['PHA'], bins=phabins, right=False)

    # Create a pulse profile for each bin and save it as a .csv
    for bin_range, group in df.groupby('phabin', observed=True):
        binnedprofile = binphases(group['PHASES'].values, nbrBins=nbrpulsebins)
        pulseProfile = {
            'ppBins': binnedprofile['ppBins'],
            'ppBinsRange': binnedprofile['ppBinsRange'],
            'countRate': binnedprofile['ctsBins'] / (exposure / nbrpulsebins),  # Approximate rate
            'countRateErr': binnedprofile['ctsBinsErr'] / (exposure / nbrpulsebins)  # Approximate rate
        }

        pulseProfile_df = pd.DataFrame(pulseProfile)
        pulseProfile_df.to_csv(os.path.join(outdirectory, 'pulseprofile_{}_{}-{}_pha.csv'.format(srcname,
                                                                                                 int(bin_range.left),
                                                                                                 int(bin_range.right))),
                               index=False)
        del group
    df.drop(columns='phabin', inplace=True)
    return None


def gbmpulsarpipe(poshistfile, detectors, date, phalow, phahigh, srcra, srcdec, parfile, inputdir,
                  usecrimp=False, n_workers=1, threads_per_worker=8):
    """
    Run the full GBM pulsar pipeline over one or more detectors.
    Returns a pandas DataFrame of combined results.
    """
    # Determine detectors list
    if detectors.lower() == 'all':
        nai_det_list, _ = definegbmdetectors()
    else:
        nai_det_list = [d.strip() for d in detectors.split(',')]

    # Read and prepare barycenter times
    poshist_alltimes_df = readposhistforttebarycenter(poshistfile)
    if not usecrimp:
        # Preparing things for phase calculation with PINT
        min_mjd = Time(poshist_alltimes_df['BARY_TIME'].min(), format='fermi', scale='tdb').mjd
        max_mjd = Time(poshist_alltimes_df['BARY_TIME'].max(), format='fermi', scale='tdb').mjd
        # Build a single polycos table
        model = pint.models.get_model(parfile, allow_T2=True, allow_tcb=True)
        ptable = Polycos().generate_polycos(model, min_mjd, max_mjd, '@', 120, 10, 0)

    # Start Dask client
    dayexposure = 0
    gooddet = 0
    with Client(n_workers=n_workers, threads_per_worker=threads_per_worker) as client:
        if not usecrimp:
            ptable_future = client.scatter(ptable, broadcast=True)
        parts = []

        for det in nai_det_list:
            # GTI for this detector
            full_gti_df, _ = gtifromposhistpersource(
                poshistfile, det, srcra, srcdec,
                poshist_res=10, zentosrccutang=80,
                zentodetcutang=60, dettosrccutang=50
            )
            if full_gti_df.empty:
                logger.info(f"No GTI for detector {det}")
                continue

            # Energy filter & merge
            eneflt_dd = mergettefullday(det, date, phalow, phahigh, inputdir)
            if eneflt_dd is None:
                logger.info(f"No data for detector {det} on date {date} albeit GTI is > 0 - should not happen unless "
                            f"tte files are missing - please check corresponding files in {inputdir}")
                continue
            else:
                # Done here just in case
                dayexposure += np.sum(full_gti_df['STOP'] - full_gti_df['START'])
                gooddet += 1  # Number of "good" detectors with on-source GTI != 0

            # Stage 1: GTI filter, barycenter, MJD
            meta1 = {
                'TIME': np.float64,
                'PHA': np.uint8,
                'BARY_TIME': np.float64,
                'BARY_TIME_MJD': np.float64
            }
            stage1 = eneflt_dd.map_partitions(
                filter_bary_mjd,
                full_gti_df,
                poshist_alltimes_df,
                meta=meta1
            ).repartition(npartitions=n_workers)

            # Stage 2: phase evaluation
            if not usecrimp:
                meta2 = {**meta1, 'PHASES': np.float64}
                part_dd = stage1.map_partitions(
                    eval_phases_pint,
                    ptable_future,
                    meta=meta2
                )
                parts.append(part_dd)
            else:
                meta2 = {**meta1, 'PHASES': np.float64}
                part_dd = stage1.map_partitions(
                    eval_phases_crimp,
                    parfile,
                    meta=meta2
                )
                parts.append(part_dd)

        if not parts:
            logger.warning("Pipeline produced no output for all detectors.")
            return None

        # Concatenate and compute
        result_dd = dd.concat(parts, interleave_partitions=True)
        result_df = result_dd.compute()

    logger.info('Total exposure summed over all good detectors: {} over {} \n'
                'Total number of counts across all good detectors: {}\n'
                'Rate = {}\n'.format(dayexposure, gooddet, len(result_df), len(result_df) / dayexposure))

    return result_df, dayexposure


def main():
    parser = argparse.ArgumentParser(description="GBM pulsar pipeline runner")
    parser.add_argument("poshistfile", help="Posthist file (barycentered)", type=str)
    parser.add_argument("detectors", help="Comma-separated detectors or 'all'", type=str)
    parser.add_argument("date", help="Date of interest (YYMMDD)", type=str)
    parser.add_argument("srcra", help="Source RA", type=float)
    parser.add_argument("srcdec", help="Source Dec", type=float)
    parser.add_argument("parfile", help="Pulsar timing .par file", type=str)
    parser.add_argument("srcname", help="Pulsar name", type=str)
    parser.add_argument("inputdir", help="Input directory for continuous TTE files", type=str)
    parser.add_argument("-cr", "--usecrimp",
                        help="Flag to use crimp for phase calculation - slightly faster but only applicable to isolated "
                             "sources", type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument("-pl", "--phalow", help="Low PHA bound (we recommend to leave as is)", type=int, default=4)
    parser.add_argument("-ph", "--phahigh", help="High PHA bound (we recommend to leave as is)", type=int, default=127)
    parser.add_argument("-of", "--outputfile", help="Name of output source events and phases .parquet "
                                                    "file - only necessary in case of follow-up analysis, e.g., "
                                                    "ToA calculation (default = None)", type=str, default=None)
    parser.add_argument("-ip", "--integratedpp", help="Name of output pulse profile plot - makes sense "
                                                      "only for crab (default = None)", type=str, default=None)
    parser.add_argument("-od", "--outputdir", help="Name of output directory where to drop created files"
                                                   "(default = ./args.date)", type=str, default=None)
    parser.add_argument("-nw", "--n_workers", help="Numbers of dask workers (i.e. jobs - default = 1)",
                        type=int, default=1)
    parser.add_argument("-th", "--threads_per_worker", help="Numbers of dask threads (default = 8)",
                        type=int, default=8)
    args = parser.parse_args()

    if args.outputdir is None:
        args.outputdir = args.date

    if not os.path.exists(args.outputdir):
        # if output directory is not present - then create it
        os.makedirs(args.outputdir)

    start = time.perf_counter()
    df, dayexposure = gbmpulsarpipe(
        poshistfile=args.poshistfile,
        detectors=args.detectors,
        date=args.date,
        phalow=args.phalow,
        phahigh=args.phahigh,
        srcra=args.srcra,
        srcdec=args.srcdec,
        parfile=args.parfile,
        inputdir=args.inputdir,
        usecrimp=args.usecrimp,
        n_workers=args.n_workers,
        threads_per_worker=args.threads_per_worker
    )

    if df is None:
        sys.exit("No data to process.")

    # Save data to .parquet file
    if args.outputfile is not None:
        df.to_parquet(args.outputfile + '.parquet')

    # Create pulse profiles in several different energy bands
    buildpulseprofiles(df, dayexposure, args.outputdir, args.srcname)

    # Create a pulse profile plot of integrated energy range
    if args.integratedpp is not None:
        binnedProfile = binphases(df['PHASES'].values, nbrBins=50)
        pulseProfile = {
            'ppBins': binnedProfile['ppBins'],
            'ppBinsRange': binnedProfile['ppBinsRange'],
            'countRate': binnedProfile['ctsBins'] / (86400 / 50),  # Approximate rate
            'countRateErr': binnedProfile['ctsBinsErr'] / (86400 / 50)  # Approximate rate
        }
        plotpulseprofile(pulseProfile, outFile=args.integratedpp)

    dt = time.perf_counter() - start
    logger.info(f"Total runtime: {dt:.2f}s")

    # Renaming logger
    logfile = 'gbmpulsar_' + str(args.date) + '_' + str(args.srcra) + '_' + str(args.srcdec) + '.log'
    os.rename('gbmpulsar.log', logfile)

    # Moving created files to prefered directory
    if os.path.abspath(os.getcwd()) != os.path.abspath(args.outputdir):
        shutil.move('./' + logfile, args.outputdir)
        if args.outputfile:
            shutil.move('./' + args.outputfile + '.parquet', args.outputdir)
        if args.integratedpp:
            shutil.move('./' + args.integratedpp + '.pdf', args.outputdir)


if __name__ == '__main__':
    main()
