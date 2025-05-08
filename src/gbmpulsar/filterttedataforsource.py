"""
Filtering TTE data to maximum signal-to-noise for source
Occurs according to  zenToSrcCutAng, zenToDetCutAng, and detToSrcCutAng
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import sys
import os

from scipy.interpolate import splrep, splev

from astropy.coordinates import SkyCoord
from astropy.io import fits

import dask.dataframe as dd
from dask import delayed

from gbmpulsar.gbmeventfile import EvtFileOps

from gdt.missions.fermi.time import Time
from gdt.missions.fermi.gbm.poshist import GbmPosHist
from gdt.core.data_primitives import Gti

from gbmpulsar.gbmpulsar_logging import get_logger

# Log config
############
logger = get_logger(__name__)

sys.dont_write_bytecode = True


def mergettefullday(detector, date, enelow, enehigh, inputdir):
    portions = definegbmtteportions()

    listofallportions = []
    for portion in portions:
        tteFile_arr = glob.glob(os.path.join(inputdir, f'glg_tte_{detector}_{date}_{portion}z*'))
        # Sometimes the day is missing some chunks, not sure why
        if tteFile_arr:
            tteFile = tteFile_arr[0]
        else:
            continue
        listofallportions.append(tteFile)

    if not listofallportions:
        logger.info('No tte file found for detecotor ' + detector)
        return None
    else:
        delayed_filteredttefiles = [delayed(energyfilterttefile)(fp, enelow, enehigh) for fp in listofallportions]
        filteredttefiles_dd = dd.from_delayed(delayed_filteredttefiles)
        return filteredttefiles_dd


def gtifromposhistpersource(poshistFile, detector, srcra, srcdec, poshist_res=10, zentosrccutang=75, zentodetcutang=60,
                            dettosrccutang=50):
    # Reading poshist file
    poshist = GbmPosHist.open(poshistFile)
    # Spacecraft coordinates
    poshist_gbmstates = poshist.get_spacecraft_states()
    # Poshist GTIs - outside of SAA
    poshist_gti = np.array(
        Gti.from_boolean_mask(poshist_gbmstates['time'].value, poshist_gbmstates['good'].value).as_list())
    # Poshist strict GTIs - eliminating +/- 300 second around end/start times of SAA
    poshist_gti_strictSAA = (np.vstack((poshist_gti[:, 0] + 300, poshist_gti[:, 1] - 300))).T

    # Source location
    src_coord = SkyCoord(srcra, srcdec, frame='icrs', unit='deg')

    # Initilizing lists
    full_gti_list = []
    poshist_goodtimes_detang = []
    for gtiindex, poshist_gtiint in enumerate(poshist_gti_strictSAA):
        # creating a TIME array of resolution increments within strict SAA GTIs
        nbrsec_pergti_tmp = int((poshist_gtiint[1] - poshist_gtiint[0]) / poshist_res)
        # Sometimes gti is too small that removing +/- 300 seconds from start/end makes it negative - skip
        if nbrsec_pergti_tmp <= 0:
            continue
        else:
            poshist_gtiint_resstamp = np.unique(np.hstack((poshist_gtiint[0], poshist_gtiint[0] + np.cumsum(
                np.zeros(nbrsec_pergti_tmp) + poshist_res), poshist_gtiint[1])))

        # Moving to source-specific filtering (target-specific filtering):
        # source to Zenith angle > zentosrccutang
        # source to detector angle < dettosrccutang
        # detector to zenith angle > zentodetcutang

        # Source to zenith angle separation
        ####################################
        # 1- Fermi frame at each time-stamp
        poshist_gti_times = Time(poshist_gtiint_resstamp, format='fermi')
        gbmframe = poshist.get_spacecraft_frame()
        one_frame = gbmframe.at(poshist_gti_times)
        # 2- Earth location relative to Fermi at each time-stamp
        earth_coord_infermi = one_frame.geocenter
        # 3- Source location relative to Fermi at each time-stamp
        src_coord_infermi = src_coord.transform_to(one_frame)
        # 4- Separation
        srctoearth_sep_infermi = earth_coord_infermi.separation(src_coord_infermi).deg
        # 5- Separation should be less than zentosrccutang degrees
        zencut = srctoearth_sep_infermi > zentosrccutang

        # Detector to zenith angle separation
        ######################################
        # First let's get detector skycoord in Fermi frame (these are constants)
        # could not think of a more elegant way to do this (quickly)
        if detector == 'n0':
            det_coord_infermi = one_frame.detectors.n0.skycoord(one_frame)
        elif detector == 'n1':
            det_coord_infermi = one_frame.detectors.n1.skycoord(one_frame)
        elif detector == 'n2':
            det_coord_infermi = one_frame.detectors.n2.skycoord(one_frame)
        elif detector == 'n3':
            det_coord_infermi = one_frame.detectors.n3.skycoord(one_frame)
        elif detector == 'n4':
            det_coord_infermi = one_frame.detectors.n4.skycoord(one_frame)
        elif detector == 'n5':
            det_coord_infermi = one_frame.detectors.n5.skycoord(one_frame)
        elif detector == 'n6':
            det_coord_infermi = one_frame.detectors.n6.skycoord(one_frame)
        elif detector == 'n7':
            det_coord_infermi = one_frame.detectors.n7.skycoord(one_frame)
        elif detector == 'n8':
            det_coord_infermi = one_frame.detectors.n8.skycoord(one_frame)
        elif detector == 'n9':
            det_coord_infermi = one_frame.detectors.n9.skycoord(one_frame)
        elif detector == 'na':
            det_coord_infermi = one_frame.detectors.na.skycoord(one_frame)
        elif detector == 'nb':
            det_coord_infermi = one_frame.detectors.nb.skycoord(one_frame)
        elif detector == 'b0':
            det_coord_infermi = one_frame.detectors.b0.skycoord(one_frame)
        elif detector == 'b1':
            det_coord_infermi = one_frame.detectors.b1.skycoord(one_frame)
        else:
            raise Exception('Unknown detector {}'.format(detector))

        # Angle seperation between detector pointing and earth geocentric location
        dettoearth_sep_infermi = earth_coord_infermi.separation(det_coord_infermi).deg
        # Above separation should be less than zentodetcutang (default = 60) degrees
        zendetcut = dettoearth_sep_infermi > zentodetcutang

        # Source to Detector angle separation
        ######################################
        dettosource_sep_infermi = src_coord_infermi.separation(det_coord_infermi).deg
        dettosrccut = dettosource_sep_infermi < dettosrccutang

        # Merging cuts
        ##############
        allcuts = np.array([all(tup) for tup in zip(zencut, zendetcut, dettosrccut)])

        # Retrieving detector to source angle for each good poshist timestamp
        good_det_src_ang = dettosource_sep_infermi[allcuts]

        # Creating GTIs from all these cuts
        ###################################
        # Find indices where the boolean value is 1
        indices = np.where(allcuts == 1)[0]

        if indices.size == 0:
            # No timestamp within inital GTI satisfied good observing time
            continue
        else:
            # Group contiguous indices together.
            # np.diff(indices) computes the difference between adjacent indices.
            # Wherever the difference is not 1, we have a break between consecutive 1's.
            groups = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)
            # For each group, extract the start time and the stop time.
            intervals = [(poshist_gtiint_resstamp[group[0]], poshist_gtiint_resstamp[group[-1]]) for group in groups]

        gti_df = pd.DataFrame(intervals, columns=['START', 'STOP'])
        full_gti_list.append(gti_df)

        # Creating the dataframe of good poshist timestamps and corresponding detector to source angle
        poshist_goodtimes_detang_pergti = pd.DataFrame({'POSHIST_TIMES': poshist_gtiint_resstamp[allcuts],
                                                        'DET_SRC_ANG': good_det_src_ang})
        poshist_goodtimes_detang.append(poshist_goodtimes_detang_pergti)

    if not full_gti_list:
        full_gti_df = pd.DataFrame()
        poshist_goodtimes_detang = pd.DataFrame()
        return full_gti_df, poshist_goodtimes_detang
    else:
        # GTI
        full_gti_df = pd.concat(full_gti_list, ignore_index=True)
        full_gti_df['DET'] = detector  # Adding detector for bookkeeping
        # All good time-stamps and corresponding detector angles - if i ever decide to add weights
        poshist_goodtimes_detang = pd.concat(poshist_goodtimes_detang, ignore_index=True)
        poshist_goodtimes_detang['DET'] = detector  # Adding detector for bookkeeping
        # Logging important info
        logger.info("Calculated Good Time Intervals for:\n"
                    "PosHist file: {}\n "
                    "Detector: {}\n"
                    "Source RA: {}\n"
                    "Source DEC: {}\n"
                    "GTI:\n{}\n"
                    "Total exposure {}\n".format(poshistFile, detector, srcra, srcdec, full_gti_df,
                                               np.sum(full_gti_df['STOP']-full_gti_df['START'])))

        return full_gti_df, poshist_goodtimes_detang


def energyfilterttefile(tteevtfile, enelow, enehigh):
    """
    Simple wrapper function to filtering energy of a tte event file
    :param tteevtfile: name of tte event file
    :type tteevtfile: str
    :param enelow: low energy (PHA) for filtering
    :type enelow: int
    :param enehigh: high energy (PHA) for filtering
    :type enehigh: int
    :return: filtered dataframe of TIME and PHA values
    :rtype: pandas.DataFrame
    """
    evt_ops = EvtFileOps(tteevtfile)
    return evt_ops.filtenergy(phalow=enelow, phahigh=enehigh)


def timefilterttedf(tte_df, gti):
    """
    Simple wrapper function to filter time of a tte event file according to a GTI built from gtifromposhistpersource
    :param tte_df: TTE TIME, PI, DET dataframe
    :type tte_df: pandas.DataFrame
    :param gti: START and STOP dataframe
    :type gti: pandas.DataFrame
    :return tte_df: time-filtered TIME, PI, DET dataframe
    :rtype: pandas.DataFrame
    """
    mask = pd.Series(False, index=tte_df.index)
    for _, row in gti.iterrows():
        mask |= tte_df['TIME'].between(row['START'], row['STOP'])
    return tte_df.loc[mask].copy()


def barycenter_ttetimes(poshist_alltimes_df, ttetimes):
    # Create a spline fit to the barycentered data
    tck = splrep(poshist_alltimes_df['TIME'], poshist_alltimes_df['BARY_TIME'], s=0, k=3)
    # Interpolate the ttetimes
    ttetimes_bary = splev(ttetimes.values, tck)
    return pd.Series(ttetimes_bary, index=ttetimes.index)


def readposhistforttebarycenter(poshistfile):
    hdulist = fits.open(poshistfile)
    poshist_alltimes_df = pd.DataFrame(np.vstack((hdulist[1].data.field('SCLK_UTC'),
                                                  hdulist[1].data.field('BARY_TIME'))).T,
                                       columns=['TIME', 'BARY_TIME'])
    return poshist_alltimes_df


def definegbmtteportions():
    """
    Define the GBM detectors
    :return: fermigbmdetectors
    :rtype: numpy array
    """
    gbmtteportions = (["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13",
                       "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"])
    return gbmtteportions


def definegbmdetectors():
    """
    Define the GBM detectors
    :return: fermigbmdetectors
    :rtype: numpy array
    """
    gbmnaidetectors = (["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "na", "nb"])
    gbmbgodetectors = (["b0", "b1"])
    return gbmnaidetectors, gbmbgodetectors
