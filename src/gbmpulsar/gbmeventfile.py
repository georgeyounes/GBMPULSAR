"""
gbmeventfile.py is a module to perform simple operations on a Fermi GBM event file,
reading useful header words, filter for energy etc.
This is a simpler version of the eventfile.py module in CRIMP
"""

import sys

import numpy
import numpy as np
import pandas
import pandas as pd

from astropy.io import fits

from gbmpulsar.gbmpulsar_logging import get_logger

# Log config
############
logger = get_logger(__name__)

sys.dont_write_bytecode = True


class EvtFileOps:
    """
        A class to operate on a fits event file

        Attributes
        ----------
        evtFile : str
            name of the fits event file

        Methods
        -------
        readEF(): reads essential keywords from an event file
        readGTI(): reads GTI table from an event
        filtEneEF(eneLow, eneHigh): filters the event list accoring to energy (in keV)
        filttime(gti): filters the event list accoring to GTI
        """

    def __init__(self, evtFile: str):
        """
        Constructs the necessary attribute

        :param evtFile: name of the fits event file
        :type evtFile: str
        """
        self.evtFile = evtFile

    #################################################################
    def readEF(self):  # Reading fits event file from X-ray satellites
        """
        Reads essential keywords from an event file
        :return: evtFileKeyWords - dictionary of essential keywords
        :rtype: dict
        """
        hdulist = fits.open(self.evtFile)

        # Reading some essential keywords
        TELESCOPE = hdulist['EVENTS'].header['TELESCOP']
        INSTRUME = hdulist['EVENTS'].header['INSTRUME']
        TSTART = hdulist['EVENTS'].header['TSTART']
        TSTOP = hdulist['EVENTS'].header['TSTOP']
        TIMESYS = hdulist['EVENTS'].header['TIMESYS']
        DATEOBS = hdulist['EVENTS'].header['DATE-OBS']
        MJDREF = hdulist['EVENTS'].header['MJDREFI'] + hdulist['EVENTS'].header['MJDREFF']
        DATEEND = hdulist['EVENTS'].header['DATE-END']
        DETNAME = hdulist['EVENTS'].header['DETNAM']
        DATATYPE = hdulist['EVENTS'].header['DATATYPE']

        evtFileKeyWords = {'TELESCOPE': TELESCOPE, 'INSTRUME': INSTRUME, 'TSTART': TSTART,
                           'TSTOP': TSTOP, 'TIMESYS': TIMESYS, 'DATEOBS': DATEOBS,
                           'MJDREF': MJDREF, 'DATEEND': DATEEND, 'DETNAME': DETNAME, 'DATATYPE': DATATYPE}

        return evtFileKeyWords

    #################################################################
    def readGTI(self):  # Reading fits event file GTI lists
        """
        Reads GTI table from an event
        :return: gtiList
        :rtype: numpy.ndarray
        """
        # Reading EF for some necessary keywords
        evtfilekeywords = self.readEF()
        DATATYPE = evtfilekeywords["DATATYPE"]

        hdulist = fits.open(self.evtFile)
        GTIdata = hdulist["GTI"].data
        ST_GTI = GTIdata.field("START")
        ET_GTI = GTIdata.field("STOP")
        gtiList = (np.vstack((ST_GTI, ET_GTI))).T

        df_gti = pd.DataFrame(gtiList, columns=['ST_GTI', 'ET_GTI'])

        if DATATYPE == "TTE":
            logger.warning("Default GTI of GBM TTE file is simply start and end time of event file")

        return evtfilekeywords, df_gti

    ################################################################################
    def filtenergy(self, phalow: float, phahigh: float):  # Filtering event file according to energy
        """
        Filters the event list according to pulse-height (in detector unit, integer 0 to 128)
        :param phalow: low pulse-height cutoff
        :type phalow: float
        :param phahigh: high pulse-height cutoff
        :type phahigh: float
        :return: dataTP_eneFlt, pandas dataframe of TIME and PHA, filtered for pulse-height
        :rtype: pandas.DataFrame
        """
        # Reading columns TIME and PI (pulse-invariant - proxy for photon energy) from binary table
        hdulist = fits.open(self.evtFile)
        tbdata = hdulist['EVENTS'].data

        dataTP = pd.DataFrame(np.vstack((tbdata.field('TIME'), tbdata.field('PHA'))).T, columns=['TIME', 'PHA'])
        dataTP_eneFlt = dataTP.loc[((dataTP['PHA'] >= phalow) & (dataTP['PHA'] <= phahigh))]

        # Close the file handle
        hdulist.close()

        return dataTP_eneFlt

    ################################################################################
    def filttime(self, gti):  # Filtering event file according to energy
        """
        Filters the event list according to GTI
        :param gti: dataframe with START and STOP columns
        :type gti: pandas.DataFrame
        :return: dataTP_timFlt, pandas dataframe of TIME and PI, filtered for GTI
        :rtype: pandas.DataFrame
        """
        # Reading columns TIME and PI (pulse-invariant - proxy for photon energy) from binary table
        hdulist = fits.open(self.evtFile)
        tbdata = hdulist['EVENTS'].data

        dataTP = pd.DataFrame(np.vstack((tbdata.field('TIME'), tbdata.field('PHA'))).T, columns=['TIME', 'PHA'])

        mask = pd.Series(False, index=dataTP.index)
        for _, row in gti.iterrows():
            mask |= dataTP['TIME'].between(row['START'], row['STOP'])
        dataTP_timFlt = dataTP.loc[mask]

        # Close the file handle
        hdulist.close()

        return dataTP_timFlt
