"""
Module to calculate phases of GBM tte times - relies on PINT
"""
# Suppress PINT extended-precision RuntimeWarnings globally
import warnings
warnings.filterwarnings(
    "ignore",
    message="This platform does not support extended precision floating-point"
)
import pint.models
import pint.residuals
import pint.polycos as polycos
pint.logging.setup(level="WARNING")

import sys
import pandas as pd
import numpy as np

from gdt.missions.fermi.time import Time

from gbmpulsar.gbmpulsar_logging import get_logger

# Log config
############
logger = get_logger(__name__)

sys.dont_write_bytecode = True


class CalcttePhases:
    """
        A class to operate on an event file dataframe containing at least TIME, PHA

        Attributes
        ----------
        tteevent_df : pandas.DataFrame
            name of the [TIME, PHA] pandas dataframe

        Methods
        -------
        calcttephases(): reads essential keywords from an event file
        fermittetomjd(): filters the event list accoring to energy (in keV)
        """

    def __init__(self, tteevent_df: pd.DataFrame):
        """
        Constructs the necessary attribute

        :param tteevent_df: name of the event file dataframe (create with EvtFileOps.filtenergy in gbmeventfile module)
        :type tteevent_df: pandas.DataFrame
        """
        self.tteevent_df = tteevent_df

    def calcttephases(self, parfile, allow_T2=True, allow_tcb=True):
        logger.info("Using PINT photonphase (with polycos) to get tte pulse phases")

        # Polycos parameters
        segLength = 120  # in minutes
        ncoeff = 10
        obsfreq = 0

        tteeventwithmjd_df = self.fermittetomjd()
        ttemjds = tteeventwithmjd_df['BARY_TIME_MJD'].to_numpy()

        if ttemjds.size == 0:
            # No rows ⇒ no phases
            return np.empty(0, dtype=float)

        minmjd, maxmjd = ttemjds.min(), ttemjds.max()

        # Read in model
        modelin = pint.models.get_model(parfile, allow_T2=allow_T2, allow_tcb=allow_tcb)

        telescope_n = "@"
        p = polycos.Polycos()
        ptable = p.generate_polycos(
            modelin, minmjd, maxmjd, telescope_n, segLength, ncoeff, obsfreq
        )

        phases = ptable.eval_phase(ttemjds)
        phases[phases < 0] += 1.0

        return phases

    def fermittetomjd(self):
        self.tteevent_df.loc[:, 'BARY_TIME_MJD'] = Time(self.tteevent_df['BARY_TIME'], format='fermi', scale='tdb').mjd
        return self.tteevent_df

