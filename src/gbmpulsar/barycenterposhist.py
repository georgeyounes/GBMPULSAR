"""
A module that barycenters GBM poshist data - note that the GBM team does not provide a barycentering tool,
for this reason, we use the LAT provided one, gtbary. gtbary is much faster than a python based implementation
(e.g. astropy), and has been extensively tested in the past (it is also based on the HEASOFT tools axbary and barycorr).

We barycenter a poshist file rather than a GBM TTE file because the latter is millions of time-tagged events (rows),
compared to approximately 86400 timestamps (daily poshist file at the 1 second resolution).

The steps to achive the above are a little unconventional:
- Swap the time column of a placeholder LAT event file with the timestamps from a GBM daily poshist file
- Barycenter correct these times using gtbary
- Add the resulting barycentered timestamps as a BARY_TIME column to the poshist .fits file

This module provides a script to achieve the above, called "barycenterposhist" which the user can run from the command
line.

Barycentering a TTE datafile is done elsewhere by reading-in the non-barycentered and the barycentered poshist
timestamps, creating a spline fit to all data points, and using a linear interpolation between pairs of data points to
estimate the barycenter correction at each TTE event timestamp

Warnings: you must be in the fermitools (fermi LAT data analysis software) environment to be able to run "gtbary" and
by extension "barycenterposhist". fermitools does not talk well with astro-gdt-fermi (the fermi GBM data analysis
software). Hence, the recommendation for the time being is to barycenter whatever poshist you need separately, and then
run the gbmpulsar sofware on that data, i.e., do not try to wrap the two together in a nice little function or script,
unless you know what you are doing.
"""

import os
import argparse

from astropy.io import fits
from astropy.table import Table, Column


def barycenterposhist(poshitfile, latplaceholder, latspacecraftfile, srcra, srcdec):
    # Read time column from poshistfile
    hdulist_poshist = fits.open(poshitfile, mode='update')
    poshist_hdrEvents = hdulist_poshist['GLAST POS HIST'].header
    poshist_times = hdulist_poshist['GLAST POS HIST'].data.field("SCLK_UTC")

    # Copy LAT placeholder file to a tmp file in the current working directory
    latplaceholder_tmp = 'tmp_' + latplaceholder
    os.system('cp ' + latplaceholder + ' ./' + latplaceholder_tmp)

    # Open the latplaceholder file in update mode
    hdulEF_LAT = fits.open(latplaceholder_tmp, mode='update')
    hdrEvents = hdulEF_LAT['EVENTS'].header

    # Reading EVENT table, and replacing the TIME column and adding counts columns
    evtTable = Table.read(latplaceholder_tmp, format='fits', hdu='EVENTS')

    # Replacing TIME column
    evtTable['TIME'] = poshist_times[:]

    # Updating TSTART and TSTOP and adding exposure to header
    hdrEvents['EXPOSURE'] = (poshist_times[-1] - poshist_times[0], 'Full day exposure')
    hdrEvents['TSTART'] = (poshist_times[0], 'START of POSHIST file in MET')
    hdrEvents['TSTOP'] = (poshist_times[-1], 'END of POSHIST file in MET')

    # Updating event file with new EVENTS table and header keywords
    newhdulEF = fits.BinTableHDU(data=evtTable, header=hdrEvents, name='EVENTS')
    fits.update(latplaceholder_tmp, newhdulEF.data, newhdulEF.header, 'EVENTS')

    # Barycenter data with gtbary and add phase column according to a timing solution
    outputbaryfile = 'barycentered_' + poshitfile
    command = ('gtbary evfile=' + latplaceholder_tmp + ' scfile=' + latspacecraftfile + ' outfile=' + outputbaryfile +
               ' ra=' + str(srcra) + ' dec=' + str(srcdec))
    os.system(command)

    # Read barycentered time column from outputbaryfile and add it to the poshist file
    hdulbaryposhist = fits.open(outputbaryfile)
    barycenteredposhisttimes = hdulbaryposhist['EVENTS'].data.field("TIME")

    #####################################
    # Creating an astropy Table that corresponds to the EVENTS table
    baryTable = Table.read(poshitfile, format='fits', hdu='GLAST POS HIST')
    barytime = Column(name='BARY_TIME', data=barycenteredposhisttimes, dtype=float)
    baryTable.add_column(barytime)

    # Updating event file
    newhdulEFPH = fits.BinTableHDU(data=baryTable, header=poshist_hdrEvents, name='GLAST POS HIST')
    fits.update(poshitfile, newhdulEFPH.data, newhdulEFPH.header, 1)  # cannot specify name 'GLAST POS HIST'

    # remove secondary files
    os.system("rm -rf " + outputbaryfile + ' ' + latplaceholder_tmp)

    return poshitfile


def main():
    parser = argparse.ArgumentParser(description='barycenter correct TIME column in GBM poshist file')
    parser.add_argument("poshitfile", help="Posthist file name", type=str)
    parser.add_argument("srcra", help="Source right ascension", type=float)
    parser.add_argument("srcdec", help="Source declination", type=float)
    parser.add_argument("latplaceholder", help="LAT placeholder event file", type=str)
    parser.add_argument("latspacecraftfile", help="LAT spacecraft file", type=str)
    args = parser.parse_args()

    barycenterposhist(args.poshitfile, args.latplaceholder, args.latspacecraftfile, args.srcra, args.srcdec)

    return


if __name__ == "__main__":
    main()
