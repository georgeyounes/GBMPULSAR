"""
Module to work with daily gbm .csv pulse profiles. The latter are built with gbmpulsarpipe module.
The main functions are:
1- mergecsvpulseprofiles(): merge the CSV pulse profiles over specific dates and channels (energy range)
2- waterfallpulseprofile(): merge data over specific channel, but bin over time - produce a waterfall plot

The second option exists as a command line tool: "gbmmergepulseprofiles -h" for more information
"""
import os

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import argparse

from gbmpulsar.gbmpulsar_logging import get_logger

# Log config
############
logger = get_logger(__name__)


def channel_to_range(channel_indices):
    # Define bin edges
    energybins = [8, 10, 30, 50, 100, 200, 300, 500, 1000]  # keV (bins are [start, end))
    phabins = [4, 7, 22, 33, 51, 73, 86, 104, 128]

    # Define valid channel indices
    accepted_indices = list(range(len(phabins) - 1))  # [0, 1, ..., 7]

    # Handle 'all' option
    if channel_indices == "all":
        channel_indices = accepted_indices
    elif isinstance(channel_indices, list):
        channel_indices = [int(i) for i in channel_indices]

    # Validate channel indices
    if not all(i in accepted_indices for i in channel_indices):
        logger.error(f"Invalid channel index found. Accepted values are: {accepted_indices}")
        raise ValueError(f"Invalid channel index found. Accepted values are: {accepted_indices}")

    # Return user-selected pharange, energyrange, channel_indices, and full phabins
    pharange = [f"{phabins[i]}-{phabins[i + 1]}" for i in channel_indices]
    energyrange = [f"{energybins[i]}-{energybins[i + 1]}" for i in channel_indices]

    return pharange, energyrange, channel_indices, phabins


def parse_date_folder(foldername):
    return datetime.strptime(foldername, "%y%m%d")


def get_date_range(start_date, end_date):
    start = datetime.strptime(start_date, "%Y-%m-%d")
    end = datetime.strptime(end_date, "%Y-%m-%d")
    return [(start + timedelta(days=i)).strftime("%y%m%d") for i in range((end - start).days + 1)]


def mergecsvpulseprofiles(root_dir, start_date, end_date, channels):
    date_folders = get_date_range(start_date, end_date)
    _, _, channel_ranges, _ = channel_to_range(channels)

    daily_profiles = []
    for date_folder in date_folders:
        folder_path = Path(root_dir) / date_folder
        if not folder_path.exists():
            continue

        dfs = []
        for f in folder_path.glob("pulseprofile_*_pha.csv"):
            if any(f"{rng}_pha.csv" in f.name for rng in channel_ranges):
                df = pd.read_csv(f)
                dfs.append(df)

        if not dfs:
            continue

        # Combine for the day
        df_day = dfs[0][["ppBins", "ppBinsRange"]].copy()
        count_rates = [df["countRate"] for df in dfs]
        count_errs = [df["countRateErr"] ** 2 for df in dfs]

        df_day["countRate"] = np.sum(count_rates, axis=0)
        df_day["countRateErr"] = np.sqrt(np.sum(count_errs, axis=0))

        daily_profiles.append(df_day)

    if not daily_profiles:
        raise ValueError("No valid profiles found in the specified date range and channels.")

    # Average over days
    sum_df = daily_profiles[0][["ppBins", "ppBinsRange"]].copy()
    count_rates_all = np.array([df["countRate"].values for df in daily_profiles])
    count_errs_all = np.array([df["countRateErr"].values for df in daily_profiles])

    sum_df["countRate"] = np.mean(count_rates_all, axis=0)
    sum_df["countRateErr"] = np.sqrt(np.sum(count_errs_all ** 2, axis=0)) / len(date_folders)  # latter is number of days

    # Also return as dict
    profile_dict = sum_df.to_dict(orient="series")

    return sum_df, profile_dict


def rebin_pulseprofile(df, num_bins):
    if num_bins >= len(df):
        logger.info(f"Number of bins in pulse profile must be less than number of original bins (= {len(df)}) - "
                    f"setting the number of bins to original binning")
        num_bins = len(df)

    # create new binning
    bin_edges = (np.linspace(0, 1, num_bins + 1))  # + (1 / num_bins) / 2
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_widths = np.diff(bin_edges)

    binned_countRate = []
    binned_countRateErr = []

    for i in range(num_bins):
        in_bin = (df["ppBins"] >= bin_edges[i]) & (df["ppBins"] < bin_edges[i + 1])
        bin_cr = df.loc[in_bin, "countRate"].sum() / in_bin.sum()
        bin_err = np.sqrt((df.loc[in_bin, "countRateErr"] ** 2).sum()) / in_bin.sum()

        binned_countRate.append(bin_cr)
        binned_countRateErr.append(bin_err)

    rebinned_df = pd.DataFrame({
        "ppBins": bin_centers,
        "ppBinsRange": bin_widths,
        "countRate": binned_countRate,
        "countRateErr": binned_countRateErr
    })

    return rebinned_df


def rebin_pulseprofiles_over_time(daily_profiles, daily_labels, time_binsize):
    if len(daily_profiles) < time_binsize:
        raise ValueError("Time bin size is larger than the number of days available.")

    num_bins = len(daily_profiles) // time_binsize
    binned_profiles = []
    binned_labels = []

    for i in range(num_bins):
        start = i * time_binsize
        end = start + time_binsize
        avg_profile = np.mean(daily_profiles[start:end], axis=0)
        avg_date = daily_labels[start:end].mean()
        binned_profiles.append(avg_profile)
        binned_labels.append(avg_date)

    return np.array(binned_profiles), [d.strftime("%Y-%m-%d") for d in binned_labels]


def waterfallpulseprofile(root_dir, start_date, end_date, channels="all", pp_rebin=None, time_binsize=None,
                          dispplot=False, plotname=None):
    # Number of days to merge through
    totnbrdays = len(get_date_range(start_date, end_date))

    # Get all date folders
    date_folders = get_date_range(start_date, end_date)
    _, _, channel_ranges, _ = channel_to_range(channels)

    daily_profiles = []
    daily_labels = []

    for date_folder in date_folders:
        folder_path = Path(root_dir) / date_folder
        if not folder_path.exists():
            continue

        dfs = []
        for f in folder_path.glob("pulseprofile_*_pha.csv"):
            if any(f"{rng}_pha.csv" in f.name for rng in channel_ranges):
                df = pd.read_csv(f)
                dfs.append(df)

        if not dfs:
            continue

        # Combine channels for the day
        df_day = dfs[0][["ppBins", "ppBinsRange"]].copy()
        count_rates = [df["countRate"] for df in dfs]
        count_errs = [df["countRateErr"] ** 2 for df in dfs]

        df_day["countRate"] = np.sum(count_rates, axis=0)
        df_day["countRateErr"] = np.sqrt(np.sum(count_errs, axis=0))

        if pp_rebin:
            df_day = rebin_pulseprofile(df_day, pp_rebin)

        daily_profiles.append(df_day["countRate"].values)
        daily_labels.append(datetime.strptime(date_folder, "%y%m%d"))

    if not daily_profiles:
        raise ValueError("No valid profiles found in the specified date range and channels.")

    daily_profiles = np.array(daily_profiles)
    daily_labels = np.array(daily_labels)

    # Rebin in time
    if time_binsize:
        data, y_labels = rebin_pulseprofiles_over_time(daily_profiles, daily_labels, time_binsize)
    else:
        data = daily_profiles
        y_labels = [d.strftime("%Y-%m-%d") for d in daily_labels]

    x_vals = df_day["ppBins"].values

    # Merge all profiles for the top panel
    df_merged, _ = mergecsvpulseprofiles(root_dir, start_date, end_date, channels)
    if pp_rebin:
        df_merged = rebin_pulseprofile(df_merged, pp_rebin)

    plot_waterfall_with_merged_profile(df_merged=df_merged, data=data, x_vals=x_vals, y_labels=y_labels,
                                       time_binsize=time_binsize, plotname=plotname, dispplot=dispplot)

    return data, x_vals, y_labels  # Waterfall data


def plot_waterfall_with_merged_profile(df_merged, data, x_vals, y_labels, time_binsize=None, plotname=None,
                                       dispplot=False):
    """
    Plot a two-panel figure:
      - Top: Integrated pulse profile (from mergecsvpulseprofiles)
      - Bottom: Waterfall plot over time
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), height_ratios=[1, 3], sharex=True)

    # Top panel: Integrated pulse profile
    #ax1.plot(df_merged["ppBins"], df_merged["countRate"], color="black", lw=1.5)
    ax1.step(df_merged["ppBins"], df_merged["countRate"], 'k+-', where='mid')
    ax1.errorbar(df_merged["ppBins"], df_merged["countRate"], yerr=df_merged["countRateErr"], fmt='+k', markersize=1)

    ax1.set_ylabel("Integrated\ncountRate")
    ax1.grid(True)

    # Bottom panel: Waterfall
    y_vals = np.arange(len(data))
    ax2.imshow(data, aspect='auto', origin='lower',
               extent=[x_vals[0], x_vals[-1], 0, len(y_vals)],
               cmap='viridis')
    ax2.set_yticks(np.arange(len(y_labels)) + 0.5)
    ax2.set_yticklabels(y_labels)
    ax2.set_xlabel("Phase (cycles)")
    ax2.set_ylabel("Date" if not time_binsize else f"{time_binsize}-day bins")
    # fig.colorbar(c, ax=ax2, label="countRate")

    fig.tight_layout()

    if plotname:
        plt.savefig(plotname + '.pdf', dpi=200, bbox_inches='tight', format='pdf')
        logger.info(f"Waterfall plot saved to {plotname}")
    if dispplot:
        plt.show()
    else:
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="Process and visualize GBM daily pulse profiles")
    parser.add_argument("start_date", type=str, help="Start date in YYYY-MM-DD format.")
    parser.add_argument("end_date", type=str, help="End date in YYYY-MM-DD format.")

    parser.add_argument("-rd", "--root_dir", type=str, help="Root directory with date-named folders "
                                                            "(default = current directory)", default=os.getcwd())
    parser.add_argument("-ch", "--channels", nargs="+", help="List of channel indices (0–7) or 'all' "
                                                             "(default)", default="all")
    parser.add_argument("-pr", "--pp_rebin", type=int, help="Number of bins to rebin pulse profiles into "
                                                            "(default = none)", default=None)
    parser.add_argument("-tb", "--time_binsize", type=int, help="Number of days to average together in "
                                                                "time (default = none)", default=None)
    parser.add_argument("-pn", "--plotname", type=str, help="Filename to save the waterfall plot "
                                                            "(default = none)", default=None)
    parser.add_argument("-dp", "--dispplot", action="store_true", help="Display the plot interactively")

    args = parser.parse_args()

    # Handle 'all' channels properly
    channels = args.channels
    if channels == "all":
        channels = "all"
    else:
        channels = [int(c) for c in channels]

    # Call the waterfall function
    try:
        waterfallpulseprofile(root_dir=args.root_dir, start_date=args.start_date, end_date=args.end_date,
                              channels=channels, pp_rebin=args.pp_rebin, time_binsize=args.time_binsize,
                              dispplot=args.dispplot, plotname=args.plotname)

    except Exception as e:
        logger.error(f"Error: {e}")
        raise


if __name__ == "__main__":
    main()
