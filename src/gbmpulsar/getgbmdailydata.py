#!/usr/bin/env python3
"""
This code downloads GBM CTTE (TTE), CSPEC, CTIME, or poshist daily data given a full date
Be aware CTTE data is quite large and requires plenty of disk space
"""
from ftplib import FTP_TLS
import argparse
import os
import fnmatch


class FermiGBMDownloader:
    """
    A class to download Fermi GBM daily data via FTPS.

    Attributes:
        host (str): FTP server hostname.
        base_dir (str): Base directory for GBM daily data on the FTP server.
        ftp (ftplib.FTP_TLS): FTP connection object.
        output_dir (str): Local directory to save downloaded files.
    """

    def __init__(self,
                 host="heasarc.gsfc.nasa.gov",
                 base_dir="/FTP/fermi/data/gbm/daily",
                 output_dir=None):
        self.host = host
        self.base_dir = base_dir
        # Set output directory (default: current working directory)
        self.output_dir = output_dir or os.getcwd()
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        self.ftp = None

    def connect(self):
        """Connect to the FTP server using FTP over TLS (FTPS)."""
        try:
            print(f"Connecting to FTP server: {self.host} ...")
            self.ftp = FTP_TLS(self.host)
            self.ftp.login()  # anonymous login
            self.ftp.prot_p()  # secure data connection
            print("Connected successfully using TLS!")
        except Exception as e:
            print(f"Error connecting to FTP server: {e}")
            exit(1)

    def disconnect(self):
        """Disconnect from the FTP server."""
        if self.ftp:
            self.ftp.quit()
            print("Disconnected from FTP server.")

    def download_files(self, fullDate, dataType="ctime", det="all", portion="all"):
        """
        Change to the remote directory for the given date, list files,
        filter by the pattern, skip those already present, and download the rest.
        """
        remote_dir = f"{self.base_dir}/{fullDate}/current"
        try:
            print(f"Changing to remote directory: {remote_dir}")
            self.ftp.cwd(remote_dir)
        except Exception as e:
            print(f"Error accessing remote directory {remote_dir}: {e}")
            return

        pattern = build_pattern(fullDate, dataType, det, portion)
        print(f"Using file pattern: {pattern}")

        try:
            files = self.ftp.nlst()
        except Exception as e:
            print(f"Error listing files in directory: {e}")
            return

        # Filter by pattern
        matching = [f for f in files if fnmatch.fnmatch(f, pattern)]
        if not matching:
            print("No files found matching the given criteria.")
            return

        # Skip files already in output_dir
        to_download = [f for f in matching if not os.path.exists(os.path.join(self.output_dir, f))]
        skipped = set(matching) - set(to_download)
        if skipped:
            print(f"Skipping {len(skipped)} file(s) already present: {sorted(skipped)}")
        if not to_download:
            print("All matching files already exist locally. Nothing to do.")
            return

        print(f"Downloading {len(to_download)} new file(s): {to_download}")
        for fname in to_download:
            self.download_file(fname)

    def download_file(self, filename):
        """
        Download a single file from the current remote directory,
        unless it already exists.
        """
        local_filename = os.path.join(self.output_dir, filename)
        if os.path.exists(local_filename):
            print(f"[SKIP] {filename} already exists locally.")
            return

        print(f"Downloading file: {filename}")
        try:
            with open(local_filename, 'wb') as f:
                self.ftp.retrbinary(f"RETR {filename}", f.write)
            print(f"Downloaded {filename} successfully!")
        except Exception as e:
            print(f"Error downloading {filename}: {e}")


def build_obsid(fullDate):
    """
    Convert a date string from "YYYY/MM/DD" to "YYMMDD".
    """
    date_str = "".join(fullDate.split("/"))
    if len(date_str) < 6:
        raise ValueError("Invalid date format, expected YYYY/MM/DD")
    return date_str[2:]


def build_pattern(fullDate, dataType, det, portion):
    """
    Build the filename search pattern based on user parameters.
    """
    obs_id = build_obsid(fullDate)
    dataType = dataType.lower()
    det = det.lower()
    portion = portion.lower()

    if dataType == 'tte':
        detector_part = "*" if det == 'all' else det
        if portion == 'all':
            pattern = f"glg_{dataType}_{detector_part}_{obs_id}*"
        else:
            pattern = f"glg_{dataType}_{detector_part}_{obs_id}_{portion}*"
    elif dataType == 'poshist':
        pattern = f"glg_{dataType}_all_{obs_id}*"
    else:
        detector_part = "*" if det == 'all' else det
        pattern = f"glg_{dataType}_{detector_part}_{obs_id}*"
    return pattern


def main():
    parser = argparse.ArgumentParser(
        description="Download Fermi GBM daily data (CTIME, CSPEC, TTE, or POSHIST) using FTP."
    )
    parser.add_argument("fullDate",
                        help="Date in the format YYYY/MM/DD, e.g., 2018/03/01",
                        type=str)
    parser.add_argument("-t", "--dataType",
                        help="Data type to download [ctime, cspec, tte, poshist] (default: ctime)",
                        type=str, default="ctime")
    parser.add_argument("-d", "--det",
                        help="Detector [all, b0, b1, n0–n9, na, nb] (default: all)",
                        type=str, default="all")
    parser.add_argument("-p", "--portion",
                        help="For tte data: which hourly portion 00–23 (default: all)",
                        type=str, default="all")
    parser.add_argument("-o", "--outdir",
                        help="Directory to save downloaded files (default: current directory)",
                        type=str,
                        default=os.getcwd())
    args = parser.parse_args()

    downloader = FermiGBMDownloader(output_dir=args.outdir)
    downloader.connect()
    downloader.download_files(args.fullDate, args.dataType, args.det, args.portion)
    downloader.disconnect()


if __name__ == '__main__':
    main()
