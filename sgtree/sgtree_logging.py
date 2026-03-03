import os
import sys
import glob
import datetime

from sgtree.config import Config


def write_logfile(cfg: Config, timings: dict):
    """Write runtime statistics to logfile."""
    logfile = os.path.join(
        cfg.outdir,
        "logfile_" + cfg.start_time[:-7].replace(" ", "_").replace("-", "_") + ".txt",
    )
    sep = "=" * 80

    try:
        with open(logfile, "w") as f:
            f.write(f"{cfg.outdir}\n{sep}\n")
            f.write(f"Sg_Tree v.2\nstart time: {cfg.start_time}\n{sep}\n")
            f.write(f"Genomes database {cfg.genomedir} contains {cfg.genome_count} genomes\n{sep}\n")
            f.write(f"Marker database {cfg.modeldir} contains {cfg.model_file_count} models\n{sep}\n")
            f.write(f"{sys.argv}\n{sep}\n")

            for step_name, (timestamp, runtime) in timings.items():
                f.write(f"{timestamp}- {step_name}")
                f.write(f"\n{step_name} done - runtime: {runtime} seconds\n{sep}\n")

            f.write(f"Sgtree endtime: \n{datetime.datetime.now()}\n")
            f.write(f"{cfg.genome_count} genomes\n{cfg.model_file_count} models\n")
            f.write(f"{cfg.start_time}\n{datetime.datetime.now()}\n{sep}\n")
    except Exception:
        print("Error writing logfile", sys.exc_info())


def append_logfile(cfg: Config, timings: dict):
    """Append marker selection timings to existing logfile."""
    logfile = os.path.join(
        cfg.outdir,
        "logfile_" + cfg.start_time[:-7].replace(" ", "_").replace("-", "_") + ".txt",
    )
    sep = "=" * 80

    try:
        with open(logfile, "a") as f:
            for step_name, (timestamp, runtime) in timings.items():
                f.write(f"{timestamp}- {step_name}")
                f.write(f"\n{step_name} done - runtime: {runtime} seconds\n{sep}\n")

            f.write(f"Sgtree start, endtime (with marker selection): \n")
            f.write(f"{cfg.start_time}\n{datetime.datetime.now()}\n")
    except Exception:
        print("Error writing logfile", sys.exc_info())
