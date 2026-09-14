"""Write human-readable run settings and stage timings beside pipeline results."""

import logging
import sys
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

from sgtree._version import DISPLAY_VERSION
from sgtree.config import Config

logger = logging.getLogger("sgtree")


def write_logfile(cfg: Config, timings: Mapping[str, tuple[datetime, float]]) -> None:
    """Write runtime statistics to logfile."""
    logfile = _logfile_path(cfg)
    sep = "=" * 80

    try:
        with logfile.open("w", encoding="utf-8") as f:
            f.write(f"{cfg.outdir}\n{sep}\n")
            f.write(f"{DISPLAY_VERSION}\nstart time: {cfg.start_time}\n{sep}\n")
            f.write(
                f"Genomes database {cfg.genomedir} contains "
                f"{cfg.genome_count} genomes\n{sep}\n"
            )
            f.write(
                f"Marker database {cfg.modeldir} contains "
                f"{cfg.model_file_count} models\n{sep}\n"
            )
            f.write(f"{sys.argv}\n{sep}\n")

            for step_name, (timestamp, runtime) in timings.items():
                f.write(f"{timestamp}- {step_name}")
                f.write(f"\n{step_name} done - runtime: {runtime} seconds\n{sep}\n")

            end_time = datetime.now(timezone.utc).astimezone().replace(tzinfo=None)
            f.write(f"Sgtree endtime: \n{end_time}\n")
            f.write(f"{cfg.genome_count} genomes\n{cfg.model_file_count} models\n")
            f.write(f"{cfg.start_time}\n{end_time}\n{sep}\n")
    except OSError:
        logger.exception("Error writing logfile")
        print("Error writing logfile", sys.exc_info())


def append_logfile(cfg: Config, timings: Mapping[str, tuple[datetime, float]]) -> None:
    """Append marker selection timings to existing logfile."""
    logfile = _logfile_path(cfg)
    sep = "=" * 80

    try:
        with logfile.open("a", encoding="utf-8") as f:
            for step_name, (timestamp, runtime) in timings.items():
                f.write(f"{timestamp}- {step_name}")
                f.write(f"\n{step_name} done - runtime: {runtime} seconds\n{sep}\n")

            f.write("Sgtree start, endtime (with marker selection): \n")
            end_time = datetime.now(timezone.utc).astimezone().replace(tzinfo=None)
            f.write(f"{cfg.start_time}\n{end_time}\n")
    except OSError:
        logger.exception("Error writing logfile")
        print("Error writing logfile", sys.exc_info())


def _logfile_path(cfg: Config) -> Path:
    timestamp = cfg.start_time.split(".")[0].replace(" ", "_").replace("-", "_")
    return Path(cfg.outdir) / f"logfile_{timestamp}.txt"
