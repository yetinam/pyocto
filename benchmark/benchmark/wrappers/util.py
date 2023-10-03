import logging
import os
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

logger = logging.getLogger("benchmark")


@contextmanager
def set_cwd(path: Path):
    """
    Sets the cwd within the context and resets it afterwards

    Adapted from https://dev.to/teckert/changing-directory-with-a-python-context-manager-2bj8
    """
    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


class StationMapper:
    """
    Many tools only allow fixed-width station names that might be incompatible with some deployments.
    This mapper recodes station names into fixed names and back.
    """

    def __init__(self, stations: pd.DataFrame) -> None:
        self._build_station_dict(stations)

    def _build_station_dict(self, stations: pd.DataFrame) -> None:
        self._station_dict = {
            station: f"S{i:05d}" for i, station in enumerate(stations["id"].unique())
        }
        self._inv_station_dict = {v: k for k, v in self._station_dict.items()}
        if len(self._station_dict) > 100000:
            logger.warning("This module is not designed to handle > 100,000 stations.")

    def translate_station(self, station: str, inv: bool = False) -> str:
        if inv:
            return self._inv_station_dict[station]
        else:
            if station not in self._station_dict:
                i = 1 + max(int(x[1:]) for x in self._station_dict.values())
                self._station_dict[station] = f"S{i:05d}"
                self._inv_station_dict[f"S{i:05d}"] = station
            return self._station_dict[station]
