from typing import Any

import numpy as np
import pandas as pd

import pyocto

from .base import AbstractAssociator


class PyOctoAssociator(AbstractAssociator):
    def __init__(self, config: dict[str, Any], **kwargs):
        self.associator = pyocto.OctoAssociator(**config)

        super().__init__(**kwargs)

    def get_events(
        self, picks: pd.DataFrame, stations: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        if len(picks) == 0:
            return pd.DataFrame(), pd.DataFrame()

        picks_org = picks.copy(deep=True)
        picks_org["idx"] = np.arange(len(picks_org))

        picks = picks.copy(deep=True)
        stations = stations.copy(deep=True)
        picks["time"] = picks["time"].apply(lambda x: x.timestamp())
        stations = self.transform_station_coords(stations)
        stations.rename(
            columns={"x(km)": "x", "y(km)": "y", "z(km)": "z"}, inplace=True
        )
        stations["z"] = -stations["z"]

        catalog, assignments = self.associator.associate(picks, stations)

        if len(catalog) > 0:
            catalog.rename(
                columns={"x": "x(km)", "y": "y(km)", "z": "z(km)"}, inplace=True
            )
            catalog = self.transform_event_coords(catalog)
            catalog.sort_values("time", inplace=True)

        assignments = pd.merge(
            assignments, picks_org, left_on="pick_idx", right_on="idx", validate="m:1"
        )
        assignments.drop(columns="idx", inplace=True)

        return catalog, assignments
