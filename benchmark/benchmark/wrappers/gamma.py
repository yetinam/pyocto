from typing import Any

import pandas as pd
from gamma.utils import association as gamma_association

from .base import AbstractAssociator


class GammaAssociator(AbstractAssociator):
    def __init__(self, config: dict[str, Any], **kwargs) -> None:
        super().__init__(**kwargs)

        self.config = self._adjust_config(config)

    @staticmethod
    def _adjust_config(config: dict) -> dict:
        if "bfgs_bounds" not in config:
            config["bfgs_bounds"] = [
                (config["x(km)"][0] - 1, config["x(km)"][1] + 1),
                (config["y(km)"][0] - 1, config["y(km)"][1] + 1),
                (config["z(km)"][0] - 1, config["z(km)"][1] + 1),
                (None, None),
            ]

        return config

    def get_events(
        self, picks: pd.DataFrame, stations: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        if len(picks) == 0:
            return pd.DataFrame(), pd.DataFrame()

        picks_org = picks.copy(deep=True)
        picks.rename(
            columns={
                "time": "timestamp",
                "phase": "type",
                "probability": "prob",
                "station": "id",
            },
            inplace=True,
        )

        if "prob" not in picks.columns:
            picks["prob"] = 1

        if "amp" not in picks.columns:
            picks["amp"] = 0

        picks["type"] = picks["type"].apply(lambda x: x.lower())

        if any(f"{c}(km)" not in stations for c in "xyz"):
            stations = self.transform_station_coords(stations)

        # noinspection PyCallingNonCallable
        catalog, assignments = gamma_association(
            picks,
            stations,
            self.config,
            event_idx0=0,
            method=self.config["method"],
        )

        catalog = pd.DataFrame(catalog)

        catalog.rename(
            columns={
                "event_index": "idx",
                "sigma_time": "gamma_sigma_time",
                "sigma_amp": "gamma_sigma_amp",
                "cov_time_amp": "gamma_cov_time_amp",
            },
            inplace=True,
        )
        picks_org["idx"] = picks_org.index

        assignments = pd.DataFrame(
            assignments, columns=["pick_idx", "event_idx", "probability_gamma"]
        )
        merged = pd.merge(assignments, picks_org, left_on="pick_idx", right_on="idx")
        merged.drop(columns="idx", inplace=True)  # Remove redundant column

        if len(catalog) > 0:
            catalog = self.transform_event_coords(catalog)
            catalog.sort_values("time", inplace=True)

        return catalog, merged
