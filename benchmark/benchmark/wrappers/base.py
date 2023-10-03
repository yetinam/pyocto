from abc import ABC, abstractmethod
from typing import Optional

import pandas as pd
from pyproj import CRS, Transformer


class AbstractAssociator(ABC):
    def __init__(self, crs: Optional[CRS] = None, global_crs: Optional[CRS] = None):
        if global_crs is None:
            self._global_crs = CRS.from_epsg(4326)
        else:
            self._global_crs = global_crs

        self.crs = crs

    def setup(self, stations: pd.DataFrame) -> None:
        pass

    @abstractmethod
    def get_events(
        self, picks: pd.DataFrame, stations: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        pass

    @property
    def crs(self):
        return self._crs

    @crs.setter
    def crs(self, value):
        self._crs = value
        if value is None:
            self.transformer = None
            self.inv_transformer = None
        else:
            self.transformer = Transformer.from_crs(self._global_crs, self._crs)
            self.inv_transformer = Transformer.from_crs(self._crs, self._global_crs)

    def transform_station_coords(self, stations: pd.DataFrame) -> pd.DataFrame:
        stations["x(km)"] = stations.apply(
            lambda x: self.transformer.transform(x["latitude"], x["longitude"])[0]
            / 1e3,
            axis=1,
        )
        stations["y(km)"] = stations.apply(
            lambda x: self.transformer.transform(x["latitude"], x["longitude"])[1]
            / 1e3,
            axis=1,
        )
        stations["z(km)"] = stations["elevation"] / 1e3

        return stations

    def transform_event_coords(self, catalog: pd.DataFrame) -> pd.DataFrame:
        if len(catalog) == 0:
            return catalog

        catalog["latitude"] = catalog.apply(
            lambda x: self.inv_transformer.transform(
                x["x(km)"] * 1e3, x["y(km)"] * 1e3
            )[0],
            axis=1,
        )
        catalog["longitude"] = catalog.apply(
            lambda x: self.inv_transformer.transform(
                x["x(km)"] * 1e3, x["y(km)"] * 1e3
            )[1],
            axis=1,
        )
        catalog["depth"] = catalog["z(km)"]
        catalog.drop(columns=["x(km)", "y(km)", "z(km)"], inplace=True)

        return catalog
