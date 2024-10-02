import datetime
import logging
import os
import struct
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional, Tuple, Union

import numpy as np
import pandas as pd
import pyproj
from pyproj import CRS, Transformer

from . import _core as backend

logger = logging.getLogger("pyocto")


class VelocityModel(ABC):
    """
    An abstract velocity model class.
    Only used to define the interface and provide type hints.
    """

    @abstractmethod
    def to_cpp(self, stations: list[backend.Station]) -> backend.VelocityModel:
        """
        Converts the Python representation of the object to a cpp representation.

        :param stations: A list of stations as cpp objects
        :return: The cpp representation of the object
        """
        pass


class VelocityModel0D(VelocityModel):
    """
    A homogeneous velocity model

    :param p_velocity: P wave velocity in km/s
    :param s_velocity: S wave velocity in km/s
    :param tolerance: Velocity model tolerance in s
    :param association_cutoff_distance: Only use stations up to this distance for space-partitioning association
    :param location_cutoff_distance: Only use stations up to this distance for location.
                                     Only such picks will be included in the assignments of the associator.
    """

    def __init__(
        self,
        p_velocity: float,
        s_velocity: float,
        tolerance: float,
        association_cutoff_distance: float = None,
        location_cutoff_distance: float = None,
    ):
        self.p_velocity = p_velocity
        self.s_velocity = s_velocity
        self.tolerance = tolerance
        self.association_cutoff_distance = association_cutoff_distance
        self.location_cutoff_distance = location_cutoff_distance

    def to_cpp(self, stations: list[backend.Station]) -> backend.VelocityModel:
        # Do the conversion to float explicitly to avoid crashes
        association_cutoff_distance = self.association_cutoff_distance
        if association_cutoff_distance is None:
            association_cutoff_distance = 1e9  # Infinity

        location_cutoff_distance = self.location_cutoff_distance
        if location_cutoff_distance is None:
            location_cutoff_distance = 1e9  # Infinity

        model = backend.VelocityModel0D(
            float(self.p_velocity),
            float(self.s_velocity),
            float(self.tolerance),
            float(association_cutoff_distance),
            float(location_cutoff_distance),
        )
        for station in stations:
            model.add_station(station)
        return model


class VelocityModel1D(VelocityModel):
    """
    A 1D layered velocity model. PyOcto uses a binary representation of the travel-time tables.
    To create this representation, please use :py:func:`create_model`. This step only needs
    to be executed once.

    :param path: Path to the travel-time table
    :param tolerance: Velocity model tolerance in s
    :param association_cutoff_distance: Only use stations up to this distance for space-partitioning association
    :param location_cutoff_distance: Only use stations up to this distance for location.
                                     Only such picks will be included in the assignments of the associator.
    :param surface_p_velocity: P wave velocity used for elevation correction in km/s
    :param surface_s_velocity: S wave velocity used for elevation correction in km/s

    .. warning::

        The VelocityModel1D does currently not allow search spaces above the surface, i.e., negative values of z.
    """

    def __init__(
        self,
        path: Union[str, Path],
        tolerance: float,
        association_cutoff_distance: float = None,
        location_cutoff_distance: float = None,
        surface_p_velocity: float = None,
        surface_s_velocity: float = None,
    ):
        self.path = path
        self.tolerance = tolerance
        self.association_cutoff_distance = association_cutoff_distance
        self.location_cutoff_distance = location_cutoff_distance
        self.surface_p_velocity = surface_p_velocity
        self.surface_s_velocity = surface_s_velocity

    def to_cpp(self, stations: list[backend.Station]) -> backend.VelocityModel:
        model = backend.VelocityModel1D(str(self.path))
        model.tolerance = self.tolerance
        if self.surface_p_velocity is not None:
            model.surface_p_velocity = self.surface_p_velocity
        if self.surface_s_velocity is not None:
            model.surface_s_velocity = self.surface_s_velocity
        if self.association_cutoff_distance is not None:
            model.association_cutoff_distance = self.association_cutoff_distance
        if self.location_cutoff_distance is not None:
            model.location_cutoff_distance = self.location_cutoff_distance
        for station in stations:
            model.add_station(station)
        return model

    @staticmethod
    def create_model(
        model: pd.DataFrame,
        delta: float,
        xdist: float,
        zdist: float,
        path: Union[str, Path],
    ) -> None:
        """
        Create a velocity model for PyOcto from a data frame

        :param model: DataFrame with columns depth, vp, vs
        :param delta: Grid spacing in kilometer
        :param xdist: Maximum distance in horizontal direction in km
        :param zdist: Maximum distance in vertical direction in km
        :param path: Output path
        """
        try:
            from pyrocko.modelling import eikonal
        except ImportError:
            raise ImportError(
                "Creating a velocity model requires pyrocko. "
                "Please install the package and try again."
            )

        nx = int(xdist / delta)
        nz = int(zdist / delta)

        p_speeds = np.ones((nx, nz))
        p_times = -np.ones_like(p_speeds)
        p_times[0, 0] = 0

        s_speeds = np.ones((nx, nz))
        s_times = -np.ones_like(s_speeds)
        s_times[0, 0] = 0

        for i in range(len(model) - 1):
            depth1 = model.iloc[i]["depth"]
            depth2 = model.iloc[i + 1]["depth"]
            vp1 = model.iloc[i]["vp"]
            vp2 = model.iloc[i + 1]["vp"]
            vs1 = model.iloc[i]["vs"]
            vs2 = model.iloc[i + 1]["vs"]

            depth1, depth2, vp1, vp2, vs1, vs2 = VelocityModel1D._adjust_depth_velocity(
                depth1, depth2, vp1, vp2, vs1, vs2, zdist
            )

            p_speeds[:, int(depth1 / delta) : int(depth2 / delta)] = np.linspace(
                vp1, vp2, int(depth2 / delta) - int(depth1 / delta)
            )
            s_speeds[:, int(depth1 / delta) : int(depth2 / delta)] = np.linspace(
                vs1, vs2, int(depth2 / delta) - int(depth1 / delta)
            )

        last = model.iloc[-1]
        p_speeds[:, int(last["depth"] / delta) :] = last["vp"]
        s_speeds[:, int(last["depth"] / delta) :] = last["vs"]

        eikonal.eikonal_solver_fmm_cartesian(p_speeds, p_times, delta)
        eikonal.eikonal_solver_fmm_cartesian(s_speeds, s_times, delta)

        with open(path, "wb") as f:
            f.write(struct.pack("iid", nx, nz, delta))
            f.write(p_times.tobytes())
            f.write(s_times.tobytes())

    @staticmethod
    def _adjust_depth_velocity(
        depth1: float,
        depth2: float,
        vp1: float,
        vp2: float,
        vs1: float,
        vs2: float,
        zdist: float,
    ) -> Tuple[float, float, float, float, float, float]:
        """
        Clip the interval [depth1, depth2] to lie fully within [0, zdist].
        Interpolates vp1 and vp2 linearly in between.

        :param depth1:
        :param depth2:
        :param vp1:
        :param vp2:
        :param vs1:
        :param vs2:
        :param zdist:
        :return:
        """
        slope_p = (vp2 - vp1) / (depth2 - depth1)
        depth_ref = depth1
        v_p_ref = vp1

        def func_vp(depth: float) -> float:
            return v_p_ref + (depth - depth_ref) * slope_p

        slope_s = (vs2 - vs1) / (depth2 - depth1)
        v_s_ref = vs1

        def func_vs(depth: float) -> float:
            return v_s_ref + (depth - depth_ref) * slope_s

        depth1 = min(max(depth1, 0), zdist)
        depth2 = min(max(depth2, 0), zdist)

        return (
            depth1,
            depth2,
            func_vp(depth1),
            func_vp(depth2),
            func_vs(depth1),
            func_vs(depth2),
        )


class OctoAssociator:
    """
    The OctoAssociator is the main class of PyOcto. An instance of this associator
    describes the configuration of the algorithm. To start the actual association,
    use the :py:func:`associate` function. You can also use one of the alternative
    interfaces :py:func:`associate_gamma`, :py:func:`associate_real`, or
    :py:func:`associate_seisbench`.

    In addition to the core functionality, this class offers convenience functions
    related to association. If you want to define your search area using latitude
    and longitude with an automatic coordinate projection, use :py:func:`from_area`.
    If you want to identify an appropriate projection for you stations, use
    :py:func:`get_crs`. For coordinate projections, :py:func:`transform_stations`
    and :py:func:`transform_events` are useful helper functions. To convert obspy
    inventory object to station data frames for PyOcto, use :py:func:`inventory_to_df`.

    As the PyOcto output locations are only preliminary, there is a helper function to
    convert the outputs to the NonLinLoc input format. Check out :py:func:`to_nonlinloc`.

    This documentation explains all parameters from a technical perspective.
    For a more user-centric view on how to set appropriate parameters,
    check out this :ref:`guide on parameter choices<parameters>`.

    :param xlim: Limit of the search space in kilometers in the x direction.
    :param ylim: Limit of the search space in kilometers in the y direction.
    :param zlim: Limit of the search space in kilometers in depth direction.
                 Negative values indicated above surface locations, positive values below surface.
    :param velocity_model: The velocity model.
    :param time_before: The overlap between consecutive time slices.
    :param min_node_size: Minimum node size for association. If a node becomes smaller, the event
                          creation process with localisation and pick refinement is triggered.
    :param min_node_size_location: Minimum node size, i.e., precision regarding discretisation,
                                   for the location algorithm. Usually smaller than `min_node_size`.
    :param pick_match_tolerance: Maximum difference between the predicted travel time and the observed time
                                 for associating a pick to an origin in the refinement step.
    :param min_interevent_time: Minimum time required between two events.
    :param exponential_edt: Exponentiate the individual term in the EDT loss.
                            This will make the loss surface more spiky, leading to better locations.
                            However, to accurately find minima, the `location_split_depth` and
                            `location_split_return` need to be increased, leading to higher computational
                            cost.
    :param edt_pick_std: Standard deviation for the EDT loss. Only relevant if exponential_edt is enabled.
    :param max_pick_overlap: Maximum number of picks shared between two events.
                             Note that overlaps are only possible at the intersection of different time blocks.
    :param n_picks: Minumum required number picks for an event.
    :param n_p_picks: Minumum required number P picks for an event.
    :param n_s_picks: Minumum required number S picks for an event.
    :param n_p_and_s_picks: Minumum required number of stations that have both P and S pick for an event.
    :param refinement_iterations: The number of localisation and pick matching iterations.
    :param time_slicing: The size of each time block.
    :param node_log_interval: If the value is larger than zero, each thread prints every time the number of nodes
                              explored so far is divisible by the value.
    :param queue_memory_protection_dfs_size: Maximum size of the priority queue per thread.
                                             If the queue is full, all further nodes are explored using
                                             depth-first search.
    :param location_split_depth: Search depth for location splits.
    :param location_split_return: Part of search depth for location splits that is not descended but only used to
                                  evenly sample the space. Always needs to be smaller than `location_split_depth`.
    :param min_pick_fraction: A distance based pick criterion.
                              If for a station, less than this fraction of closer stations have at least one pick,
                              the station picks are discarded.
    :param second_pass_overwrites: If not None, PyOcto will perform a second pass of the association procedure.
                                   In this second pass, only picks not associated in the first round will be used.
                                   The results of both passes are concatenated. The overwrites define all parameters
                                   that should be different in the second pass. A typical use case for this feature
                                   might be to associate events with many P picks but few S picks. This is, for example,
                                   a common scenario for large events picked with ML pickers. The waveform saturates
                                   and in effect the ML picker fails to identify the S picks. This case could be caught
                                   by setting a high ``n_picks`` and a low ``n_s_picks`` and a low ``n_p_and_s_picks``.
                                   At the same time, these settings might not be advisable for a larger processing
                                   because it will lead to many missed small events with lower numbers of P picks.
    :param n_threads: The number of threads to use.
                      By default, the number of threads will be set to the number of available cores.
    :param velocity_model_location: The velocity model for location.
                                    If not set, the same model as for association will be used.
    :param crs: The coordinate reference system. Required for all helper functions for coordinate transformation.
    """

    def __init__(
        self,
        xlim: tuple[float, float],
        ylim: tuple[float, float],
        zlim: tuple[float, float],
        velocity_model: VelocityModel,
        time_before: float,  # Should match travel time through the network
        min_node_size: float = 10.0,
        min_node_size_location: int = 1.5,
        pick_match_tolerance: float = 1.5,
        min_interevent_time: float = 3.0,
        exponential_edt: float = False,
        edt_pick_std: float = 1.0,
        max_pick_overlap: int = 4,
        n_picks: int = 10,
        n_p_picks: int = 3,
        n_s_picks: int = 3,
        n_p_and_s_picks: int = 3,
        refinement_iterations: int = 3,
        time_slicing: float = 1200.0,
        node_log_interval: int = 0,  # No logs
        queue_memory_protection_dfs_size: int = 500000,
        location_split_depth: int = 6,
        location_split_return: int = 4,
        min_pick_fraction: float = 0.25,
        second_pass_overwrites: Optional[dict[str, Any]] = None,
        n_threads: Optional[int] = None,  # default to number of available cores
        velocity_model_location: Optional[
            VelocityModel
        ] = None,  # defaults to velocity model
        crs: Optional[CRS] = None,
    ):
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim

        self.n_threads = n_threads
        self.node_log_interval = node_log_interval
        self.queue_memory_protection_dfs_size = queue_memory_protection_dfs_size
        self.min_pick_fraction = min_pick_fraction
        self.location_split_depth = location_split_depth
        self.location_split_return = location_split_return
        self.time_slicing = time_slicing
        self.refinement_iterations = refinement_iterations
        self.n_p_and_s_picks = n_p_and_s_picks
        self.n_s_picks = n_s_picks
        self.n_p_picks = n_p_picks
        self.n_picks = n_picks
        self.max_pick_overlap = max_pick_overlap
        self.min_interevent_time = min_interevent_time
        self.exponential_edt = exponential_edt
        self.edt_pick_std = edt_pick_std
        self.pick_match_tolerance = pick_match_tolerance
        self.min_node_size_location = min_node_size_location
        self.min_node_size = min_node_size
        self.time_before = time_before
        self.second_pass_overwrites = second_pass_overwrites

        self.velocity_model_association = velocity_model
        if velocity_model_location is None:
            self.velocity_model_location = velocity_model
        else:
            self.velocity_model_location = velocity_model_location

        self.crs = crs
        self._cached_pointers = (
            {}
        )  # References that need to be kept in memory to avoid automatic garbage collection

        self._os_check()
        self._pick_count_warnings()

    def _pick_count_warnings(self):
        if self.n_picks < self.n_p_picks + self.n_s_picks:
            logger.warning(
                f"The required number of picks per event ({self.n_picks}) is lower than the sum of the "
                f"required P ({self.n_p_picks}) and S ({self.n_s_picks}). The effective number of picks "
                f"required will be {self.n_p_picks + self.n_s_picks}."
            )

        if self.n_picks < 2 * self.n_p_and_s_picks:
            logger.warning(
                f"The required number of picks per event ({self.n_picks}) is lower than twice the number "
                f"of stations with both P and S pick ({self.n_p_and_s_picks}). The effective number of "
                f"picks required will be {2 * self.n_p_and_s_picks}."
            )

        if self.n_p_picks < self.n_p_and_s_picks:
            logger.warning(
                f"The required number of P picks per event ({self.n_p_picks}) is lower than the number "
                f"of stations with both P and S pick ({self.n_p_and_s_picks}). The effective number of "
                f"P picks required will be {self.n_p_and_s_picks}."
            )

        if self.n_s_picks < self.n_p_and_s_picks:
            logger.warning(
                f"The required number of S picks per event ({self.n_s_picks}) is lower than the number "
                f"of stations with both P and S pick ({self.n_p_and_s_picks}). The effective number of "
                f"S picks required will be {self.n_p_and_s_picks}."
            )

    def _os_check(self) -> None:
        if self.n_threads is None or self.n_threads != -1:
            if os.name == "nt":
                logger.warning(
                    "PyOcto does not support multi-threading on Windows. "
                    "Only one thread will be used for association. "
                    "Set n_threads=1 to suppress this warning."
                )

    @property
    def crs(self) -> Optional[pyproj.CRS]:
        """
        Get and set the local coordinate reference system if defined.

        :return: The local coordinate reference system.
        """
        return self._crs

    @crs.setter
    def crs(self, val):
        if val is None:
            self._crs = None
            self._transformer = None
            self._global_crs = None
        else:
            self._meter_model = False
            for axis_info in val.axis_info:
                assert axis_info.unit_code in [
                    "9001",
                    "9036",
                ], "CRS units need to be meter or kilometer."
                if axis_info.unit_code == "9001":
                    self._meter_model = True

            self._crs = val
            self._global_crs = CRS.from_epsg(4326)  # WSG-84
            self._transformer = Transformer.from_crs(self._global_crs, self._crs)

    def _build_config(
        self, stations: list[backend.Station], second_pass: bool = False
    ) -> backend.OctoTreeConfig:
        """
        Build the config as cpp object

        :param stations:
        :param second_pass: If true, uses overwrites for the second pass.
        :return:
        """
        config = backend.OctoTreeConfig()

        config.xlim = self.xlim
        config.ylim = self.ylim
        config.zlim = self.zlim

        config.node_log_interval = self.node_log_interval
        config.queue_memory_protection_dfs_size = self.queue_memory_protection_dfs_size
        config.time_slicing = self.time_slicing
        config.refinement_iterations = self.refinement_iterations
        config.n_p_and_s_picks = self.n_p_and_s_picks
        config.n_s_picks = self.n_s_picks
        config.n_p_picks = self.n_p_picks
        config.n_picks = self.n_picks
        config.max_pick_overlap = self.max_pick_overlap
        config.min_interevent_time = self.min_interevent_time
        config.pick_match_tolerance = self.pick_match_tolerance
        config.min_node_size_location = self.min_node_size_location
        config.min_node_size = self.min_node_size
        config.time_before = self.time_before
        config.location_split_depth = self.location_split_depth
        config.location_split_return = self.location_split_return
        config.min_pick_fraction = self.min_pick_fraction
        config.exponential_edt = self.exponential_edt
        config.edt_pick_std = self.edt_pick_std

        if self.n_threads is not None:
            config.n_threads = self.n_threads

        self._cached_pointers[
            "velocity_model_association"
        ] = self.velocity_model_association.to_cpp(stations)
        self._cached_pointers[
            "velocity_model_location"
        ] = self.velocity_model_location.to_cpp(stations)

        config.velocity_model_association = self._cached_pointers[
            "velocity_model_association"
        ]
        config.velocity_model_location = self._cached_pointers[
            "velocity_model_location"
        ]

        if second_pass:
            for key, value in self.second_pass_overwrites.items():
                if key not in dir(config):
                    raise KeyError(
                        f"Invalid overwrite '{key}'. All overwrites need to be valid config parameters."
                    )
                config.__setattr__(key, value)

        return config

    @staticmethod
    def _convert_picks(picks: pd.DataFrame) -> list[backend.Pick]:
        """
        Converts the picks from a data frame to a list of cpp objects

        :param picks:
        :return:
        """
        return [
            backend.Pick(i, row["time"], row["station"], row["phase"])
            for i, (_, row) in enumerate(picks.iterrows())
        ]

    @staticmethod
    def _convert_stations(stations: pd.DataFrame) -> list[backend.Station]:
        """
        Converts the stations from a data frame to a list of cpp objects

        :param stations:
        :return:
        """
        stations = stations.copy()

        for phase in "ps":
            # Add fake station terms
            if f"{phase}_residual" not in stations.columns:
                stations[f"{phase}_residual"] = 0.0
            # Replace missing station terms with 0
            stations[f"{phase}_residual"] = np.where(
                np.isnan(stations[f"{phase}_residual"]),
                0.0,
                stations[f"{phase}_residual"],
            )

        return [
            backend.Station(
                row["id"],
                row["x"],
                row["y"],
                row["z"],
                row["p_residual"],
                row["s_residual"],
            )
            for _, row in stations.iterrows()
        ]

    @staticmethod
    def _parse_events(
        events_cpp: list[backend.Event],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Translates the events list from cpp objects to data frames for events and assignments

        :param events_cpp:
        :return:
        """
        events = pd.DataFrame(
            [
                {
                    "idx": i,
                    "time": event.time,
                    "x": event.x,
                    "y": event.y,
                    "z": event.z,
                    "picks": len(event.picks),
                }
                for i, event in enumerate(events_cpp)
            ]
        )

        assignments = []
        for i, event in enumerate(events_cpp):
            assignments += [
                (i, pick.idx, residual)
                for pick, residual in zip(event.picks, event.residuals)
            ]
        assignments = pd.DataFrame(
            assignments, columns=["event_idx", "pick_idx", "residual"]
        )
        return events, assignments

    def associate(
        self, picks: pd.DataFrame, stations: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run the PyOcto associator.
        For details on the data formats, see :ref:`the data format description<data_formats>`.

        :param picks: The picks in PyOcto format
        :param stations: The stations in PyOcto format
        :return: Two dataframes.
                 The first contains the events.
                 The second one the assignment of picks to events.
        """
        if any(c not in stations.columns for c in "xyz"):
            raise ValueError(
                f"Colums 'x', 'y', and 'z' expected in station dataframe but not present. "
                f"The associator needs a local coordinate projection to run."
                f"Use the transform_stations function to add it or add it manually."
            )

        picks_cpp = self._convert_picks(picks)
        stations_cpp = self._convert_stations(stations)

        config = self._build_config(stations_cpp)

        associator = backend.OctoAssociator(config)

        events_cpp = associator.associate(picks_cpp)

        self._cached_pointers = {}  # Delete cached objects and allow garbage collection

        events, assignments = self._parse_events(events_cpp)
        picks_org = picks.copy(deep=True)
        picks_org["idx"] = np.arange(len(picks_org))

        if self.second_pass_overwrites is not None:
            # Perform second association pass with the unused picks
            picks_sub = picks_org[~picks_org["idx"].isin(assignments["pick_idx"])].copy(
                deep=True
            )  # Only unused picks
            picks_sub["idx2"] = np.arange(len(picks_sub))
            picks_cpp = self._convert_picks(picks_sub)

            config = self._build_config(stations_cpp, second_pass=True)

            associator = backend.OctoAssociator(config)

            events_cpp = associator.associate(picks_cpp)
            self._cached_pointers = (
                {}
            )  # Delete cached objects and allow garbage collection
            events2, assignments2 = self._parse_events(events_cpp)

            if len(events2) > 0:  # Skip this if there is anyhow no event
                # Update event_idx to be non-intersecting with original output
                if len(events) > 0:
                    event_idx_correction = events["idx"].max() + 1
                    events2["idx"] += event_idx_correction
                    assignments2["event_idx"] += event_idx_correction

                # Map pick idx to the orginal one (idx) instead of the subset one (idx2)
                assignments2 = pd.merge(
                    assignments2,
                    picks_sub,
                    left_on="pick_idx",
                    right_on="idx2",
                    validate="m:1",
                )
                assignments2.drop(columns="pick_idx", inplace=True)
                assignments2.rename(columns={"idx": "pick_idx"}, inplace=True)
                assignments2 = assignments2[
                    ["event_idx", "pick_idx", "residual"]
                ].copy()

                events = pd.concat([events, events2])
                assignments = pd.concat([assignments, assignments2])

                events.sort_values("time", inplace=True)
                events.reset_index(drop=True, inplace=True)
                assignments.reset_index(drop=True, inplace=True)

        assignments = pd.merge(
            assignments, picks_org, left_on="pick_idx", right_on="idx", validate="m:1"
        )
        assignments.drop(columns="idx", inplace=True)

        return events, assignments

    # Note: Missing type annotation to avoid SeisBench import
    def associate_seisbench(
        self, picks, stations: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run associator on a list of SeisBench picks.

        :param picks: A list of picks as output by SeisBench.
                      This list can be a SeisBench `PickList` instance.
        :param stations: Stations in PyOcto format
        :return: The outputs from `associate`
        """
        pick_df = []
        for p in picks:
            pick_df.append(
                {
                    "station": p.trace_id,
                    "time": p.peak_time.timestamp,
                    "probability": p.peak_value,
                    "phase": p.phase,
                }
            )

        pick_df = pd.DataFrame(pick_df)

        return self.associate(pick_df, stations)

    def associate_gamma(
        self, picks: pd.DataFrame, stations: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run associator on the GaMMA input format

        :param picks: Picks in GaMMA format
        :param stations: Stations in GaMMA format
        :return: The outputs from `associate`
        """
        picks = picks.copy()
        stations = stations.copy()

        stations.rename(
            columns={"x(km)": "x", "y(km)": "y", "z(km)": "z"},
            inplace=True,
        )
        stations["z"] = -stations["z"]

        picks.rename(
            columns={"timestamp": "time", "type": "phase", "id": "station"},
            inplace=True,
        )

        picks["time"] = picks["time"].apply(lambda x: x.timestamp())
        picks["phase"] = picks["phase"].apply(str.upper)

        return self.associate(picks, stations)

    def associate_real(
        self, pick_path: Union[str, Path], station_path: Union[str, Path]
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run associator on the REAL input format

        :param pick_path: Path of the directory containing the pick files for REAL
        :param station_path: Path of the station file for REAL
        :return: The outputs from `associate`
        """
        picks = self._real_parse_picks(Path(pick_path))
        stations = self._real_parse_stations(Path(station_path))

        stations = self.transform_stations(stations)

        return self.associate(picks, stations)

    @staticmethod
    def _real_parse_stations(station_path: Path) -> pd.DataFrame:
        """
        Parse stations from the format of the REAL associator and convert them to PyOcto format

        :param station_path:
        :return:
        """
        stations = pd.read_csv(
            station_path,
            names=[
                "longitude",
                "latitude",
                "network",
                "station",
                "channel",
                "elevation",
            ],
            sep=r"\s+",
        )
        stations[
            "elevation"
        ] *= 1e3  # To meters, will be transformed back in transform stations
        stations["id"] = stations.apply(
            lambda x: f"{x['network']}.{x['station']}.", axis=1
        )

        return stations

    @staticmethod
    def _real_parse_picks(pick_path: Path) -> pd.DataFrame:
        """
        Parse picks from the format of the REAL associator and convert them to PyOcto format

        :param pick_path:
        :return:
        """
        picks = []
        for pick_file in pick_path.iterdir():
            if not pick_file.name.endswith(".txt"):
                continue
            parts = pick_file.name.split(".")
            if len(parts) != 4:
                continue
            net, sta, phase, _ = parts

            pick_block = pd.read_csv(
                pick_file, names=["time", "confidence", "amplitude"], sep=r"\s+"
            )
            pick_block["phase"] = phase.upper()
            pick_block["station"] = f"{net}.{sta}."

            picks.append(pick_block)

        return pd.concat(picks)

    # Note: Missing type annotation to avoid obspy import
    def inventory_to_df(self, inventory) -> pd.DataFrame:
        """
        Convert an obspy inventory to a dataframe. Applies the coordinate projection if set.

        :param inventory: An obspy inventory object
        :return: A data frame with the station that can be input to :py:func:`associate`
        """
        station_df = []
        for network in inventory:
            for station in network:
                locs = [channel.location_code for channel in station] + [""]

                for loc in set(locs):
                    station_df.append(
                        {
                            "id": f"{network.code}.{station.code}.{loc}",
                            "longitude": station.longitude,
                            "latitude": station.latitude,
                            "elevation": station.elevation,
                        }
                    )

        return self.transform_stations(pd.DataFrame(station_df))

    @staticmethod
    def get_crs(stations: pd.DataFrame, warning_limit_deg: float = 15.0) -> CRS:
        """
        Get a transverse Mercator projection coordinate reference system centered in the middle
        of the station distribution.

        :param stations: A data frame with coordinates in latitude and longitude
        :param warning_limit_deg: If the along-axis distance between too stations is higher than this value,
                                  a warning is printed.
        :return: A coordinate reference system
        """
        max_lat = stations["latitude"].max()
        min_lat = stations["latitude"].min()

        min_lon = stations["longitude"].min()
        max_lon = stations["longitude"].max()

        if (max_lon - min_lon) > 180:
            logger.warning("Projection spanning longitude discontinuity")
            min_lon = stations[stations["longitude"] > 0]["longitude"].min()
            max_lon = stations[stations["longitude"] < 0]["longitude"].max() + 360

        lat0 = (min_lat + max_lat) / 2
        lon0 = (min_lon + max_lon) / 2

        if lon0 > 180:
            lon0 -= 360

        if max(max_lon - min_lon, max_lat - min_lat) > warning_limit_deg:
            logger.warning(
                f"Stations span more than {warning_limit_deg} degrees, "
                f"the coordinate projection might be inaccurate."
            )

        return CRS(f"+proj=tmerc +lat_0={lat0} +lon_0={lon0} +units=km")

    @classmethod
    def from_area(
        cls,
        lat: tuple[float, float],
        lon: tuple[float, float],
        zlim: tuple[float, float],
        velocity_model: VelocityModel,
        time_before: float,
        **kwargs,
    ):
        """
        Create an associator instance based on a bounding box in latitude, longitude and depth.

        :param lat: Minimum and maximum latitude of study area in degrees
        :param lon: Minimum and maximum longitude of study area in degrees
        :param zlim: Minimum and maximum depth of study area in km
        :param velocity_model: see the class constructor
        :param time_before: see the class constructor
        :param kwargs: passed to class constructor
        :return: an instance of OctoAssociator
        """
        fake_stations = pd.DataFrame(
            {"latitude": [lat[0], lat[1]], "longitude": [lon[0], lon[1]]}
        )
        crs = cls.get_crs(fake_stations)
        global_crs = CRS.from_epsg(4326)  # WSG-84
        transformer = Transformer.from_crs(global_crs, crs)

        left, bottom, right, top = transformer.transform_bounds(
            lat[0], lon[0], lat[1], lon[1]
        )

        xlim = (left, right)
        ylim = (bottom, top)

        return cls(xlim, ylim, zlim, velocity_model, time_before, crs=crs, **kwargs)

    def transform_stations(self, stations: pd.DataFrame) -> pd.DataFrame:
        """
        Project stations from cartesian coordinates into a local coordinate system.
        Requires the `crs` attribute to be set.
        Note that the original data frame is modified in-place.

        :param events: A data frame with the stations containing the latitude, longitude
                       and elevation columns. Elevation needs to be provided in meters.
                       Note that the transform flips the sign convention. `elevation` is in
                       meters above zero, `z` is in kilometers below zero.
        :return: A dataframe with additional x, y and z.
        """
        if self.crs is None:
            raise ValueError(
                "CRS is not set. Use 'associator.crs = associator.get_crs(stations)' to "
                "automatically set a CRS. Alternatively, set a CRS manually."
            )

        factor = 1
        if self._meter_model:
            factor = 1000

        stations["x"] = stations.apply(
            lambda x: self._transformer.transform(x["latitude"], x["longitude"])[0]
            / factor,
            axis=1,
        )
        stations["y"] = stations.apply(
            lambda x: self._transformer.transform(x["latitude"], x["longitude"])[1]
            / factor,
            axis=1,
        )
        if "elevation" in stations.columns:
            stations["z"] = -stations["elevation"] / 1e3

        return stations

    def transform_events(self, events: pd.DataFrame) -> pd.DataFrame:
        """
        Project event coordinates from local coordinate system to global coordinate system.
        Requires the `crs` attribute to be set.
        Note that the original data frame is modified in-place.

        :param events: A data frame with the events as output by :py:func:`associate`
        :return: A dataframe with additional latitude, longitude and depth columns.
        """
        if self.crs is None:
            raise ValueError(
                "CRS is not set. Most likely, the associator was not used to transform the "
                "station coordinates. Please manually transform the event coordinates."
            )

        if len(events) == 0:
            return events

        factor = 1
        if self._meter_model:
            factor = 1000

        events["latitude"] = events.apply(
            lambda x: self._transformer.transform(
                x["x"] * factor, x["y"] * factor, direction="INVERSE"
            )[0],
            axis=1,
        )
        events["longitude"] = events.apply(
            lambda x: self._transformer.transform(
                x["x"] * factor, x["y"] * factor, direction="INVERSE"
            )[1],
            axis=1,
        )
        events["depth"] = events["z"]

        return events

    @staticmethod
    def to_nonlinloc(
        assignments: pd.DataFrame, path: Union[str, Path], pick_std: float = 0.05
    ) -> None:
        """
        Write the outputs to the .obs format that can be parsed by NonLinLoc

        :param assignments: The assignments as output by PyOcto
        :param path: Output path for the observations
        :param pick_std: Gaussian uncertainty of the picks in seconds.
                         Currently, does not support individual uncertainties per pick.
        """
        with open(path, "w") as f:
            for event_idx, event_catalog in assignments.groupby("event_idx"):
                f.write(f"PUBLIC_ID E{event_idx:08d}\n")
                for _, pick in event_catalog.iterrows():
                    # Example: GRX    ?    ?    ? P      U 19940217 2216   44.9200 GAU  2.00e-02 -1.00e+00 -1.00e+00 -1.00e+00
                    if isinstance(pick["time"], datetime.datetime):
                        time = pick["time"]
                    else:
                        time = datetime.datetime.fromtimestamp(pick["time"])
                    phase = pick["phase"].upper()
                    station = pick["station"]
                    daystr = time.strftime("%Y%m%d")
                    hourmin = time.strftime("%H%M")
                    second = time.strftime("%S.%f")[:-2]

                    f.write(
                        f"{station:6s} ?    ?    ? {phase}      ? {daystr} {hourmin}   {second} "
                        f"GAU  {pick_std:.2e} -1.00e+00 -1.00e+00 -1.00e+00\n"
                    )
                f.write("\n")
