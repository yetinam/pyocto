import datetime
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .util import StationMapper, logger, set_cwd

try:
    from obspy.taup import TauPyModel
    from obspy.taup.taup_create import build_taup_model
except ImportError:  # Running without obspy, no 1D models available
    TauPyModel = None
    build_taup_model = None
    logger.warning(
        "Did not find obspy installation. 1D velocity models will not be supported."
    )

from .base import AbstractAssociator

FAKE_NETWORK_CODE = "XX"


class REALAssociator(AbstractAssociator):
    def __init__(
        self,
        *,
        search_range_horizontal_deg: float,
        search_range_vertical_km: float,
        grid_size_horizontal_deg: float,
        grid_size_vertical_km: float,
        event_separation_sec: float,
        vp: float,
        vs: float,
        p_picks: int,
        s_picks: int,
        total_picks: int,
        p_and_s_picks: int,
        max_residual_std: float,
        min_p_to_s_separation: float,
        nrt: float = 1.5,  # Not sure what this does.
        drt: float = 0.5,  # Not sure what this does
        nxd: float = 1.0,  # Distance criterion, 1.0 deactivates it
        tolerance_multiplier: float = 4.0,
        shallow_vp: float = np.nan,
        shallow_vs: float = np.nan,
        elevation_correction: bool = False,
        max_azimuthal_gap: float = 360,
        max_distance_pick_deg: float = 180,
        latref0: Optional[float] = None,
        lonref0: Optional[float] = None,
        velocity_model: Optional[str] = None,  # A model that can be read by taupy
        tt_range_horizontal_deg: Optional[float] = None,
        tt_grid_size_horizontal_deg: Optional[float] = None,
        tt_grid_size_vertical_km: Optional[float] = None,
        real_prefix: str = "",
        **kwargs,
    ) -> None:
        self.search_range_horizontal_deg = search_range_horizontal_deg
        self.search_range_vertical_km = search_range_vertical_km
        self.grid_size_horizontal_deg = grid_size_horizontal_deg
        self.grid_size_vertical_km = grid_size_vertical_km
        self.event_separation_sec = event_separation_sec
        self.vp = vp
        self.vs = vs
        self.p_picks = p_picks
        self.s_picks = s_picks
        self.total_picks = total_picks
        self.p_and_s_picks = p_and_s_picks
        self.max_residual_std = max_residual_std
        self.min_p_to_s_separation = min_p_to_s_separation
        self.nrt = nrt
        self.drt = drt
        self.nxd = nxd
        self.tolerance_multiplier = tolerance_multiplier
        self.shallow_vp = shallow_vp
        self.shallow_vs = shallow_vs
        self.elevation_correction = elevation_correction
        self.max_azimuthal_gap = max_azimuthal_gap
        self.max_distance_pick_deg = max_distance_pick_deg
        self.latref0 = latref0
        self.lonref0 = lonref0
        self.velocity_model = velocity_model
        self.tt_range_horizontal_deg = tt_range_horizontal_deg
        self.tt_range_vertical_km = self.search_range_vertical_km
        self.tt_grid_size_horizontal_deg = tt_grid_size_horizontal_deg
        self.tt_grid_size_vertical_km = tt_grid_size_vertical_km
        self._real_prefix = real_prefix

        self.latitude_center = None
        self._working_directory = None
        self._station_mapper = None

        self._verify_real()
        super().__init__(**kwargs)

    def _verify_real(self):
        p = subprocess.run(
            [self._real_prefix + "REAL"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if p.returncode != 255:
            raise ImportError(
                "Missing or incorrectly installed dependency REAL. "
                "Installation instructions at https://github.com/Dal-mzhang/REAL/ . "
                "Run REAL command in your shell for details."
            )

    def setup(self, stations: pd.DataFrame) -> None:
        self._station_mapper = StationMapper(stations)
        self._working_directory = Path(tempfile.TemporaryDirectory().name)

        self.latitude_center = (
            stations["latitude"].min() + stations["latitude"].max()
        ) / 2

        if self.velocity_model is None:
            pass
        else:
            self._calc_tt_tables()

    def _calc_tt_tables(self) -> bool:
        if TauPyModel is None:
            raise ImportError(
                "1D models require obspy. Please install obspy and rerun."
            )

        # Code inspired by https://github.com/Dal-mzhang/REAL/blob/master/demo_syn/tt_db/taup_tt.py
        if str(self.velocity_model).endswith(".nd"):
            build_taup_model(
                self.velocity_model, self._working_directory, verbose=False
            )
            self.velocity_model = self._working_directory / (
                self.velocity_model.split("/")[-1][:-3] + ".npz"
            )

        model = TauPyModel(self.velocity_model)

        table = ""
        for dep in np.arange(
            0,
            self.tt_range_vertical_km + self.tt_grid_size_vertical_km,
            self.tt_grid_size_vertical_km,
        ):
            for dist in np.arange(
                0,
                self.tt_range_horizontal_deg + self.tt_grid_size_horizontal_deg,
                self.tt_grid_size_horizontal_deg,
            ):
                arrivals = model.get_travel_times(
                    source_depth_in_km=dep,
                    distance_in_degree=dist,
                    phase_list=["P", "p", "S", "s"],
                )

                phases = {}
                for arr in arrivals:
                    if arr.name.lower() not in phases:
                        ray_param = arr.ray_param * 2 * np.pi / 360
                        slowness = -(ray_param / 111.19) / np.tan(
                            arr.takeoff_angle * np.pi / 180
                        )
                        phases[arr.name.lower()] = (
                            arr.time,
                            ray_param,
                            slowness,
                            arr.name,
                        )

                table += f"{dist} {dep}"
                for vp, vs in zip(phases["p"], phases["s"]):
                    table += f" {vp} {vs}"
                table += "\n"

        with open(self._working_directory / "tt_db.txt", "w") as f:
            f.write(table)

        return True

    def get_events(
        self, picks: pd.DataFrame, stations: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        # [REAL can process seismic picks recorded in one day (or a few days but only up to 31 days, e.g.,
        # 2016/10/01 – 2016/10/31 but not eligible for 2016/10/02 – 2016/11/01). All picks are relative to
        # ZERO of the day (e.g., 60.00 corresponds to 2016/10/14 00:01:00.00 and 86460 corresponds to
        # 2016/10/15 00:01:00.00)]
        if "probability" not in picks.columns:
            picks["probability"] = 1

        picks["segment"] = picks["time"].apply(lambda x: x.strftime("%Y/%m"))

        catalog = []
        assignments = []

        for i, (_, segment_picks) in enumerate(picks.groupby("segment")):
            t0 = (
                segment_picks["time"]
                .min()
                .replace(hour=0, minute=0, second=0, microsecond=0)
            )

            tmp_dir = self._working_directory
            tmp_dir.mkdir(exist_ok=True, parents=True)

            self.write_stations(tmp_dir, stations)
            self.write_picks(tmp_dir, segment_picks, t0)

            with set_cwd(
                tmp_dir
            ):  # some paths are hard coded, so we need to set them manually
                p = subprocess.run(
                    [self._real_prefix + "REAL"] + self.get_real_args(t0),
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                if p.returncode != 0:
                    raise ValueError(f"REAL exited with code {p}")

            seg_catalog, seg_assignments = self.parse_outputs(tmp_dir, t0)
            catalog.append(seg_catalog)
            assignments.append(seg_assignments)

        catalog = pd.concat(catalog)
        assignments = pd.concat(assignments)

        pick_lookup = picks.copy()
        pick_lookup["pick_idx"] = pick_lookup.index
        pick_lookup.drop(columns="event")
        pick_lookup["timestamp"] = pick_lookup["time"].apply(lambda x: x.timestamp())
        pick_lookup["timestamp"] = np.round(pick_lookup["timestamp"], 4)
        assignments["timestamp"] = assignments["time"].apply(lambda x: x.timestamp())

        assignments = pd.merge(
            assignments, pick_lookup, on=["station", "phase", "timestamp"]
        )
        assignments.drop(columns="timestamp", inplace=True)

        return catalog, assignments

    def get_real_args(self, t0: datetime.datetime) -> list[str]:
        # Variable names identical to REAL flag names
        D = t0.strftime("%Y/%m/%d") + f"/{self.latitude_center:.2f}"

        R = (
            f"{self.search_range_horizontal_deg:.2f}/{self.search_range_vertical_km:.1f}/"
            f"{self.grid_size_horizontal_deg:.3f}/{self.grid_size_vertical_km:.2f}/"
            f"{self.event_separation_sec:.1f}/{self.max_azimuthal_gap:.1f}/"
            f"{self.max_distance_pick_deg:.2f}"
        )

        if self.latref0 is not None and self.lonref0 is not None:
            R += f"/{self.latref0:.2f}/{self.lonref0:.2f}"

        V = f"{self.vp}/{self.vs}/{self.shallow_vp}/{self.shallow_vs}/{int(self.elevation_correction)}"

        S = (
            f"{self.p_picks}/{self.s_picks}/{self.total_picks}/{self.p_and_s_picks}/{self.max_residual_std}/"
            f"{self.min_p_to_s_separation}/{self.nrt}/{self.drt}/{self.nxd}/{self.tolerance_multiplier}/0"
        )

        flags = ["-D" + D, "-R" + R, "-V" + V, "-S" + S]
        args = ["stations.dat", "picks"]
        if self.velocity_model is not None:
            G = (
                f"{self.tt_range_horizontal_deg}/{self.tt_range_vertical_km}/"
                f"{self.tt_grid_size_horizontal_deg}/{self.tt_grid_size_vertical_km}"
            )
            flags.extend(["-G" + G])
            args.append("tt_db.txt")

        return flags + args

    def parse_outputs(
        self, tmp_dir: Path, t0: datetime.datetime
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        catalog = []
        assignments = []
        with open(tmp_dir / "phase_sel.txt", "r") as f:
            event_idx = -1
            for line in f:
                if not line.strip():
                    continue

                parts = line.strip().split()
                if len(parts) == 17:
                    # Event line
                    # num, year, mon, day, time (hh:mm:ss), origin time (relative to ZERO, sec), residual (sec), lat.,
                    # lon., dep., mag., mag var (uncertainty), number of P picks, number of S picks, total number of
                    # picks, number of stations with both P and S, station gap
                    event_idx += 1
                    catalog.append(
                        {
                            "idx": event_idx,
                            "time": t0 + datetime.timedelta(seconds=float(parts[5])),
                            "real_residual": float(parts[6]),
                            "latitude": float(parts[7]),
                            "longitude": float(parts[8]),
                            "depth": float(parts[9]),
                            "number_p_picks": int(parts[12]),
                            "number_s_picks": int(parts[13]),
                            "number_picks": int(parts[14]),
                            "number_p_and_s_picks": int(parts[15]),
                            "real_station_gap": float(parts[16]),
                        }
                    )

                    pass
                else:
                    # Pick line
                    # network, station, phase name, absolute travetime (relative to ZERO, sec), traveltime
                    # relative to event origin time (sec), phase amplitude in millimeter,
                    # individual phase residual (sec), weight, azimuth
                    assignments.append(
                        {
                            "event_idx": event_idx,
                            "station": self._station_mapper.translate_station(
                                parts[1], inv=True
                            ),
                            "phase": parts[2],
                            "time": t0 + datetime.timedelta(seconds=float(parts[3])),
                            "real_residual": float(parts[6]),
                            "real_weight": float(parts[7]),
                        }
                    )

        return pd.DataFrame(catalog), pd.DataFrame(assignments)

    def write_stations(self, tmp_dir: Path, stations: pd.DataFrame) -> None:
        with open(tmp_dir / "stations.dat", "w") as f:
            for _, station in stations.iterrows():
                # lon., lat., network, station, component, elevation (km)
                # Fake network and channel to provide fixed station id
                mapped_station_code = self._station_mapper.translate_station(
                    station["id"]
                )
                f.write(
                    f"{station['longitude']} {station['latitude']} {FAKE_NETWORK_CODE} "
                    f"{mapped_station_code} HHZ {station['elevation'] / 1e3}\n"
                )

    def write_picks(
        self, tmp_dir: Path, picks: pd.DataFrame, t0: datetime.datetime
    ) -> None:
        pick_dir = tmp_dir / "picks"
        pick_dir.mkdir(exist_ok=True)

        groups = defaultdict(list)
        t0_timestamp = t0.timestamp()

        for _, pick in picks.iterrows():
            mapped_station_code = self._station_mapper.translate_station(
                pick["station"]
            )
            group = (
                f"{FAKE_NETWORK_CODE}.{mapped_station_code}.{pick['phase'].upper()}.txt"
            )
            # arrivaltime (sec), stalta_ratio or phase_probability, amplitude_in_millimeter
            groups[group].append(
                f"{pick['time'].timestamp() - t0_timestamp} {pick['probability']} 0.0\n"
            )

        for group, group_picks in groups.items():
            with open(pick_dir / group, "w") as f:
                f.write("".join(group_picks))
