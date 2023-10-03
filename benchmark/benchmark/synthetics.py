import datetime
import logging
from pathlib import Path
from typing import Tuple

import hydra
import numpy as np
import pandas as pd
from obspy.geodetics import gps2dist_azimuth
from pyrocko.modelling import eikonal
from tqdm import tqdm

logger = logging.getLogger("benchmark")


@hydra.main(
    version_base=None, config_path="../synth_configs", config_name="synthetics.yaml"
)
def main(cfg):
    if "log_level" in cfg:
        logger.setLevel(cfg["log_level"])

    logger.debug("Instantiating config")
    cfg = hydra.utils.instantiate(cfg)

    logger.debug("Creating TT helper")
    helper = TTHelper(**cfg.data.tt_conf)

    logger.debug("Generating events")
    catalog, picks = generate_events(
        **cfg.base,
        stations=cfg.data.stations,
        candidates=cfg.data.candidates,
        helper=helper,
    )

    name = build_name(cfg)

    path = Path(cfg.output_base) / name
    logger.debug(f"Writing output to: {path}")
    path.mkdir(exist_ok=True, parents=True)
    catalog.to_parquet(path / "catalog", index=False)
    picks.to_parquet(path / "picks", index=False)
    cfg.data.stations.to_parquet(path / "stations", index=False)


def build_name(cfg) -> str:
    return f"{cfg.data.name}_{cfg.base.n_events}_{cfg.base.noise_factor:.1f}"


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def detection_probability(x, r, w, p):
    return (1 - sigmoid((x - r) / w)) * (1 - p)


def correlated_bernoulli(rho, p, n=1):
    x = np.random.rand(n) < p
    u = np.random.rand(n)
    y = np.where(x, u < rho * (1 - p) + p, u < p - rho * p)
    return x, y


class TTHelper:
    def __init__(
        self,
        model: pd.DataFrame,
        delta: float = 0.2,
        xdist: int = 1500,
        zdist: int = 300,
        tt_var_rel: float = 0.01,
        tt_var_abs: float = 0.4,
    ):
        self.model = model
        self.delta = delta
        self.xdist = xdist
        self.zdist = zdist
        self.tt_var_rel = tt_var_rel
        self.tt_var_abs = tt_var_abs

        self.p_times, self.s_times = self.get_tt_grids()

    def get_tt_grids(self) -> Tuple[np.ndarray, np.ndarray]:
        delta = self.delta
        model = self.model

        nx = int(self.xdist / delta)
        nz = int(self.zdist / delta)

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

        return p_times, s_times

    def get_tt(self, dist: float, depth: float, elev: float, phase: str):
        if phase == "P":
            t = self.p_times[int(dist / self.delta), int(depth / self.delta)]
            t += elev / 6.1
        elif phase == "S":
            t = self.s_times[int(dist / self.delta), int(depth / self.delta)]
            t += elev / 3.6
        else:
            raise ValueError(f"Unknown phase: {phase}")

        std = max(t * self.tt_var_rel, self.tt_var_abs)
        t += np.random.randn() * std
        return t


def generate_event(
    event_id: int,
    event_candidates: pd.DataFrame,
    stations: pd.DataFrame,
    tt_helper: TTHelper,
    duration: float,
    A: float = 120,
    B: float = 80,
):
    idx = np.random.randint(len(event_candidates))
    event = event_candidates.iloc[idx].to_dict()
    event["time"] = np.random.rand() * duration  # Seconds in one day
    event["MA"] = np.random.choice(event_candidates["MA"])
    event["r"] = A * event["MA"] + B
    event["event"] = event_id

    picks = []

    n_p_and_s = 0

    for _, station in stations.iterrows():
        dist = (
            gps2dist_azimuth(
                event["LAT"], event["LON"], station["latitude"], station["longitude"]
            )[0]
            / 1e3
        )
        hypo = np.sqrt(dist**2 + (event["DEPTH"] + station["elevation"] / 1e3) ** 2)

        p = detection_probability(
            hypo, r=event["r"], w=max(0.1 * event["r"], 30), p=0.2
        )
        has_p, has_s = correlated_bernoulli(0.5, p, n=1)

        if has_p:
            picks.append(
                (
                    station["id"],
                    "P",
                    event["time"]
                    + tt_helper.get_tt(
                        dist, event["DEPTH"], station["elevation"] / 1e3, "P"
                    ),
                )
            )
        if has_s:
            picks.append(
                (
                    station["id"],
                    "S",
                    event["time"]
                    + tt_helper.get_tt(
                        dist, event["DEPTH"], station["elevation"] / 1e3, "S"
                    ),
                )
            )

        if has_p and has_s:
            n_p_and_s += 1

    if n_p_and_s < 4 or len(picks) < 10:
        # Event is not detectable
        return generate_event(
            event_id, event_candidates, stations, tt_helper, duration, A, B
        )

    picks = pd.DataFrame(picks, columns=["station", "phase", "time"])
    picks["event"] = event_id

    return event, picks


def generate_false_picks(n: int, stations: pd.DataFrame, duration: float):
    return pd.DataFrame(
        {
            "station": np.random.choice(stations["id"], n),
            "phase": np.random.choice(["P", "S"], n),
            "time": np.random.rand(n) * duration,
            "event": -1,
        }
    )


def generate_events(
    n_events: int,
    candidates: pd.DataFrame,
    stations: pd.DataFrame,
    helper: TTHelper,
    noise_factor: float,
    duration: float = 24 * 60 * 60,
    **kwargs,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    events = []
    picks = []

    for i in tqdm(range(n_events)):
        event, ev_picks = generate_event(
            i, candidates, stations, helper, duration=duration, **kwargs
        )
        events.append(event)
        picks.append(ev_picks)

    events = pd.DataFrame(events)
    picks = pd.concat(picks)

    picks = pd.concat(
        [
            picks,
            generate_false_picks(int(noise_factor * len(picks)), stations, duration),
        ]
    )

    events.rename(
        columns={"LAT": "latitude", "LON": "longitude", "DEPTH": "depth", "MA": "mag"},
        inplace=True,
    )
    picks.sort_values("time", inplace=True)
    picks["time"] = picks["time"].apply(datetime.datetime.fromtimestamp)
    events["time"] = events["time"].apply(datetime.datetime.fromtimestamp)

    return events, picks


def random_station(
    n: int, xlim: tuple[float, float], ylim: tuple[float, float]
) -> pd.DataFrame:
    stations = {
        "id": [f"S{i:04d}" for i in range(n)],
        "latitude": np.random.rand(n) * (xlim[1] - xlim[0]) + xlim[0],
        "longitude": np.random.rand(n) * (ylim[1] - ylim[0]) + ylim[0],
        "elevation": np.zeros(n),
    }

    return pd.DataFrame(stations)


def grid_stations(
    nx: int, ny: int, xlim: tuple[float, float], ylim: tuple[float, float]
) -> pd.DataFrame:
    x, y = np.meshgrid(np.linspace(*xlim, nx), np.linspace(*ylim, ny))
    n = nx * ny
    stations = {
        "id": [f"S{i:04d}" for i in range(n)],
        "latitude": x.reshape(-1),
        "longitude": y.reshape(-1),
        "elevation": np.zeros(n),
    }

    return pd.DataFrame(stations)


def random_events(
    n: int,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    zlim: tuple[float, float],
    a: float = -0.5,
    b: float = 1,
) -> pd.DataFrame:
    x = np.linspace(a, 9, 10000)

    N = 10 ** (-b * x)  # inverse CDF
    nn_pdf = N[1:] - N[:-1]  # non-normed pdf
    nn_cdf = np.cumsum(nn_pdf)  # non-normed cdf
    cdf = nn_cdf / nn_cdf[-1]
    cdf = np.concatenate([np.zeros(1), cdf])

    def get_magnitude():
        return np.max(x[np.random.rand() > cdf])

    events = {
        "LAT": np.random.rand(n) * (xlim[1] - xlim[0]) + xlim[0],
        "LON": np.random.rand(n) * (ylim[1] - ylim[0]) + ylim[0],
        "DEPTH": np.random.rand(n) * (zlim[1] - zlim[0]) + zlim[0],
        "MA": np.array([get_magnitude() for _ in range(n)]),
    }

    return pd.DataFrame(events)


if __name__ == "__main__":
    main()
