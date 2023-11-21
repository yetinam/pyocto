import argparse
from collections import defaultdict
from pathlib import Path
from queue import Queue

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from cmcrameri import cm

textwidth = 7  # Arbitrary value, but this helps to orient everything to the same size
colwidth = 3.4

paper_figure_path = Path("/home/munchmej/publications/pyocto-paper/graphics")
paper_table_path = Path("/home/munchmej/publications/pyocto-paper/tables")
base_config = {}

name_lookup = {
    "gamma": "GaMMA",
    "real": "REAL",
    "real_1d": "REAL1D",
    "pyocto": "PyOcto",
    "pyocto_1d": "PyOcto1D",
    "pyocto_cutoff": "PyOcto",
    "pyocto_1d_cutoff": "PyOcto1D",
}
iquique_suffix = "_0.05"

if not paper_figure_path.is_dir():
    print("Did not find figure path, logging locally.")
    paper_figure_path = Path("paper/figures")
    paper_figure_path.mkdir(exist_ok=True, parents=True)

if not paper_table_path.is_dir():
    print("Did not find table path, logging locally.")
    paper_table_path = Path("paper/tables")
    paper_table_path.mkdir(exist_ok=True, parents=True)


def decimal_formatter(i: int) -> str:
    s = str(i)
    o = []
    for j, c in enumerate(reversed(s)):
        if j > 0 and j % 3 == 0:
            o.append(",")
        o.append(c)
    return "".join(reversed(o))


class FunctionManager:
    def __init__(self):
        self.function_dict = {}
        self.markers = defaultdict(list)

    def register(self, function):
        self.function_dict[function.__name__] = function
        return function

    def register_with_markers(self, markers):
        def register(func):
            for marker in markers:
                self.markers[marker] += [func]
            return self.register(func)

        return register

    def run(self, function):
        if function in self.function_dict:
            self.function_dict[function]()
        else:
            raise ValueError(f"Function {function} unknown")

    def run_marker(self, marker):
        if marker == "all":
            for func in self.function_dict.values():
                print(f"Running {func.__name__}")
                func()
        else:
            for func in self.markers[marker]:
                print(f"Running {func.__name__}")
                func()


manager = FunctionManager()


def init_style():
    sns.set(font_scale=0.6)
    sns.set_style("ticks")
    sns.set_palette("colorblind")  # Choose colorblind friendly color palette
    base_config["cmap"] = cm.batlow
    mpl.rcParams["figure.dpi"] = 300.0
    mpl.rcParams["figure.facecolor"] = "#00000000"  # Transparent background
    mpl.rcParams["lines.linewidth"] = 1.0  # Thinner lines


def get_benchmark_results() -> pd.DataFrame:
    results = pd.read_csv("pred/results.csv")

    def truncate(x):
        if isinstance(x, float):
            return ""
        else:
            return x.split(".")[-1]

    def add_1D_tag(row):
        tag = ""
        if (
            "associator.config.velocity_model._target_" in row
            and row["associator.config.velocity_model._target_"]
            == "pyocto.VelocityModel1D"
        ):
            tag = "1D"
        if "associator.tt_grid_size_horizontal_deg" in row and not np.isnan(
            row["associator.tt_grid_size_horizontal_deg"]
        ):
            tag = "1D"
        return row["associator"] + tag

    def add_exp_tag(row):
        if row["comment"] == "pyocto_cutoff":
            return row["associator"] + "V1"
        return row["associator"]

    if "comment" not in results.columns:
        results["comment"] = ""

    results["associator"] = results["associator._target_"].apply(truncate)
    results["data.path"] = results["data.path"].apply(lambda x: x.split("/")[-1])
    results["exp"] = results["data.path"].apply(lambda x: x.split("_")[0])
    results["events"] = results["data.path"].apply(lambda x: int(x.split("_")[1]))
    results["noise"] = results["data.path"].apply(lambda x: float(x.split("_")[2]))

    results["associator"] = results.apply(add_1D_tag, axis=1)
    # results["associator"] = results.apply(add_exp_tag, axis=1)
    results["associator"] = results["associator"].apply(
        lambda x: x.replace("Associator", "")
    )

    mask = (~results["associator"].isin(["PyOcto", "PyOcto1D"])) | (
        results["comment"] == "pyocto_cutoff"
    )
    results = results[mask].copy()

    results.sort_values(
        ["exp", "events", "noise", "associator", "comment"], inplace=True
    )

    return results[
        [
            "exp",
            "events",
            "noise",
            "associator",
            "comment",
            "precision",
            "recall",
            "f1",
            "missing_picks_per_event",
            "additional_picks_per_event",
            "runtime",
            "pred_path",
        ]
    ]


@manager.register_with_markers(["plot", "paper"])
def octotree_schematic_empty():
    octotree_schematic(empty=True)


@manager.register_with_markers(["plot", "paper"])
def octotree_schematic(empty=False):
    np.random.seed(5005)
    fig = plt.figure(figsize=(textwidth, 0.4 * textwidth))
    ax = fig.add_subplot(111)

    p_vel = 7
    s_vel = 4
    global_xlim = (0, 40)
    global_ylim = (-5, 75)
    stations = np.arange(0, 71, 10)

    ax.plot(0 * stations, stations, "k^", clip_on=False)
    ax.set_xlim(*global_xlim)
    ax.set_ylim(*global_ylim)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Distance [km]")

    p_picks = []
    s_picks = []

    events = [(4, 23), (24, 56)]
    for t0, x0 in events:
        t_p = t0 + np.abs(stations - x0) / p_vel
        t_s = t0 + np.abs(stations - x0) / s_vel

        p_picks += [(t, x) for t, x in zip(t_p, stations)]
        s_picks += [(t, x) for t, x in zip(t_s, stations)]

        if not empty:
            lw = 0.5
            ax.plot([t0, t0 + 100], [x0, x0 + 100 * p_vel], "k-", lw=lw)
            ax.plot([t0, t0 + 100], [x0, x0 - 100 * p_vel], "k-", lw=lw)
            ax.plot([t0, t0 + 100], [x0, x0 + 100 * s_vel], "k--", lw=lw)
            ax.plot([t0, t0 + 100], [x0, x0 - 100 * s_vel], "k--", lw=lw)
            ax.plot(t0, x0, "r*")

    p_noise = np.random.rand(17) * 40
    s_noise = np.random.rand(17) * 40
    p_picks += [
        (t, x) for t, x in zip(p_noise, np.random.choice(stations, len(p_noise), True))
    ]
    s_picks += [
        (t, x) for t, x in zip(s_noise, np.random.choice(stations, len(s_noise), True))
    ]

    p_picks = np.array(p_picks)
    s_picks = np.array(s_picks)

    if not empty:
        blocks = Queue()
        blocks.put((global_xlim, global_ylim))

        def min_max_dist(p, p0, p1):
            return (
                np.where(
                    (p0 < p) & (p < p1), 0, np.minimum(np.abs(p - p0), np.abs(p - p1))
                ),
                np.maximum(np.abs(p - p0), np.abs(p - p1)),
            )

        def matching_picks(xlim, ylim, picks, vel):
            ymin, ymax = min_max_dist(picks[:, 1], *ylim)
            tmin = xlim[0] + ymin / vel
            tmax = xlim[1] + ymax / vel

            return (tmin <= picks[:, 0]) & (picks[:, 0] < tmax)

        cmap = plt.get_cmap("Blues")

        def get_color(xlim, ylim, n_picks):
            area = (xlim[1] - xlim[0]) * (ylim[1] - ylim[0])
            return cmap(0.3 * n_picks / area)

        def split2(xlim, ylim):
            if xlim[1] - xlim[0] > (ylim[1] - ylim[0]) / s_vel:
                m = (xlim[0] + xlim[1]) / 2
                return ((xlim[0], m), ylim), ((m, xlim[1]), ylim)
            else:
                m = (ylim[0] + ylim[1]) / 2
                return (xlim, (ylim[0], m)), (xlim, (m, ylim[1]))

        while not blocks.empty():
            xlim, ylim = blocks.get()
            if xlim == (20.0, 25.0) and ylim == (55.0, 65.0):
                pass
            n_p_picks = np.sum(matching_picks(xlim, ylim, p_picks, p_vel))
            n_s_picks = np.sum(matching_picks(xlim, ylim, s_picks, s_vel))

            if n_p_picks < 6 or n_s_picks < 6:
                continue
            # print(xlim, ylim, n_p_picks, n_s_picks)

            patch = patches.Rectangle(
                (xlim[0], ylim[0]),
                xlim[1] - xlim[0],
                ylim[1] - ylim[0],
                facecolor=get_color(xlim, ylim, n_p_picks + n_s_picks),
                linewidth=0.5,
                edgecolor="k",
            )
            ax.add_patch(patch)

            if max(xlim[1] - xlim[0], (ylim[1] - ylim[0]) / s_vel) > 2:
                block1, block2 = split2(xlim, ylim)
                blocks.put(block1)
                blocks.put(block2)

    ax.plot(p_picks[:, 0], p_picks[:, 1], "+", c="C1", label="P")
    ax.plot(s_picks[:, 0], s_picks[:, 1], "+", c="C2", label="S")
    ax.legend()

    if empty:
        fig.savefig(
            paper_figure_path / "octotree_schematic_empty.png", bbox_inches="tight"
        )
    else:
        fig.savefig(paper_figure_path / "octotree_schematic.png", bbox_inches="tight")


@manager.register_with_markers(["table", "paper"])
def full_benchmark_table():
    results = get_benchmark_results()
    results["exp"][results["exp"] == "chile"] = "subduction"
    results.sort_values(
        ["exp", "events", "noise", "associator", "comment"], inplace=True
    )
    results = results[
        [
            "exp",
            "events",
            "noise",
            "associator",
            "precision",
            "recall",
            "f1",
            "missing_picks_per_event",
            "additional_picks_per_event",
            "runtime",
        ]
    ]
    results.rename(
        columns={
            "exp": "Scenario",
            "events": "Events",
            "noise": "Noise",
            "associator": "Associator",
            "precision": "Precision",
            "recall": "Recall",
            "f1": "F1",
            "missing_picks_per_event": "Missing ppe",
            "additional_picks_per_event": "Additional ppe",
            "runtime": "Run time [s]",
        },
        inplace=True,
    )
    results.to_latex(
        paper_table_path / "full_benchmark.tex",
        index=False,
        float_format="%.2f",
        longtable=True,
        label="tab:full_benchmark",
        caption="Full results of synthetic benchmark. We abbreviate \\textit{picks per event} as \\textit{ppe}.",
    )


@manager.register_with_markers(["table", "paper"])
def scenario_tables():
    base = Path("synthetcs")
    stats = []
    for path in base.iterdir():
        if not path.is_dir() or path.name == "base":
            continue
        else:
            exp, events, noise = path.name.split("_")
            events = int(events)
            picks = pd.read_parquet(path / "picks")
            stations = pd.read_parquet(path / "stations")
            event_picks = picks[picks["event"] != -1]
            noise_picks = picks[picks["event"] == -1]

            path_stats = {
                "exp": exp,
                "Events": events,
                "Noise": noise,
                "Event picks": len(event_picks),
                "Noise picks": len(noise_picks),
                "Total picks": len(picks),
                "Picks per event": len(event_picks) / events,
                "Picks per station": len(picks) / len(stations),
            }
            stats.append(path_stats)

    stats = pd.DataFrame(stats)

    desc_lookup = {
        "shallow": "shallow seismicity",
        "chile": "subduction",
    }

    for exp, sub in stats.groupby("exp"):
        sub = sub.sort_values(["Events", "Noise"])
        table_path = paper_table_path / f"statistics_{exp}.tex"
        sub.drop(columns=["exp"]).to_latex(
            table_path,
            index=False,
            float_format="%.2f",
            longtable=False,
            label=f"tab:statistics_{exp}",
            caption=f"Dataset statistics for the {desc_lookup[exp]} scenario. "
            f"We do not differentiate between P and S picks as both are generated in almost equal number. "
            f"The picks per station include the noise picks.",
            formatters={
                "Events": decimal_formatter,
                "Event picks": decimal_formatter,
                "Noise picks": decimal_formatter,
                "Total picks": decimal_formatter,
            },
        )

        with open(table_path, "r") as f:
            cont = f.read()
        with open(table_path, "w") as f:
            f.write(cont.replace("{table}", "{table*}"))


@manager.register_with_markers(["plot", "paper"])
def benchmark_results():
    results = get_benchmark_results()
    for exp, sub in results.groupby("exp"):
        fig = plt.figure(figsize=(textwidth, 1.37 * textwidth))
        axs = fig.subplots(5, 1, sharex=True, gridspec_kw={"hspace": 0.05})

        axs[0].set_ylabel("Precision")
        axs[1].set_ylabel("Recall")
        axs[2].set_ylabel("F1 score")
        axs[3].set_ylabel("Missing/additional\npicks per event")
        axs[4].set_ylabel("Run time [s]")

        events = sub["events"].unique()
        noises = sub["noise"].unique()
        associators = sub["associator"].unique()

        wi = 0.1
        wj = len(associators) * wi + wi
        wk = len(noises) * wj + wi

        ticks = []
        labels = []

        x_pos = [0.1, 0.1, 0.1, 0, 3]

        for k, event in enumerate(events):
            for j, noise in enumerate(noises):
                ticks.append(j * wj + k * wk + (len(associators) // 2) * wi)
                labels.append(f"{event} Events\n{noise:.1f} Noise")
                for i, associator in enumerate(associators):
                    px = i * wi + j * wj + k * wk
                    row = sub[
                        (sub["events"] == event)
                        & (sub["noise"] == noise)
                        & (sub["associator"] == associator)
                    ]
                    assert (
                        len(row) <= 1
                    ), f"Found more than one entry for {exp} {event} {noise} {associator}"

                    if len(row) == 0:
                        for ax, y in zip(axs, x_pos):
                            ax.text(
                                px,
                                y,
                                "X",
                                color="gray",
                                weight="bold",
                                va="center",
                                ha="center",
                            )
                    else:
                        row = row.iloc[0]

                    axs[0].bar(px, row["precision"], width=wi, color=f"C{i}")
                    axs[1].bar(px, row["recall"], width=wi, color=f"C{i}")
                    axs[2].bar(px, row["f1"], width=wi, color=f"C{i}")
                    axs[3].bar(
                        px, row["additional_picks_per_event"], width=wi, color=f"C{i}"
                    )
                    axs[3].bar(
                        px, -row["missing_picks_per_event"], width=wi, color=f"C{i}"
                    )
                    axs[4].bar(px, row["runtime"], width=wi, color=f"C{i}")

        for ax in axs[:3]:
            ax.set_ylim(0, 1.03)

        for i, associator in enumerate(associators):
            axs[-1].bar(-1, 1, color=f"C{i}", label=associator)
        axs[-1].legend(loc="upper left")

        axs[-1].set_yscale("log")
        axs[-1].set_xticks(ticks)
        axs[-1].set_xticklabels(labels, rotation=90)
        axs[-1].set_xlim(ticks[0] - 3 * wi, ticks[-1] + 3 * wi)
        fig.savefig(paper_figure_path / f"results_{exp}.png", bbox_inches="tight")


@manager.register_with_markers(["plot", "paper"])
def iquique_sections():
    fig = plt.figure(figsize=(textwidth, 1.37 * textwidth))
    exps = ["pyocto_cutoff", "pyocto_1d_cutoff", "real", "real_1d", "gamma"]
    base_path = Path("/home/munchmej/code/ml-catalog/catalogs")

    axs = fig.subplots(
        len(exps), 1, sharex=True, sharey=True, gridspec_kw={"hspace": 0.05}
    )

    for i, exp in enumerate(exps):
        ax = axs[i]
        events = pd.read_csv(base_path / f"benchmark_{exp}{iquique_suffix}/events.csv")
        ax.scatter(events["longitude"], events["depth"], s=1, c="k")
        ax.set_ylabel(name_lookup[exp] + "\nDepth [km]")

    axs[-1].set_xlim(-71.3, -68.1)
    axs[-1].set_ylim(250, -10)
    axs[-1].set_xlabel("Longitude [$^\circ$]")

    fig.savefig(paper_figure_path / f"iquique_sections.png", bbox_inches="tight")


@manager.register_with_markers(["plot", "paper"])
def iquique_section_pyocto():
    fig = plt.figure(figsize=(0.5 * textwidth, 0.3 * textwidth))
    exp = "pyocto_1d_cutoff"
    base_path = Path("/home/munchmej/code/ml-catalog/catalogs")

    ax = fig.add_subplot(111)

    events = pd.read_csv(base_path / f"benchmark_{exp}{iquique_suffix}/events.csv")
    ax.scatter(events["longitude"], events["depth"], s=0.7, c="k", lw=0)
    ax.set_ylabel(name_lookup[exp] + "\nDepth [km]")

    ax.set_xlim(-71.3, -68.1)
    ax.set_ylim(250, -10)
    ax.set_xlabel("Longitude [$^\circ$]")

    fig.savefig(paper_figure_path / f"iquique_section_pyocto.png", bbox_inches="tight")


@manager.register_with_markers(["plot", "paper"])
def iquique_maps():
    fig = plt.figure(figsize=(textwidth, 1.2 * textwidth))
    exps = ["pyocto_cutoff", "pyocto_1d_cutoff", "real", "real_1d", "gamma"]
    base_path = Path("/home/munchmej/code/ml-catalog/catalogs")

    axs = []
    for i, _ in enumerate(exps):
        axs.append(fig.add_subplot(23 * 10 + i + 1, projection=ccrs.Mercator()))

    for i, exp in enumerate(exps):
        ax = axs[i]
        ax.coastlines(zorder=51)
        events = pd.read_csv(base_path / f"benchmark_{exp}{iquique_suffix}/events.csv")
        cb = ax.scatter(
            events["longitude"],
            events["latitude"],
            s=1,
            c=events["depth"],
            transform=ccrs.PlateCarree(),
            cmap=cm.batlow,
            vmin=0,
            vmax=250,
            zorder=100,
        )
        ax.set_extent([-71.4, -68, -25, -18], crs=ccrs.PlateCarree())

        ax.set_title(name_lookup[exp])

        gl = ax.gridlines(
            draw_labels=True,
            dms=True,
            x_inline=False,
            y_inline=False,
            zorder=50,
            linewidth=0.5,
        )
        gl.right_labels = False
        gl.bottom_labels = False

    pos = axs[-1].get_position().corners()
    left = axs[2].get_position().corners()[0, 0]
    cax = fig.add_axes([left, pos[0, 1], 0.10, (pos[1, 1] - pos[0, 1])])
    fig.colorbar(cb, cax=cax, label="Depth [km]")
    cax.invert_yaxis()

    fig.savefig(paper_figure_path / f"iquique_maps.png", bbox_inches="tight")


@manager.register_with_markers(["plot", "paper"])
def iquique_rates():
    fig = plt.figure(figsize=(colwidth, 0.8 * textwidth))
    exps = ["gamma", "pyocto_cutoff", "pyocto_1d_cutoff", "real", "real_1d"]
    base_path = Path("/home/munchmej/code/ml-catalog/catalogs")

    axs = fig.subplots(3, 1, sharex=True, gridspec_kw={"hspace": 0.06})

    days = None
    for exp in exps:
        events = pd.read_csv(base_path / f"benchmark_{exp}{iquique_suffix}/events.csv")
        if days is None:
            days = sorted(events["group"].unique())

        stats = {}
        for day, day_events in events.groupby("group"):
            stats[day] = (
                len(day_events),
                day_events["number_p_picks"].mean(),
                day_events["number_s_picks"].mean(),
            )
        stats = np.array([stats[day] for day in days])
        x_days = np.arange(len(days)) + 0.5
        axs[0].step(x_days, stats[:, 0], label=name_lookup[exp])
        axs[1].step(x_days, stats[:, 1])
        axs[2].step(x_days, stats[:, 2])

    for ax in axs:
        ax.axvline(17.5, c="k", linestyle="--", lw=1.5)
        ax.axvline(1.5, c="k", linestyle="--", lw=1.0)

    w = 5
    axs[-1].set_xticks(np.arange(len(days))[1::w])
    labels = [day[:10] for day in days]
    axs[-1].set_xticklabels(labels[1::w], rotation=90)

    axs[0].set_ylim(0, 1050)
    axs[1].set_ylim(8, 12)
    axs[2].set_ylim(5.5, 8.5)

    axs[-1].set_xlim(0.5, 30.5)

    axs[0].set_ylabel("# Events per day")
    axs[1].set_ylabel("P picks per event")
    axs[2].set_ylabel("S picks per event")
    axs[0].legend(ncol=3, bbox_to_anchor=(0.5, 1.03), loc="lower center")

    fig.savefig(paper_figure_path / f"iquique_rates.png", bbox_inches="tight")


def load_picks(path: Path):
    picks = []
    for file in (path / "SeisBenchPicker/picks").iterdir():
        picks.append(pd.read_parquet(file))
    return pd.concat(picks)


@manager.register_with_markers(["table", "paper"])
def iquique_table():
    exps = ["gamma", "pyocto", "pyocto_1d", "real", "real_1d"]
    base_path = Path("/home/munchmej/code/ml-catalog/catalogs")

    run_times = {
        "pyocto": 21,
        "pyocto_1d": 30,
        "real": 1487,
        "real_1d": 1557,
        "gamma": 1021,
    }

    stats = []
    for exp in exps:
        exp_path = base_path / f"benchmark_{exp}{iquique_suffix}"
        events = pd.read_csv(exp_path / "events.csv")
        picks = load_picks(exp_path)
        stats.append(
            {
                "Associator": name_lookup[exp],
                "Events": len(events),
                "Ppe": events["number_picks"].sum() / len(events),
                "P ppe": events["number_p_picks"].sum() / len(events),
                "S ppe": events["number_s_picks"].sum() / len(events),
                "Associated": np.sum(events["number_picks"]) / len(picks),
                "P associated": np.sum(events["number_p_picks"])
                / np.sum(picks["phase"] == "P"),
                "S associated": np.sum(events["number_s_picks"])
                / np.sum(picks["phase"] == "S"),
                "Total picks": len(picks),
                "Time [s]": run_times[exp],
            }
        )

    stats = pd.DataFrame(stats)

    table_path = paper_table_path / "iquique.tex"
    stats.to_latex(
        table_path,
        index=False,
        float_format="%.2f",
        formatters={
            "Events": decimal_formatter,
            "Total picks": decimal_formatter,
        },
        longtable=False,
        label=f"tab:iquique",
        caption="Catalog statistics for the Iquique sequence catalog with different associators. "
        "The table shows the number of events, picks per event, the fraction of "
        "associated picks among all picks, and the total number of picks. "
        "We abbreviate \\textit{picks per event} as \\textit{ppe}. "
        "Times refer to average run times per day of data.",
    )

    with open(table_path, "r") as f:
        cont = f.read()
    with open(table_path, "w") as f:
        f.write(cont.replace("{table}", "{table*}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--function", type=str, required=False)
    parser.add_argument("--marker", type=str, required=False)
    args = parser.parse_args()

    init_style()

    if args.function is not None:
        functions = args.function.split(",")
        for function in functions:
            manager.run(function)

    if args.marker is not None:
        markers = args.marker.split(",")
        for marker in markers:
            manager.run_marker(marker)
