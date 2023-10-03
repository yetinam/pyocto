import datetime
import logging
import time
from pathlib import Path
from typing import Union

import hydra
import numpy as np
import pandas as pd
from omegaconf import OmegaConf

from .wrappers import AbstractAssociator

logger = logging.getLogger("benchmark")


@hydra.main(version_base=None, config_path="../configs", config_name="eval.yaml")
def main(cfg):
    if "log_level" in cfg:
        logger.setLevel(cfg["log_level"])

    output_path = get_output_path(cfg.base_path)
    logger.debug(f"Output path is {output_path}")
    with open(output_path / "config.yaml", "w") as f:
        OmegaConf.save(cfg, f)

    logger.debug("Setting up")
    cfg = hydra.utils.instantiate(cfg)

    associator: AbstractAssociator = cfg.associator

    associator.setup(cfg.data["stations"])

    logger.debug(f"Starting association")
    t0 = time.time()
    catalog, assignments = associator.get_events(
        cfg.data["picks"], cfg.data["stations"]
    )
    t1 = time.time()

    logger.debug(f"Writing outputs")
    catalog.to_parquet(output_path / "catalog", index=False)
    assignments.to_parquet(output_path / "assignments", index=False)

    logger.debug(f"Calculating metrics")
    metrics = get_metrics(cfg.data["catalog"], cfg.data["picks"], catalog, assignments)
    metrics["runtime"] = t1 - t0
    metrics = {
        k: float(v) for k, v in metrics.items()
    }  # OmegaConf crashes on np.float64
    metrics = OmegaConf.create(metrics)
    with open(output_path / "metrics.yaml", "w") as f:
        OmegaConf.save(metrics, f)


def load_data(path: str) -> dict[str, pd.DataFrame]:
    path = Path(path)
    return {
        "catalog": pd.read_parquet(path / "catalog"),
        "picks": pd.read_parquet(path / "picks"),
        "stations": pd.read_parquet(path / "stations"),
    }


def get_output_path(base: Union[Path, str]) -> Path:
    base = Path(base)
    base.mkdir(exist_ok=True, parents=True)
    for i in range(100000):
        output_path = base / str(i)
        try:
            output_path.mkdir(exist_ok=False)
            break
        except FileExistsError:
            pass
    else:
        raise ValueError("Did not find free output path")

    return output_path


def get_metrics(
    true_catalog: pd.DataFrame,
    true_assignments: pd.DataFrame,
    pred_catalog: pd.DataFrame,
    pred_assignments: pd.DataFrame,
    threshold: float = 0.6,
) -> dict:
    matches = []
    for new_idx, new_df in pred_assignments.groupby("event_idx"):
        for old_idx, old_df in true_assignments.groupby("event"):
            new_picks = set(new_df["pick_idx"])
            old_picks = set(old_df.index)
            intersect = new_picks & old_picks

            score = len(intersect) / max(len(new_picks), len(old_picks))

            if score >= threshold:
                matches.append(
                    (
                        new_idx,
                        old_idx,
                        score,
                        len(new_picks),
                        len(old_picks),
                        len(intersect),
                    )
                )
                break  # at most one match per event

    matches = pd.DataFrame(
        matches,
        columns=["new_idx", "old_idx", "score", "new_picks", "old_picks", "intersect"],
    )

    tp = len(matches)
    fn = len(true_catalog) - len(matches)
    fp = len(pred_catalog) - len(matches)

    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1 = 2 * (prec * rec) / (prec + rec)

    true_picks_per_event = np.sum(matches["old_picks"]) / len(matches)
    pred_picks_per_event = np.sum(matches["new_picks"]) / len(matches)
    additional_picks_per_event = np.sum(
        matches["new_picks"] - matches["intersect"]
    ) / len(matches)
    missing_picks_per_event = np.sum(matches["old_picks"] - matches["intersect"]) / len(
        matches
    )

    return {
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "true_picks_per_event": true_picks_per_event,
        "pred_picks_per_event": pred_picks_per_event,
        "additional_picks_per_event": additional_picks_per_event,
        "missing_picks_per_event": missing_picks_per_event,
    }


if __name__ == "__main__":
    main()
