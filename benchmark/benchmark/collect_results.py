from pathlib import Path

import pandas as pd
from omegaconf import OmegaConf


def main():
    base_path = Path("pred")

    configs = []
    for exp_path in base_path.iterdir():
        if not exp_path.is_dir():
            continue
        try:
            exp_cfg = OmegaConf.load(exp_path / "config.yaml")
            metrics = OmegaConf.load(exp_path / "metrics.yaml")
        except FileNotFoundError:
            continue

        exp_cfg.update(metrics)
        exp_cfg["pred_path"] = exp_path.name
        configs.append(exp_cfg)

    results = reduce_configs(configs)
    results.to_csv(base_path / "results.csv", index=False)


def reduce_configs(configs: list[OmegaConf]) -> pd.DataFrame:
    configs = pd.json_normalize([OmegaConf.to_container(c) for c in configs], sep=".")

    if len(configs) == 1:
        return configs

    del_list = []
    for col in configs.columns:
        if col in [
            "precision",
            "recall",
            "f1",
            "missing_picks_per_event",
            "additional_picks_per_event",
            "runtime",
        ]:
            continue

        # Convert columns to str to allow unifying
        if configs[col].dtype == "O":
            configs[col] = configs[col].astype("str")

        if len(configs[col].unique()) == 1:
            del_list.append(col)

    configs.drop(columns=del_list, inplace=True)

    return configs


if __name__ == "__main__":
    main()
