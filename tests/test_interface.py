import logging
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from pyproj import CRS

import pyocto


def test_get_crs():
    stations = pd.DataFrame(
        {
            "latitude": [5, 7, 8, 9, 10],
            "longitude": [-20, -10, -15, -17.5, -19.2],
        }
    )
    crs = pyocto.OctoAssociator.get_crs(stations)
    crs_dict = crs.to_dict()

    assert crs_dict["lat_0"] == 7.5
    assert crs_dict["lon_0"] == -15.0

    # Longitude discontinuity - left side
    stations = pd.DataFrame(
        {
            "latitude": [-20, -10, -15, -19.2],
            "longitude": [178, 178.1, -179, -179.5],
        }
    )
    crs = pyocto.OctoAssociator.get_crs(stations)
    crs_dict = crs.to_dict()

    assert crs_dict["lat_0"] == -15.0
    assert crs_dict["lon_0"] == 179.5

    # Longitude discontinuity - right side
    stations = pd.DataFrame(
        {
            "latitude": [-20, -10, -15, -19.2],
            "longitude": [-178, -178.1, 179, 179.5],
        }
    )
    crs = pyocto.OctoAssociator.get_crs(stations)
    crs_dict = crs.to_dict()

    assert crs_dict["lat_0"] == -15.0
    assert crs_dict["lon_0"] == -179.5


def test_get_crs_warning(caplog):
    caplog.clear()
    stations = pd.DataFrame(
        {
            "latitude": [5, 10],
            "longitude": [-20, -10],
        }
    )
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator.get_crs(stations, warning_limit_deg=15.0)
    assert "inaccurate" not in caplog.text

    caplog.clear()
    stations = pd.DataFrame(
        {
            "latitude": [5, 25],
            "longitude": [-20, -10],
        }
    )
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator.get_crs(stations, warning_limit_deg=15.0)
    assert "inaccurate" in caplog.text

    caplog.clear()
    stations = pd.DataFrame(
        {
            "latitude": [5, 10],
            "longitude": [-20, 10],
        }
    )
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator.get_crs(stations, warning_limit_deg=15.0)
    assert "inaccurate" in caplog.text


def test_unit_assertion():
    velocity_model = pyocto.VelocityModel0D(7.0, 4.0, 2.0)
    associator = pyocto.OctoAssociator((0, 1), (0, 1), (0, 1), velocity_model, 300)

    associator.crs = CRS.from_proj4("+proj=tmerc +units=km")
    assert not associator._meter_model

    associator.crs = CRS.from_proj4("+proj=tmerc +units=m")
    assert associator._meter_model

    with pytest.raises(AssertionError):
        associator.crs = CRS.from_proj4("+proj=tmerc +units=mm")


@pytest.mark.parametrize(
    "velocity_model_1d",
    [False, True],
)
def test_associate(tmp_path, velocity_model_1d):
    if velocity_model_1d:
        model_path = tmp_path / "model"
        layers = pd.read_csv("tests/data/graeber.csv")
        pyocto.VelocityModel1D.create_model(layers, 1, 400, 250, model_path)
        velocity_model = pyocto.VelocityModel1D(model_path, 2.0)
    else:
        velocity_model = pyocto.VelocityModel0D(7.0, 4.0, 2.0)

    stations = pd.read_parquet("tests/data/stations")
    picks = pd.read_parquet("tests/data/picks")
    picks["time"] = picks["time"].apply(lambda x: x.timestamp())

    crs = CRS.from_epsg(9155)

    associator = pyocto.OctoAssociator(
        xlim=(250.0, 600.0),
        ylim=(7200.0, 8000.0),
        zlim=(0.0, 250.0),
        time_before=300.0,
        velocity_model=velocity_model,
        n_picks=10,
        n_p_picks=2,
        n_s_picks=2,
        n_p_and_s_picks=4,
        pick_match_tolerance=2.0,
        crs=crs,
    )

    associator.transform_stations(stations)

    events, assignments = associator.associate(picks, stations)

    # There are 50 events contained so this is a fairly easy condition.
    assert len(events) > 25
    assert len(assignments) > 500


@pytest.mark.parametrize(
    "velocity_model_1d",
    [False, True],
)
def test_associate_cutoff(tmp_path, velocity_model_1d):
    if velocity_model_1d:
        model_path = tmp_path / "model"
        layers = pd.read_csv("tests/data/graeber.csv")
        pyocto.VelocityModel1D.create_model(layers, 1, 400, 250, model_path)
        velocity_model = pyocto.VelocityModel1D(model_path, 2.0)
    else:
        velocity_model = pyocto.VelocityModel0D(7.0, 4.0, 2.0)

    stations = pd.read_parquet("tests/data/stations")
    picks = pd.read_parquet("tests/data/picks")
    picks["time"] = picks["time"].apply(lambda x: x.timestamp())

    crs = CRS.from_epsg(9155)

    associator = pyocto.OctoAssociator(
        xlim=(250.0, 600.0),
        ylim=(7200.0, 8000.0),
        zlim=(0.0, 250.0),
        time_before=300.0,
        velocity_model=velocity_model,
        n_picks=10,
        n_p_picks=2,
        n_s_picks=2,
        n_p_and_s_picks=4,
        pick_match_tolerance=2.0,
        crs=crs,
    )

    associator.transform_stations(stations)

    velocity_model.association_cutoff_distance = 10
    events, assignments = associator.associate(picks, stations)

    # There are 50 events contained so this is a fairly easy condition.
    assert len(events) == 0
    assert len(assignments) == 0


def test_create_model(tmp_path):
    model = pd.DataFrame(
        {
            "depth": [0, 10, 50, 100],
            "vp": [4, 5, 6, 7],
            "vs": [3, 4, 5, 6],
        }
    )

    # Regular call
    pyocto.VelocityModel1D.create_model(model, 2, 200, 100, tmp_path / "test")

    # zdist shallower than deepest element
    pyocto.VelocityModel1D.create_model(model, 2, 200, 20, tmp_path / "test")

    # Negative entry
    model = pd.DataFrame(
        {
            "depth": [-5, 10, 50, 100],
            "vp": [4, 5, 6, 7],
            "vs": [3, 4, 5, 6],
        }
    )
    pyocto.VelocityModel1D.create_model(model, 2, 200, 100, tmp_path / "test")


def test_adjust_depth_velocity():
    # No cutoff
    depth1, depth2, vp1, vp2, vs1, vs2 = pyocto.VelocityModel1D._adjust_depth_velocity(
        5, 15, 1, 2, 3, 4, 20
    )
    assert depth1 == 5
    assert depth2 == 15
    assert vp1 == 1
    assert vp2 == 2
    assert vs1 == 3
    assert vs2 == 4

    # Deep cutoff
    depth1, depth2, vp1, vp2, vs1, vs2 = pyocto.VelocityModel1D._adjust_depth_velocity(
        5, 15, 1, 2, 3, 4, 10
    )
    assert depth1 == 5
    assert depth2 == 10
    assert vp1 == 1
    assert vp2 == 1.5
    assert vs1 == 3
    assert vs2 == 3.5

    # Shallow cutoff
    depth1, depth2, vp1, vp2, vs1, vs2 = pyocto.VelocityModel1D._adjust_depth_velocity(
        -5, 5, 1, 2, 3, 4, 10
    )
    assert depth1 == 0
    assert depth2 == 5
    assert vp1 == 1.5
    assert vp2 == 2
    assert vs1 == 3.5
    assert vs2 == 4


def test_windows_warning(caplog):
    with patch("os.name", "posix"):
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            pyocto.OctoAssociator(
                xlim=(250.0, 600.0),
                ylim=(7200.0, 8000.0),
                zlim=(0.0, 250.0),
                time_before=300.0,
                velocity_model=None,
            )
        assert "on Windows" not in caplog.text

    with patch("os.name", "nt"):
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            pyocto.OctoAssociator(
                xlim=(250.0, 600.0),
                ylim=(7200.0, 8000.0),
                zlim=(0.0, 250.0),
                time_before=300.0,
                velocity_model=None,
            )
        assert "on Windows" in caplog.text


def test_pick_count_warnings(caplog):
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator(
            xlim=(250.0, 600.0),
            ylim=(7200.0, 8000.0),
            zlim=(0.0, 250.0),
            time_before=300.0,
            velocity_model=None,
            n_picks=10,
            n_p_picks=5,
            n_s_picks=5,
            n_p_and_s_picks=5,
        )
    assert "The required number" not in caplog.text

    caplog.clear()
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator(
            xlim=(250.0, 600.0),
            ylim=(7200.0, 8000.0),
            zlim=(0.0, 250.0),
            time_before=300.0,
            velocity_model=None,
            n_picks=7,
            n_p_picks=5,
            n_s_picks=5,
            n_p_and_s_picks=3,
        )
    assert "The required number of picks per event " in caplog.text

    caplog.clear()
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator(
            xlim=(250.0, 600.0),
            ylim=(7200.0, 8000.0),
            zlim=(0.0, 250.0),
            time_before=300.0,
            velocity_model=None,
            n_picks=7,
            n_p_picks=4,
            n_s_picks=3,
            n_p_and_s_picks=4,
        )
    assert "The required number of picks per event " in caplog.text

    caplog.clear()
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator(
            xlim=(250.0, 600.0),
            ylim=(7200.0, 8000.0),
            zlim=(0.0, 250.0),
            time_before=300.0,
            velocity_model=None,
            n_picks=20,
            n_p_picks=3,
            n_s_picks=4,
            n_p_and_s_picks=4,
        )
    assert "The required number of P picks per event" in caplog.text

    caplog.clear()
    with caplog.at_level(logging.WARNING):
        pyocto.OctoAssociator(
            xlim=(250.0, 600.0),
            ylim=(7200.0, 8000.0),
            zlim=(0.0, 250.0),
            time_before=300.0,
            velocity_model=None,
            n_picks=20,
            n_p_picks=4,
            n_s_picks=3,
            n_p_and_s_picks=4,
        )
    assert "The required number of S picks per event" in caplog.text


@pytest.mark.parametrize(
    "crs",
    [
        CRS(f"+proj=tmerc +lat_0={-22} +lon_0={-70} +units=m"),
        CRS(f"+proj=tmerc +lat_0={-22} +lon_0={-70} +units=km"),
    ],
)
def test_coordinate_transform_invtransform(crs):
    # Test coordinate transforms back and forth with m and km scale
    n = 20
    lat0 = -22
    lon0 = -70
    stations = pd.DataFrame(
        {
            "id": [f"S{i:04d}" for i in range(n)],
            "latitude": lat0 + np.random.randn(n),
            "longitude": lon0 + np.random.randn(n),
            "elevation": 1000 * np.random.random(n),
        }
    )
    associator = pyocto.OctoAssociator(
        xlim=(250.0, 600.0),
        ylim=(7200.0, 8000.0),
        zlim=(0.0, 250.0),
        time_before=300.0,
        velocity_model=None,
        crs=crs,
    )

    associator.transform_stations(stations)
    stations2 = stations.copy()
    stations2.drop(columns=["latitude", "longitude"], inplace=True)

    associator.transform_events(stations2)
    assert np.allclose(stations["latitude"], stations2["latitude"])
    assert np.allclose(stations["longitude"], stations2["longitude"])
