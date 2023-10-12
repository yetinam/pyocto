//
// Created by munchmej on 04/07/23.
//

#include "OctoAssociator.h"
#include "VelocityModel.h"
#include "types.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace octoassociator {

PYBIND11_MODULE(_core, m) {
  m.doc() = "pybind11 example plugin"; // optional module docstring

  py::class_<OctoTreeConfig>(m, "OctoTreeConfig")
      .def(py::init())
      .def_readwrite("min_node_size", &OctoTreeConfig::min_node_size)
      .def_readwrite("min_node_size_location",
                     &OctoTreeConfig::min_node_size_location)
      .def_readwrite("pick_match_tolerance",
                     &OctoTreeConfig::pick_match_tolerance)
      .def_readwrite("min_interevent_time",
                     &OctoTreeConfig::min_interevent_time)
      .def_readwrite("edt_pick_std", &OctoTreeConfig::edt_pick_std)
      .def_readwrite("exponential_edt", &OctoTreeConfig::exponential_edt)
      .def_readwrite("max_pick_overlap", &OctoTreeConfig::max_pick_overlap)
      .def_readwrite("n_picks", &OctoTreeConfig::n_picks)
      .def_readwrite("n_p_picks", &OctoTreeConfig::n_p_picks)
      .def_readwrite("n_s_picks", &OctoTreeConfig::n_s_picks)
      .def_readwrite("n_p_and_s_picks", &OctoTreeConfig::n_p_and_s_picks)
      .def_readwrite("n_threads", &OctoTreeConfig::n_threads)
      .def_readwrite("refinement_iterations",
                     &OctoTreeConfig::refinement_iterations)
      .def_readwrite("node_log_interval", &OctoTreeConfig::node_log_interval)
      .def_readwrite("queue_memory_protection_dfs_size",
                     &OctoTreeConfig::queue_memory_protection_dfs_size)
      .def_readwrite("time_before", &OctoTreeConfig::time_before)
      .def_readwrite("time_slicing", &OctoTreeConfig::time_slicing)
      .def_readwrite("location_split_depth",
                     &OctoTreeConfig::location_split_depth)
      .def_readwrite("location_split_return",
                     &OctoTreeConfig::location_split_return)
      .def_readwrite("min_pick_fraction", &OctoTreeConfig::min_pick_fraction)
      .def_readwrite("xlim", &OctoTreeConfig::xlim)
      .def_readwrite("ylim", &OctoTreeConfig::ylim)
      .def_readwrite("zlim", &OctoTreeConfig::zlim)
      .def_readwrite("velocity_model_association",
                     &OctoTreeConfig::velocity_model_association)
      .def_readwrite("velocity_model_location",
                     &OctoTreeConfig::velocity_model_location);

  auto velocity_model = py::class_<VelocityModel>(m, "VelocityModel");
  velocity_model.def("add_station", &VelocityModel::add_station);

  py::class_<VelocityModel0D>(m, "VelocityModel0D", velocity_model)
      .def(py::init<double, double, double, double, double>())
      .def_readwrite("p_velocity", &VelocityModel0D::p_velocity)
      .def_readwrite("s_velocity", &VelocityModel0D::s_velocity)
      .def_readwrite("association_cutoff_distance",
                     &VelocityModel0D::association_cutoff_distance)
      .def_readwrite("location_cutoff_distance",
                     &VelocityModel0D::location_cutoff_distance)
      .def_readwrite("tolerance", &VelocityModel0D::tolerance);

  py::class_<VelocityModel1D>(m, "VelocityModel1D", velocity_model)
      .def(py::init<char *>())
      .def_readwrite("surface_p_velocity", &VelocityModel1D::surface_p_velocity)
      .def_readwrite("surface_s_velocity", &VelocityModel1D::surface_s_velocity)
      .def_readwrite("association_cutoff_distance",
                     &VelocityModel1D::association_cutoff_distance)
      .def_readwrite("location_cutoff_distance",
                     &VelocityModel1D::location_cutoff_distance)
      .def_readwrite("tolerance", &VelocityModel1D::tolerance);

  py::class_<Pick>(m, "Pick")
      .def(py::init<int, double, std::string, char>())
      .def_readonly("idx", &Pick::idx)
      .def_readonly("time", &Pick::time)
      .def_readonly("station", &Pick::station)
      .def_readonly("phase", &Pick::phase);

  py::class_<Event>(m, "Event")
      .def_readonly("time", &Event::time)
      .def_readonly("x", &Event::x)
      .def_readonly("y", &Event::y)
      .def_readonly("z", &Event::z)
      .def_readonly("picks", &Event::picks)
      .def_readonly("residuals", &Event::residuals);

  py::class_<Station>(m, "Station")
      .def(py::init<std::string, double, double, double, double, double>())
      .def_readonly("id", &Station::id)
      .def_readonly("x", &Station::x)
      .def_readonly("y", &Station::y)
      .def_readonly("z", &Station::z)
      .def_readonly("p_residual", &Station::p_residual)
      .def_readonly("s_residual", &Station::s_residual);

  py::class_<OctoAssociator>(m, "OctoAssociator")
      .def(py::init<OctoTreeConfig *>())
      .def("associate", &OctoAssociator::associate,
           py::return_value_policy::copy);
}

} // namespace octoassociator
