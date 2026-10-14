//
// Created by munchmej on 17/06/23.
//

#include <algorithm>
#include <cmath>

#include "VelocityModel.h"

std::pair<double, double>
octoassociator::VelocityModel::min_max_dist(double x, double wx, double px) {
  double min_dist, max_dist;

  if (x - wx / 2 <= px and px <= x + wx / 2) {
    min_dist = 0; // Point in interval
  } else {
    min_dist = std::min(std::abs(x - wx / 2 - px), std::abs(x + wx / 2 - px));
  }

  max_dist = std::max(std::abs(x - wx / 2 - px), std::abs(x + wx / 2 - px));

  return {min_dist, max_dist};
}

namespace octoassociator {

std::vector<unsigned int>
VelocityModel::get_active_stations(const std::vector<Pick *> &picks) {
  std::vector<unsigned int> active_stations;
  for (const auto &pick : picks) {
    if (pick->station < stations.size()) {
      active_stations.push_back(pick->station);
    }
  }
  std::sort(active_stations.begin(), active_stations.end());
  active_stations.erase(
      std::unique(active_stations.begin(), active_stations.end()),
      active_stations.end());
  return active_stations;
}

bool VelocityModel0D::contains(const Volume &volume, const Pick *pick) {
  if (pick->station >= stations.size()) {
    printf("Warning: Station %u not found\n", pick->station);
    return false;
  }
  auto &station = stations[pick->station];
  double min_dist, max_dist, centroid_dist;
  distance_volume_pick(volume, min_dist, max_dist, centroid_dist, station);

  if (std::isnan(max_dist) or std::isnan(centroid_dist))
    return false;

  if (min_dist > association_cutoff_distance)
    return false; // Fail if outside cutoff distance

  max_dist = std::min(max_dist,
                      association_cutoff_distance); // Truncate maximum distance

  double velocity = get_phase_velocity(pick->phase);

  // First and last arrival time from a pick in the volume
  double t0 = volume.t - volume.wt / 2 + station.residual(pick->phase) -
              tolerance + min_dist / velocity;
  double t1 = volume.t + volume.wt / 2 + station.residual(pick->phase) +
              tolerance + max_dist / velocity;

  return (t0 < pick->time and pick->time < t1);
}

double VelocityModel0D::get_phase_velocity(char phase) const {
  if (phase == 'P') {
    return p_velocity;
  } else {
    return s_velocity;
  }
}

void VelocityModel0D::distance_volume_pick(const Volume &volume,
                                           double &min_dist, double &max_dist,
                                           double &centroid_dist,
                                           const Station &station) {
  auto x_dist = min_max_dist(volume.x, volume.wx, station.x);
  auto y_dist = min_max_dist(volume.y, volume.wy, station.y);
  auto z_dist = min_max_dist(volume.z, volume.wz, station.z);

  min_dist =
      sqrt(pow(x_dist.first, 2) + pow(y_dist.first, 2) + pow(z_dist.first, 2));
  max_dist = sqrt(pow(x_dist.second, 2) + pow(y_dist.second, 2) +
                  pow(z_dist.second, 2));
  centroid_dist =
      sqrt(pow(volume.x - station.x, 2) + pow(volume.y - station.y, 2) +
           pow(volume.z - station.z, 2));
}

double VelocityModel0D::travel_time(const Volume &volume, unsigned int station,
                                    char phase) {
  double velocity = get_phase_velocity(phase);
  if (station >= stations.size()) {
    printf("Warning: Station %u not found\n", station);
    return NAN;
  }
  auto &station_obj = stations[station];

  double dist =
      sqrt(pow(volume.x - station_obj.x, 2) + pow(volume.y - station_obj.y, 2) +
           pow(volume.z - station_obj.z, 2));

  if (dist > location_cutoff_distance) {
    return NAN;
  } else {
    return dist / velocity + station_obj.residual(phase);
  }
}

VelocityModel1D::VelocityModel1D(char *path) {
  auto f = fopen(path, "rb");
  fread(&nx, 4, 1, f);
  fread(&nz, 4, 1, f);
  fread(&delta, 8, 1, f);

  p_times = new double[nx * nz];
  s_times = new double[nx * nz];
  for (int i = 0; i < nx * nz; i++) {
    fread(p_times + i, 8, 1, f);
  }
  for (int i = 0; i < nx * nz; i++) {
    fread(s_times + i, 8, 1, f);
  }
  fclose(f);

  precalculate_extrema();
}

double VelocityModel1D::travel_time(const Volume &volume, unsigned int station,
                                    char phase) {
  return travel_time_for_station(volume, station, phase);
}

double VelocityModel1D::travel_time_for_station(const Volume &volume,
                                                unsigned int station,
                                                char phase) {
  if (station >= stations.size()) {
    printf("Warning: Station %u not found\n", station);
    return NAN;
  }
  const auto &station_obj = stations[station];

  double dx =
      sqrt(pow(volume.x - station_obj.x, 2) + pow(volume.y - station_obj.y, 2));
  double tt = travel_time(dx, volume.z, phase);
  double elev_correction = elevation_correction(phase, station_obj);

  return tt + elev_correction + station_obj.residual(phase);
}

double VelocityModel1D::travel_time(double x, double z, char phase) {
  if (x > location_cutoff_distance)
    return NAN;

  double *times;
  if (phase == 'P')
    times = p_times;
  else
    times = s_times;

  int ix = floor(x / delta), iz = floor(z / delta);
  double ax = x / delta - ix, az = z / delta - iz;

  if (ix >= nx - 1)
    return NAN;
  if (iz >= nz - 1)
    return NAN;

  return times[index(ix, iz)] * (1 - ax) * (1 - az) +
         times[index(ix + 1, iz)] * ax * (1 - az) +
         times[index(ix, iz + 1)] * (1 - ax) * az +
         times[index(ix + 1, iz + 1)] * ax * az;
}

double VelocityModel1D::elevation_correction(char phase,
                                             const Station &station) const {
  if (phase == 'P')
    return -station.z / surface_p_velocity;
  else
    return -station.z / surface_s_velocity;
}

bool VelocityModel1D::contains(const Volume &volume, const Pick *pick) {
  return contains(volume, pick, pick->station);
}

bool VelocityModel1D::contains(const Volume &volume, const Pick *pick,
                               unsigned int station) {
  // The function is only exact up to grid points
  // Calculate min-max distance 2D
  if (station >= stations.size()) {
    printf("Warning: Station %u not found\n", station);
    return false;
  }
  auto &station_obj = stations[station];

  auto x_dist = min_max_dist(volume.x, volume.wx, station_obj.x);
  auto y_dist = min_max_dist(volume.y, volume.wy, station_obj.y);

  auto min_dist = sqrt(pow(x_dist.first, 2) + pow(y_dist.first, 2));
  auto max_dist = sqrt(pow(x_dist.second, 2) + pow(y_dist.second, 2));

  if (min_dist > association_cutoff_distance)
    return false;
  max_dist = std::min(association_cutoff_distance, max_dist);

  // Determine indices
  int ix_min = floor(min_dist / delta); // Take closer one (underestimate min)
  int ix_max = ceil(max_dist / delta);  // Take farther one (overestimate max)
  if (ix_min >= nx)
    return false;
  if (ix_max >= nx)
    ix_max = nx - 1;

  // Find local maxima/minima
  std::vector<int> local_minima, local_maxima;
  double *times;
  if (pick->phase == 'P') {
    local_minima = local_minima_p[ix_min];
    local_maxima = local_maxima_p[ix_max];
    times = p_times;
  } else {
    local_minima = local_minima_s[ix_min];
    local_maxima = local_maxima_s[ix_max];
    times = s_times;
  }

  int iz_min = floor(volume.z_min() / delta);
  int iz_max = ceil(volume.z_max() / delta);

  if (iz_min < 0)
    iz_min = 0;
  if (iz_min >= nz)
    return false;
  if (iz_max >= nz)
    iz_max = nz - 1;

  // Candidates for extrema are the end points of the depth interval and all
  // extrema in between
  std::vector<double> min_candidates = {times[index(ix_min, iz_min)],
                                        times[index(ix_min, iz_max)]};
  std::vector<double> max_candidates = {times[index(ix_max, iz_min)],
                                        times[index(ix_max, iz_max)]};

  for (auto z_idx : local_minima) {
    if (volume.z_min() <= z_idx * delta and z_idx * delta <= volume.z_max())
      min_candidates.push_back(times[index(ix_min, z_idx)]);
  }
  for (auto z_idx : local_maxima) {
    if (volume.z_min() <= z_idx * delta and z_idx * delta <= volume.z_max())
      max_candidates.push_back(times[index(ix_max, z_idx)]);
  }

  auto min_tt = *std::min_element(min_candidates.begin(), min_candidates.end());
  auto max_tt = *std::max_element(max_candidates.begin(), max_candidates.end());

  double elev_correction = elevation_correction(pick->phase, station_obj);

  auto t0 = volume.t_min() + min_tt + elev_correction +
            station_obj.residual(pick->phase) - tolerance;
  auto t1 = volume.t_max() + max_tt + elev_correction +
            station_obj.residual(pick->phase) + tolerance;

  return (t0 < pick->time and pick->time < t1);
}

void VelocityModel1D::precalculate_extrema() {
  local_maxima_p.resize(nx);
  local_minima_p.resize(nx);
  local_maxima_s.resize(nx);
  local_minima_s.resize(nx);
  for (int x_idx = 0; x_idx < nx; ++x_idx) {
    for (int z_idx = 1; z_idx < nz - 1; ++z_idx) {
      if (p_times[index(x_idx, z_idx - 1)] <= p_times[index(x_idx, z_idx)] and
          p_times[index(x_idx, z_idx + 1)] <= p_times[index(x_idx, z_idx)])
        local_maxima_p[x_idx].push_back(z_idx);
      if (s_times[index(x_idx, z_idx - 1)] <= s_times[index(x_idx, z_idx)] and
          s_times[index(x_idx, z_idx + 1)] <= s_times[index(x_idx, z_idx)])
        local_maxima_s[x_idx].push_back(z_idx);
      if (p_times[index(x_idx, z_idx - 1)] >= p_times[index(x_idx, z_idx)] and
          p_times[index(x_idx, z_idx + 1)] >= p_times[index(x_idx, z_idx)])
        local_minima_p[x_idx].push_back(z_idx);
      if (s_times[index(x_idx, z_idx - 1)] >= s_times[index(x_idx, z_idx)] and
          s_times[index(x_idx, z_idx + 1)] >= s_times[index(x_idx, z_idx)])
        local_minima_s[x_idx].push_back(z_idx);
    }
  }
}

StationSpecificVelocityModel1D::StationSpecificVelocityModel1D(char *path)
    : model_dir(path) {
  if (!fs::exists(model_dir) || !fs::is_directory(model_dir)) {
    printf("Error: %s does not exist or is not a directory\n",
           model_dir.c_str());
    throw std::runtime_error(
        "Failed to initialize StationSpecificVelocityModel1D: path should be a "
        "directory containing binary travel time tables");
  }
  auto padding_info_file = model_dir / "n_padding";
  if (!fs::exists(padding_info_file) ||
      !fs::is_regular_file(padding_info_file)) {
    throw std::runtime_error(
        "File \"n_padding\" specifying the number of padding nodes is missing. "
        "Please rebuild the model using create_model().");
  }

  auto f = fopen(padding_info_file.string().c_str(), "rb");
  fread(&n_padding, 4, 1, f);
  fclose(f);
}

void StationSpecificVelocityModel1D::add_station(const Station &station) {
  auto model_file_path = model_dir / (station.id + std::string(".pyocto"));
  if (!fs::exists(model_file_path) || !fs::is_regular_file(model_file_path)) {
    std::string err_msg = std::string("Travel time table for station") +
                          station.id + std::string("not found");
    throw std::runtime_error(err_msg);
  }

  // using a pointer instead of a local variable, to avoid double free of
  // ptimes/stimes
  VelocityModel1D *model = new VelocityModel1D(model_file_path.string().data());
  stations.push_back(station);
  models.push_back(model);
  Station station_no_z = Station(station.id, station.x, station.y, 0.0,
                                 station.p_residual, station.s_residual);
  models.back()->add_station(station_no_z);
  models.back()->tolerance = tolerance;
  models.back()->association_cutoff_distance = association_cutoff_distance;
  models.back()->surface_p_velocity = 1e10;
  models.back()->surface_s_velocity = 1e10;
}

bool StationSpecificVelocityModel1D::contains(const Volume &volume,
                                              const Pick *pick) {
  if (pick->station >= stations.size()) {
    printf("Warning: Station %u not found\n", pick->station);
    return false;
  }
  auto &station = stations[pick->station];
  auto *current_model = models[pick->station];

  // Note that we shouldn't use (current_model->stations)[...].z here, because
  // we have set 0 elevation for stations in each VelocityModel1D* member to
  // disable elevation correction
  double z_shift = -station.z + n_padding * current_model->get_delta();
  const Volume volume_new_z =
      Volume(volume.x, volume.y, volume.z + z_shift, volume.t, volume.wx,
             volume.wy, volume.wz, volume.wt);
  // Each per-station model holds a single station, stored at local index 0
  return current_model->contains(volume_new_z, pick, 0);
}

double StationSpecificVelocityModel1D::travel_time(const Volume &volume,
                                                   unsigned int station,
                                                   char phase) {
  if (station >= stations.size()) {
    printf("Warning: Station %u not found\n", station);
    return NAN;
  }

  const auto &station_obj = stations[station];
  auto *current_model = models[station];

  double z_shift = -station_obj.z + n_padding * current_model->get_delta();
  const Volume volume_new_z =
      Volume(volume.x, volume.y, volume.z + z_shift, volume.t, volume.wx,
             volume.wy, volume.wz, volume.wt);
  // Each per-station model holds a single station, stored at local index 0
  return current_model->travel_time_for_station(volume_new_z, 0, phase);
}

} // namespace octoassociator
