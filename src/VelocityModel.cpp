//
// Created by munchmej on 17/06/23.
//

#include <cmath>
#include <set>

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

std::set<Station, StationCompare>
VelocityModel::get_active_stations(const std::vector<Pick *> &picks) {
  std::set<Station, StationCompare> active_stations;
  for (const auto &pick : picks) {
    if (stations.find(pick->station) != stations.end()) {
      auto station = stations[pick->station];
      active_stations.insert(station);
    }
  }
  return active_stations;
}

bool VelocityModel0D::contains(const Volume &volume, const Pick *pick) {
  auto station_it = stations.find(pick->station);
  if (station_it == stations.end()) {
    printf("Warning: Station %s not found\n", pick->station.c_str());
    return false;
  }
  auto &station = station_it->second;
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

double VelocityModel0D::travel_time(const Volume &volume,
                                    const std::string &station, char phase) {
  double velocity = get_phase_velocity(phase);
  auto station_it = stations.find(station);
  if (station_it == stations.end()) {
    printf("Warning: Station %s not found\n", station.c_str());
    return NAN;
  }
  auto &station_obj = station_it->second;

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

double VelocityModel1D::travel_time(const Volume &volume,
                                    const std::string &station, char phase) {
  auto station_it = stations.find(station);
  if (station_it == stations.end()) {
    printf("Warning: Station %s not found\n", station.c_str());
    return NAN;
  }
  const auto &station_obj = station_it->second;

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
  // The function is only exact up to grid points
  // Calculate min-max distance 2D
  auto station_it = stations.find(pick->station);
  if (station_it == stations.end()) {
    printf("Warning: Station %s not found\n", pick->station.c_str());
    return false;
  }
  auto &station = station_it->second;

  auto x_dist = min_max_dist(volume.x, volume.wx, station.x);
  auto y_dist = min_max_dist(volume.y, volume.wy, station.y);

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

  for (auto z : local_minima) {
    if (volume.z_min() <= z and z <= volume.z_max())
      min_candidates.push_back(times[index(ix_min, z)]);
  }
  for (auto z : local_maxima) {
    if (volume.z_min() <= z and z <= volume.z_max())
      max_candidates.push_back(times[index(ix_max, z)]);
  }

  auto min_tt = *std::min_element(min_candidates.begin(), min_candidates.end());
  auto max_tt = *std::max_element(max_candidates.begin(), max_candidates.end());

  double elev_correction = elevation_correction(pick->phase, station);

  auto t0 = volume.t_min() + min_tt + elev_correction +
            station.residual(pick->phase) - tolerance;
  auto t1 = volume.t_max() + max_tt + elev_correction +
            station.residual(pick->phase) + tolerance;

  return (t0 < pick->time and pick->time < t1);
}

void VelocityModel1D::precalculate_extrema() {
  local_maxima_p.resize(nx);
  local_minima_p.resize(nx);
  local_maxima_s.resize(nx);
  local_minima_s.resize(nx);
  for (int x = 0; x < nx; ++x) {
    for (int z = 1; z < nz - 1; ++z) {
      if (p_times[index(x, z - 1)] <= p_times[index(x, z)] and
          p_times[index(x, z + 1)] <= p_times[index(x, z)])
        local_maxima_p[x].push_back(z);
      if (s_times[index(x, z - 1)] <= s_times[index(x, z)] and
          s_times[index(x, z + 1)] <= s_times[index(x, z)])
        local_maxima_s[x].push_back(z);
      if (p_times[index(x, z - 1)] >= p_times[index(x, z)] and
          p_times[index(x, z + 1)] >= p_times[index(x, z)])
        local_minima_p[x].push_back(z);
      if (s_times[index(x, z - 1)] >= s_times[index(x, z)] and
          s_times[index(x, z + 1)] >= s_times[index(x, z)])
        local_minima_s[x].push_back(z);
    }
  }
}

} // namespace octoassociator
