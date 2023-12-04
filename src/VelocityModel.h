//
// Created by munchmej on 17/06/23.
//

#ifndef OCTO_ASSOCIATOR_SRC_VELOCITYMODEL_H_
#define OCTO_ASSOCIATOR_SRC_VELOCITYMODEL_H_

#include <ciso646> // For windows compatibility
#include <map>
#include <stdlib.h>
#include <utility>
#include <vector>

#include "types.h"

namespace octoassociator {

class VelocityModel {
public:
  virtual bool contains(const Volume &volume, const Pick *pick) = 0;
  virtual double travel_time(const Volume &volume, const std::string &station,
                             char phase) = 0;
  double travel_time(const Volume &volume, const Pick *pick) {
    return travel_time(volume, pick->station, pick->phase);
  }; // Travel time from centroid to pick
  void add_station(const Station &station) { stations[station.id] = station; }
  std::map<std::string, Station> stations;
  std::set<Station, StationCompare>
  get_active_stations(const std::vector<Pick *> &picks);
  /**
   * Determines the minimum and maximum distance of a point px to the interval
   * [x - wx / 2, x + wx / 2]
   */
  static std::pair<double, double> min_max_dist(double x, double wx, double px);
};

class VelocityModel0D : public VelocityModel {
public:
  double p_velocity;
  double s_velocity;
  double tolerance = 1; // Tolerance in seconds
  double association_cutoff_distance =
      1e9; // Ignore all travel times outside this range in km from contains
           // queries
  double location_cutoff_distance = 1e9; // Ignore all travel times outside this
                                         // range in km for travel time queries

  VelocityModel0D(double p_velocity, double s_velocity, double tolerance,
                  double association_cutoff_distance,
                  double location_cutoff_distance)
      : p_velocity(p_velocity), s_velocity(s_velocity), tolerance(tolerance),
        association_cutoff_distance(association_cutoff_distance),
        location_cutoff_distance(location_cutoff_distance){};
  bool contains(const Volume &volume, const Pick *pick) override;
  double travel_time(const Volume &volume, const std::string &station,
                     char phase) override;
  void distance_volume_pick(const Volume &volume, double &min_dist,
                            double &max_dist, double &centroid_dist,
                            const Station &station);
  double get_phase_velocity(char phase) const;
};

class VelocityModel1D : public VelocityModel {
public:
  double tolerance = 1.;
  // Both distances hypocentral
  double association_cutoff_distance =
      1e9; // Ignore all travel times outside this range in km from contains
           // queries
  double location_cutoff_distance = 1e9; // Ignore all travel times outside this
                                         // range in km for travel time queries

  double surface_p_velocity = 5.3;
  double surface_s_velocity = 3.1;
  explicit VelocityModel1D(char *path);
  ~VelocityModel1D() {
    delete[] p_times;
    delete[] s_times;
  }
  bool contains(const Volume &volume, const Pick *pick) override;
  double travel_time(const Volume &volume, const std::string &station,
                     char phase) override;
  double travel_time(double x, double z, char phase);

private:
  int nx{}, nz{};  // Number of nodes in each direction
  double delta{};  // Spacing between nodes
  double *p_times; // 2D array of travel times
  double *s_times; // 2D array of travel times

  double elevation_correction(char phase, const Station &station) const;
  int index(int x, int z) const {
    return x * nz + z;
  }; // Helper function for indexing dynamic 2D array

  void precalculate_extrema();
  std::vector<std::vector<int>> local_minima_p;
  std::vector<std::vector<int>> local_minima_s;
  std::vector<std::vector<int>> local_maxima_p;
  std::vector<std::vector<int>> local_maxima_s;
};

} // namespace octoassociator

#endif // OCTO_ASSOCIATOR_SRC_VELOCITYMODEL_H_
