//
// Created by munchmej on 15/06/23.
//

#ifndef OCTO_ASSOCIATOR_TYPES_H
#define OCTO_ASSOCIATOR_TYPES_H

#include <algorithm>
#include <ciso646> // For windows compatibility
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>

namespace octoassociator {

struct Pick {
  unsigned int idx;
  double time;
  std::string station;
  char phase;
  Pick(int idx, double time, std::string station, char phase)
      : idx(idx), time(time), station(station), phase(phase){};

  bool operator<(const Pick &b) const {
    return time < b.time or (time == b.time and station < b.station) or
           (time == b.time and station == b.station and phase < b.phase);
  };
  char *to_string(char *output) const {
    sprintf(output, "%s\t%c\t%f", station.c_str(), phase, time);
    return output;
  };
};

struct Station {
  std::string id;
  double x, y, z;
  // Additive, i.e., corrected times are tt_pred + correction.
  // Same convention as NonLinLoc
  double p_residual = 0., s_residual = 0.;
  Station() : x(0), y(0), z(0){};
  Station(std::string id, double x, double y, double z, double p_residual,
          double s_residual)
      : id(id), x(x), y(y), z(z), p_residual(p_residual),
        s_residual(s_residual){};
  double residual(char phase) const {
    if (phase == 'P')
      return p_residual;
    return s_residual;
  };
};

class StationCompare {
public:
  bool operator()(const Station &a, const Station &b) const {
    return a.id < b.id;
  }
};

class Volume {
  // Example: element goes from x - wx / 2 to x + wx / 2
public:
  double x, y, z, t;     // Location in space and time
  double wx, wy, wz, wt; // Width in each axis
  double maxw() const { return std::max({wx, wy, wz, wt * 5}); };
  double maxw_wo_time() const { return std::max({wx, wy, wz}); };
  std::pair<Volume, Volume> split2() const;
  std::vector<Volume> split4() const;
  std::vector<Volume> multisplit(unsigned int splits, double min_size);
  std::vector<std::pair<Volume, std::vector<Volume>>>
  grouped_multisplit(unsigned int splits, unsigned int grouping_depth,
                     double min_size);

  double x_min() const { return x - wx / 2; };
  double x_max() const { return x + wx / 2; };
  double y_min() const { return y - wy / 2; };
  double y_max() const { return y + wy / 2; };
  double z_min() const { return z - wz / 2; };
  double z_max() const { return z + wz / 2; };
  double t_min() const { return t - wt / 2; };
  double t_max() const { return t + wt / 2; };

  Volume(double x, double y, double z, double t, double wx, double wy,
         double wz, double wt)
      : x(x), y(y), z(z), t(t), wx(wx), wy(wy), wz(wz), wt(wt){};
  explicit Volume(Volume *source)
      : x(source->x), y(source->y), z(source->z), t(source->t), wx(source->wx),
        wy(source->wy), wz(source->wz), wt(source->wt){};
};

struct Event {
  double time;
  double x;
  double y;
  double z;
  std::vector<Pick *> picks;
  std::vector<double> residuals;

  bool operator<(const Event &b) const {
    return time < b.time or (time == b.time and x < b.x) or
           (time == b.time and x == b.x and y < b.y) or
           (time == b.time and x == b.x and y == b.y and z < b.z);
  };

  Event() = delete;
  Event(Volume volume, const std::vector<Pick *> &picks,
        const std::vector<double> &residuals)
      : x(volume.x), y(volume.y), z(volume.z), time(volume.t), picks(picks),
        residuals(residuals){};
  Event(double time, double x, double y, double z,
        const std::vector<Pick *> &picks, const std::vector<double> &residuals)
      : time(time), x(x), y(y), z(z), picks(picks), residuals(residuals){};
};

} // namespace octoassociator

#endif // OCTO_ASSOCIATOR_TYPES_H
