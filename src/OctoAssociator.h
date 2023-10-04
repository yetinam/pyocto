//
// Created by munchmej on 15/06/23.
//

#ifndef OCTO_ASSOCIATOR_OCTOASSOCIATOR_H
#define OCTO_ASSOCIATOR_OCTOASSOCIATOR_H

#define NODE_CONTINUE 0
#define NODE_BREAK_NO_EVENT 1
#define NODE_BREAK_EVENT 2

#define QUEUE_DFS 0
#define QUEUE_PRIORITY 1

// All Windows version are caught with this
// As Windows does not have pthreads, currently only a single thread
// implementation is available.
#ifdef _WIN32
#define SINGLE_THREAD
#endif

#include <algorithm>
#include <ciso646> // For windows compatibility
#include <mutex>
#include <queue>
#include <set>
#include <stdlib.h>
#include <string>
#ifndef SINGLE_THREAD
#include <thread>
#endif
#include <vector>

#include "VelocityModel.h"
#include "types.h"

namespace octoassociator {

struct OctoTreeConfig {
public:
  OctoTreeConfig() = default;
  explicit OctoTreeConfig(VelocityModel *velocity_model_association,
                          VelocityModel *velocity_model_location)
      : xlim(25, 700), ylim(6300, 7450), zlim(0, 250),
        velocity_model_association(velocity_model_association),
        velocity_model_location(velocity_model_location){};
  double dfs_node_size = 0.; // Trigger DFS for nodes smaller than this limit -
                             // Not publicly available
  double min_node_size = 10.;
  double min_node_size_location = 1.5;
  double pick_match_tolerance = 1.5;
  double min_interevent_time = 3.0; // In seconds
  bool exponential_edt = false;     // If false, skip exponential step
  double edt_pick_std = 1.0;        // Standard deviation of picks for EDT loss
  int max_pick_overlap =
      4; // Maximum number of picks allowed to overlap between two events
  unsigned int n_picks = 10;
  unsigned int n_p_picks = 5;
  unsigned int n_s_picks = 5;
  unsigned int n_p_and_s_picks = 3;
#ifdef SINGLE_THREAD
  unsigned int n_threads = 1;
#else
  unsigned int n_threads = std::thread::hardware_concurrency();
#endif
  unsigned int refinement_iterations = 3;
  unsigned int node_log_interval = 0;
  unsigned int location_split_depth = 6;
  unsigned int location_split_return = 4;
  double min_pick_fraction = 0.25;
  unsigned int pick_cache_hash_size = 10000;
  unsigned int queue_memory_protection_dfs_size =
      500000; // Start DFS if queue gets too big to avoid memory overflow
  double time_before = 300; // Maximum travel time from event to station
  double time_slicing = 1200;
  std::pair<double, double> xlim = {0., 0.}, ylim = {0., 0.}, zlim = {0., 0.};
  VelocityModel *velocity_model_association = nullptr;
  VelocityModel *velocity_model_location = nullptr;

  Volume root_volume() const;
};

class PickContainer {
  // A pick container that counts statistics about the picks
public:
  std::vector<Pick *> picks;
  unsigned long n_picks() const { return picks.size(); };
  unsigned long n_p_stations() const { return p_stations.size(); };
  unsigned long n_s_stations() const { return s_stations.size(); };
  unsigned long n_p_and_s_stations() const;
  void push_back(Pick *pick);
  PickContainer() = default;
  explicit PickContainer(const std::vector<Pick *> &picks);

private:
  std::set<std::string> p_stations;
  std::set<std::string> s_stations;
};

class PickGrouper {
public:
  std::map<std::string, std::vector<Pick *>> pick_groups;
  explicit PickGrouper(const std::vector<Pick *> &picks);
  Pick *pick_at_time(const std::string &station, char phase, double time,
                     double tolerance);
};

class PickGroupCache {
public:
  explicit PickGroupCache(unsigned int pick_hash_size)
      : pick_hash_size(pick_hash_size){};
  bool contains_group(const std::vector<Pick *> &picks);
  void insert_group(const std::vector<Pick *> &picks);
  void reset() { cache = {}; };

private:
  unsigned int hash_picks(const std::vector<unsigned int> &picks)
      const; // Calculate a hash for a set of picks
  std::map<unsigned int, std::vector<unsigned int>> cache;
  unsigned int pick_hash_size;
  std::vector<unsigned int>
  picks_to_idx(const std::vector<Pick *> &picks) const;
};

class Node {
public:
  Volume volume;
  double &x = volume.x, &y = volume.y, &z = volume.z, &t = volume.t;
  double &wx = volume.wx, &wy = volume.wy, &wz = volume.wz, &wt = volume.wt;
  PickContainer picks;
  VelocityModel *velocity_model;

  Node(double x, double y, double z, double t, double wx, double wy, double wz,
       double wt, VelocityModel *velocity_model)
      : volume(x, y, z, t, wx, wy, wz, wt), velocity_model(velocity_model){};
  Node(double x, double y, double z, double t, double wx, double wy, double wz,
       double wt, std::vector<Pick *> &picks, VelocityModel *velocity_model)
      : volume(x, y, z, t, wx, wy, wz, wt), picks(picks),
        velocity_model(velocity_model){};

  void split16();
  void split2();
  void push_nodes(Node *node1, Node *node2);
  double maxw() const { return volume.maxw(); };

  std::vector<Node *> children;
};

class NodeCompare {
public:
  bool operator()(Node *a, Node *b) {
    return a->picks.n_picks() < b->picks.n_picks();
  }
};

class OctoAssociator {
public:
  explicit OctoAssociator(OctoTreeConfig *config) : config(config){};
  std::vector<Event> associate(std::vector<Pick *> &picks);

private:
  OctoTreeConfig *config;
  void filter_distance_criterion(
      Node *node, const std::set<Station, StationCompare> &active_stations);
  std::vector<Pick *> filter_distance_criterion(
      const Volume &volume, const std::vector<Pick *> &picks,
      const std::set<Station, StationCompare> &active_stations);
  void filter_used_picks(Node *node, const std::set<Pick *> &used_picks) const;
  std::vector<Event> deduplicate_events(const std::vector<Event> &events);
  int node_status(const Node *node);
  int queue_choice(const Node *node, unsigned int priority_queue_size);
  bool picks_sufficient(const PickContainer &picks);
  unsigned long process_tree(std::vector<Event> &events,
                             std::deque<Node *> &base_nodes, int thread_id,
                             unsigned int &thread_count);
  Volume locate_edt(Volume &volume, const std::vector<Pick *> &picks);
  double edt_loss(Volume &volume, const std::vector<Pick *> &picks,
                  bool print_residuals);
  // Returns true if an event has been created
  bool create_event(std::vector<Event> &events, std::set<Pick *> &used_picks,
                    PickGrouper &grouped_picks,
                    PickGroupCache &pick_group_cache,
                    const std::set<Station, StationCompare> &active_stations,
                    Node *node);
  std::pair<double, double> min_max_hypo_dist(const Volume &volume,
                                              const Station &station);
};
} // namespace octoassociator

#endif // OCTO_ASSOCIATOR_OCTOASSOCIATOR_H
