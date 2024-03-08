//
// Created by munchmej on 15/06/23.
//

#include "OctoAssociator.h"
#include <algorithm>
#include <cmath>
#include <csignal>
#include <pybind11/pybind11.h>
#include <queue>
#ifndef SINGLE_THREAD
#include <thread>
#endif
#include <vector>

namespace py = pybind11;

namespace octoassociator {

#ifndef SINGLE_THREAD
std::mutex events_mutex;
std::mutex node_mutex;
std::mutex thread_count_mutex;
#endif
volatile sig_atomic_t exitRequested = 0;

unsigned long bin_search(double time, const std::vector<Pick *> &picks) {
  // Returns index of the last pick smaller than or equal to the given time
  unsigned long p0 = 0, p1 = picks.size();
  while (p1 - p0 > 1) {
    unsigned long m = (p0 + p1) / 2;
    if (picks[m]->time > time)
      p1 = m;
    else
      p0 = m;
  }
  return p0;
}

std::vector<Event> OctoAssociator::associate(std::vector<Pick *> &picks) {
  exitRequested = 0; // Reset flag
  std::sort(picks.begin(), picks.end(),
            [](const Pick *a, const Pick *b) { return a->time < b->time; });

  std::set<std::string> missing_stations;
  std::vector<Pick *> filtered_picks;
  for (const auto &pick : picks) {
    if (config->velocity_model_location->stations.find(pick->station) !=
            config->velocity_model_location->stations.end() and
        config->velocity_model_association->stations.find(pick->station) !=
            config->velocity_model_association->stations.end()) {
      filtered_picks.push_back(pick);
    } else {
      missing_stations.insert(pick->station);
    }
  }

  if (not missing_stations.empty()) {
    printf("Found stations with picks but without station information:");
    for (const auto &station : missing_stations)
      printf(" %s", station.c_str());
    printf("\n");
  }

  double min_t = filtered_picks[0]->time;
  double max_t = filtered_picks[filtered_picks.size() - 1]->time;

  PickGrouper grouped_picks = PickGrouper(filtered_picks);

  std::vector<Event> events = {};
  auto active_nodes = std::deque<Node *>();

  while (min_t < max_t + config->time_before) {
    double t0 = min_t - config->time_before;
    double t1 = min_t + config->time_slicing;

    unsigned long p0 = bin_search(t0, filtered_picks),
                  p1 = bin_search(t1, filtered_picks);

    std::vector<Pick *> slice_picks = std::vector<Pick *>(
        filtered_picks.begin() + p0, filtered_picks.begin() + p1 + 1);

    if (slice_picks.size() >= config->n_picks) {
      // printf("Node %f %f %ld\n", t0, t1, slice_picks.size());

      auto *root =
          new Node((config->xlim.first + config->xlim.second) / 2,
                   (config->ylim.first + config->ylim.second) / 2,
                   (config->zlim.first + config->zlim.second) / 2,
                   (t0 + t1) / 2, config->xlim.second - config->xlim.first,
                   config->ylim.second - config->ylim.first,
                   config->zlim.second - config->zlim.first, (t1 - t0),
                   slice_picks, config->velocity_model_association);

      active_nodes.push_back(root);
    }
    min_t += config->time_slicing;
  }

  std::set<Pick *> used_picks;

#ifdef SINGLE_THREAD
  unsigned int thread_count = config->n_threads;
  py::gil_scoped_release release;
  process_tree(events, active_nodes, 0, thread_count);
  py::gil_scoped_acquire acquire;
#else
  std::vector<std::thread> threads;
  unsigned int thread_count = config->n_threads;
  threads.reserve(config->n_threads);
  for (int i = 0; i < config->n_threads; i++) {
    threads.emplace_back(&OctoAssociator::process_tree, this, std::ref(events),
                         std::ref(active_nodes), i, std::ref(thread_count));
  }

  py::gil_scoped_release release;
  while (thread_count > 0) {
    // Avoid busy wait
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    py::gil_scoped_acquire acquire;
    if (PyErr_CheckSignals() != 0) {
      exitRequested = 1;
      break;
    }
    py::gil_scoped_release release2;
  }
  for (auto &thread : threads)
    thread.join();
  py::gil_scoped_acquire acquire;
#endif

  if (exitRequested)
    throw py::error_already_set(); // Throw error if interrupt was caught

  std::sort(events.begin(), events.end());

  std::vector<Event> dedup_events;
  {
    py::gil_scoped_release release_dedup;
    dedup_events = deduplicate_events(events);
    py::gil_scoped_acquire acquire_dedup;
  }

  return dedup_events;
}

std::vector<Event>
OctoAssociator::deduplicate_events(const std::vector<Event> &events) {
  std::map<unsigned int, std::vector<unsigned int>>
      overlaps; // Maps picks to events containing them
  unsigned int event_counter = 0;
  for (const auto &event : events) {
    for (auto pick : event.picks) {
      if (overlaps.find(pick->idx) != overlaps.end()) {
        overlaps[pick->idx].push_back(event_counter);
      } else {
        overlaps[pick->idx] = {event_counter};
      }
    }
    event_counter++;
  }

  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
      intersection_counter;
  for (const auto &pick_events : overlaps) {
    for (int i = 0; i < pick_events.second.size(); i++) {
      for (int j = 0; j < i; j++) {
        std::pair<unsigned int, unsigned int> pair = {pick_events.second[i],
                                                      pick_events.second[j]};
        if (intersection_counter.find(pair) == intersection_counter.end())
          intersection_counter[pair] = 0;
        intersection_counter[pair]++;
      }
    }
  }

  std::set<unsigned int> delete_list;
  for (const auto &conflict : intersection_counter) {
    if (conflict.second <= config->max_pick_overlap)
      continue;
    auto event1 = conflict.first.first, event2 = conflict.first.second;
    if (events[event1].picks.size() > events[event2].picks.size() or
        (events[event1].picks.size() == events[event2].picks.size() and
         events[event1].time < events[event2].time)) {
      delete_list.insert(event2);
    } else {
      delete_list.insert(event1);
    }
  }

  std::vector<Event> new_events;
  int duplicates = 0;
  for (unsigned int i = 0; i < events.size(); i++) {
    if (delete_list.find(i) != delete_list.end()) {
      duplicates++;
      continue;
    }

    const auto &event = events[i];
    if (new_events.empty() or
        event.time - new_events.back().time > config->min_interevent_time) {
      new_events.push_back(event);
    } else {
      duplicates++;
      if (event.picks.size() > new_events.back().picks.size()) {
        new_events.pop_back();
        new_events.push_back(event);
      } // else ignore new events
    }
  }

  // if (duplicates > 0)
  //   printf("Removed %d duplicate events\n", duplicates);

  return new_events;
}

unsigned long OctoAssociator::process_tree(std::vector<Event> &events,
                                           std::deque<Node *> &base_nodes,
                                           int thread_id,
                                           unsigned int &thread_count) {
  // printf("Started thread %d\n", thread_id);
  std::set<Pick *> used_picks;
  std::set<Station, StationCompare> active_stations;
  std::deque<Node *> dfs_nodes;
  std::priority_queue<Node *, std::vector<Node *>, NodeCompare> active_nodes;

  PickGrouper grouped_picks = PickGrouper({});
  PickGroupCache pick_group_cache(config->pick_cache_hash_size);
  unsigned long node_count = 0;
  Node *node;
  while (true) {
    // printf("Thread %d - Loop start\n", thread_id);
#ifdef SINGLE_THREAD
    // As this check is not performed in the main loop anymore, it needs to be
    // done here
    if (node_count % 1000 == 0) {
      py::gil_scoped_acquire acquire;
      if (PyErr_CheckSignals() != 0) {
        exitRequested = 1;
      }
    }
#endif

    if (exitRequested)
      break; // Break with Ctrl-C

    if (not dfs_nodes.empty()) { // DFS always has preference
      node = dfs_nodes.front();
      dfs_nodes.pop_front();
    } else if (active_nodes.empty()) { // Fetch new root note
#ifndef SINGLE_THREAD
      std::lock_guard<std::mutex> node_lock(node_mutex);
#endif
      if (base_nodes.empty())
        break; // Everything is already processed
      node = base_nodes.front();
      base_nodes.pop_front();
      // printf("Thread %d - got new base node\n", thread_id);
      // Reset processing status
      used_picks = {};
      grouped_picks = PickGrouper(node->picks.picks);
      pick_group_cache.reset();
      active_stations = config->velocity_model_association->get_active_stations(
          node->picks.picks);
    } else { // Use local node
      node = active_nodes.top();
      active_nodes.pop();
    }
    // printf("Node size %lu\n", node->picks.n_picks());

    // printf("Thread %d - Got node\n", thread_id);
    ++node_count;
    if (config->node_log_interval > 0 and
        node_count % config->node_log_interval == 0)
      printf("Thread %d - node %lu\n", thread_id, node_count);

    // printf("Processing node
    // %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n%lld\n",
    //        node->x, node->y, node->z, node->t, node->wx, node->wy, node->wz,
    //        node->wt, (long long int) node);

    // Filter node only to picks that have not been used in events so far
    filter_used_picks(node, used_picks);
    // Filter picks failing the distance criterion
    filter_distance_criterion(node, active_stations);

    // if (pick_group_cache.contains_group(node->picks.picks)) continue;

    // printf("Thread %d - Checking node status\n", thread_id);
    switch (node_status(node)) {
    case NODE_CONTINUE:
      node->split2();
      for (auto *child : node->children) {
        int target_queue = queue_choice(node, active_nodes.size());
        if (target_queue == QUEUE_PRIORITY) {
          active_nodes.push(child);
        } else {
          dfs_nodes.push_front(child);
        }
      }
      break;
    case NODE_BREAK_EVENT: {
      create_event(events, used_picks, grouped_picks, pick_group_cache,
                   active_stations, node);
      break;
    }
      // case NODE_BREAK_NO_EVENT: noop
    }
    // printf("Thread %d - Trying to delete node\n", thread_id);
    delete node;
    // printf("Thread %d - Loop complete\n", thread_id);
  }

  {
#ifndef SINGLE_THREAD
    std::lock_guard<std::mutex> thread_count_lock(thread_count_mutex);
#endif
    thread_count--;
    // printf("Thread %d exiting. New count %u.\n", thread_id, thread_count);
  }

  return node_count;
}

bool OctoAssociator::create_event(
    std::vector<Event> &events, std::set<Pick *> &used_picks,
    PickGrouper &grouped_picks, PickGroupCache &pick_group_cache,
    const std::set<Station, StationCompare> &active_stations, Node *node) {
  // printf("Creating event (%lu picks)\n", node->picks.n_picks());
  if (pick_group_cache.contains_group(node->picks.picks))
    return false;
  pick_group_cache.insert_group(node->picks.picks);
  auto root_volume = config->root_volume();
  auto volume = locate_edt(
      root_volume, node->picks.picks); // Search in full volume for source

  std::vector<Pick *> new_picks;
  for (int i = 0; i < config->refinement_iterations; i++) {
    new_picks = {};
    for (const auto &station_pair : grouped_picks.pick_groups) {
      for (char phase : {'P', 'S'}) {
        auto tt = config->velocity_model_location->travel_time(
            volume, station_pair.first, phase);
        if (std::isnan(tt))
          continue;

        auto pick =
            grouped_picks.pick_at_time(station_pair.first, phase, volume.t + tt,
                                       config->pick_match_tolerance);
        if (pick != nullptr and used_picks.find(pick) == used_picks.end()) {
          new_picks.push_back(pick);
        }
      }
    }

    // Filter picks not matching the distance criterion
    new_picks = filter_distance_criterion(volume, new_picks, active_stations);

    // Event is not valid anymore
    if (not picks_sufficient(PickContainer(new_picks)))
      return false;

    root_volume = config->root_volume();
    volume = locate_edt(root_volume, new_picks);

    if (std::isnan(volume.t))
      // No location found
      return false;
  }

  std::vector<double> residuals;
  for (Pick *pick : new_picks) {
    double tt = config->velocity_model_location->travel_time(volume, pick);
    residuals.push_back(pick->time - (volume.t + tt));
  }

  {
#ifndef SINGLE_THREAD
    std::lock_guard<std::mutex> event_lock(events_mutex);
#endif
    events.emplace_back(volume, new_picks, residuals);
  }

  for (auto pick : new_picks)
    used_picks.insert(pick);
  // printf("Created event (%lu picks, %lu events)\n", node->picks.n_picks(),
  // events.size());
  return true;
}

int OctoAssociator::node_status(const Node *node) {
  if (not picks_sufficient(node->picks)) {
    return NODE_BREAK_NO_EVENT;
  } else if (node->maxw() <= config->min_node_size) {
    return NODE_BREAK_EVENT;
  } else {
    return NODE_CONTINUE;
  }
}

int OctoAssociator::queue_choice(const Node *node,
                                 unsigned int priority_queue_size) {
  if (node->maxw() > config->dfs_node_size and
      priority_queue_size < config->queue_memory_protection_dfs_size) {
    return QUEUE_PRIORITY;
  } else {
    return QUEUE_DFS;
  }
}

bool OctoAssociator::picks_sufficient(const PickContainer &picks) {
  return (picks.n_picks() >= config->n_picks and
          picks.n_p_stations() >= config->n_p_picks and
          picks.n_s_stations() >= config->n_s_picks and
          picks.n_p_and_s_stations() >= config->n_p_and_s_picks);
}

Volume OctoAssociator::locate_edt(Volume &volume,
                                  const std::vector<Pick *> &picks) {
  /* printf("Volume %.2f %.2f %.2f %.2f %.2f %.2f\n",
         volume.x, volume.y, volume.z,
         volume.wx, volume.wy, volume.wz); */
  if (volume.maxw_wo_time() < config->min_node_size_location) {
    edt_loss(volume, picks, false); // Calculate time for volume
    return volume;
  } else {
    auto children_groups = volume.grouped_multisplit(
        config->location_split_depth, config->location_split_return,
        config->min_node_size_location);
    unsigned int min_i = 0;
    double min_val = INFINITY;
    for (unsigned int i = 0; i < children_groups.size(); i++) {
      for (auto &child : children_groups[i].second) {
        double volume_loss = edt_loss(child, picks, false);
        if (volume_loss < min_val) {
          min_i = i;
          min_val = volume_loss;
        }
      }
    }
    return locate_edt(children_groups[min_i].first, picks);
    /*if (edt_loss(children.first, picks, false) <
        edt_loss(children.second, picks, false)) {
      return locate_edt(children.first, picks);
    } else {
      return locate_edt(children.second, picks);
    }*/
  }
}

double OctoAssociator::edt_loss(Volume &volume,
                                const std::vector<Pick *> &picks,
                                bool print_residuals) {
  std::vector<double> tts;
  std::vector<double> arrivals;
  double origin_time = 0;

  for (auto pick : picks) {
    auto tt_tmp = config->velocity_model_location->travel_time(volume, pick);
    if (std::isnan(tt_tmp))
      continue;

    tts.push_back(tt_tmp);
    arrivals.push_back(pick->time);
    origin_time += (arrivals.back() - tts.back());
  }
  origin_time /= arrivals.size();
  volume.t = origin_time;

  if (arrivals.empty())
    // Not a single arrival within range
    return INFINITY;

  if (print_residuals) {
    // printf("Event - Number of residuals: %lu\n", picks.size());
    for (auto pick : picks) {
      auto tt = config->velocity_model_location->travel_time(volume, pick);
      auto diff = pick->time - (origin_time + tt);
      // printf("%s %c %f\n", pick->station.c_str(), pick->phase, diff);
    }
  }

  double edt = 0;
  double edt_exp = 0;
  double double_sigma_squared = 2 * config->edt_pick_std * config->edt_pick_std;
  for (int i = 0; i < arrivals.size(); i++) {
    for (int j = 0; j < i; j++) {
      edt += pow((tts[i] - tts[j]) - (arrivals[i] - arrivals[j]), 2);
      edt_exp += exp(-pow((tts[i] - tts[j]) - (arrivals[i] - arrivals[j]), 2) /
                     double_sigma_squared);
    }
  }
  edt_exp = -edt_exp; // Flip sign so target is minimizing (even though it's
                      // actually a density)
  edt_exp /=
      (arrivals.size() * (arrivals.size() - 1)) / 2; // Norm by number of picks
  edt /=
      (arrivals.size() * (arrivals.size() - 1)) / 2; // Norm by number of picks
  if (config->exponential_edt) {
    return edt_exp;
  } else {
    return edt;
  }
}

struct DistanceQueueEntry {
  DistanceQueueEntry(const Station &station, double dist, bool has_pick)
      : station_id(station.id), dist(dist), has_pick(has_pick){};
  std::string station_id;
  double dist;
  bool has_pick;
};

bool operator<(const DistanceQueueEntry &a, const DistanceQueueEntry &b) {
  return a.dist > b.dist; // Invert to get ascending order
}

std::vector<Pick *> OctoAssociator::filter_distance_criterion(
    const Volume &volume, const std::vector<Pick *> &picks,
    const std::set<Station, StationCompare> &active_stations) {
  if (config->min_pick_fraction == 0)
    return picks; // Skip filtering calculations
  auto node_stations =
      config->velocity_model_association->get_active_stations(picks);

  // Criterion is met globally, so it does not make sense to filter
  if ((double)node_stations.size() / (double)active_stations.size() >
      config->min_pick_fraction)
    return picks;

  // Priority queue with <station_id, dist, has_pick> sorted by distance
  std::priority_queue<DistanceQueueEntry> entries;

  for (const auto &station : active_stations) {
    auto min_max_dist = min_max_hypo_dist(volume, station);
    if (node_stations.find(station) ==
        node_stations.end()) { // Station has no picks
      entries.emplace(station, min_max_dist.second, false);
    } else { // Station has picks
      entries.emplace(station, min_max_dist.first, true);
    }
  }

  std::set<std::string> delete_stations;
  double pick_count = 0, total_count = 0;
  while (not entries.empty()) {
    auto entry = entries.top();
    entries.pop();
    total_count++;
    if (entry.has_pick) {
      pick_count++;
      if (pick_count / total_count < config->min_pick_fraction)
        delete_stations.insert(entry.station_id);
    }
    if (pick_count / total_count >= config->min_pick_fraction)
      delete_stations.clear();
  }

  std::vector<Pick *> new_picks;
  for (auto pick : picks) {
    if (delete_stations.find(pick->station) == delete_stations.end()) {
      new_picks.push_back(pick);
    }
  }
  return new_picks;
}

void OctoAssociator::filter_distance_criterion(
    Node *node, const std::set<Station, StationCompare> &active_stations) {
  auto new_picks = filter_distance_criterion(node->volume, node->picks.picks,
                                             active_stations);
  node->picks = PickContainer(new_picks);
}

void OctoAssociator::filter_used_picks(
    Node *node, const std::set<Pick *> &used_picks) const {
  PickContainer new_picks;
  for (auto pick : node->picks.picks) {
    if (used_picks.find(pick) == used_picks.end()) {
      new_picks.push_back(pick);
    }
  }
  node->picks = new_picks;
}

std::pair<double, double>
OctoAssociator::min_max_hypo_dist(const Volume &volume,
                                  const Station &station) {
  auto x_dist = config->velocity_model_association->min_max_dist(
      volume.x, volume.wx, station.x);
  auto y_dist = config->velocity_model_association->min_max_dist(
      volume.y, volume.wy, station.y);

  double min_dist = sqrt(pow(x_dist.first, 2) + pow(y_dist.first, 2));
  double max_dist = sqrt(pow(x_dist.second, 2) + pow(y_dist.second, 2));

  return {min_dist, max_dist};
}

void Node::split16() {
  for (double ox : {-1., 1.}) {
    for (double oy : {-1., 1.}) {
      for (double oz : {-1., 1.}) {
        for (double ot : {-1., 1.}) {
          auto *node = new Node(x + ox * wx / 4, y + oy * wy / 4,
                                z + oz * wz / 4, t + ot * wt / 4, wx / 2,
                                wy / 2, wz / 2, wt / 2, velocity_model);
          for (const auto pick : picks.picks) {
            if (velocity_model->contains(volume, pick))
              node->picks.push_back(pick);
          }
          children.push_back(node);
        }
      }
    }
  }
}

void Node::split2() {
  Node *node1, *node2;
  if (maxw() == wx) {
    node1 = new Node(x - wx / 4, y, z, t, wx / 2, wy, wz, wt, velocity_model);
    node2 = new Node(x + wx / 4, y, z, t, wx / 2, wy, wz, wt, velocity_model);
  } else if (maxw() == wy) {
    node1 = new Node(x, y - wy / 4, z, t, wx, wy / 2, wz, wt, velocity_model);
    node2 = new Node(x, y + wy / 4, z, t, wx, wy / 2, wz, wt, velocity_model);
  } else if (maxw() == wz) {
    node1 = new Node(x, y, z - wz / 4, t, wx, wy, wz / 2, wt, velocity_model);
    node2 = new Node(x, y, z + wz / 4, t, wx, wy, wz / 2, wt, velocity_model);
  } else {
    node1 = new Node(x, y, z, t - wt / 4, wx, wy, wz, wt / 2, velocity_model);
    node2 = new Node(x, y, z, t + wt / 4, wx, wy, wz, wt / 2, velocity_model);
  }

  for (const auto pick : picks.picks) {
    if (velocity_model->contains(node1->volume, pick))
      node1->picks.push_back(pick);
    if (velocity_model->contains(node2->volume, pick))
      node2->picks.push_back(pick);
  }
  push_nodes(node1, node2);
}

void Node::push_nodes(Node *node1, Node *node2) {
  // This function is unfortunately broken. It is possible that node1 and node2
  // have the same picks, node1 has the lower residual but not actually a
  // hypocenter supporting all picks while node 2 has one.

  // Push in the order that the node with more picks is processed first
  // Note that the order is reversed afterwards by the push front call
  if (node1->picks.picks.size() > node2->picks.picks.size()) {
    children.push_back(node2);
    children.push_back(node1);
  } else {
    children.push_back(node1);
    children.push_back(node2);
  }
}

void PickContainer::push_back(Pick *pick) {
  if (pick->phase == 'P') {
    p_stations.insert(pick->station);
  } else if (pick->phase == 'S') {
    s_stations.insert(pick->station);
  } else {
    printf("Unknown phase '%c'\n", pick->phase);
  }
  picks.push_back(pick);
}

unsigned long PickContainer::n_p_and_s_stations() const {
  std::set<std::string> intersect;
  set_intersection(p_stations.begin(), p_stations.end(), s_stations.begin(),
                   s_stations.end(),
                   std::inserter(intersect, intersect.begin()));
  return intersect.size();
}

PickContainer::PickContainer(const std::vector<Pick *> &picks) {
  for (auto *pick : picks) {
    push_back(pick);
  }
}

Volume OctoTreeConfig::root_volume() const {
  return {(xlim.first + xlim.second) / 2, (ylim.first + ylim.second) / 2,
          (zlim.first + zlim.second) / 2, 0,
          xlim.second - xlim.first,       ylim.second - ylim.first,
          zlim.second - zlim.first,       0};
}

unsigned int
PickGroupCache::hash_picks(const std::vector<unsigned int> &pick_idx) const {
  // Use some primes as factors
  std::vector<unsigned int> factors = {2,  3,  5,  7,  11, 13, 17, 19, 23,
                                       29, 31, 37, 41, 43, 47, 53, 59, 61,
                                       67, 71, 73, 79, 83, 89, 97};

  unsigned int hash = 0;
  for (unsigned int i = 0; i < pick_idx.size(); i++) {
    hash += pick_idx[i] * factors[i % factors.size()];
  }
  return hash % pick_hash_size;
}

bool PickGroupCache::contains_group(const std::vector<Pick *> &picks) {
  std::vector<unsigned int> pick_idx = picks_to_idx(picks);
  auto hash = hash_picks(pick_idx);

  if (cache.find(hash) == cache.end())
    return false;

  return not(cache[hash] == pick_idx);
}

std::vector<unsigned int>
PickGroupCache::picks_to_idx(const std::vector<Pick *> &picks) const {
  std::vector<unsigned int> pick_idx;
  pick_idx.reserve(picks.size());
  for (auto pick : picks)
    pick_idx.push_back(pick->idx);
  std::sort(pick_idx.begin(), pick_idx.end());
  return pick_idx;
}

void PickGroupCache::insert_group(const std::vector<Pick *> &picks) {
  std::vector<unsigned int> pick_idx = picks_to_idx(picks);
  auto hash = hash_picks(pick_idx);

  cache[hash] = pick_idx;
}

PickGrouper::PickGrouper(const std::vector<Pick *> &picks) {
  for (Pick *pick : picks) {
    auto group = pick_groups.find(pick->station);
    if (group == pick_groups.end()) {
      pick_groups[pick->station] = {pick};
    } else {
      group->second.push_back(pick);
    }
  }

  for (auto &group : pick_groups) {
    std::sort(group.second.begin(), group.second.end(),
              [](const Pick *a, const Pick *b) { return a->time < b->time; });
  }
}

Pick *PickGrouper::pick_at_time(const std::string &station, char phase,
                                double time, double tolerance) {
  // Returns the chronologically first pick in the interval (time - tolerance,
  // time + tolerance) that fits the phase
  auto group = pick_groups.find(station);
  if (group == pick_groups.end())
    return nullptr;

  const auto &group_picks = group->second;
  auto p0 = bin_search(time + tolerance, group_picks);

  if (group_picks[p0]->time < time - tolerance)
    return nullptr; // Before first arrival

  // Move p0 to the left until the time is too early
  while (p0 > 0 and group_picks[p0 - 1]->time > time - tolerance) {
    p0--;
    if (p0 == 0)
      break;
  }

  // Mote to the right until the end of the sentence, exiting the tolerance or
  // the first match
  for (; p0 < group_picks.size(); p0++) {
    if (group_picks[p0]->time >= time + tolerance)
      break;
    if (group_picks[p0]->phase == phase)
      return group_picks[p0];
  }
  return nullptr;
}

} // namespace octoassociator
