#include "types.h"

namespace octoassociator {

std::pair<Volume, Volume> Volume::split2() const {
  if (maxw_wo_time() == wx) {
    return {Volume(x - wx / 4, y, z, t, wx / 2, wy, wz, wt),
            Volume(x + wx / 4, y, z, t, wx / 2, wy, wz, wt)};
  } else if (maxw_wo_time() == wy) {
    return {Volume(x, y - wy / 4, z, t, wx, wy / 2, wz, wt),
            Volume(x, y + wy / 4, z, t, wx, wy / 2, wz, wt)};
  } else {
    return {Volume(x, y, z - wz / 4, t, wx, wy, wz / 2, wt),
            Volume(x, y, z + wz / 4, t, wx, wy, wz / 2, wt)};
  }
}

std::vector<Volume> Volume::split4() const {
  auto sub = split2();
  auto parts1 = sub.first.split2();
  auto parts2 = sub.second.split2();
  return {parts1.first, parts1.second, parts2.first, parts2.second};
}

std::vector<Volume> Volume::multisplit(unsigned int splits, double min_size) {
  if (maxw_wo_time() < min_size) {
    std::vector<Volume> result;
    result.push_back(*this);
    return result;
  }
  auto base_split = split2();
  if (splits == 1) {
    return {base_split.first, base_split.second};
  } else {
    auto part1 = base_split.first.multisplit(splits - 1, min_size);
    auto part2 = base_split.second.multisplit(splits - 1, min_size);

    part1.insert(part1.end(), part2.begin(), part2.end());
    return part1;
  }
}

std::vector<std::pair<Volume, std::vector<Volume>>>
Volume::grouped_multisplit(unsigned int splits, unsigned int grouping_depth,
                           double min_size) {
  std::vector<std::pair<Volume, std::vector<Volume>>> output;
  if (splits == grouping_depth) {
    output.emplace_back(Volume(this), multisplit(splits, min_size));
  } else {
    if (maxw_wo_time() < min_size) {
      std::vector<Volume> children = {Volume(this)};
      output.emplace_back(Volume(this), children);
    } else {
      auto base_split = split2();
      auto part1 = base_split.first.grouped_multisplit(
          splits - 1, grouping_depth, min_size);
      auto part2 = base_split.second.grouped_multisplit(
          splits - 1, grouping_depth, min_size);

      part1.insert(part1.end(), part2.begin(), part2.end());
      return part1;
    }
  }
  return output;
}

} // namespace octoassociator
