.. _parameters:

How to set the parameters
=========================

As most associators, PyOcto has numerous parameters and it can be hard to navigate the configuration jungle.
Therefore, here is a quick guide on the parameters structured by observed behaviour.
For every behaviour, we provide a list of suggestions on how to mitigate the issue and potential drawbacks.
As always, recommended behaviour is to change one parameter at a time and observe the impact.
For a definition of each parameter, check out the :ref:`associator page<associator>`.

Setting the study zone
----------------------
The parameters `xlim`, `ylim`, and `zlim` define your study zone.
As an easy rule, the larger the study zone the longer the run time
and the more memory will be needed.
In addition, you will need to set a velocity model.
We'd recommend starting with a homogeneous model and then switching
to a 1D model.
The last mandatory parameter is the `time_before`, the overlap between subsequent time slices.
This value should usually be roughly the maximum S travel time between any station and
the boundary of the study zone.
Setting this value too high will cause increased runtime but a safe margin is recommended.
Typical values for regional deployments (50 to 2000 km diameter) are between 60 and 300 s.

PyOcto requires a local coordinate transform and local coordinates for the limits.
If you instead want to just provide the latitude, longitude and depth limits,
you can use the :py:meth:`pyocto.associator.OctoAssociator.from_area` method to create your associator.
It will automatically select a reasonable projection and set the correct limits.
Afterwards, to transform the stations,
the :py:meth:`pyocto.associator.OctoAssociator.transform_stations` function might help.

PyOcto runs for too long
----------------------------

* Increase the values of `n_picks`, `n_p_picks`, `n_s_picks`, or `n_p_and_s_picks`.
  PyOcto uses early stopping, so any node that does not contain enough picks any more will be discarded immediately.
  Higher thresholds therefore lead to better run times at the cost of fewer events.
  Note that events with few picks are anyhow often false positive detections.
* Increase the `min_node_size`. In general this should be adapted to your station spacing.
  If stations are far apart, the value should be larger and vice versa.
  Typical values for regional deployments (50 to 2000 km diameter) are 5 to 25 kilometers.
* Increase `min_node_size_location`. This will make your location estimates less accurate but will improve runtime.
  Note that less accurate location estimates might also lead to more missing picks for the events
  as they are used for matching picks to events.
  In general, `min_node_size_location` should always be smaller than `min_node_size`.
* Reduce `location_split_depth` or `location_split_return`. A lot of time for PyOcto is spend in estimating locations.
  By reducing `location_split_depth` fewer cells are tested in each split, leading to quicker convergence.
  By reducing `location_split_return`, fewer points per cell are sampled, leading to quicker decent and thereby quicker
  convergence as well. On the other hand, lower values for either of these parameters increase the susceptibility to
  local minima and might impact the location quality. Remember, `location_split_depth` always needs to be
  bigger than `location_split_return`.
* Set (or reduce) the `association_cutoff_distance`. This distance defines picks up to which distance are included
  in the space-partitioning search. Especially for large regions, this often makes sense as usually earthquakes are
  substantially easier to detect at nearby stations. Note that further away stations will still be used for the
  localisation process (and therefore be annotated in the assignments output). To limit these assignments, use the
  `location_cutoff_distance` parameter. This will usually only have little impact on the runtime.
* Increase `min_pick_fraction`. This value cuts off picks at stations with higher distance if closer stations did not
  pick the event. In contrast to the `association_cutoff_distance`, this is a relative parameter, i.e., it only removes
  far away stations if closer stations exists. This behaviour is especially helpful for deployments with large numbers
  of stations that are likely to produce spurious associations from far away station sets. Reasonable values are
  between 0 and 0.7. (Maybe even higher, depending on the noise characteristics and quality of your stations.)
  Too high values will lead to missed events because good picks are discarded.
  Setting `min_pick_fraction` to zero will disable the distance check completely and thereby improve runtime.
  This is especially useful for smaller deployments.
* **Minor:** Decrease `refinement_iterations`. Each event is located multiple times. After each iteration, PyOcto associates
  all picks consistent with the location to the event. Lowering the value will usually improve runtimes slightly
  but might also lead to fewer associated picks.
* **Carefully:** Reduce the `tolerance` in your velocity model.
  Picks will need to be more exactly within the nodes if the value is lower.
  On the other hand, too low tolerance will discard events if they are right on the boundary between boxes.

PyOcto runs out of memory
-------------------------
* Reduce `n_threads`. Memory usage is roughly linear with the number of threads, so lowering this number will lower
  the memory used. However, if your memory usage is dominated by one time segment, this might not help to much.
  Note that the lower memory usage is bought at the cost of a higher run time.
* Reduce the `time_slicing` so that every slice needs less memory. However, this will also lead to more overlap,
  i.e., more picks being processed multiple times. This will increase run times.
* Increase `n_picks`, `n_p_picks`, `n_s_picks`, `n_p_and_s_picks`, `min_node_size_location` or `min_pick_fraction`
  or decrease `tolerance`. All of this lead to fewer explored nodes and, therefore, lower memory usage.
  See above for the advantages and caveats of each of these parameters.
* Reduce `queue_memory_protection_dfs_size`. This is the most direct way to control memory usage, as it effectively
  caps the memory used. Lower values will lead to lower memory consumption at the cost of higher runtimes. In most
  cases, there should be now substantial impacts on the resulting catalogs.

PyOcto only finds few events
----------------------------
* Decrease your requirements (`n_picks`, `n_p_picks`, `n_s_picks`, `n_p_and_s_picks`, `min_pick_fraction`) or increase
  your tolerance (`tolerance` in velocity models, `pick_match_tolerance`). Note that excessive values might lead to
  high numbers of false positives.
* Use a better velocity model. Especially in large and deep study areas, 1D models usually find more events.
* **Rarely:** Decrease `min_node_size` or `min_node_size_location`. If these values are too high, PyOcto might fail
  to locate events or associate picks to them. This will usually only be relevant if your study area is very small
  or you increased the values substantially above the default values.
* Make sure there are actually more events in your data. Maybe you need a more sensitive picker?

Event locations show stripes
-----------------------------
* Increase `location_split_depth` or `location_split_return`. This will make locations more accurate but at the cost
  of longer runtimes. See above for details on the two parameters.
* Switch to a better velocity model. Good 1D velocity models usually show fewer stripes.
* Accept it and use a dedicated tool for relocalisation to get better locations. PyOcto is first and foremost an
  associator. Dedicated tools like NonLinLoc or HypoDD will without question produce more accurate locations and
  in addition come with well-calibrated uncertainty estimates and many more helpful features. PyOcto offers an
  interface to write picks directly in a format for NonLinLoc.

PyOcto does not use all my CPU cores
------------------------------------

* You're running on Windows? Unfortunately, PyOcto does not support multi-threading on Windows as the OS does not
  support POSIX threads. For larger computations, we suggest switching to a different operating system. If you are
  an experienced C++ developer and would like to add multi-threading support for Windows, we'd happily receive a PR.
* By default, PyOcto will use all available cores. Check if you accidentally overwrote the `n_threads` variable to a
  lower values. It is generally not recommend to set `n_threads` to a higher number than the available cores as this
  will increase memory usage without significant runtime decreases.
* If PyOcto starts running on all cores but later only uses a few or even one core, this is because some time blocks
  take longer to associate. You can reduce the `time_slicing` to distribute the load more evenly. However, this will
  also lead to more overlap, i.e., more picks being processed multiple times.

PyOcto uses all my CPU cores and I don't want that
--------------------------------------------------
* Set `n_threads` to the number of cores you want to use.

PyOcto crashes or does something else unexpected
------------------------------------------------
Open an issue on the `PyOcto Github repository <https://github.com/yetinam/pyocto>`_.
