.. _data_formats:

Data formats
=========================

PyOcto relies on `Pandas DataFrames <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`_ for all inputs and outputs.
This allows easy input from and serialisation to common formats such as csv.
In addition, PyOcto offers compatibility interfaces for SeisBench, Obspy, GaMMA, REAL, and NonLinLoc.
These interfaces are documented on the :ref:`associator documentation<associator>`.

Pick format
-----------

The pick dataframe has one pick per row and the following columns:

* **station:** a string identifier for the station
* **phase:** the phase label, either an uppercase `P` or an uppercase `S`
* **time:** the time as a timestamp, i.e., a float in seconds since a common time base.
  In general, unix timestamps are recommended but all other time bases work just as well.

There is no implementation limit to the number of picks provided to PyOcto.
There is also no limit to the duration spanned by the picks, i.e., picks for several days can be provided at once.
However, PyOcto does not write intermediate results, so it might be advisable to chunk the input by time in case
of large datasets.

Station format
--------------

The station dataframe has one station per row and the following columns:

* **id:** a string identifier for the station.
  This needs to match the identifier in the `station` column of the picks.
* **x**, **y**, **z:** the coordinates of the station in km in a local coordinate projection.
  The z axis is oriented downwards, i.e., stations above ground have negative sign, below the ground positive.
* *(Optional)* **p_residual** and **s_residual:** Travel time residuals for the P and S wave.
  If not present or NaN, these values are set to 0. We implement additive station terms, i.e.,
  the station term is added to the predicted travel term from the velocity model.
  This is the same sign convention as used by NonLinLoc.

Output format
-------------

The output of PyOcto consists of two dataframes, the events and the assignments.
The assignment contains all picks that have been associated to an event together with the respective event index.
To merge the two frames, join them on the `idx` (events) and `event_idx` (assignments) columns.
We recommend the `Pandas merge method <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html>`_,
i.e., `pd.merge(events, assignments, left_on="idx", right_on="event_idx")`.

The event dataframe contains on event per row the following columns:

* **idx:** an event identifier. Event identifier are unique.
* **time:** the source time of the event in seconds to the time base of the picks, i.e., usually a unix timestamp.
* **x**, **y**, **z:** the coordinates of the origin in km in a local coordinate projection
* **picks:** the number of merged picks

The assignments dataframe contains one pick per row and the following columns:

* **event_idx:** the event identifier
* **pick_idx:** the pick identifier. Pick identifier are assigned integers starting at zero in the order the picks
  have been provided.
* **residual:** the travel time residual from the localisation procedure in seconds. The residuals use the same
  sign convention as the station terms in the station dataframe, i.e., the predicted time plus the residual equals
  the observed travel time.

In addition, the assignments dataframe contains all columns that were contained in the original pick dataframe.
This includes additional columns that have not been used.
