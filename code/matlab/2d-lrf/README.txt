/******************************************************************************
 * Copyright (C) 2012 by Jerome Maye                                          *
 * jerome.maye@gmail.com                                                      *
 *                                                                            *
 * This program is free software; you can redistribute it and/or modify       *
 * it under the terms of the Lesser GNU General Public License as published by*
 * the Free Software Foundation; either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * Lesser GNU General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the Lesser GNU General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 ******************************************************************************/

This folder contains some MATLAB scripts related to a 2d robot equipped with
a Laser Range Finder (LRF). Two extrinsic calibration models are considered for
the LRF, namely with a single offset parameter, and with 2 offset parameters
and one angle. For each of the model, we have:

- calibration
- localization,
- Simultaneous Localization And Mapping (SLAM),
- SLAM and calibration

For each of the above problems, we have implemented:

- Extended Kalman Filter (EKF)
- non-linear least square with or without regularization
- iterative non-linear least square with batch selection

The SLAM problem considered here is landmark-based and assume known data
association.

There is real dataset from Tim Barfoot that can be used for testing the
algorithms in the case of a single offset parameter. Two simulation scripts are
as well provided for playing with the parameters.
