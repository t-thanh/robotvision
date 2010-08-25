/**
 * @author Hauke Strasdat
 *
 * Copyright (C) 2010 Hauke Strasdat
 *                    Imperial College London
 *
 * linear_camera.h is part of RobotVision.
 *
 * RobotVision is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or any later version.
 *
 * RobotVision is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and the GNU Lesser General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RV_LINEAR_CAMERA_H
#define RV_LINEAR_CAMERA_H

#include "abstract_camera.h"


namespace RobotVision {
  class LinearCamera : public AbstractCamera
  {
  public:
    static const int NUM_PARAMS = 4;

    // Construct Identity Camera
    LinearCamera();

    // Contract from Linear portion of cam
    LinearCamera(const AbstractCamera& cam);

    // Construct Linear camera given intrinsics matrix
    LinearCamera(const TooN::Matrix<3,3> K,
                 const CVD::ImageRef & size);

    // Construct Linear camera given parameters
    LinearCamera(const TooN::Vector<2> & focal_length,
                 const TooN::Vector<2> & principle_point,
                 const CVD::ImageRef & size);

    const CVD::ImageRef& size() const
    {
      return image_size;
    }

    /// Principle Point
    const TooN::Vector<2>& P() const
    {
      return principle_point;
    }

    /// Focal Length
    const TooN::Vector<2>& F() const
    {
      return focal_length;
    }

    /// Focal Length
    const TooN::Vector<2>& Finv() const
    {
      return inv_focal_length;
    }

    const TooN::Matrix<3,3>& K() const
    {
      return intrinsics;
    }

    const TooN::Matrix<3,3>& Kinv() const
    {
      return inv_intrinsics;
    }

    const TooN::Vector<> params() const
    {
      return fl_and_pp;
    }

    TooN::Vector<2> map(const TooN::Vector<2>& camframe) const;
    TooN::Vector<2> unmap(const TooN::Vector<2>& imframe) const;

    TooN::Matrix<2> jacobian(const TooN::Vector<2>& camframe) const;
    TooN::Matrix<2> inv_jacobian(const TooN::Vector<2>& imframe) const;

    TooN::Matrix<2> jacobian() const;
    TooN::Matrix<2> inv_jacobian() const;

  private:
    void initialise();

    //Trick: Do computation once, and store ambiguous represenations
    TooN::Vector<> fl_and_pp; //
    TooN::Vector<2> focal_length;
    TooN::Vector<2> principle_point;
    CVD::ImageRef image_size;
    TooN::Matrix<3> intrinsics;
    TooN::Matrix<3> inv_intrinsics;
    TooN::Vector<2> inv_focal_length;

  };
}
#endif // RV_LINEAR_CAMERA_H
