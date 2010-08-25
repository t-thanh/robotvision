/**
 * @author  Steven Lovegrove, Hauke Strasdat
 *
 * Copyright (C) 2010  Steven Lovegrove, Hauke Strasdat
 *                     Imperial College London
 *
 * abstract_camera.h is part of RobotVision.
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

#ifndef RV_ABSTRACT_CAMERA_H
#define RV_ABSTRACT_CAMERA_H

#include <TooN/TooN.h>
#include <cvd/image_ref.h>

namespace RobotVision {

  class AbstractCamera
  {
  public:
    /// Image Size
    virtual const CVD::ImageRef& size() const = 0;

    /// Principle Point
    virtual const TooN::Vector<2>& P() const = 0;

    /// Focal Length
    virtual const TooN::Vector<2>& F() const = 0;

    virtual const TooN::Matrix<3,3>& K() const = 0;
    virtual const TooN::Matrix<3,3>& Kinv() const = 0;


    //applies camera intrinsics including distortion
    virtual TooN::Vector<2> map(const TooN::Vector<2>& camframe) const = 0;

    //undo camera intrinsics including undistortion
    virtual TooN::Vector<2> unmap(const TooN::Vector<2>& imframe) const = 0;


    virtual TooN::Matrix<2> jacobian(const TooN::Vector<2>& camframe) const = 0;


    virtual const TooN::Vector<> params() const = 0;

    // Project p_cam in to the camera, applying intrinsics and distortion
    inline TooN::Vector<2> projectmap(const TooN::Vector<3>& p_cam) const
    {
      return map( TooN::project(p_cam) );
    }

    inline int width() const { return size().x; }
    inline int height() const { return size().y; }
    inline int pixel_area() const { return size().x * size().y; }
    inline double aspect() const { return (double)size().x / (double)size().y; }
  };

}

#endif // RV_ABSTRACT_CAMERA_H
