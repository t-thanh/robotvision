/**
 * @author Hauke Strasdat
 *
 * Copyright (C) 2010 Hauke Strasdat
 *                    Imperial College London
 *
 * linear_camera.cpp is part of RobotVision.
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


#include "linear_camera.h"

#include "../maths_utils.h"

using namespace TooN;

namespace RobotVision{

  LinearCamera::LinearCamera()
    : fl_and_pp(4),
    focal_length(makeVector(1,1)),
    principle_point(Zeros),
    image_size(-1,-1) ,
    intrinsics(Identity),
    inv_intrinsics(Identity(3))
  {
    fl_and_pp.slice<0,2>() = focal_length;
    fl_and_pp.slice<2,2>() = principle_point;
  }

  LinearCamera::LinearCamera(const AbstractCamera& cam)
    :fl_and_pp(4),
    focal_length(cam.F()),
    principle_point(cam.P()),
    image_size(cam.size()) ,
    intrinsics(Identity),
    inv_intrinsics(Identity(3))
  {
    initialise();
  }

  LinearCamera::LinearCamera(const TooN::Matrix<3,3> K,
               const CVD::ImageRef & size)
    :fl_and_pp(4),
    focal_length(makeVector( K[0][0]/K[2][2], K[1][1]/K[2][2])),
    principle_point(makeVector( K[0][2]/K[2][2], K[1][2]/K[2][2])),
    image_size(size),
    intrinsics(Identity)
  {
    initialise();
  }

  LinearCamera::LinearCamera(const Vector<2> & focal_length,
                             const Vector<2> & principle_point,
                             const CVD::ImageRef & size)
                               :fl_and_pp(4),
                               focal_length(focal_length),
                               principle_point(principle_point),
                               image_size(size) ,
                               intrinsics(Identity),
                               inv_intrinsics(Identity(3))
  {
    initialise();
  }

  void LinearCamera::initialise()
  {
    fl_and_pp.slice<0,2>() = focal_length;
    fl_and_pp.slice<2,2>() = principle_point;
    inv_focal_length = makeVector(1./focal_length[0], 1./focal_length[1]);
    intrinsics(0,0) = focal_length[0];
    intrinsics(1,1) = focal_length[1];
    intrinsics(0,2) = principle_point[0];
    intrinsics(1,2) = principle_point[1];
    inv_intrinsics(0,0) = inv_focal_length[0];
    inv_intrinsics(1,1) = inv_focal_length[1];
    inv_intrinsics(0,2) = -principle_point[0]*inv_focal_length[0];
    inv_intrinsics(1,2) = -principle_point[1]*inv_focal_length[1];
  }


  TooN::Vector<2> LinearCamera::map(const TooN::Vector<2>& point) const
  {
    return element_product(focal_length,point) + principle_point;
  }

  TooN::Vector<2> LinearCamera::unmap(const TooN::Vector<2>& dist_point) const
  {
    return element_product(inv_focal_length,(dist_point-principle_point));
  }

  TooN::Matrix<2>
      LinearCamera::jacobian(const TooN::Vector<2>& camframe) const
  {
    return jacobian();
  }

  TooN::Matrix<2>
      LinearCamera::inv_jacobian(const TooN::Vector<2>& imframe) const
  {
    return inv_jacobian();
  }

  TooN::Matrix<2> LinearCamera::jacobian() const
  {
    return intrinsics.slice<0,0,2,2>();
  }

  TooN::Matrix<2> LinearCamera::inv_jacobian() const
  {
    return inv_intrinsics.slice<0,0,2,2>();
  }

}
