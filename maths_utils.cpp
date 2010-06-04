/**
 * @author Hauke Strasdat, Steven Lovegrove
 *
 * Copyright (C) 2010  Hauke Strasdat, Steven Lovegrove
 *                     Imperial College London
 *
 * maths_utils.cpp is part of RobotVision.
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

#include <TooN/helpers.h>
#include <TooN/SVD.h>
#include <TooN/LU.h>

#include "maths_utils.h"
#include "sim3.h"
#include <assert.h>
#include <limits>

using namespace std;
using namespace TooN;


namespace RobotVision
{
  TooN::Vector<3> deltaR(const TooN::Matrix<3> & R)
  {
    return TooN::makeVector(R(2,1)-R(1,2),R(0,2)-R(2,0),R(1,0)-R(0,1));
  }

  Vector<3> trans2center(const SE3<double>& pose){
    return pose.inverse().get_translation();
  }

  Vector<3> trans2center(const Sim3<double>& pose){
    return pose.inverse().get_translation();
  }

  TooN::SO3<> Rotation(double yaw, double pitch, double roll )
  {
    const double cyaw = cos(yaw);
    const double syaw = sin(yaw);
    const double croll = cos(roll);
    const double sroll = sin(roll);
    const double cpitch = cos(pitch);
    const double spitch = sin(pitch);
    const Matrix<3,3> R = Data(
        cyaw * croll,
        syaw*spitch - cyaw*sroll*cpitch,
        cyaw*sroll*spitch + syaw*cpitch,
        sroll,
        croll*cpitch,
        -croll*spitch,
        -syaw*croll,
        syaw*sroll*cpitch + cyaw*spitch,
        -syaw*sroll*spitch + cyaw*cpitch
        );
    return SO3<>(R);
  }

  double norm1(const Vector<> & v)
  {
    double sum = 0;
    for (int i=0; i<v.size(); i++)
    {
      sum += fabs(v[i]);
    }
    return sqrt(sum);
  }

  double norm_max(const Vector<> & v)
  {
    double max = -1;
    for (int i=0; i<v.size(); i++)
    {
      double abs = fabs(v[i]);
      if(abs>max){
        max = abs;
      }
    }
    return max;
  }

  double norm_max(const Vector<> & v, int & idx)
  {
    double max = -1;
    idx = -1;
    for (int i=0; i<v.size(); i++)
    {
      double abs = fabs(v[i]);
      if(abs>max){
        max = abs;
        idx = i;
      }
    }
    return max;
  }

  double depth(const SE3<double> & pose, const Vector<3> & XYZ)
  {
    Matrix<3,4> P;
    P.slice(0,0,3,3) = pose.get_rotation().get_matrix();
    P.T()[3] = pose.get_translation();
    Matrix<3> M = pose.get_rotation().get_matrix();
    SVD<> SVD_M(M);
    Vector<4> p3 = P[2];
    Vector<3> m3 = M[2];
    double w = p3*unproject(XYZ);
    return (sign(SVD_M.determinant())*w) / (norm(m3));
  }

  SE3<> CameraPose2(const Vector<3>& t_w,
                    const Vector<3>& lookat_w,
                    const Vector<3>& up_w )
  {
    TooN::Matrix<3,3> R_cw;
    R_cw.T()[2] = unit((lookat_w - t_w));
    R_cw.T()[0] = up_w ^ R_cw.T()[2];
    R_cw.T()[1] = R_cw.T()[2] ^ R_cw.T()[0];

    return SE3<>( SO3<>(R_cw.T()), -R_cw.T()*t_w);
  }

}


