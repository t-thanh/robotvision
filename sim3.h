/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * sim3.h is part of RobotVision.
 *
 * RobotVision is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or any later version.
 *
 * RobotVision is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR Layout PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and the GNU Lesser General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef RV_SIMILARITY_H
#define RV_SIMILARITY_H

#include "maths_utils.h"

#include <TooN/se3.h>
#include <TooN/LU.h>

#define RV_PROPER_SIM3 1


namespace RobotVision
{

  /**
   * This class implements the Lie group and the corresponding Lie
   * algebra of 3D similarity transformationns Sim3 as described in:
   *
   * > H. Strasdat, J.M.M. Montiel, A.J. Davison:
   *   "Scale Drift-Aware Large Scale Monocular SLAM",
   *   Proc. of Robotics: Science and Systems (RSS),
   *   Zaragoza, Spain, 2010.
   *   http://www.roboticsproceedings.org/rss06/p10.html <
   *
   * The over-all design of this class is inspired by an older
   * version (April 2009) of the SE3 class of the TooN library
   * (which was written by Tom Drummond and others).
   * See: http://mi.eng.cam.ac.uk/~er258/cvd/toon.html
   */
  template <typename Precision = double> class Sim3
                                 {
                                 public:

    inline Sim3() : my_translation(TooN::Zeros), my_scale(1) {}

    template <int Size, typename Layout>
        Sim3(const TooN::SO3<> & R,
             const TooN::Vector<Size, Precision, Layout>& t,
             const Precision s)
               : my_rotation(R), my_translation(t), my_scale(s) {}

    template <int Size, typename Layout>
        Sim3(const TooN::Matrix<3> & R,
             const TooN::Vector<Size, Precision, Layout>& t,
             const Precision s)
               : my_rotation(TooN::SO3<>(R)), my_translation(t), my_scale(s) {}

    template <int Size, typename Layout>
        Sim3(const TooN::Vector<Size, Precision, Layout> & v)
    {
      *this = exp(v);
    }

    inline TooN::SO3<Precision>& get_rotation()
    {
      return my_rotation;
    }

    inline const TooN::SO3<Precision>& get_rotation() const
    {
      return my_rotation;
    }

    inline TooN::Vector<3, Precision>& get_translation()
    {
      return my_translation;
    }

    inline const TooN::Vector<3, Precision>& get_translation() const
    {
      return my_translation;
    }

    inline Precision& get_scale()
    {
      return my_scale;
    }

    inline const Precision& get_scale() const
    {
      return my_scale;
    }

    template <int Size, typename Layout>
        static inline Sim3
        exp(const TooN::Vector<Size, Precision, Layout>& vect)
    {
#ifdef RV_PROPER_SIM3
      TooN::Vector<3,Precision> upsilon = vect.slice(0,3);
      TooN::Vector<3,Precision> omega = vect.slice(3,3);
      Precision sigma = vect[6];
      Precision eps = 0.00001;

      Precision theta = norm(omega);
      TooN::Matrix<3,3,Precision> Omega = RobotVision::skew(omega);
      TooN::Matrix<3,3,Precision> Omega2 = Omega*Omega;
      TooN::Matrix<3,3,Precision> R;

      TooN::Matrix<3,3,Precision> I = TooN::Identity;;
      Precision s = std::exp(sigma);

      Precision A,B,C;
      if (fabs(sigma)<eps)
      {
        C = 1;
        if (theta<eps)
        {
          A = 1./2.;
          B = 1./6.;
          R = (I + Omega + Omega*Omega);
        }
        else
        {
          A = (1-cos(theta))/Po2(theta);
          B = (theta-sin(theta))/Po3(theta);
          R = I + sin(theta)/theta *Omega + (1-cos(theta))/(Po2(theta))*Omega2;
        }
      }
      else
      {
        C=(s-1)/sigma;
        if (theta<eps)
        {
          A = ((sigma-1)*s+1)/Po2(sigma);
          B= ((0.5*Po2(sigma)-sigma+1)*s)/Po3(sigma);
          R = (I + Omega + Omega2);
        }
        else
        {
          R = I + sin(theta)/theta *Omega + (1-cos(theta))/(Po2(theta))*Omega2;



          Precision a=s*sin(theta);
          Precision b=s*cos(theta);
          Precision c=Po2(theta)+Po2(sigma);
          A = (a*sigma+ (1-b)*theta)/(theta*c);
          B = (C-((b-1)*sigma+a*theta)/(c))*1./(Po2(theta));
        }
      }

      TooN::Matrix<3,3,Precision> W = A*Omega + B*Omega2 + C*I;
      TooN::Vector <3,Precision> t = W*upsilon;
      return Sim3(R, t, s);
#else
      TooN::SE3<Precision> se3 (vect.slice(0,6));
      return Sim3(se3.get_rotation(), se3.get_translation(), std::exp(vect[6]));
#endif
    }

    static inline TooN::Vector<7,Precision> ln(const Sim3& sim3);

    inline TooN::Vector<7,Precision> ln() const
    {
      return Sim3::ln(*this);
    }

    inline Sim3 inverse() const
    {
      const TooN::SO3<Precision> Rinv = get_rotation().inverse();
      return Sim3(Rinv, -(1./my_scale)*(Rinv*my_translation), 1./my_scale);
    }

    inline Sim3& operator *=(const Sim3& sim3)
                            {
      get_translation()
          += get_scale()*(get_rotation() * sim3.get_translation());
      get_rotation() *= sim3.get_rotation();
      get_scale() *= sim3.my_scale;
      return *this;
    }

    inline Sim3 operator *(const Sim3& sim3) const
    {
      return Sim3(get_rotation()*sim3.get_rotation(),
                  get_scale()*(get_rotation()*sim3.get_translation())
                  + get_translation(),
                  get_scale()*sim3.get_scale());
    }

  protected:
    TooN::SO3<Precision> my_rotation;
    TooN::Vector<3, Precision> my_translation;
    Precision my_scale;
  };

  template <typename Precision>
      inline std::ostream& operator <<(std::ostream& out_str,
                                       const Sim3<Precision>& sim3)
  {
    for(int i=0; i<3; ++i)
    {
      out_str << sim3.get_rotation().get_matrix()[i]
          << sim3.get_translation()[i] << std::endl;
    }
    out_str << sim3.get_scale() << std::endl;
    return out_str;
  }

  template <typename Precision>
      inline TooN::Vector<7, Precision>
      Sim3<Precision>::ln(const Sim3<Precision>& sim3)
  {
#ifdef RV_PROPER_SIM3
    TooN::Vector<7, Precision> res;
    Precision s = sim3.my_scale;
    Precision sigma = log(s);

    TooN::Vector<3,Precision> t = sim3.my_translation;
    TooN::Matrix<3,3,Precision> R = sim3.my_rotation.get_matrix();
    Precision d =  0.5*(R(0,0)+R(1,1)+R(2,2)-1);

    TooN::Vector<3,Precision> omega;
    TooN::Vector<3,Precision> upsilon;
    TooN::Matrix<3,3,Precision> Omega;

    Precision eps = 0.00001;
    TooN::Matrix<3,3,Precision> I = TooN::Identity;;

    Precision A,B,C;
    if (fabs(sigma)<eps)
    {
      C = 1;
      if (d>1-eps)
      {
        omega=0.5*RobotVision::deltaR(R);
        Omega = RobotVision::skew(omega);
        A = 1./2.;
        B = 1./6.;

      }
      else
      {
        Precision theta = acos(d);
        omega = theta/(2*sqrt(1-d*d))*RobotVision::deltaR(R);
        Omega = RobotVision::skew(omega);
        A = (1-cos(theta))/Po2(theta);
        B = (theta-sin(theta))/Po3(theta);
      }
    }
    else
    {
      C=(s-1)/sigma;
      if (d>1-eps)
      {
        omega=0.5*RobotVision::deltaR(R);
        Omega = RobotVision::skew(omega);
        A = ((sigma-1)*s+1)/Po2(sigma);
        B = ((0.5*Po2(sigma)-sigma+1)*s)/Po3(sigma);
      }
      else
      {
        Precision theta = acos(d);
        omega = theta/(2*sqrt(1-d*d))*RobotVision::deltaR(R);
        Omega = RobotVision::skew(omega);
        Precision a=s*sin(theta);
        Precision b=s*cos(theta);
        Precision c=Po2(theta)+Po2(sigma);
        A = (a*sigma+ (1-b)*theta)/(theta*c);
        B = (C-((b-1)*sigma+a*theta)/(c))*1./(Po2(theta));
      }
    }

    TooN::Matrix<3,3,Precision> W = A*Omega + B*Omega*Omega + C*I;
    TooN::LU<3,Precision> LU_W(W);
    upsilon = LU_W.backsub(t);
    res.slice(0,3) = upsilon;
    res.slice(3,3) = omega;
    res[6] = sigma;
    return res;
#else
    TooN::Vector<7, Precision> res;
    res.slice(0,6) = TooN::SE3<Precision>(sim3.get_rotation(),
                                          sim3.get_translation()).ln();
    res[6] = log(sim3.my_scale);
    return res;
#endif

  }
}


#endif // RV_SIMILARITY_H
