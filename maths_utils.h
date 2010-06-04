/**
 * @author Hauke Strasdat, Steven Lovegrove
 *
 * Copyright (C) 2010  Hauke Strasdat, Steven Lovegrove
 *                     Imperial College London
 *
 * maths_utils.h is part of RobotVision.
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


#ifndef MATHSUTILS_H
#define MATHSUTILS_H

#include <math.h>
#include <set>
#include <assert.h>
#include <TooN/TooN.h>

#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <TooN/SVD.h>
#include <cvd/image_io.h>

#include "./Camera/abstract_camera.h"



// Bodge C99 isfinite macro for windows
#ifdef WIN32
#define isfinite(x) (_finite(x))
#endif


namespace RobotVision
{
  template <typename Precision>
      class Sim3;


  template <typename T>
      inline T Po2(const T& value)
  {
    return value * value;
  }

  template <typename T>
      inline T Po3(const T& value)
  {
    return Po2(value) * value;
  }

  template <typename T>
      T Po4(const T& value)
  {
    return Po2(Po2(value));
  }

  template <typename T>
      T Po5(const T& value)
  {
    return Po4(value) * value;
  }

  template <typename T>
      T Po6(const T& value)
  {
    return Po4(value) * Po2(value);
  }


  template <typename T>
      T Po7(const T& value)
  {
    return Po6(value) * value;
  }


  template <typename T>
      T Po8(const T& value)
  {
    return Po2(Po4(value));
  }

  template <typename T>
      inline int sign(const T& value)
  {
    return value < 0 ? -1 : 1;
  }

  template<typename Precision, typename WTF>
  TooN::Matrix<3,3,Precision> skew(const TooN::Vector<3,Precision,WTF> & v)
  {
    TooN::Matrix<3,3,Precision> S = TooN::Zeros;
    S(0,1) = -v[2];
    S(0,2) = v[1];
    S(1,2) = -v[0];

    S(1,0) = v[2];
    S(2,0) = -v[1];
    S(2,1) = v[0];

    return S;
  }

  TooN::Vector<3> deltaR(const TooN::Matrix<3> & R);

  TooN::Vector<3> trans2center(const TooN::SE3<double>& pose);
  TooN::Vector<3> trans2center(const RobotVision::Sim3<double>& pose);

  TooN::SO3<> Rotation(double yaw, double pitch, double roll );

  template<int Size, typename Precision, typename Base>
  inline TooN::Vector<Size,Precision,Base> element_product
      (
          const TooN::Vector<Size,Precision,Base>& v1,
          const TooN::Vector<Size,Precision,Base>& v2
          )
  {
    TooN::Vector<Size,Precision,Base> r;
    for( int i=0; i < Size; ++i )
      r[i] = v1[i] * v2[i];
    return r;
  }

  template<int K, int L,int M, int N, typename Precision>
  TooN::Matrix<TooN::Dynamic,TooN::Dynamic,Precision>
      kron(const TooN::Matrix<K,L,Precision> & A,
           const TooN::Matrix<M,N,Precision> & B)
  {
    int k = A.num_rows();
    int l = A.num_cols();
    int m = B.num_rows();
    int n = B.num_cols();
    int q = k*m;
    int r = l*n;
    TooN::Matrix<TooN::Dynamic,TooN::Dynamic,Precision> C(q,r);
    for (int i=0; i<k; ++i)
    {
      for (int j=0; j<l; ++j)
      {
        C.slice(i*m,j*n,m,n) = A(i,j)*B;
      }
    }
    return C;
  }



  template<int K, int L, typename Precision>
  TooN::Vector<TooN::Dynamic,Precision>
      vec(const TooN::Matrix<K,L,Precision> & M)
  {
    uint rows = M.num_rows();
    uint cols = M.num_cols();
    TooN::Vector<TooN::Dynamic,Precision> v(rows*cols);

    for (uint i=0; i<cols;++i)
    {
      v.slice(i*rows,rows) = M.T()[i];
    }

    return v;
  }



  template<int K,  typename Precision, typename Order>
  TooN::Matrix<TooN::Dynamic,TooN::Dynamic,Precision>
      unvec(const TooN::Vector<K,Precision,Order> & v,
            int num_rows,
            int num_cols)
  {
    assert(v.size() == num_rows*num_cols);

    TooN::Matrix<TooN::Dynamic,TooN::Dynamic,Precision>  M(num_rows,num_cols);

    for (int i=0; i<num_cols;++i)
    {
      M.T()[i] = v.slice(i*num_rows,num_rows);
    }

    return M;

  }

  double norm1(const TooN::Vector<> & v);

  double norm_max(const TooN::Vector<> & v);

  double norm_max(const TooN::Vector<> & v, int & idx);

  double depth(const TooN::SE3<double> & pose, const TooN::Vector<3> & XYZ);

  template <typename T>
      T median(const std::set<T> & m_set)
  {
    int size = m_set.size();
    typename std::set<T>::const_iterator it = m_set.begin();
    if (m_set.size()%2==1)
    {
      for (int i=0; i<size/2; ++i)
      {
        it++;
      }
      return *it;
    }

    for (int i=0; i<size/2-1; ++i)
    {
      it++;
    }
    T val = *it;
    it++;
    val += *it;
    return 0.5*val;
  }


  TooN::SE3<> CameraPose2(const TooN::Vector<3>& t_w,
                          const TooN::Vector<3>& lookat_w,
                          const TooN::Vector<3>& up_w );

}




#endif //FRAMEUTILS_H
