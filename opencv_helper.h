/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * opencv_helper.h is part of RobotVision.
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

#ifndef RV_OPENCV_HELPER_H
#define RV_OPENCV_HELPER_H

#include <opencv/cv.h>
#include <cvd/image_io.h>

namespace RobotVision
{

  void cvd2opencv(const CVD::SubImage<CVD::byte> &cvd_Im,
                  IplImage & openCv_Im);
  void cvd2opencv(const CVD::SubImage<CVD::Rgb<CVD::byte> > &cvd_Im,
                  IplImage & openCv_Im);
  void cvd2opencv(const CVD::SubImage<CVD::Rgba<CVD::byte> > &cvd_Im,
                  IplImage & openCv_Im);

  template<int R, int C>
  inline cv::Mat intmatrix2cv_mat(TooN::Matrix<R,C,int> & M){
    return cv::Mat(M.num_rows(),M.num_cols(),CV_32S,(M.my_data));
  }

  template<int R, int C>
  inline cv::Mat matrix2cv_mat(TooN::Matrix<R,C,float> & M){
    return cv::Mat(M.num_rows(),M.num_cols(),CV_32FC1,(M.my_data));
  }

  template<int R, int C>
  inline cv::Mat matrix2cv_mat(TooN::Matrix<R,C,double> & M){
    return cv::Mat(M.num_rows(),M.num_cols(),CV_64FC1,(M.my_data));
  }

  template<int R, int C>
  inline CvMat matrix2cvc_mat(TooN::Matrix<R,C,double> & M){
    cv::Mat cv_mat = cv::Mat(M.num_rows(),M.num_cols(),CV_64FC1,(M.my_data));
    return (CvMat)cv_mat;
  }

  template<int R, int C>
  inline CvMat matrix2cvc_mat(TooN::Matrix<R,C,float> & M){
    cv::Mat cv_mat = cv::Mat(M.num_rows(),M.num_cols(),CV_32FC1,(M.my_data));
    return (CvMat)cv_mat;
  }

  template<int R, int C>
  inline CvMat intmatrix2cvc_mat(TooN::Matrix<R,C,int> & M){
    cv::Mat cv_mat = cv::Mat(M.num_rows(),M.num_cols(),CV_32S,(M.my_data));
    return (CvMat)cv_mat;
  }


}

#endif // RV_OPENCV_HELPER_H
