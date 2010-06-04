/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * image.h is part of RobotVision.
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


#ifndef RV_IMAGE_H
#define RV_IMAGE_H

#ifdef RV_OPENCV_SUPPORT
#include <opencv/cv.h>
#endif

#include <cvd/image_io.h>

namespace RobotVision
{

  class MonoImage : public CVD::Image<CVD::byte>
  {
  public:
    MonoImage();
    MonoImage(const MonoImage& copy);
    MonoImage(const CVD::ImageRef & size);
#ifdef RV_OPENCV_SUPPORT
    cv::Mat cv_mat;
    IplImage cvc_image;
    CvMat * cvc_mat_pt;
  private:
    CvMat cvc_mat_stub;
#endif
  };


  class ColorImage : public CVD::Image<CVD::Rgb<CVD::byte> >
  {
  public:
    ColorImage();
    ColorImage(const ColorImage& copy);
    ColorImage(const CVD::ImageRef & size);
#ifdef RV_OPENCV_SUPPORT
    cv::Mat cv_mat;
    IplImage cvc_image;
    CvMat * cvc_mat_pt;
  private:
    CvMat cvc_mat_stub;
#endif
  };



  class RgbaImage : public CVD::Image<CVD::Rgba<CVD::byte> >
  {
  public:
    RgbaImage();
    RgbaImage(const RgbaImage& copy);
    RgbaImage(const CVD::ImageRef & size);
#ifdef RV_OPENCV_SUPPORT
    cv::Mat cv_mat;
    IplImage cvc_image;
    CvMat * cvc_mat_pt;
  private:
    CvMat cvc_mat_stub;
#endif
  };

}

#endif // RV_IMAGE_H
