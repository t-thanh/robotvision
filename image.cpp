/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * image.cpp is part of RobotVision.
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

#include "image.h"

#ifdef RV_OPENCV_SUPPORT
#include "opencv_helper.h"
#endif

using namespace CVD;
using namespace RobotVision;

MonoImage::MonoImage() : Image<byte>(ImageRef(1,1))
{
#ifdef RV_OPENCV_SUPPORT
  cvd2opencv(*this, cvc_image);
  cv_mat = cv::Mat(&cvc_image,false);
  cvc_mat_pt = cvGetMat(&cvc_image,&cvc_mat_stub);
#endif
}


MonoImage::MonoImage(const MonoImage& copy) : Image<byte>(copy)
{
#ifdef RV_OPENCV_SUPPORT
  cvd2opencv(*this, cvc_image);
  cv_mat = cv::Mat(&cvc_image,false);
  cvc_mat_pt = cvGetMat(&cvc_image,&cvc_mat_stub);
#endif
}



MonoImage::MonoImage(const ImageRef & size) : Image<byte>(size)
{
#ifdef RV_OPENCV_SUPPORT
  cvd2opencv(*this, cvc_image);
  cv_mat = cv::Mat(&cvc_image,false);
  cvc_mat_pt = cvGetMat(&cvc_image,&cvc_mat_stub);
#endif
}





ColorImage::ColorImage() : Image<Rgb<byte> >()
{

}


ColorImage::ColorImage(const ColorImage& copy) : Image<Rgb<byte> >(copy)
{
#ifdef RV_OPENCV_SUPPORT
  cvd2opencv(*this, cvc_image);
  cv_mat = cv::Mat(&cvc_image,false);

  cvc_mat_pt = cvGetMat(&cvc_image,&cvc_mat_stub);
#endif
}



ColorImage::ColorImage(const ImageRef & size) : Image<Rgb<byte> >(size)
{
#ifdef RV_OPENCV_SUPPORT
  cvd2opencv(*this, cvc_image);
  cv_mat = cv::Mat(&cvc_image,false);
  cvc_mat_pt = cvGetMat(&cvc_image,&cvc_mat_stub);
#endif
}




RgbaImage::RgbaImage() : Image<Rgba<byte> >()
{

}


RgbaImage::RgbaImage(const RgbaImage& copy) : Image<Rgba<byte> >(copy)
{
#ifdef RV_OPENCV_SUPPORT
  cvd2opencv(*this, cvc_image);
  cv_mat = cv::Mat(&cvc_image,false);

  cvc_mat_pt = cvGetMat(&cvc_image,&cvc_mat_stub);
#endif
}



RgbaImage::RgbaImage(const ImageRef & size) : Image<Rgba<byte> >(size)
{
#ifdef RV_OPENCV_SUPPORT
  cvd2opencv(*this, cvc_image);
  cv_mat = cv::Mat(&cvc_image,false);
  cvc_mat_pt = cvGetMat(&cvc_image,&cvc_mat_stub);
#endif
}




