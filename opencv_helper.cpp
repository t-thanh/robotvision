/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * opencv_helper.cpp is part of RobotVision.
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

#ifdef RV_OPENCV_SUPPORT

#include "opencv_helper.h"

using namespace CVD;

namespace RobotVision
{

  void cvd2opencv(const SubImage<byte> &cvd_Im, IplImage & openCv_Im)
  {
    int width = cvd_Im.size().x;
    int height = cvd_Im.size().y;
    int channels = 1;

    openCv_Im.nSize = sizeof(IplImage);
    openCv_Im.ID = 0;
    openCv_Im.nChannels = channels;
    openCv_Im.depth = sizeof(byte)*8;
    openCv_Im.dataOrder = 0;
    openCv_Im.origin = 0;
    openCv_Im.width = width;
    openCv_Im.height = height;
    openCv_Im.roi = NULL;
    openCv_Im.maskROI = NULL;
    openCv_Im.imageId = NULL;
    openCv_Im.tileInfo = NULL;
    openCv_Im.imageSize = width*height*channels;
    openCv_Im.imageData = (char *)(cvd_Im.data());
    openCv_Im.widthStep = width*channels;
    openCv_Im.imageDataOrigin = NULL;
  }

  void cvd2opencv(const SubImage<Rgb<byte> > &cvd_Im, IplImage & openCv_Im)
  {
    int width = cvd_Im.size().x;
    int height = cvd_Im.size().y;
    int channels = 3;

    openCv_Im.nSize = sizeof(IplImage);
    openCv_Im.ID = 0;
    openCv_Im.nChannels = channels;
    openCv_Im.depth = sizeof(byte)*8;
    openCv_Im.dataOrder = 0;
    openCv_Im.origin = 0;
    openCv_Im.width = width;
    openCv_Im.height = height;
    openCv_Im.roi = NULL;
    openCv_Im.maskROI = NULL;
    openCv_Im.imageId = NULL;
    openCv_Im.tileInfo = NULL;
    openCv_Im.imageSize = width*height*channels;
    openCv_Im.imageData = (char *)(cvd_Im.data());
    openCv_Im.widthStep = width*channels;
    openCv_Im.imageDataOrigin = NULL;
  }

  void cvd2opencv(const SubImage<Rgba<byte> > &cvd_Im, IplImage & openCv_Im)
  {
    int width = cvd_Im.size().x;
    int height = cvd_Im.size().y;
    int channels = 4;

    openCv_Im.nSize = sizeof(IplImage);
    openCv_Im.ID = 0;
    openCv_Im.nChannels = channels;
    openCv_Im.depth = sizeof(byte)*8;
    openCv_Im.dataOrder = 0;
    openCv_Im.origin = 0;
    openCv_Im.width = width;
    openCv_Im.height = height;
    openCv_Im.roi = NULL;
    openCv_Im.maskROI = NULL;
    openCv_Im.imageId = NULL;
    openCv_Im.tileInfo = NULL;
    openCv_Im.imageSize = width*height*channels;
    openCv_Im.imageData = (char *)(cvd_Im.data());
    openCv_Im.widthStep = width*channels;
    openCv_Im.imageDataOrigin = NULL;
  }

}

#endif
