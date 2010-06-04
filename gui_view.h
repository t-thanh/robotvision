/**
 * @author  Hauke Strasdat, Steven Lovegrove
 *
 * Copyright (C) 2010  Hauke Strasdat, Steven Lovegrove
 *                     Imperial College London
 *
 * gui_view.h is part of RobotVision.
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


#ifndef RV_GUI_VIEW_H
#define RV_GUI_VIEW_H

#include <TooN/TooN.h>
#include <TooN/se3.h>


#include <cvd/image_io.h>

#include "gui_window.h"
#include "Camera/linear_camera.h"
#include "quadtree.h"
#include "sim3.h"
#include "abstract_view.h"

namespace RobotVision
{
  class GuiView : public AbstractView
  {
  public:
    GuiView()
    {}

    GuiView(const GuiView & view);
    GuiView(const CVD::ImageRef & ir,
            const RobotVision::LinearCamera & cam_pars);
    GuiView(const CVD::ImageRef & ir,
            const RobotVision::LinearCamera & cam_pars,
            const TooN::Vector<3>& axis_angle,
            const TooN::Vector<3>& trans );

    void activate2D();
    void activate3D();
    void activate3D(const TooN::SE3<double> & T_cw);

      protected:
    RobotVision::LinearCamera cam_pars;

  };


}


#endif // RV_GUI_VIEW_H
