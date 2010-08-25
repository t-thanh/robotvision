  /**
 * @author  Hauke Strasdat, Steven Lovegrove
 *
 * Copyright (C) 2010  Hauke Strasdat, Steven Lovegrove,
 *                     Imperial College London
 *
 * gui_view.cpp is part of RobotVision.
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

#include <GL/gl.h>
#include <GL/glext.h>

#include "gui_view.h"

using namespace CVD;
using namespace RobotVision;
using namespace TooN;
using namespace std;


GuiView::GuiView(const ImageRef & ir,
                 const RobotVision::LinearCamera & cam_pars)
  : AbstractView(ir), cam_pars(cam_pars)
{

}

GuiView::GuiView(const ImageRef & ir,
                 const RobotVision::LinearCamera & cam_pars,
                 const Vector<3>& axis_angle,
                 const Vector<3>& trans)
                   : AbstractView(ir,axis_angle, trans), cam_pars(cam_pars)
{

}

GuiView::GuiView(const GuiView & view)
  : AbstractView(view)
{
}

void GuiView::activate2D()
{
  glShadeModel(GL_FLAT);
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
  glLineWidth (1);

  activate();

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  ::glOrtho (0, pixel_size[0], pixel_size[1], 0, 0, 1);

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity();

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
}

void GuiView::activate3D()
{
  activate3D(T_cw_);
}

void GuiView::activate3D(const SE3<double> & T_wc)
{
  glShadeModel(GL_FLAT);

  activate();

  glMatrixMode(GL_PROJECTION);

  double fu = cam_pars.F()[0];
  double fv = -cam_pars.F()[1];
  double px = cam_pars.P()[0]+0.5;
  double py = cam_pars.P()[1]+0.5;

  glLoadIdentity();


  GLdouble left, right, bottom, top;
  double nr = 0.1;
  double fr = 1000;

  left = -nr * px / fu;
  top = nr * py / fv;
  right = nr * ( pixel_size[0] - px ) / fu;
  bottom = -nr * ( pixel_size[1] - py ) / fv;

  ::glFrustum( left, right, bottom, top, nr, fr );

  Vector<3> aa = T_wc.get_rotation().ln();
  double angle = norm(aa);

  glScaled(1,1,-1);

  glTranslatef(T_wc.get_translation()[0],T_wc.get_translation()[1],
               T_wc.get_translation()[2]);
  glRotated(angle* 180.0 / M_PI, aa[0], aa[1], aa[2]);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
}
