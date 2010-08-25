/**
 * @author Hauke Strasdat, Steven Lovegrove, Andrew J. Davison
 *
 * Copyright (C) 2010 Hauke Strasdat, Steven Lovegrove, Andrew J. Davison
 *                    Imperial College London
 *
 * abstract_view.cpp is part of RobotVision.
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


#include "abstract_view.h"

#include <iostream>

#include <TooN/LU.h>
#include <TooN/SymEigen.h>

#include <cvd/gl_helpers.h>

#include "maths_utils.h"

using namespace CVD;
using namespace RobotVision;
using namespace TooN;
using namespace std;

void AbstractView::updatePose(const TooN::SE3<>& T_nc )
{
  T_cw_ = T_nc * T_cw_;
}

void AbstractView::updatePose(const TooN::SO3<>& R_nc )
{
  T_cw_ = SE3<>(R_nc,makeVector(0,0,0)) * T_cw_;
}

void AbstractView::setPose(const TooN::SE3<>& T_wc )
{
  T_cw_ = T_wc.inverse();
}

TooN::SE3<> AbstractView::getPose()
{
  return T_cw_.inverse();
}

AbstractView::AbstractView(const ImageRef & size)
  : Viewport(size), handler(NULL), T_cw_(Identity(3), makeVector(0,0,3))

{
  qtree = new QuadTree<int>(RobotVision::Rectangle(0,0,size.x,size.y),1);
}

AbstractView::AbstractView(const ImageRef & size,
                           const Vector<3>& axis_angle,
                           const Vector<3>& trans)
                             : Viewport(size),handler(NULL),
                             T_cw_(SO3<>(axis_angle),trans)
{
  qtree = new QuadTree<int>(RobotVision::Rectangle(0,0,size.x,size.y),1);
}

AbstractView::AbstractView(const AbstractView & view)
  :  Viewport(view.pixel_size), handler(NULL)
{
  qtree
      = new QuadTree<int>(Rectangle(0,0,
                                    view.pixel_size[0],view.pixel_size[1]),1);
}

CVD::GLWindow::EventHandler* AbstractView::get_handler()
{
  return handler;
}

void AbstractView::set_handler( CVD::GLWindow::EventHandler* handler )
{
  this->handler = handler;
}

GuiWindow* AbstractView::parent()
{
  return win;
}


bool AbstractView::isNavigating()
{
  return win && win->get_active_view() == this;
}

TooN::Vector<3> AbstractView::intersectScene( const CVD::ImageRef& ip )
{
  activate();
  glEnable(GL_DEPTH_TEST);

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  static const GLdouble modelMatrix[] = {
    1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
  };
  const GLdouble projMatrix[] = {
    glK[0][0], glK[1][0], glK[2][0], glK[3][0],
    glK[0][1], glK[1][1], glK[2][1], glK[3][1],
    glK[0][2], glK[1][2], glK[2][2], glK[3][2],
    glK[0][3], glK[1][3], glK[2][3], glK[3][3]
  };

  const GLfloat winx = viewport[0]+((float)ip.x / pixel_size[0])*viewport[2];
  const GLfloat winy = viewport[1]+((float)ip.y / pixel_size[1])*viewport[3];
  GLfloat winz;
  glReadPixels((GLint)winx,(GLint)winy,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&winz);

  if( 0 <= winz && winz < 1.0 )
  {
    GLdouble x,y,z;
    gluUnProject(winx,winy,winz,modelMatrix,projMatrix,viewport,&x,&y,&z);
    return makeVector(x,y,z);
  }else{
    return makeVector(0,0,0);
  }
}

void AbstractView::drawTexture2D(const SubImage<float> & img)
{
  glEnable(GL_TEXTURE_2D);
  GLuint texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(img);

  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glBindTexture(GL_TEXTURE_2D, texture);

  float x1 = 0;
  float y1  = 0;
  float x2 = pixel_size[0];
  float y2  = pixel_size[1];

  glBegin(GL_QUADS);

  glTexCoord2f(0.0f, 0.0f); glVertex2i(0, 0);
  glTexCoord2f(0.0f, 1.0f); glVertex2i(0, y2-y1);
  glTexCoord2f(1.0f, 1.0f); glVertex2i(x2-x1, y2-y1);
  glTexCoord2f(1.0f, 0.0f); glVertex2i(x2-x1, 0);

  glEnd();

  glDeleteTextures(1,&texture);
  glDisable(GL_TEXTURE_2D);
}

void AbstractView::drawTexture2D(const SubImage<CVD::byte> & img)
{
  glEnable(GL_TEXTURE_2D);
  GLuint texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(img);

  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glBindTexture(GL_TEXTURE_2D, texture);

  float x1 = 0;
  float y1  = 0;
  float x2 = pixel_size[0];
  float y2  = pixel_size[1];

  glBegin(GL_QUADS);

  glTexCoord2f(0.0f, 0.0f); glVertex2i(0, 0);
  glTexCoord2f(0.0f, 1.0f); glVertex2i(0, y2-y1);
  glTexCoord2f(1.0f, 1.0f); glVertex2i(x2-x1, y2-y1);
  glTexCoord2f(1.0f, 0.0f); glVertex2i(x2-x1, 0);

  glEnd();

  glDeleteTextures(1,&texture);
  glDisable(GL_TEXTURE_2D);
}

void AbstractView::drawTexture2D(const SubImage<Rgb<CVD::byte> > & img,
                                 const ImageRef & offset){
  glEnable(GL_TEXTURE_2D);
  GLuint texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(img);

  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glBindTexture(GL_TEXTURE_2D, texture);

  glBegin(GL_QUADS);

  glTexCoord2f(0.0f, 0.0f);
  glVertex2i(offset.x, offset.y);

  glTexCoord2f(0.0f, 1.0f);
  glVertex2i(offset.x, offset.y + img.size().y);

  glTexCoord2f(1.0f, 1.0f);
  glVertex2i(offset.x+img.size().x, offset.y + img.size().y);

  glTexCoord2f(1.0f, 0.0f);
  glVertex2i(offset.x+img.size().x, offset.y);

  glEnd();
  glDeleteTextures(1,&texture);
  glDisable(GL_TEXTURE_2D);
}


void AbstractView::drawTexture2D(const SubImage<Rgba<CVD::byte> > & img,
                                 const ImageRef & offset)
{
  glEnable(GL_TEXTURE_2D);
  GLuint texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(img);

  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glBindTexture(GL_TEXTURE_2D, texture);

  glBegin(GL_QUADS);

  glTexCoord2f(0.0f, 0.0f);
  glVertex2i(offset.x, offset.y);

  glTexCoord2f(0.0f, 1.0f);
  glVertex2i(offset.x, offset.y + img.size().y);

  glTexCoord2f(1.0f, 1.0f);
  glVertex2i(offset.x+img.size().x, offset.y + img.size().y);

  glTexCoord2f(1.0f, 0.0f);
  glVertex2i(offset.x+img.size().x, offset.y);

  glEnd();
  glDeleteTextures(1,&texture);
  glDisable(GL_TEXTURE_2D);
}




void AbstractView::drawTexture2D(const SubImage<CVD::byte> & img,
                                 const ImageRef & offset)
{
  glEnable(GL_TEXTURE_2D);
  GLuint texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(img);

  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glBindTexture(GL_TEXTURE_2D, texture);

  glBegin(GL_QUADS);

  glTexCoord2f(0.0f, 0.0f);
  glVertex2i(offset.x, offset.y);

  glTexCoord2f(0.0f, 1.0f);
  glVertex2i(offset.x, offset.y + img.size().y);

  glTexCoord2f(1.0f, 1.0f);
  glVertex2i(offset.x+img.size().x, offset.y + img.size().y);

  glTexCoord2f(1.0f, 0.0f);
  glVertex2i(offset.x+img.size().x, offset.y);

  glEnd();
  glDeleteTextures(1,&texture);

  glDisable(GL_TEXTURE_2D);
}

void AbstractView::drawTexture2D(const SubImage<Rgb<CVD::byte> > & img)
{
  glEnable(GL_TEXTURE_2D);
  GLuint texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(img);

  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  glBindTexture(GL_TEXTURE_2D, texture);

  float x1 = 0;
  float y1  = 0;
  float x2 = pixel_size[0];
  float y2  = pixel_size[1];

  glBegin(GL_QUADS);

  glTexCoord2f(0.0f, 0.0f); glVertex2i(0, 0);
  glTexCoord2f(0.0f, 1.0f); glVertex2i(0, y2-y1);
  glTexCoord2f(1.0f, 1.0f); glVertex2i(x2-x1, y2-y1);
  glTexCoord2f(1.0f, 0.0f); glVertex2i(x2-x1, 0);
  glEnd();
  glDeleteTextures(1,&texture);

  glDisable(GL_TEXTURE_2D);
}




Vector<2> AbstractView::pixel2opengl(const Vector<2> & v)
{
  //In the image space, the center of the first pixel is denoted by (0,0)
  //whereas the position of the top left corner
  //of the first pixel is actually (-0.5,-0.5).

  return makeVector(v[0]+0.5,v[1]+0.5);
}


void AbstractView::drawBall3D(const TooN::Vector<3> & p1, double radius)
{
  glPushMatrix();
  glTranslate(p1);
  gluSphere((GLUquadric* )win->quad, radius,4,4);
  glPopMatrix();
}

void AbstractView::drawPose3D(const TooN::SE3<double> & pose)
{
  glPushMatrix();
  Vector<3> center = RobotVision::trans2center(pose);

  glTranslate(center);
  Vector<3> aa = pose.inverse().get_rotation().ln();
  double angle = norm(aa);
  if ( angle != 0.0 ) {
    glRotatef( angle * 180.0 / M_PI, aa[0], aa[1], aa[2]);
  }
  drawLine3D(TooN::makeVector(0,0,0), TooN::makeVector(0.05, 0, 0));
  drawLine3D(TooN::makeVector(0,0,0), TooN::makeVector(0, 0.05, 0));
  drawLine3D(TooN::makeVector(0,0,0), TooN::makeVector(0, 0, 0.05));
  drawLine3D(TooN::makeVector(0.025,0.025,0),
             TooN::makeVector(-0.025, 0.025, 0));
  drawLine3D(TooN::makeVector(-0.025,-0.025,0),
             TooN::makeVector(-0.025, 0.025, 0));
  drawLine3D(TooN::makeVector(-0.025,-0.025,0),
             TooN::makeVector(0.025, -0.025, 0));
  drawLine3D(TooN::makeVector(0.025,0.025,0),
             TooN::makeVector(0.025, -0.025, 0));
  glPopMatrix();
}


void AbstractView::drawSimilarity3D(const RobotVision::Sim3<double> & pose)
{

  glPushMatrix();
  Vector<3> center = RobotVision::trans2center(pose);

  glTranslate(center);
  Vector<3> aa = pose.inverse().get_rotation().ln();
  double angle = norm(aa);
  if ( angle != 0.0 ) {
    glRotatef( angle * 180.0 / M_PI, aa[0], aa[1], aa[2]);
  }

  drawLine3D(TooN::makeVector(0,0,0), TooN::makeVector(0.05, 0, 0));

  drawLine3D(TooN::makeVector(0,0,0), TooN::makeVector(0, 0.05, 0));

  drawLine3D(TooN::makeVector(0,0,0), TooN::makeVector(0, 0, 0.05));


  double size = 0.025;//*1./pose.get_scale();

  drawLine3D(TooN::makeVector(size,size,0), TooN::makeVector(-size, size, 0));
  drawLine3D(TooN::makeVector(-size,-size,0), TooN::makeVector(-size, size, 0));
  drawLine3D(TooN::makeVector(-size,-size,0), TooN::makeVector(size, -size, 0));
  drawLine3D(TooN::makeVector(size,size,0), TooN::makeVector(size, -size, 0));

  glPopMatrix();

}



void AbstractView::drawCircle2D(const Vector<2> & p,
                                double inner_radius,
                                double outer_radius)
{
  Vector<2> v = pixel2opengl(p);
  glPushMatrix();
  glTranslatef(v[0],v[1],0);
  gluDisk((GLUquadric* )win->quad, inner_radius, outer_radius, 20, 2 );
  glPopMatrix();
}

void AbstractView::drawCovariance2D(const Vector<2> & mu,
                                    const Matrix<2>& Sigma,
                                    double number_of_sigma,
                                    double ring_thickness)
{
  Vector<2> v = pixel2opengl(mu);
  glPushMatrix();
  glTranslatef(v[0],v[1],0);

  SymEigen<2> e_system(Sigma);
  Matrix<2> e_vectors = e_system.get_evectors();
  Vector<2> e_values = e_system.get_evalues();

  // Eigenvectors become rotation matrix
  // Turn into 4*4 transformation

  //NOT SURE WHY THIS IS NECESASSY (THE MINUS SIGNS).
  //PROBABLY BECAUSE OF THE DIFFERENCE BETWEEN
  //LEFT-HAND (Image) AND RIGHT-HAND (OpenGL) COORDINATE SYSTEMS.
  GLfloat Varray[ 16 ] = {e_vectors( 0, 0 ), -e_vectors( 1, 0 ), 0.0, 0.0,
                          -e_vectors( 0, 1 ), e_vectors( 1, 1 ), 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 1.0};
  glMultMatrixf( Varray );

  double factor_x = sqrt( e_values[0] ) * number_of_sigma;
  double factor_y = sqrt( e_values[1] ) * number_of_sigma;

  glScalef( factor_x, factor_y, 1);

  double inner_size = ( factor_x - ring_thickness ) / factor_x;
  if ( inner_size < 0.0 )
    inner_size = 0.0;

  gluDisk((GLUquadric* )win->quad , inner_size, 1.0, 20,20);
  glPopMatrix();
}

void AbstractView::drawLine3D(const Vector<3> & p1, const Vector<3> & p2)
{
  glBegin(GL_LINES);
  glVertex3f(p1[0], p1[1], p1[2]);
  glVertex3f(p2[0], p2[1], p2[2]);
  glEnd();
}

void AbstractView::drawLine2D(const Vector<2> & p1, const Vector<2> & p2)
{
  Vector<2> v1 = pixel2opengl(p1);
  Vector<2> v2 = pixel2opengl(p2);

  glBegin(GL_LINES);
  glVertex2f(v1[0], v1[1]);
  glVertex2f(v2[0], v2[1]);
  glEnd();
}

void AbstractView::drawBox2D(const RobotVision::Rectangle & r)
{
  Vector<2> v1 = pixel2opengl(makeVector(r.x1,r.y1));
  Vector<2> v2 = pixel2opengl(makeVector(r.x2,r.y2));

  glBegin(GL_LINE_LOOP);
  glVertex2f(v1[0], v1[1]);
  glVertex2f(v1[0], v2[1]);
  glVertex2f(v2[0], v2[1]);
  glVertex2f(v2[0], v1[1]);
  glEnd();
}

void AbstractView::activate()
{
  const float x = bounding_box.x1*win->size().x;
  const float y = bounding_box.y1*win->size().y;
  const float w = (bounding_box.x2-bounding_box.x1)*win->size().x;
  const float h = (bounding_box.y2-bounding_box.y1)*win->size().y;
  glViewport (x, y, w, h);
}

void AbstractView::init(GuiWindow * win, const RobotVision::Rectangle & box){
  win->connectView(*this,box);
}

void AbstractView::drawCovariance3D(const Vector<3> & trans,
                                    const Matrix<3> & pose_unc,
                                    double number_of_sigma){
  //ToDo: Check whether this is doing the right thing!
  //Should be alright now!
  glPushMatrix();
  glTranslatef(trans[0], trans[1], trans[2]);

  SymEigen<3> e_system(pose_unc);
  Matrix<3> & V = e_system.get_evectors();
  GLfloat Varray[ 16 ] = {V( 0, 0 ), V( 0, 1 ), V( 0, 2 ), 0.0,
                          V( 1, 0 ), V( 1, 1 ), V( 1, 2 ), 0.0,
                          V( 2, 0 ), V( 2, 1 ), V( 2, 2 ), 0.0,
                          0.0, 0.0, 0.0, 1.0};
  glMultMatrixf( Varray );
  Vector<3> v = e_system.get_evalues();
  glScalef( sqrt( v[0] ) * number_of_sigma,
            sqrt( v[1] ) * number_of_sigma,
            sqrt( v[2] ) * number_of_sigma );
  gluSphere((GLUquadric* )win->quad, 1, 20, 20 );
  glPopMatrix();


}


void DoomViewHandler::on_event(GLWindow& win, int event_id)
{
  if( event_id == GuiWindow::EVENT_VIEW_ACTIVATED )
  {
    on_activated(win);
  }else if( event_id == GuiWindow::EVENT_VIEW_DEACTIVATED )
  {
    on_deactivated(win);
  }
}

void DoomViewHandler::on_activated(CVD::GLWindow& /*win*/)
{
}

void DoomViewHandler::on_deactivated(CVD::GLWindow& /*win*/)
{
  enabled = false;
}

void DoomViewHandler::on_key_down(CVD::GLWindow& /*win*/, int key)
{
  if( enabled) {
    switch(key){
    case 119:
      view->updatePose(SE3<>(Identity,makeVector(0,0,-0.03) ) );
      break;
    case 115:
      view->updatePose(SE3<>(Identity,makeVector(0,0,+0.03) ) );
      break;
    case 97:
      view->updatePose(SE3<>(Identity,makeVector(-0.03,0,0) ) );
      break;
    case 100:
      view->updatePose(SE3<>(Identity,makeVector(+0.03,0,0) ) );
      break;
    case 114:
      view->updatePose(SE3<>(Identity,makeVector(0,+0.03,0) ) );
      break;
    case 102:
      view->updatePose(SE3<>(Identity,makeVector(0,-0.03,0) ) );
      break;
    case 113:
      view->updatePose(Rotation(0,0,-0.05));
      break;
    case 101:
      view->updatePose(Rotation(0,0,+0.05));
      break;
    }
  }
}

void DoomViewHandler::on_mouse_down(CVD::GLWindow& /*win*/,
                                    CVD::ImageRef pos,
                                    int /*state*/,
                                    int /*button*/)
{
  enabled = !enabled;
  last_pos = pos;
}

void DoomViewHandler::on_mouse_up(CVD::GLWindow& /*win*/,
                                  CVD::ImageRef /*where*/,
                                  int /*state*/,
                                  int /*button*/)
{
}

void DoomViewHandler::on_mouse_move(CVD::GLWindow& /*win*/,
                                    CVD::ImageRef pos,
                                    int /*state*/)
{
  if( enabled ) {
    view->updatePose( Rotation(
        (pos.x-last_pos.x)*0.01,
        (pos.y-last_pos.y)*0.01,
        0
        ));
  }
  last_pos = pos;
}

