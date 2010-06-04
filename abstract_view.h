/**
 * @author Hauke Strasdat, Steven Lovegrove, Andrew J. Davison
 *
 * Copyright (C) 2010 Hauke Strasdat, Steven Lovegrove, Andrew J. Davison
 *                    Imperial College London
 *
 * abstract_view.h is part of RobotVision.
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


#ifndef RV_ABSTRACT_VIEW_H
#define RV_ABSTRACT_VIEW_H


#include <TooN/TooN.h>
#include <TooN/se3.h>

#include <cvd/glwindow.h>
#include <cvd/image_io.h>

#include "gui_window.h"
#include "quadtree.h"
#include "sim3.h"

namespace RobotVision
{
  class GuiWindow;

  class AbstractView
  {
    friend class GuiWindow;
  public:
    AbstractView()
    {}

    ~AbstractView(){
      if (qtree!=0)
        delete qtree;
    }

    AbstractView(const AbstractView & view);
    AbstractView(const CVD::ImageRef & size);
    AbstractView(const CVD::ImageRef & size,
                 const TooN::Vector<3>& axis_angle,
                 const TooN::Vector<3>& trans );

    // Attach view to GuiWindow
    void init(GuiWindow * win, const Rectangle & bounding_box);

    // Return default event handler
    CVD::GLWindow::EventHandler* get_handler();

    // Set default event handler
    void set_handler( CVD::GLWindow::EventHandler* handler );

    // Return parent GuiWindow
    GuiWindow* parent();

    // Return true iff the view is currently active
    bool isNavigating();

    // Return intersection of image point ip with the rendered scene
    // in frame of reference of the view's camera
    TooN::Vector<3> intersectScene( const CVD::ImageRef& ip );

    // Active Viewport
    void activate();

    virtual void activate2D() = 0;
    virtual void activate3D() = 0;
    virtual void activate3D(const TooN::SE3<double> & pose) = 0;

    void drawTexture2D(const CVD::SubImage<CVD::byte> & img);
    void drawTexture2D(const CVD::SubImage<CVD::Rgb<CVD::byte> > & img);
    void drawTexture2D(const CVD::SubImage<CVD::byte> & img,
                       const CVD::ImageRef & ir);
    void drawTexture2D(const CVD::SubImage<CVD::Rgb<CVD::byte> > & img,
                       const CVD::ImageRef & ir);

    void drawTexture2D(const CVD::SubImage<CVD::Rgba<CVD::byte> > & img);
    void drawTexture2D(const CVD::SubImage<CVD::Rgba<CVD::byte> > & img,
                       const CVD::ImageRef & ir);

    void drawCircle2D(const TooN::Vector<2> & p1,
                      double inner_radius,
                      double outer_radius);
    void drawLine2D(const TooN::Vector<2> & p1,
                    const TooN::Vector<2> & p2);
    void drawBox2D(const Rectangle & r);
    void drawCovariance2D(const TooN::Vector<2> & mu,
                          const TooN::Matrix<2>& Sigma,
                          double number_of_sigma,
                          double ring_thickness );

    void drawLine3D(const TooN::Vector<3> & p1, const TooN::Vector<3> & p2);
    void drawBall3D(const TooN::Vector<3> & p1, double radius);
    void drawPose3D(const TooN::SE3<double> & pose);
    void drawSimilarity3D(const RobotVision::Sim3<double> & pose);
    void drawCovariance3D(const TooN::Vector<3> & trans,
                          const TooN::Matrix<3> & pose_unc,
                          double number_of_sigma);

    // Update pose given delta transformation Current to New.
    void updatePose(const TooN::SE3<>& T_nc );
    void updatePose(const TooN::SO3<>& R_nc );

    // Set pose of camera, camera to world
    void setPose(const TooN::SE3<>& T_wc );

    // Return pose T_wc of virtual camera
    TooN::SE3<> getPose();

    inline CVD::ImageRef size() const {
      return pixel_size;
    }

    inline CVD::ImageRef viewSize() const {
      return CVD::ImageRef(bounding_box.Width()*win->size().x,
                           bounding_box.Height()*win->size().y );
    }

  protected:
    CVD::GLWindow::EventHandler* handler;
    CVD::ImageRef pixel_size;
    TooN::Matrix<4> glK;
    TooN::SE3<double>  T_cw_;
    TooN::Vector<2> pixel2opengl(const TooN::Vector<2> & v);
    GuiWindow * win;
    int id;
    Rectangle bounding_box;
    TooN::QuadTree<int>  * qtree;
  };

  class DoomViewHandler : public CVD::GLWindow::EventHandler
  {
  public:
    DoomViewHandler( AbstractView* view ) : view(view), enabled(false) {}
    void on_key_down(CVD::GLWindow& /*win*/, int /*key*/);
    void on_mouse_down(CVD::GLWindow& /*win*/,
                       CVD::ImageRef /*where*/,
                       int /*state*/,
                       int /*button*/);
    void on_mouse_up(CVD::GLWindow& /*win*/,
                     CVD::ImageRef /*where*/,
                     int /*state*/,
                     int /*button*/);
    void on_mouse_move(CVD::GLWindow& /*win*/,
                       CVD::ImageRef /*where*/,
                       int /*state*/);
    void on_event(CVD::GLWindow& /*win*/,
                  int /*event_id*/);
    void on_activated(CVD::GLWindow& /*win*/);
    void on_deactivated(CVD::GLWindow& /*win*/);
  protected:
    AbstractView* view;
    CVD::ImageRef last_pos;
    bool enabled;
  };


}
#endif // RV_ABSTRACT_VIEW_H
