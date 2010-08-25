/**
 * @author  Hauke Strasdat, Steven Lovegrove
 *
 * Copyright (C) 2010  Hauke Strasdat, Steven Lovegrove
 *                     Imperial College London
 *
 * gui_window.h is part of RobotVision.
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

#ifndef RV_GUI_WINDOW_H
#define RV_GUI_WINDOW_H

#include <vector>
#include <list>
#include <set>

#include <cvd/glwindow.h>

#include "rectangle.h"
#include "quadtree.h"

namespace RobotVision
{
  // Forward declarations
  class GuiWindow;

  class Viewport
  {
    friend class GuiWindow;

  public:
    Viewport() {}
    Viewport(const CVD::ImageRef& size )
      : pixel_size(TooN::makeVector(size.x,size.y)) { }
    Viewport(const TooN::Vector<2>& size )
      : pixel_size(size) { }

    virtual CVD::GLWindow::EventHandler* get_handler() = 0;

    inline Rectangle& BoundingBox() { return bounding_box; }
    inline TooN::Vector<2>& PixelSize() { return pixel_size; }

    inline float Left() { return bounding_box.x1; }
    inline float Top() { return bounding_box.y1; }
    inline float Width() { return bounding_box.Width(); }
    inline float Height() { return bounding_box.Height(); }
  protected:
    GuiWindow* win;
    TooN::Vector<2> pixel_size;
    Rectangle bounding_box;
  };

  class GuiWindow : public CVD::GLWindow, public CVD::GLWindow::EventHandler
  {
    friend class Viewport;
    friend class ViewEventHandler;

  public:
    const static int EVENT_VIEW_ACTIVATED = 2134;
    const static int EVENT_VIEW_DEACTIVATED = 2135;

    GuiWindow(const CVD::ImageRef & ir,
              const CVD::GLWindow* sharedContext = NULL);
    ~GuiWindow();

    void nextFrame(CVD::GLWindow::EventHandler& handler);

    void handle_events_default();

    bool closed();

    void on_key_down(CVD::GLWindow& /*win*/, int /*key*/);
    void on_key_up(CVD::GLWindow& /*win*/, int /*key*/);
    void on_mouse_move(CVD::GLWindow& /*win*/,
                       CVD::ImageRef /*where*/,
                       int /*state*/);
    void on_mouse_down(CVD::GLWindow& /*win*/,
                       CVD::ImageRef /*where*/,
                       int /*state*/,
                       int /*button*/);
    void on_mouse_up(CVD::GLWindow& /*win*/,
                     CVD::ImageRef /*where*/,
                     int /*state*/, int /*button*/);
    void on_resize(CVD::GLWindow& /*win*/,
                   CVD::ImageRef /*size*/);
    void on_event(CVD::GLWindow& /*win*/,
                  int /*event*/);

    void set_active_view(Viewport* view);

    Viewport* get_active_view();

    // Connect view to this window within region box
    void connectView(Viewport & view, const Rectangle & box);

    std::set<int> active_point_set;

    void * quad;

  protected:

    // Return view that contains window point wp
    Viewport* get_view( const CVD::ImageRef& wp, CVD::ImageRef* vp );

    // Return coordinates relative to view
    CVD::ImageRef view_coords( const Viewport* view,
                               const CVD::ImageRef wp ) const;

//    int num_id;
    std::list<Viewport*> viewpt_list;
    CVD::ImageRef active_pos;
    Viewport * active_view;
    int mode;
  };
}

#endif // RV_GUI_WINDOW_H
