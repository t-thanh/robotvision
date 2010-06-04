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
  class AbstractView;

  class GuiWindow : public CVD::GLWindow, public CVD::GLWindow::EventHandler
  {
    friend class AbstractView;
    friend class ViewEventHandler;

  public:
    const static int EVENT_VIEW_ACTIVATED = 2134;
    const static int EVENT_VIEW_DEACTIVATED = 2135;

    GuiWindow(const CVD::ImageRef & ir,
              const CVD::GLWindow* sharedContext = NULL);
    ~GuiWindow();

    void nextFrame(CVD::GLWindow::EventHandler& handler);

    void handle_events_default();

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

    void set_active_view(AbstractView* view);

    std::set<int> active_point_set;

  private:
    // Connect view to this window within region box
    void connectView(AbstractView & view, const Rectangle & box);

    // Return view that contains window point wp
    AbstractView* get_view( const CVD::ImageRef& wp, CVD::ImageRef* vp );

    // Return coordinates relative to view
    CVD::ImageRef view_coords( const AbstractView* view,
                               const CVD::ImageRef wp ) const;

    int num_id;
    void * quad;
    std::list<AbstractView*> viewpt_list;
    CVD::ImageRef active_pos;
    AbstractView * active_view;
    int mode;

  };
}

#endif // RV_GUI_WINDOW_H
