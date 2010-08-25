/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * stopwatch.h is part of RobotVision.
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

#ifndef RV_STOP_WATCH_H
#define RV_STOP_WATCH_H

#include <cassert>

#include <cvd/timer.h>

namespace RobotVision
{
  class StopWatch{
  public:
    StopWatch()
    {
      running = false;
      time = 0;
    }

    inline void start()
    {
      assert(!running);
      timer.reset();
      running = true;
    }

    inline void stop()
    {
      assert(running);
      time = timer.get_time();
      running = false;
    }

    inline double readCurrentTime()
    {
      assert(running);
      return timer.get_time();
    }

    inline double getStoppedTime()
    {
      assert(!running);
      return time;
    }

    inline void reset(){
      time = 0;
    }

  private:
    CVD::cvd_timer timer;
    double time;
    bool running;
  };
}

#endif //RV_STOP_WATCH_H
