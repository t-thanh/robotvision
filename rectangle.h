/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * rectangle.h is part of RobotVision.
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


#ifndef RV_RECTANGLE_H
#define RV_RECTANGLE_H

#include <cvd/image_io.h>

#include "maths_utils.h"

namespace RobotVision
{

  class Rectangle
  {
  public:
    Rectangle(){}

    Rectangle(const Rectangle & r)
    {
      this->x1 = r.x1;
      this->x2 = r.x2;
      this->y2 = r.y2;
      this->y1 = r.y1;
    }

    Rectangle(double x1, double y1, double x2, double y2)
    {
      this->x1 = x1;
      this->x2 = x2;
      this->y2 = y2;
      this->y1 = y1;
    }

    Rectangle(const TooN::Vector<2> & tl, const TooN::Vector<2> & br)
    {
      this->x1 = tl[0];
      this->x2 = br[0];
      this->y2 = br[1];
      this->y1 = tl[1];
    }

    double x1, y1, x2, y2;

    double Width() const
    {
      return std::max(0.,x2-x1);
    }

    double Height() const
    {
      return std::max(0.,y2-y1);
    }

    bool intersectsWith(const Rectangle & other) const
    {
      if(y2 <= other.y1)	return false;
      if(y1    >= other.y2)	return false;
      if(x2  <= other.x1)	return false;
      if(x1   >= other.x2)	return false;
      return true;
    }

    bool contains(const Rectangle & other) const
    {
      if(  y1  > other.y1)   return false;
      if(  x1 > other.x1)  return false;
      if( x2 < other.x2) return false;
      if(y2 < other.y2) return false;
      return true;
    }

    bool contains(const TooN::Vector<2> & v) const
    {
      if(  y1  > v[1])   return false;
      if(  x1 > v[0])  return false;
      if( x2 < v[0]) return false;
      if(y2 < v[1]) return false;
      return true;
    }
  };


  template <class T>
      class RectObj  : public Rectangle
  {
  public:
    RectObj<T> () {}

    RectObj<T> (const Rectangle & r, const T & t)
    {
      this->x1 = r.x1;
      this->x2 = r.x2;
      this->y2 = r.y2;
      this->y1 = r.y1;
      this->t = t;
    }

    RectObj<T> (const TooN::Vector<2> & tl,
                const TooN::Vector<2> & br,
                const T & t)
    {
      this->x1 = tl[0];
      this->x2 = br[0];
      this->y2 = br[1];
      this->y1 = tl[1];
      this->t = t;
    }

    RectObj<T> (double x1, double y1, double x2, double y2, const T & t)
    {
      this->x1 = x1;
      this->x2 = x2;
      this->y2 = y2;
      this->y1 = y1;
      this->t = t;
    }

    T t;
  };


  class IRectangle{
  public:
    IRectangle(){}

    IRectangle(const IRectangle & r)
    {
      this->x1 = r.x1;
      this->x2 = r.x2;
      this->y2 = r.y2;
      this->y1 = r.y1;
    }

    IRectangle(int x1, int y1, int x2, int y2)
    {
      this->x1 = x1;
      this->x2 = x2;
      this->y2 = y2;
      this->y1 = y1;
    }

    IRectangle(const CVD::ImageRef & tl, const CVD::ImageRef & br)
    {
      this->x1 = tl.x;
      this->x2 = br.x;
      this->y2 = br.y;
      this->y1 = tl.y;
    }

    int x1, y1, x2, y2;

    int Width()
    {
      return std::max(0,x2-x1);
    }

    int Height(){
      return std::max(0,y2-y1);
    }

    bool IntersectsWith(const IRectangle & other)
    {
      if(y2 <other.y1)	return false;
      if(y1    > other.y2)	return false;
      if(x2  < other.x1)	return false;
      if(x1   > other.x2)	return false;
      return true;
    }

    bool Contains(const IRectangle & other){
      if(  y1  >= other.y1)   return false;
      if(  x1 >= other.x1)  return false;
      if( x2 <= other.x2) return false;
      if(y2 <= other.y2) return false;
      return true;
    }
  };

  template <class T>
      class IRectObj  : public IRectangle
  {
  public:
    IRectObj<T> () {}

    IRectObj<T> (const IRectangle & r, const T & t)
    {
      this->x1 = r.x1;
      this->x2 = r.x2;
      this->y2 = r.y2;
      this->y1 = r.y1;
      this->t = t;
    }

    IRectObj<T> (const CVD::ImageRef & tl,
                 const CVD::ImageRef & br,
                 const T & t)
    {
      this->x1 = tl.x;
      this->x2 = br.x;
      this->y2 = br.y;
      this->y1 = tl.y;
      this->t = t;
    }

    IRectObj<T> (int x1, int y1, int x2, int y2, const T & t)
    {
      this->x1 = x1;
      this->x2 = x2;
      this->y2 = y2;
      this->y1 = y1;
      this->t = t;
    }

    T t;
  };

}

#endif // RV_RECTANGLE_H
