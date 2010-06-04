/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * quadtree.h is part of RobotVision.
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

#ifndef RV_QUADTREE_H
#define RV_QUADTREE_H

#include <list>
#include <queue>
#include <stack>
#include <cvd/random.h>

#include "maths_utils.h"
#include "rectangle.h"

namespace TooN
{

  template <typename T> class QuadTreeElement
  {
  public:
    QuadTreeElement() : dfs_visited(false)
    {
    }

    QuadTreeElement (const Vector<2> & pos,
                     const T & t): pos(pos), content(t), dfs_visited(false)
    {
    }

    Vector<2> pos;
    T content;
    bool dfs_visited;
  };

  template <typename T> class QuadTreeNode
  {
  public:

    class QuadTree;
    friend class QuadTree;

    QuadTreeNode(double delta,
                 const RobotVision::Rectangle & bbox)
                   : bbox(bbox), delta(delta), empty(true)
    {
    }

    bool deleteAtPos(const Vector<2> & pos)
    {
      if (children.size() == 0)
      {
        if(empty)
        {
          return false;
        }
        else
        {
          if(norm(elem.pos-pos) < delta)
          {
            empty = true;
            return true;
          }
          return false;
        }
      }
      else
      {
        assert(children.size()==4);
        typename std::list<QuadTreeNode<T> >::iterator it =  children.begin();


        if (it->bbox.contains(pos)){
          it->deleteAtPos(pos);
        }
        else
        {
          ++it;
          if (it->bbox.contains(pos)){
            it->deleteAtPos(pos);
          }
          else
          {
            ++it;
            if (it->bbox.contains(pos)){
              it->deleteAtPos(pos);
            }
            else
            {
              ++it;
              if (it->bbox.contains(pos)){
                it->deleteAtPos(pos);
              }
            }
          }
        }
        assert(false);
        return false;
      }

    }


    bool insert(const QuadTreeElement<T> & new_elem)
    {
      assert(bbox.contains(new_elem.pos));
      assert(isnan(new_elem.pos)==false);
      if (children.size() == 0)
      {
        if(empty)
        {
          elem = new_elem;
          empty = false;
        }
        else{


          if (norm(elem.pos-new_elem.pos)<delta){
            return false;
          }

          double x_diff = bbox.x2 - bbox.x1;
          double y_diff = bbox.y2 - bbox.y1;

          double x0 = bbox.x1;
          double x1 = bbox.x1 + x_diff*0.5;
          double x2 = bbox.x2;
          double y0 = bbox.y1;
          double y1 = bbox.y1 + y_diff*0.5;
          double y2 = bbox.y2;

          QuadTreeNode xy(delta,RobotVision::Rectangle(x0,y0,x1,y1));
          QuadTreeNode xY(delta,RobotVision::Rectangle(x0,y1,x1,y2));
          QuadTreeNode Xy(delta,RobotVision::Rectangle(x1,y0,x2,y1));
          QuadTreeNode XY(delta,RobotVision::Rectangle(x1,y1,x2,y2));

          insertMacro(xy,xY,Xy,XY,elem,bbox);
          insertMacro(xy,xY,Xy,XY,new_elem,bbox);


          children.push_back(xy);
          children.push_back(xY);
          children.push_back(Xy);
          children.push_back(XY);

        }
      }

      else{
        assert(children.size()==4);
        typename std::list<QuadTreeNode<T> >::iterator it =  children.begin();
        typename std::list<QuadTreeNode<T> >::iterator xy = it;
        ++it;
        typename std::list<QuadTreeNode<T> >::iterator xY = it;
        ++it;
        typename std::list<QuadTreeNode<T> >::iterator Xy = it;
        ++it;
        typename std::list<QuadTreeNode<T> >::iterator XY = it;

        insertMacro(*xy,*xY,*Xy,*XY,new_elem,bbox);
      }
      return true;

    }

    void query(const RobotVision::Rectangle & win,
               std::list<QuadTreeElement<T> > & content_list) const
    {

      if (children.size() == 0){
        if(empty == false){
          if (win.contains(elem.pos)){
            content_list.push_back(elem);
          }
        }
      }
      else{
        assert(children.size()==4);
        typename std::list<QuadTreeNode<T> >::const_iterator it
            =  children.begin();

        if (it->bbox.intersectsWith(win)){
          it->query(win,content_list);
        }
        ++it;

        if (it->bbox.intersectsWith(win)){
          it->query(win,content_list);
        }
        ++it;

        if (it->bbox.intersectsWith(win)){
          it->query(win,content_list);
        }
        ++it;

        if (it->bbox.intersectsWith(win)){
          it->query(win,content_list);
        }
      }
    }


    //This method is faster than query, if you only wanna know
    //whether there is any element within 'win'.
    bool isEmpty(const RobotVision::Rectangle & win) const
    {
      if (children.size() == 0)
      {
        if(empty)
        {
          return true;
        }
        else
        {
          if (win.contains(elem.pos))
          {
            return false;
          }
          return true;
        }
      }
      else
      {
        assert(children.size()==4);
        typename std::list<QuadTreeNode>::const_iterator it =  children.begin();

        if (it->bbox.intersectsWith(win))
        {
          if (!it->isEmpty(win))
            return false;
        }
        ++it;
        if (it->bbox.intersectsWith(win))
        {
          if (!it->isEmpty(win))
            return false;
        }
        ++it;
        if (it->bbox.intersectsWith(win))
        {
          if (!it->isEmpty(win))
            return false;
        }
        ++it;
        if (it->bbox.intersectsWith(win))
        {
          if (!it->isEmpty(win))
            return false;
        }

      }
      return true;
    }




    //Mainly for debugging:
    //Shows the whole tree -- all nodes (quads) and leafs (elems)
    void traverse(std::list<QuadTreeElement<T> > & elem_list,
                  std::list<RobotVision::Rectangle> & quad_list) const
    {
      quad_list.push_back(bbox);

      if (children.size() == 0)
      {
        if(empty == false)
        {
          elem_list.push_back(elem);
        }
      }
      else
      {
        assert(children.size()==4);
        typename std::list<QuadTreeNode>::const_iterator it =  children.begin();

        it->traverse(elem_list, quad_list);
        ++it;

        it->traverse(elem_list, quad_list);
        ++it;

        it->traverse(elem_list, quad_list);
        ++it;

        it->traverse(elem_list, quad_list);
      }
    }


    static void insertMacro(QuadTreeNode & xy,
                            QuadTreeNode & xY,
                            QuadTreeNode & Xy,
                            QuadTreeNode & XY,
                            const QuadTreeElement<T> & elem,
                            const RobotVision::Rectangle & bbox)
    {
      double x_diff = bbox.x2 - bbox.x1;
      double y_diff = bbox.y2 - bbox.y1;
      double rel_x = 1-(bbox.x2 -elem.pos[0])/x_diff;
      double rel_y = 1-(bbox.y2 -elem.pos[1])/y_diff;

      assert(rel_x >= 0);
      assert(rel_y >= 0);
      assert(rel_x <= 1);
      assert(rel_y <= 1);

      if(rel_x < 0.5 && rel_y < 0.5){
        xy.insert(elem);
      }
      else if (rel_x >= 0.5 && rel_y < 0.5){
        Xy.insert(elem);
      }
      else if (rel_x < 0.5 && rel_y >= 0.5){
        xY.insert(elem);
      }
      else if (rel_x >= 0.5 && rel_y >= 0.5){
        XY.insert(elem);
      }
      else{
        assert(false);
      }
    }

    RobotVision::Rectangle bbox;
    std::list<QuadTreeNode<T> > children;
    QuadTreeElement<T> elem;
    double delta ;
    bool empty;

  };


  template <typename T> class QuadTree
  {

  public:


  public:

    QuadTree(const RobotVision::Rectangle & bbox,
             double delta) : bbox(bbox)  , root(delta, bbox)
    {
      this->delta = delta;
    }


    bool  insert(const QuadTreeElement<T> & elem)
    {
      return root.insert(elem);
    }

    bool deleteAtPos( const Vector<2> & p)
    {
      return root.deleteAtPos(p);
    }

    void query(const RobotVision::Rectangle & win,
               std::list<QuadTreeElement<T> > & l) const
    {
      root.query(win,l);
    }

    bool isEmpty(const RobotVision::Rectangle & win) const
    {
      return root.isEmpty(win);
    }

    void traverse(std::list<QuadTreeElement<T> > & elem_list,
                  std::list<RobotVision::Rectangle> & quad_list) const
    {
      return root.traverse(elem_list, quad_list);
    }

  private:

    RobotVision::Rectangle bbox;
    QuadTreeNode<T>  root;
    double delta;

  };

}



#endif // RV_QUADTREE_H
