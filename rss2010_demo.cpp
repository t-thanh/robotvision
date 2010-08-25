/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * rss2010_demo.cpp is part of RobotVision.
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

#include <vector>
#include <set>
#include <list>

#include <TooN/TooN.h>
#include <TooN/SVD.h>
#include <TooN/se3.h>

#include <cvd/gl_helpers.h>
#include <cvd/random.h>

#include "stopwatch.h"
#include "bundle_adjuster.h"
#include "maths_utils.h"
#include "graph_optimizer.h"
#include "rectangle.h"
#include "gui_view.h"
#include "image.h"


using namespace std;
using namespace CVD;
using namespace TooN;
using namespace RobotVision;

const int NUM_PYR_LEVELS = 1;

typedef std::map<int, TooN::Vector<2> >  IdObsMap;


class StructObs
{
public:
  StructObs()
  {
  }
  StructObs(const TooN::Vector<3> &p3d, const IdObsMap & frameid_obs_map) :
      p3d(p3d), frameid_obs_map(frameid_obs_map)
  {
  }
  TooN::Vector<3> p3d;
  IdObsMap frameid_obs_map;
};

typedef std::map<int, StructObs> IdStructObsMap;

class DrawItems
{
public:
  DrawItems()
  {
    mode = 0;
  }

  TooN::SE3<> first_est_pose;
  std::vector<TooN::SE3<> > update_frames;

  int num_matched_points;
  int num_active_points;
  typedef std::list<std::pair<TooN::Vector<3>,TooN::Vector<3> > > Line3DList;
  typedef std::list<std::pair<TooN::Vector<2>,TooN::Vector<2> > > LineList;
  typedef std::list< Rectangle > BoxList;
  typedef std::list<std::pair<MonoImage,CVD::ImageRef> > PatchList;

  typedef std::list<std::pair<TooN::Vector<2>, double> > CircleList;
  typedef std::list<std::pair<TooN::Vector<3>, double> > Ball3DList;
  typedef std::vector<std::pair<TooN::Vector<2>, TooN::Matrix<2> > > GaussList;

  std::list<std::pair<CVD::Rgba<float>, LineList> > line_list[NUM_PYR_LEVELS];
  std::list<std::pair<CVD::Rgba<float>, BoxList> > box_list[NUM_PYR_LEVELS];
  std::list<std::pair<CVD::Rgba<float>, PatchList> > patch_list[NUM_PYR_LEVELS];
  std::list<std::pair<CVD::Rgba<float>, CircleList> >
      circle_list[NUM_PYR_LEVELS];
  std::list<std::pair<CVD::Rgba<float>, Ball3DList> > ball3d_list;
  std::list<std::pair<CVD::Rgba<float>, Line3DList> > line3d_list;
  std::list<std::pair<CVD::Rgba<float>, GaussList> > gauss_list[NUM_PYR_LEVELS];

  std::set<int> frame_id_set;

  int mode;

  void clear(){

    update_frames.clear();
    num_matched_points = 0;
    num_active_points = 0;

    ball3d_list.clear();
    line3d_list.clear();
    for (int i=0; i<NUM_PYR_LEVELS; i++)
    {
      patch_list[i].clear();
      circle_list[i].clear();
      line_list[i].clear();
      box_list[i].clear();
      gauss_list[i].clear();
    }
    frame_id_set.clear();
  }
};


void getBaParams(const vector<SE3<> > & keyframe_vec,
                 const IdStructObsMap & globalid_structobs_map,
                 vector<int> & frame_id_vec,
                 vector<SE3<> > & ba_frame_vec,
                 vector<Vector<3> > & ba_point_vec,
                 vector<IdObs<2> > & ba_obs_vec,
                 vector<int> & globalid_vec)
{
  int point_id = 0;
  uint num_frames = frame_id_vec.size();

  for (uint i=0; i<num_frames; ++i)
  {
    ba_frame_vec.push_back(keyframe_vec[frame_id_vec[i]]);
  }

  for (IdStructObsMap::const_iterator it=globalid_structobs_map.begin();
  it!=globalid_structobs_map.end(); ++it)
  {
    vector<IdObs<2> > cand_obs_vec;

    for (uint i=0; i<num_frames; ++i)
    {
      IdObsMap::const_iterator it2
          = it->second.frameid_obs_map.find(frame_id_vec[i]);

      if (it2!=it->second.frameid_obs_map.end())
      {
        cand_obs_vec.push_back(IdObs<2>(point_id,i,it2->second));
      }
    }

    if (cand_obs_vec.size()>1)
    {
      ba_point_vec.push_back(it->second.p3d);
      for (uint i=0; i<cand_obs_vec.size(); ++i)
      {
        ba_obs_vec.push_back(cand_obs_vec[i]);
      }
      globalid_vec.push_back(it->first);
      ++point_id;
    }
  }
}


void getObsWithID(const LinearCamera & cam_pars,
                  const SE3<> & pose,
                  const vector<Vector<3> > & point_vec,
                  const set<int> & target_id_set,
                  double noise,
                  vector<IdObs<2> > & matchobs_vec,
                  vector<IdObs<2> > & newobs_vec)
{
  SE3XYZ se3xyz(cam_pars);
  for (uint i =0; i<point_vec.size(); ++i)
  {
    double d = depth(pose,point_vec[i]);

    if (d>0 && d<2)
    {
      Vector<2> obs = se3xyz.map(pose,point_vec[i])
                      + noise*makeVector(rand_g(),rand_g());
      if(obs[0]>=0 && obs[1]>=0
         && obs[0]<=cam_pars.P()[0]*2-1 && obs[1] <= cam_pars.P()[1]*2-1)
      {
        if(target_id_set.find(i)!=target_id_set.end())

          matchobs_vec.push_back(IdObs<2>(i,0,obs));
        else
          newobs_vec.push_back(IdObs<2>(i,0,obs));
      }
    }
  }
}


void getObsWith2IDs(const LinearCamera & cam_pars,
                    const SE3<> & pose,
                    const vector<Vector<3> > & point_vec,
                    const set<int> & target_id_set1,
                    const set<int> & target_id_set2,
                    double noise,

                    vector<IdObs<2> > & matchobs_vec1,
                    vector<IdObs<2> > & matchobs_vec2,
                    vector<IdObs<2> > & newobs_vec)
{
  SE3XYZ se3xyz(cam_pars);
  for (uint i =0; i<point_vec.size(); ++i)
  {
    double d = depth(pose,point_vec[i]);

    if (d>0 && d<2)
    {
      Vector<2> obs = se3xyz.map(pose,point_vec[i])
                      + noise*makeVector(rand_g(),rand_g());
      if(obs[0]>=0 && obs[1]>=0
         && obs[0]<=cam_pars.P()[0]*2-1 && obs[1] <= cam_pars.P()[1]*2-1)
      {
        if(target_id_set1.find(i)!=target_id_set1.end())

          matchobs_vec1.push_back(IdObs<2> (i,0,obs));
        else if(target_id_set2.find(i)!=target_id_set2.end())

          matchobs_vec2.push_back(IdObs<2>(i,0,obs));
        else
          newobs_vec.push_back(IdObs<2>(i,0,obs));
      }
    }
  }
}


/** Inverse depth features*/
class UvqGauss
{
public:
  UvqGauss(){}
  UvqGauss( const TooN::Vector<3> & mean,
            const TooN::Matrix<3> & inv_cov,
            const Vector<2> & obs1)
              : mean(mean), inv_cov(inv_cov), obs1(obs1), num_updates(0){}

  TooN::Vector<3> mean;
  TooN::Matrix<3> inv_cov;

  Vector<2> obs1;
  Vector<2> obs2;

  int num_updates;
  int point_id;
};




void figure2()
{
  bool show_odom = true;
  bool show_sim3 = true;
  bool show_se3 = false;

  int image_width = 640;
  int image_height = 480;
  ImageRef frame_size= ImageRef(image_width,image_height);
  Vector<2> focal_length = makeVector(195, -195);
  Vector<2> principle_point = makeVector(162, 125);

  LinearCamera cam_pars(focal_length,principle_point,frame_size);

  const int WIDTH = 512;
  const int HEIGHT = 384;

  ImageRef win_size(WIDTH*2,HEIGHT*2);
  GuiWindow glwin(win_size);
  LinearCamera cam3d(makeVector(1000,-1000),makeVector(WIDTH,HEIGHT),win_size);
  GuiView view3d(ImageRef(WIDTH,HEIGHT), cam3d,
                 makeVector(-1.4582,0.0975581,-0.458294),
                 makeVector(-2.08342,5.54298,18 ));
  view3d.init(&glwin,Rectangle(0,0,1,1));

  ifstream fp_in;

  vector<SE3<> > pose_vec;
  vector<Sim3<> > sim_vec;
  vector<Vector<3> > point_vec;
  vector<IdObs<2>  > obs_vec;

  fp_in.open("../keble_college.txt", ios::in);

  char line[256];

  double scale = 0;
  uint loop_id = -1;

  while(fp_in.getline(line,256))
  {
    static int mode = -1;
    if (strcmp(line,"Poses")==0)
    {
      mode = 0;
      continue;
    }
    else if(strcmp(line,"Scale")==0)
    {
      mode = 1;
      continue;
    }
    else if(strcmp(line,"LoopId")==0)
    {
      mode = 2;
      continue;
    }
    else if(strcmp(line,"Points")==0)
    {
      mode = 3;
      continue;
    }
    else if(strcmp(line,"Obs")==0)
    {
      mode = 4;
      continue;
    }
    else
    {
      switch (mode)
      {
      case 0:
        {
          Matrix<3> R;
          Vector<3> t;

          R(0,0) = atof(strtok(line," "));
          R(0,1) = atof(strtok(NULL," "));
          R(0,2) = atof(strtok(NULL," "));
          t[0] = atof(strtok(NULL," \n"));
          fp_in.getline(line,256);
          R(1,0) = atof(strtok(line," "));
          R(1,1) = atof(strtok(NULL," "));
          R(1,2) = atof(strtok(NULL," "));
          t[1] = atof(strtok(NULL," \n"));
          fp_in.getline(line,256);
          R(2,0) = atof(strtok(line," "));
          R(2,1) = atof(strtok(NULL," "));
          R(2,2) = atof(strtok(NULL," "));
          t[2] = atof(strtok(NULL," \n"));
          fp_in.getline(line,256);

          static int i=0;
          if (i%1==0){
            pose_vec.push_back(SE3<>(R,t));
            sim_vec.push_back(Sim3<>(R,t,1.));
          }
          ++i;
        }
        break;
      case 1:
        {
          scale = atof(strtok(line," "));
          break;
        }
      case 2:
        {
          loop_id = atof(strtok(line," "));
          break;
        }
      case 3:
        {
          Vector<3> point;
          point[0] = atof(strtok(line," "));
          point[1] = atof(strtok(NULL," "));
          point[2] = atof(strtok(NULL," "));
          point_vec.push_back(point);
        }
        break;
      case 4:
        {
          Vector<2> obs;
          double point_id =  atof(strtok(line,":, "));
          double frame_id =  atof(strtok(NULL,":, "));
          obs[0] = atof(strtok(NULL,":, "));
          obs[1] = atof(strtok(NULL,":, "));
          obs_vec.push_back(IdObs<2> (point_id,frame_id,obs));
        }
      }
    }
  }
  fp_in.close();


  int trans_id = 0;

  vector<SE3<> > cor_pose_vec;

  vector<Sim3<> > trans7_list;
  vector<Sim3<> > updated_trans7_list;

  list<Constraint<SE3<>,6> > se3_list;
  list<Constraint<Sim3<>,7> > sim3_list;
  vector<SE3<> > updated_trans6_list;
  vector<SE3<> > trans6_list;


  trans7_list.push_back(sim_vec[loop_id]);
  trans6_list.push_back(pose_vec[loop_id]);


  //compute relative constraints
  for (uint i=loop_id+1; i<pose_vec.size()-1; ++i)
  {
    trans7_list.push_back(sim_vec[i]);
    trans6_list.push_back(pose_vec[i]);


    Matrix<6,6> inf6 = TooN::Identity;
    Matrix<7,7> inf7 = TooN::Identity;


    se3_list.push_back(Constraint<SE3<>,6 >(trans_id,
                                            trans_id+1,
                                            pose_vec[i]*pose_vec[i-1].inverse(),
                                            inf6));
    sim3_list.push_back(Constraint<Sim3<>,7 >(trans_id,
                                              trans_id+1,
                                              sim_vec[i]*sim_vec[i-1].inverse(),
                                              inf7));
    ++trans_id;
  }

  int last_id = pose_vec.size()-1;

  Sim3<> loop_constraint = sim_vec[loop_id]*sim_vec[last_id].inverse();
  loop_constraint.get_scale() =  1./scale;

  Matrix<7,7> inf7 = Identity(7);
  sim3_list.push_back(Constraint<Sim3<>,7 >(trans_id,0,loop_constraint,inf7));

  list<Constraint<Sim3<>,7> >::iterator elem1 = sim3_list.begin();
  elem1->fisher_information = inf7;

  Matrix<6,6> inf6 = Identity(6);
  se3_list.push_back(
      Constraint<SE3<>,6 >(trans_id,0,
                           pose_vec[loop_id]*pose_vec[last_id].inverse(),
                           inf6));

  list<Constraint<SE3<>,6> >::iterator el1 = se3_list.begin();
  el1->fisher_information = inf6;







  SE3ConFun se3confun;

  GraphOptimizer<Sim3<>,7> opt7;
  GraphOptimizer<SE3<>,6> opt6;

  updated_trans6_list = trans6_list;

  opt6.verbose = false;

  //Perform standard 6 DoF optimisation
  opt6.optimize(updated_trans6_list,
                se3_list,
                se3confun,
                1,
                2,
                0.0000000000000000001);

  Sim3ConFun sim3confun;
  updated_trans7_list = trans7_list;
  StopWatch go_m;
  go_m.start();
  opt7.verbose = 1;
  cout << "7 DoF optimisation:" << endl;
  //Perform 7 DoF optimisation
  opt7.optimize(updated_trans7_list,
                sim3_list,
                sim3confun,
                1,
                2,
                0.0000000000000000001);
  go_m.stop();
  cout << "time in s: " << go_m.getStoppedTime() << endl<<endl;



  vector<SE3<> > cor7_pose_vec;
  for (uint i=0;i<loop_id; ++i)
  {
    cor7_pose_vec.push_back(pose_vec[i]);
  }
  // transform Sim3 --> SE3
  for (uint i_f=0; i_f<updated_trans7_list.size(); ++i_f)
  {
    cor7_pose_vec.push_back(
        SE3<>(updated_trans7_list[i_f].get_rotation(),
              updated_trans7_list[i_f].get_translation()
              /updated_trans7_list[i_f].get_scale()));
  }

  for (uint i=0;i<loop_id; ++i)
  {
    cor_pose_vec.push_back(pose_vec[i]);
  }
  for (uint i_f=0; i_f<updated_trans6_list.size(); ++i_f)
  {
    cor_pose_vec.push_back(updated_trans6_list[i_f]);
  }

  vector<Vector<3> > cor_point_vec;
  vector<Vector<3> > cor7_point_vec;

  SE3XYZ se3xyz(cam_pars);
  RobotVision::BA_SE3_XYZ ba;


  // map points into updated frames
  cor7_point_vec = point_vec;
  for (uint i=0; i<obs_vec.size(); ++i)
  {
    int frame_id = obs_vec[i].frame_id;
    int point_id = obs_vec[i].point_id;
    int id = obs_vec[i].frame_id-loop_id;
    if (id>=0)
    {
      Vector<3> rel_point = transform(sim_vec[frame_id], point_vec[point_id]);
      Vector<3> cor_point
          = transform(updated_trans7_list[id].inverse(),rel_point);


      cor7_point_vec[point_id] = cor_point;
    }

  }

  // map points into updated frames
  cor_point_vec = point_vec;

  for (uint i=0; i<obs_vec.size(); ++i)
  {
    int frame_id = obs_vec[i].frame_id;
    int point_id = obs_vec[i].point_id;
    int id = obs_vec[i].frame_id-loop_id;
    if (id>=0)
    {
      Vector<3> rel_point = transform(pose_vec[frame_id], point_vec[point_id]);
      Vector<3> cor_point
          = transform(updated_trans6_list[id].inverse(),rel_point);
      cor_point_vec[point_id] = cor_point;
    }

  }





  BA_SE3_XYZ::_TrackMap track_map;
  for (uint i=0; i<obs_vec.size();++i)
  {
    IdObs<2> & id_obs = obs_vec[i];

    BA_SE3_XYZ::_TrackMap::iterator it
        = track_map.find(id_obs.point_id);
    if (it==track_map.end())
    {
      list<IdObs<2> > obs_list;
      obs_list.push_back(id_obs);
      track_map.insert(make_pair(id_obs.point_id,obs_list));
    }
    else
    {
      it->second.push_back(id_obs);
    }
  }



  ba.calcFastStructureOnly(cor_pose_vec,
              cor_point_vec,
              se3xyz,
              track_map,
              BundleAdjusterParams(true,1,10));

  ba.verbose = true;
  StopWatch ba_m;
  ba_m.start();
  cout << "Structure-Only BA:" << endl;
  ba.calcFastStructureOnly(cor7_pose_vec,
              cor7_point_vec,
              se3xyz,
              track_map,
              BundleAdjusterParams(true,1,10));
  ba_m.stop();
  cout << "time in s: " << ba_m.getStoppedTime()  << endl << endl;;



  SE3<> cur_pose = pose_vec[pose_vec.size()-1];

  while(true)
  {
    glClearColor(1,1,1,1);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    view3d.activate3D();

    glColor3f(0,1,0);
    view3d.drawPose3D(cur_pose);

    if (show_odom)
    {
      glColor4f(0,0,0,1);

      for (uint i=0; i<point_vec.size(); ++i)
      {
        view3d.drawBall3D(point_vec[i],0.01);
      }
      for (uint i=0; i<pose_vec.size(); ++i)
      {
        view3d.drawPose3D(pose_vec[i]);

        if (i<pose_vec.size()-1)
          view3d.drawLine3D(trans2center(pose_vec[i]),
                            trans2center(pose_vec[(i+1)%pose_vec.size()]));
      }
    }

    if (show_se3)
    {
      glColor3f(0,0.5,0);
      for (uint i=0; i<updated_trans6_list.size(); ++i)

      {
        view3d.drawPose3D(updated_trans6_list[i]);

        view3d.drawLine3D(trans2center(updated_trans6_list[i]),
                          trans2center(
                              updated_trans6_list[(i+1)
                                                  %updated_trans6_list.size()])
                          );
      }

      for (uint i=0; i<cor_point_vec.size(); ++i)
      {
        view3d.drawBall3D(cor_point_vec[i],0.01);
      }
    }


    if (show_sim3)
    {
      glColor4f(1,0,0,1);

      for (uint i=0; i<updated_trans7_list.size(); ++i)
      {
        view3d.drawSimilarity3D(updated_trans7_list[i]);
        view3d.drawLine3D(trans2center(updated_trans7_list[i]),
                          trans2center(
                              updated_trans7_list[(i+1)
                                                  %updated_trans7_list.size()])
                          );
      }
      for (uint i=0; i<cor7_point_vec.size(); ++i)
      {
        view3d.drawBall3D(cor7_point_vec[i],0.01);
      }
    }
    glFlush();
    glwin.swap_buffers();
  }
  exit(0);
}

class Figure3Data
{
public:
  Figure3Data(double noise,
              double mean_rmse_se3,
              double var_rmse_se3,
              double mean_rmse_sim3,
              double var_rmse_sim3,
              double av_scale)
                : noise(noise),
                mean_rmse_se3(mean_rmse_se3),
                var_rmse_se3(var_rmse_se3),
                mean_rmse_sim3(mean_rmse_sim3),
                var_rmse_sim3(var_rmse_sim3),
                av_scale(av_scale) {}

  double noise;
  double mean_rmse_se3;
  double var_rmse_se3;
  double mean_rmse_sim3;
  double var_rmse_sim3;
  double av_scale;
};


void figure3()
{
  bool draw = true;

  bool draw_odom = true;
  bool draw_sim3 = true;
  bool draw_se3 = true;
  bool draw_ground_truth = true;

  vector<Figure3Data> fig_data;

  const int C_WIDTH = 320;
  const int C_HEIGHT = 240;

  const int WIDTH = 640;
  const int HEIGHT = 480;
  ImageRef win_size(WIDTH,HEIGHT);
  GuiWindow glwin(win_size);

  LinearCamera cam_pars(makeVector(195, -195),
                        makeVector(162 ,125),
                        ImageRef(320,240));
  LinearCamera cam3d(makeVector(1000,-1000),
                     makeVector(WIDTH/2,HEIGHT/2),
                     ImageRef(640,480));


  GuiView view3d(ImageRef(WIDTH,HEIGHT),
                 cam3d,makeVector(0,0,0),
                 makeVector(10,0,50));
  view3d.init(&glwin,Rectangle(0,0,1,1));

  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(4);


  for (double noise = 1.2; noise>=0.0; noise-=0.2)
  {

    double sum_scale=0;
    double sum_rmse_se3=0;
    double sum_rmse_sim3=0;
    double sum_scale_sq =0;
    double sum_rmse_se3_sq=0;
    double sum_rmse_sim3_sq=0;
    uint num_trials = 10;

    srand(1);

    for (uint trial=0; trial<num_trials; ++trial)
    {
      int frame_id = 2;
      bool end = false;

      vector<SE3<> > true_poses;
      vector<Vector<3> > true_points;

      IdStructObsMap globalid_structobs_map;

      Vector<2> bias = makeVector(0,0);

      int factor = 10;

      int num_frames = 72*factor ;

      int tmp=0;

      //generate circular trajectory of poses
      for (int i =0; i<num_frames; ++i)
      {
        double x = factor*cos(2*M_PI*i/num_frames) - factor;
        double y = factor*sin(2*M_PI*i/num_frames);
        true_poses.push_back(
            SE3<>(
                SO3<>(makeVector(0,0,2*M_PI*i/num_frames))
                *SO3<>(makeVector(0,M_PI_2,0)), makeVector(x,y,0)).inverse());
        tmp++;
      }
      num_frames = tmp;

      //generate points on a ring-shape distribution
      uint num_points = 5000;
      for (uint i =0; i<num_points; ++i)
      {
        double theta = rand_u();
        double r = factor+1;
        double x = r*cos(2*M_PI*theta) - factor;
        double y = r*sin(2*M_PI*theta);
        double z = rand_u()*2-1;
        true_points.push_back(makeVector(x,y,z));
      }

      vector<Vector<3> > vis_points;
      vector<Vector<3> > vis_points0;
      vector<Vector<3> > vis_points1;

      SE3<> cur_pose = true_poses[0];

      map<int, UvqGauss > new_uvq_map;
      map<int, UvqGauss > uvq_map;

      BundleAdjusterParams ba_pars_f(false,0,1,1);
      BundleAdjusterParams ba_pars_k(false,0,10);

      BA_SE3_XYZ ba;
      SE3UVQ se3uvq(cam_pars);
      SE3XYZ se3xyz(cam_pars);
      DrawItems::Line3DList line_list2;

      for (uint i =0; i<num_points; ++i)
      {
        if (depth(cur_pose,true_points[i])>0)
        {
          Vector<2> obs = se3xyz.map(true_poses[0],true_points[i]);
          if(obs[0]>=0 && obs[1]>=0 && obs[0]<=639 && obs[1] <= C_HEIGHT-1)
          {
            Matrix<3> Lambda = Zeros(3);
            Lambda(0,0) = Po2(cam_pars.F()[0]);
            Lambda(1,1) = Po2(cam_pars.F()[1]);
            Lambda(2,2) = 0;
            Vector<3> p;
            p.slice<0,2>() = cam_pars.unmap(obs);
            p[2] = 0.9;
            uvq_map.insert(make_pair(i,UvqGauss(p,Lambda,obs)));
          }
        }
      }
      cur_pose = true_poses[1];
      SE3<> ba_frame = true_poses[1]*true_poses[0].inverse();
      set<int> first_frame_point_id_set;
      map<int,int> loop_closure_map;

      for (uint i =0; i<num_points; ++i)
      {
        if (depth(cur_pose,true_points[i])>0)
        {
          Vector<2> obs = se3xyz.map(cur_pose,true_points[i]);
          if(obs[0]>=0 && obs[1]>=0
             && obs[0]<=C_WIDTH-1 && obs[1] <= C_HEIGHT-1)
          {
            map<int,UvqGauss>::iterator it = uvq_map.find(i);
            if (it!=uvq_map.end())
            {
              ba.filterSingleFeatureOnly(ba_frame,
                                         it->second.mean,
                                         it->second.inv_cov,
                                         se3uvq,
                                         obs,
                                         ba_pars_f);
              it->second.obs2 = obs;
              ++(it->second.num_updates);
              Cholesky<3> Ch(it->second.inv_cov);
              Matrix<3> Sigma = Ch.get_inverse();
              double q1 = it->second.mean[2]+sqrt(Sigma(2,2))*3;
              double q2 = max(0.0000000001,
                              it->second.mean[2]-sqrt(Sigma(2,2))*3);
              Vector<3> p1
                  = transform(true_poses.at(0).inverse(),
                              makeVector(it->second.mean[0]/q1,
                                         it->second.mean[1]/q1,
                                         1./q1)  ) ;

              Vector<3> p2
                  = transform(true_poses.at(0).inverse(),
                              makeVector(it->second.mean[0]/q2,
                                         it->second.mean[1]/q2,
                                         1./q2)  );
              line_list2.push_back(make_pair(p1,p2));

              vis_points.push_back(true_points[i]);
            }
            else
            {
              Matrix<3> Lambda = Zeros(3);
              Lambda(0,0) = Po2(cam_pars.F()[0]);
              Lambda(1,1) = Po2(cam_pars.F()[1]);
              Lambda(2,2) = 0;
              Vector<3> p;
              p.slice<0,2>() = cam_pars.unmap(obs);
              p[2] = 0.9;
              new_uvq_map.insert(make_pair(i,UvqGauss(p,Lambda,obs)));
            }
          }
        }
      }
      vector<SE3<> > ba_poses;
      vector<SE3<>  > keyframe_vec;
      ba_poses.push_back(true_poses[0]);
      ba_poses.push_back(true_poses[1]);
      vector<IdObs<2> > ba_obs;
      vector<Vector<3> > ba_points;
      vector<int> id_vec;

      SE3<>  kf1;
      SE3<>  kf2;
      kf1 = ba_poses[0];
      kf2 = ba_poses[1];

      int point_id = 0;
      for (map<int,UvqGauss>::const_iterator it=uvq_map.begin();
      it!=uvq_map.end();
      ++it)
      {
        if (it->second.num_updates<1)
          continue;
        ba_obs.push_back(IdObs<2>(point_id,0,it->second.obs1));
        ba_obs.push_back(IdObs<2>(point_id,1,it->second.obs2));

        Vector<3>  xyz
            = transform(true_poses[0].inverse(),
                        makeVector(it->second.mean[0],
                                   it->second.mean[1],1)/it->second.mean[2]);
        IdObsMap id_obs_map;
        id_obs_map.insert(make_pair(0,it->second.obs1));
        id_obs_map.insert(make_pair(1,it->second.obs2));
        globalid_structobs_map.insert(make_pair(it->first,
                                                StructObs(xyz,id_obs_map)));
        id_vec.push_back(it->first);
        ba_points.push_back(xyz);
        ++point_id;
      }

      BA_SE3_XYZ::_TrackMap track_map;
      for (uint i_obs=0; i_obs<ba_obs.size();++i_obs)
      {
        IdObs<2> & id_obs = ba_obs[i_obs];

        BA_SE3_XYZ::_TrackMap::iterator it
            = track_map.find(id_obs.point_id);
        if (it==track_map.end())
        {
          BA_SE3_XYZ::_Track obs_list;
          obs_list.push_back( id_obs);
          track_map.insert(make_pair(id_obs.point_id,obs_list));
        }
        else
        {
          it->second.push_back(id_obs);
        }
      }
      ba.calcFast(ba_poses,ba_points,se3xyz,track_map,1,ba_pars_k);
      //ba.calcFull(ba_poses,ba_points,se3xyz,ba_obs,1,0,ba_pars_k,false);

      for (uint i=0; i<id_vec.size(); ++i)
      {
        globalid_structobs_map.find(id_vec[i])->second.p3d = ba_points[i];

        first_frame_point_id_set.insert(id_vec[i]);
      }

      kf1 = ba_poses[0];
      kf2 = ba_poses[1];
      keyframe_vec.push_back(kf1);
      keyframe_vec.push_back(kf2);

      uvq_map = new_uvq_map;
      SE3<> loop_pose = kf1;
      int new_id = num_points;

      vector<Sim3<> > trans7_list;
      vector<Vector<3> > cor_point_vec;
      vector<Vector<3> > cor7_point_vec;

      vector<SE3<> > trans6_list;
      vector<SE3<> > cor7_pose_vec;

      StopWatch wait_in_the_end;
      while (!end || wait_in_the_end.readCurrentTime() < 5)
      {
        if (!end)
        {
          line_list2.clear();
          vis_points.clear();
          vector<IdObs<2> > obs_vec;
          vector<Vector<3> > point_vec;
          int point_id = 0;
          SE3<> origin = keyframe_vec[keyframe_vec.size()-1];
          SE3<>  kf;
          new_uvq_map.clear();
          int obs_count = 0;
          for (uint i =0; i<num_points; ++i)
          {
            int id =i;
            set<int>::const_iterator loop_it = first_frame_point_id_set.find(i);

            if (frame_id>10 && loop_it!=first_frame_point_id_set.end())
              //loop closure detected!
            {
              map<int,int>::iterator it = loop_closure_map.find(i);
              if (it!=loop_closure_map.end())
              {
                id = it->second;
              }
              else
              {
                id = new_id;
                ++new_id;
                loop_closure_map.insert(make_pair(i,id));
              }
            }

            Vector<2> obs = se3xyz.map(true_poses.at(frame_id),true_points[i])
                            + noise*makeVector(rand_g(),rand_g()) + bias;
            if (depth(true_poses.at(frame_id),true_points[i])>0)
            {
              if(obs[0]>=0 && obs[1]>=0
                 && obs[0]<=C_WIDTH-1 && obs[1] <= C_HEIGHT-1)
              {
                ++obs_count;
                IdStructObsMap::iterator it = globalid_structobs_map.find(id);
                if (it!=globalid_structobs_map.end())
                  // current feature is in global map?
                {
                  obs_vec.push_back(IdObs<2>(point_id,0,obs));
                  point_vec.push_back(it->second.p3d);
                  ++point_id;
                  it->second.frameid_obs_map
                      .insert(make_pair(keyframe_vec.size(),obs));

                }
                else
                {
                  map<int,UvqGauss>::iterator it = uvq_map.find(id);
                  if (it!=uvq_map.end()) //current feature is in track map
                  {
                    it->second.obs2 = obs;
                    ++(it->second.num_updates);
                    it->second.point_id = id;
                  }
                  else
                  {
                    Matrix<3> Lambda = Zeros(3);
                    Lambda(0,0) = Po2(cam_pars.F()[0]);
                    Lambda(1,1) = Po2(cam_pars.F()[1]);
                    Lambda(2,2) = 0;
                    Vector<3> p;
                    p.slice<0,2>() = cam_pars.unmap(obs);
                    p[2] = 0.9;
                    new_uvq_map.insert(make_pair(id,UvqGauss(p,Lambda,obs)));
                  }
                }
              }
            }
          }

          SE3<>  new_pose = origin;
          std::list<IdObs<2> > obs_list = std::list<IdObs<2> >(obs_vec.begin(),
                                                               obs_vec.end());

          ba.calcFastMotionOnly(new_pose,point_vec,se3xyz,obs_list,ba_pars_k);
          kf = new_pose;
          cur_pose = new_pose;
          ba_frame = cur_pose*origin.inverse();

          double scale_sum=0;
          int scale_num=0;
          for (map<int,UvqGauss>::iterator it=uvq_map.begin();
          it!=uvq_map.end();
          ++it)
          {
            if (it->second.num_updates<1)
              continue;
            ba.filterSingleFeatureOnly(ba_frame,
                                       it->second.mean,
                                       it->second.inv_cov,
                                       se3uvq,
                                       it->second.obs2,
                                       ba_pars_f);
            Vector<3>  xyz
                = transform(origin.inverse(),
                            makeVector(it->second.mean[0],
                                       it->second.mean[1],1)
                            /it->second.mean[2]);
            IdObsMap frameid_obs_map;
            frameid_obs_map
                .insert(make_pair(keyframe_vec.size()-1,it->second.obs1));

            frameid_obs_map
                .insert(make_pair(keyframe_vec.size(),it->second.obs2));

            scale_sum += 1/it->second.mean[2];
            ++scale_num;
            globalid_structobs_map
                .insert(make_pair(it->second.point_id,
                                  StructObs(xyz,frameid_obs_map)));

          }

          point_id =0;
          point_vec.clear();
          ba_obs.clear();
          ba_poses.clear();
          id_vec.clear();
          keyframe_vec.push_back(kf);
          {
            vector<int> frame_id_vec;
            uint num_frames_to_optimise = min(frame_id+1,10);
            for (uint i=0; i<num_frames_to_optimise; ++i)
              frame_id_vec
                  .push_back(keyframe_vec.size()-num_frames_to_optimise+i);

            vector<Vector<3> > ba_point_vec;
            vector<SE3<> > ba_frame_vec;
            vector<IdObs<2> > ba_obs_vec;
            vector<int> globalid_vec;

            getBaParams(keyframe_vec,
                        globalid_structobs_map,
                        frame_id_vec,
                        ba_frame_vec,
                        ba_point_vec,
                        ba_obs_vec,
                        globalid_vec);
            BA_SE3_XYZ::_TrackMap track_map;
            for (uint i_obs=0; i_obs<ba_obs_vec.size();++i_obs)
            {
              IdObs<2> & id_obs = ba_obs_vec[i_obs];

              BA_SE3_XYZ::_TrackMap::iterator it
                  = track_map.find(id_obs.point_id);
              if (it==track_map.end())
              {
                list<IdObs<2> > obs_list;
                obs_list.push_back(id_obs);
                track_map.insert(make_pair(id_obs.point_id,obs_list));
              }
              else
              {
                it->second.push_back(id_obs);
              }
            }
            ba.calcFast(ba_frame_vec,ba_point_vec,se3xyz,track_map,
                        2,
                        ba_pars_k,false);


            for (uint i=0; i<ba_frame_vec.size(); ++i)
            {
              keyframe_vec.at(frame_id_vec.at(i))= ba_frame_vec.at(i);
            }
            for (uint i=0; i<globalid_vec.size(); ++i)
            {
              globalid_structobs_map
                  .find(globalid_vec.at(i))->second.p3d = ba_point_vec.at(i);
            }
          }

          uvq_map = new_uvq_map;

          if ((uint)frame_id>=true_poses.size()-1)
          {
            set<double> scale_set;
            loop_pose = true_poses[frame_id];
            SE3<> cur_pose = keyframe_vec[keyframe_vec.size()-1];
            for (map<int,int>::iterator it = loop_closure_map.begin();
            it!=loop_closure_map.end();
            ++it)
            {
              int id0 =  it->first;
              int id1 =  it->second;

              IdStructObsMap::iterator it_id0
                  = globalid_structobs_map.find(id0);

              IdStructObsMap::iterator it_id1
                  = globalid_structobs_map.find(id1);

              if (it_id0==globalid_structobs_map.end())
              {
                exit(0);
              }

              if (it_id1==globalid_structobs_map.end())
              {
                continue;
              }
              if (it_id1->second.frameid_obs_map.find(frame_id)
                == it_id1->second.frameid_obs_map.end())
                {
                continue;
              }
              Vector<2> obs
                  = it_id1->second.frameid_obs_map.find(frame_id)->second;

              double norm0 =  norm(transform(loop_pose,true_points.at(id0)));
              double norm1 =  norm(transform(cur_pose,it_id1->second.p3d));

              it_id0->second.frameid_obs_map.insert(make_pair(frame_id,obs));

              scale_set.insert(norm1/norm0);

              vis_points0.push_back(true_points[id0]);
              vis_points1
                  .push_back(globalid_structobs_map.find(id1)->second.p3d);
            }
            SE3<> rel_pose = cur_pose * loop_pose.inverse();
            double median_scale  = median(scale_set);

            vector<Vector<3> > ba_point_vec;
            vector<SE3<> > ba_frame_vec;
            vector<IdObs<2> > ba_obs_vec;
            vector<int> globalid_vec;
            vector<int> frame_id_vec;

            for (uint i=0; i<keyframe_vec.size(); ++i)
              frame_id_vec.push_back(i);

            getBaParams(keyframe_vec,
                        globalid_structobs_map,
                        frame_id_vec,
                        ba_frame_vec,
                        ba_point_vec,
                        ba_obs_vec,
                        globalid_vec);

            list<Constraint<SE3<>,6> > se3_list;
            list<Constraint<Sim3<>,7> > sim3_list;
            list<Constraint<Sim3<>,7> > sim3_id_list;

            vector<SE3<> > se3_vec;
            vector<Sim3<> > sim3_vec;

            for (uint i=0; i<keyframe_vec.size(); ++i)
            {
              SE3<> & pose = keyframe_vec[i];
              Sim3<>  sim(pose.get_rotation(),pose.get_translation(),1.);
              se3_vec.push_back(pose);
              sim3_vec.push_back(sim);
            }

            trans7_list.push_back(sim3_vec[0]);
            trans6_list.push_back(se3_vec[0]);

            int trans_id = 0;

            for (uint i=1; i<keyframe_vec.size(); ++i)
            {
              trans7_list.push_back(sim3_vec[i]);
              trans6_list.push_back(se3_vec[i]);

              Matrix<6,6> inf6 = TooN::Identity;
              Matrix<7,7> inf7 = TooN::Identity;

              se3_list.push_back(
                  Constraint<SE3<>,6 >(trans_id,
                                       trans_id+1,
                                       se3_vec[i]*se3_vec[i-1].inverse(),
                                       inf6));

              sim3_list.push_back(
                  Constraint<Sim3<>,7 >(trans_id,
                                        trans_id+1,
                                        sim3_vec[i]*sim3_vec[i-1].inverse(),
                                        inf7));
              sim3_id_list.push_back(
                  Constraint<Sim3<>,7 >(trans_id,
                                        trans_id+1,
                                        sim3_vec[i]*sim3_vec[i-1].inverse(),
                                        Identity(7)));

              ++trans_id;
            }

            rel_pose
                =  true_poses[true_poses.size()-1]
                   *true_poses[true_poses.size()-2].inverse();

            Sim3<> loop_constraint(rel_pose.get_rotation(),
                                   rel_pose.get_translation(),
                                   1./median_scale);

            Matrix<7,7> inf7 =  Identity(7);

            sim3_list.push_back(
                Constraint<Sim3<>,7 >(trans_id,0,loop_constraint,inf7));

            sim3_id_list.push_back(
                Constraint<Sim3<>,7 >(trans_id,0,loop_constraint,Identity(7)));

            list<Constraint<Sim3<>,7> >::iterator elem1 = sim3_list.begin();
            elem1->fisher_information = inf7;

            Matrix<6,6> inf6 = Identity(6);
            se3_list.push_back(Constraint<SE3<>,6 >(trans_id,0,rel_pose,inf6));

            list<Constraint<SE3<>,6> >::iterator el1 = se3_list.begin();
            el1->fisher_information = inf6;

            list<Constraint<Sim3<>,7> >::const_iterator it = sim3_list.begin();

            GraphOptimizer<Sim3<>,7> opt7;
            GraphOptimizer<SE3<>,6> opt6;

            SO3xR3ConFun se3_confun;
            Sim3ConFun sim3_confun;

            opt7.optimize(trans7_list,
                          sim3_list,
                          sim3_confun,
                          1,
                          3,
                          0.0000000000000000001);

            opt6.optimize(trans6_list,
                          se3_list,
                          se3_confun,
                          1,
                          3,
                          0.0000000000000000001);

            for (vector<IdObs<2> >::const_iterator it=ba_obs_vec.begin();
            it!=ba_obs_vec.end();
            ++it)
            {
              int frame_id = it->frame_id;
              int point_id = it->point_id;

              if ((uint)frame_id==ba_frame_vec.size()-1)
                continue;
              Vector<3> rel_point = transform(ba_frame_vec.at(frame_id),
                                              ba_point_vec.at(point_id));
              Vector<3> cor_point
                  = transform(trans6_list.at(frame_id).inverse(),rel_point);
              cor_point_vec.push_back(cor_point);
            }

            for (vector<IdObs<2> >::const_iterator it=ba_obs_vec.begin();
            it!=ba_obs_vec.end();
            ++it)
            {
              int frame_id = it->frame_id;
              int point_id = it->point_id;
              if ((uint)frame_id==ba_frame_vec.size()-1)
                continue;
              Vector<3> rel_point
                  = transform(sim3_vec[frame_id], ba_point_vec[point_id]);

              Vector<3> cor_point
                  = transform(trans7_list[frame_id].inverse(),rel_point);

              cor7_point_vec.push_back(cor_point);
            }

            for (uint i_f=0; i_f<trans7_list.size(); ++i_f)
            {
              cor7_pose_vec.push_back(SE3<>(trans7_list[i_f].get_rotation(),
                                            trans7_list[i_f].get_translation()
                                            /trans7_list[i_f].get_scale()));
            }

            SE3CompareModScale t;
            double s=1;
            double  v1 = t.optimize(true_poses,cor7_pose_vec,s,10);

            for (uint i_f=0; i_f<cor7_pose_vec.size(); ++i_f)
            {
              cor7_pose_vec[i_f].get_translation()
                  = s*cor7_pose_vec[i_f].get_translation();
            }

            for (uint i=0; i<cor7_point_vec.size(); ++i)
            {
              cor7_point_vec[i] *= s;
            }

            s=1;
            double v2 = t.optimize(true_poses,trans6_list,s,10);
            for (uint i_f=0; i_f<trans6_list.size(); ++i_f)
            {
              trans6_list[i_f].get_translation()
                  = s*trans6_list[i_f].get_translation();
            }

            for (uint i=0; i<cor_point_vec.size(); ++i)
            {
              cor_point_vec[i] *= s;
            }

            if (median_scale<1.)
              median_scale = 1./median_scale;

            double scale_drift = median_scale-1;
            double rmse_sim3 = v1/num_frames;
            double rmse_se3 = v2/num_frames;

            sum_rmse_se3 += rmse_se3;
            sum_rmse_sim3 += rmse_sim3;
            sum_scale += scale_drift;

            sum_rmse_se3_sq += Po2(rmse_se3);
            sum_rmse_sim3_sq += Po2(rmse_sim3);
            sum_scale_sq += Po2(scale_drift);


            cout << "Trial: " << trial << "  Noise: " << noise << " ";
            cout << "Se3 error/Sim3 error/Scale drift: " <<
                rmse_se3 << "/" << rmse_sim3 << "/" <<  scale_drift << endl;

            end = true;
            wait_in_the_end.start();
          }
          ++frame_id;
        }

        if (draw)
        {
          glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
          glClearColor(1,1,1,1);
          view3d.activate3D();

          if (draw_ground_truth)
          {
            glColor3f(0.75,0.75,0.75);
            for (uint i=0; i<true_poses.size(); ++i)
            {
              view3d.drawLine3D(trans2center(true_poses[i]),
                                trans2center(true_poses[(i+1)
                                                        %true_poses.size()]));
            }
            for (uint i=0; i<true_points.size(); ++i)
            {
              view3d.drawBall3D(true_points[i],0.03);
            }
          }

          if (!end)
          {
            glColor3f(1,0,0);
            for (uint i=0; i<vis_points.size(); ++i)
            {
              view3d.drawBall3D(vis_points[i],0.04);
            }

            glColor3f(0,0,1);
            for (DrawItems::Line3DList::iterator it = line_list2.begin();
            it!=line_list2.end(); ++it)
            {
              view3d.drawLine3D(it->first,it->second);
            }

            glColor3f(0,1,0);
            view3d.drawPose3D(cur_pose);

            glColor3f(0,0,0);
            for (IdStructObsMap::iterator it=globalid_structobs_map.begin();
            it!=globalid_structobs_map.end(); ++it)
            {
              view3d.drawBall3D(it->second.p3d,0.03);
            }
            glColor3f(0,0,0);
            for (uint i=0; i<keyframe_vec.size(); ++i)
            {
              view3d.drawPose3D(keyframe_vec[i]);
            }
          }
          if (end && draw_odom)
          {
            glColor3f(1,0.5,0);


            view3d.drawPose3D(loop_pose);

            view3d.drawLine3D(trans2center(keyframe_vec[keyframe_vec.size()-1]),
                              trans2center(loop_pose));

            glColor3f(0,0,0);
            for (IdStructObsMap::iterator it=globalid_structobs_map.begin();
            it!=globalid_structobs_map.end(); ++it)
            {
              view3d.drawBall3D(it->second.p3d,0.03);
            }
            glColor3f(0,0,0);
            for (uint i=0; i<keyframe_vec.size(); ++i)
            {
              view3d.drawPose3D(keyframe_vec[i]);
            }
          }

          if (end && draw_sim3)
          {
            glColor3f(0.75,0.0,0);
            for (uint i=0; i<trans7_list.size(); ++i)
            {
              view3d.drawPose3D(cor7_pose_vec[i]);
            }
            for (uint i=0; i<cor7_point_vec.size(); ++i)
            {
              view3d.drawBall3D(cor7_point_vec[i],0.03);
            }
          }

          if (end && draw_se3)
          {
            glColor3f(0,0.75,0);
            for (uint i=0; i<trans6_list.size(); ++i)
            {
              view3d.drawPose3D(trans6_list[i]);
            }
            for (uint i=0; i<cor_point_vec.size(); ++i)
            {
              view3d.drawBall3D(cor_point_vec[i],0.03);
            }
          }
          glFlush();
          glwin.swap_buffers();
        }
      }
    }
    Figure3Data fd(noise,
                   sum_rmse_se3/num_trials,
                   sum_rmse_se3_sq/num_trials - Po2(sum_rmse_se3/num_trials),
                   sum_rmse_sim3/num_trials,
                   sum_rmse_sim3_sq/num_trials - Po2(sum_rmse_sim3/num_trials),
                   sum_scale/num_trials);
    fig_data.push_back(fd);

    cout
        << "noise | mean SE3 | var. SE3 | "
        << "mean Sim3 | var. of Sim3| Scale change "
        << endl;
    cout << fd.noise << " "
        << fd.mean_rmse_se3 << " "
        << fd.var_rmse_se3 << " "
        << fd.mean_rmse_sim3 << " "
        << fd.var_rmse_sim3 << " "
        << fd.av_scale << endl;
  }
  cout << endl
      << "noise | mean SE3 | var. SE3 | mean Sim3 | var. of Sim3| Scale change "
      << endl;
  for (uint i=0; i<fig_data.size();++i)
  {
    const Figure3Data & fd = fig_data[i];

    cout << fd.noise << " "
        << fd.mean_rmse_se3 << " "
        << fd.var_rmse_se3 << " "
        << fd.mean_rmse_sim3 << " "
        << fd.var_rmse_sim3 << " "
        << fd.av_scale << endl;
  }
  while(true)
  {
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(1,1,1,1);
    view3d.activate2D();



    double sx = WIDTH/(2*1.2);
    double sy = -HEIGHT/(2*1.6);
    Vector<2> offset = makeVector(WIDTH/4,3*HEIGHT/4);

    for (uint i=0; i<fig_data.size()-1;++i)
    {
      const Figure3Data & fd = fig_data[i];
      const Figure3Data & fd2 = fig_data[i+1];
      glColor3f(0,1,0);
      view3d.drawLine2D(makeVector(sx*fd.noise,sy*fd.mean_rmse_se3) + offset,
                        makeVector(sx*fd2.noise,sy*fd2.mean_rmse_se3) + offset);
      glColor3f(1,0,0);
      view3d.drawLine2D(makeVector(sx*fd.noise,sy*fd.mean_rmse_sim3)
                        + offset,
                        makeVector(sx*fd2.noise,sy*fd2.mean_rmse_sim3)
                        + offset);
      glColor3f(0.5,0.5,0.5);
      view3d.drawLine2D(makeVector(sx*fd.noise,sy*fd.av_scale*2) + offset,
                        makeVector(sx*fd2.noise,sy*fd2.av_scale*2) + offset);

    }
    glColor3f(0,0,0);
    view3d.drawLine2D(makeVector(0,0) + offset,
                      makeVector(0,sy*1.6) + offset);
    view3d.drawLine2D(makeVector(0,0) + offset,
                      makeVector(sx*1.2,0) + offset);
    glFlush();
    glwin.swap_buffers();
  }
}


void figure4()
{
  srand(1);

  map<int, UvqGauss > uvq_map;

  vector<Vector<3> > true_point_vec;

  vector<SE3<> > true_pose_vec;

  const int WIDTH = 640;
  const int HEIGHT = 480;
  ImageRef win_size(WIDTH*2,HEIGHT*2);
  GuiWindow glwin(win_size);
  LinearCamera cam_pars(makeVector(195, -195),makeVector(162 ,125),win_size);

  LinearCamera cam3d(makeVector(1000,-1000),
                     makeVector(WIDTH/2,HEIGHT/2),
                     win_size);
  GuiView view3d[4];

  view3d[0] = GuiView(ImageRef(WIDTH,HEIGHT),
                      cam3d,
                      makeVector(-M_PI*0.5,M_PI*0.5,M_PI*0.25),
                      makeVector(0,0,50));

  view3d[1] = GuiView(ImageRef(WIDTH,HEIGHT),
                      cam3d,makeVector(-M_PI*0.5,M_PI*0.5,M_PI*0.25),
                      makeVector(0,0,50));

  view3d[2] = GuiView(ImageRef(WIDTH,HEIGHT),
                      cam3d,makeVector(-M_PI*0.5,M_PI*0.5,M_PI*0.25),
                      makeVector(0,0,50));

  view3d[3] = GuiView(ImageRef(WIDTH,HEIGHT),
                      cam3d,
                      makeVector(-M_PI*0.5,M_PI*0.5,M_PI*0.25),
                      makeVector(0,0,50));

  view3d[2].init(&glwin,Rectangle(0,0,0.5,0.5));
  view3d[3].init(&glwin,Rectangle(0.5,0,1,0.5));
  view3d[0].init(&glwin,Rectangle(0,0.5,0.5,1));
  view3d[1].init(&glwin,Rectangle(0.5,0.5,1,1));

  uint num_true_points = 15000;

  uint num_try_poses = 50000;

  double r = 10;

  for (uint i=0; i<num_true_points; ++i)
  {

    double theta = rand_u()*M_PI*2;
    double phi = sqrt(rand_u()*M_PI);
    int sign = 1;
    if (rand_u()<0.5)
      sign = -1;

    double x = r*cos(theta)*sin(phi);
    double y = r*sin(theta)*sin(phi);
    double z = sign*r*cos(phi);
    Vector<3> xyz = makeVector(x,y,z);
    true_point_vec.push_back(xyz);
  }

  r = 11;

  for (uint i=0; i<num_try_poses; ++i)
  {
    int idx = i;
    double factor = 1;
    if (i>num_try_poses*0.5)
    {
      idx = i-num_try_poses;
      factor = -1;
    }

    double theta = 2*((idx)*M_PI*2)/(num_try_poses);
    double phi = factor*atan(theta*0.1);

    double x = r*cos(theta)*cos(phi);
    double y = r*sin(theta)*cos(phi);
    double z = -r*sin(phi);

    static Vector<3> point = makeVector(x,y,z);
    Vector<3> old_point = point;

    if (norm(makeVector(x,y,z)-old_point)>0.3)
    {
      point = makeVector(x,y,z);
      true_pose_vec.push_back(CameraPose2(point,makeVector(0,0,0),old_point));
    }

  }

  for (uint i=0; i<num_try_poses; ++i)
  {
    int idx = i;
    double factor = -1;
    if (i>num_try_poses*0.5)
    {
      idx = i-num_try_poses;
      factor = 1;
    }

    double theta = 2*((idx)*M_PI*2)/(num_try_poses);
    double phi = factor*atan(theta*0.1);

    double x = r*cos(theta)*cos(phi);
    double y = r*sin(theta)*cos(phi);
    double z = -r*sin(phi);

    static Vector<3> point = makeVector(x,y,z);
    Vector<3> old_point = point;

    if (norm(makeVector(x,y,z)-old_point)>0.3)
    {
      point = makeVector(x,y,z);
      true_pose_vec.push_back(CameraPose2(point,makeVector(0,0,0),old_point));
    }

  }


  for (uint i=0; i<num_try_poses; ++i)
  {
    int idx = i;
    double factor = 1;
    if (i>num_try_poses*0.5)
    {
      idx = i-num_try_poses;
      factor = -1;
    }

    double theta = ((idx)*M_PI*2)/(num_try_poses);
    double phi = factor*atan(theta*0.5);

    double x = r*cos(theta)*cos(phi);
    double y = r*sin(theta)*cos(phi);
    double z = -r*sin(phi);

    static Vector<3> point = makeVector(x,y,z);
    Vector<3> old_point = point;

    if (norm(makeVector(x,y,z)-old_point)>0.3)
    {
      point = makeVector(x,y,z);
      true_pose_vec.push_back(CameraPose2(point,makeVector(0,0,0),old_point));
    }

  }

  for (uint i=0; i<num_try_poses; ++i)
  {
    int idx = i;
    double factor = -1;
    if (i>num_try_poses*0.5)
    {
      idx = i-num_try_poses;
      factor = 1;
    }

    double theta = ((idx)*M_PI*2)/(num_try_poses);
    double phi = factor*atan(theta*0.5);

    double x = r*cos(theta)*cos(phi);
    double y = r*sin(theta)*cos(phi);
    double z = -r*sin(phi);

    static Vector<3> point = makeVector(x,y,z);
    Vector<3> old_point = point;

    if (norm(makeVector(x,y,z)-old_point)>0.3)
    {
      point = makeVector(x,y,z);
      true_pose_vec.push_back(CameraPose2(point,makeVector(0,0,0),old_point));
    }

  }

  vector<SE3<> > keyframe_vec;
  SE3<> cur_pose = true_pose_vec[0];

  vector<IdObs<2> > newobs_vec;
  vector<IdObs<2> > matchedtrack_obs_vec;
  vector<IdObs<2> > matchedmap_obs_vec;
  set<int> trackpoint_id_set;
  set<int> mappoint_id_set;

  double noise= 0.75;

  double rmse_sim3 = -1;
  getObsWithID(cam_pars,
               true_pose_vec[0],
               true_point_vec,
               trackpoint_id_set,
               noise,
               matchedtrack_obs_vec,
               newobs_vec);

  cur_pose = true_pose_vec[1];

  int global_num = 0;

  map<int,int> true_id_global_id_map;
  map<int,int> global_id_true_id_map;

  for(uint i=0; i<newobs_vec.size(); ++i)
  {
    const IdObs<2> & id_obs = newobs_vec[i];
    Matrix<3> Lambda = Zeros(3);
    Lambda(0,0) = Po2(cam_pars.F()[0]);
    Lambda(1,1) = Po2(cam_pars.F()[1]);
    Lambda(2,2) = 0;
    Vector<3> p;
    p.slice<0,2>() = cam_pars.unmap(id_obs.obs);
    p[2] = 0.9;
    uvq_map[id_obs.point_id] = UvqGauss(p,Lambda,id_obs.obs);
    trackpoint_id_set.insert(id_obs.point_id);
  }

  newobs_vec.clear();
  getObsWithID(cam_pars,
               cur_pose,
               true_point_vec,
               trackpoint_id_set,
               noise,
               matchedtrack_obs_vec,
               newobs_vec);

  BundleAdjusterParams ba_pars_f(false,0,1);
  BundleAdjusterParams ba_pars_k(false,0,10);
  BA_SE3_XYZ ba;
  SE3UVQ se3uvq(cam_pars);
  SE3XYZ se3xyz(cam_pars);

  DrawItems::Line3DList line_list2;

  SE3<> ba_frame = cur_pose*true_pose_vec[0].inverse();

  for(uint i=0; i<matchedtrack_obs_vec.size(); ++i)
  {
    const IdObs<2> & id_obs = matchedtrack_obs_vec[i];

    map<int, UvqGauss >::iterator uvq_it = uvq_map.find(id_obs.point_id);

    if(uvq_it!=uvq_map.end())
    {

      UvqGauss & uvq = uvq_it->second;

      ba.filterSingleFeatureOnly(ba_frame,
                                 uvq.mean,
                                 uvq.inv_cov,
                                 se3uvq,
                                 id_obs.obs,
                                 ba_pars_f);

      uvq.obs2 = id_obs.obs;
      ++(uvq.num_updates);

      Vector<3>  xyz = makeVector(uvq.mean[0],uvq.mean[1],1)/uvq.mean[2];

      Cholesky<3> Ch(uvq.inv_cov);

      Matrix<3> Sigma = Ch.get_inverse();

      double q1 = uvq.mean[2]+sqrt(Sigma(2,2))*3;
      double q2 = max(0.0000000001, uvq.mean[2]-sqrt(Sigma(2,2))*3);
      Vector<3> p1
          = transform(true_pose_vec.at(0).inverse(),
                      makeVector(uvq.mean[0]/q1, uvq.mean[1]/q1, 1./q1)  ) ;
      Vector<3> p2
          = transform(true_pose_vec.at(0).inverse(),
                      makeVector(uvq.mean[0]/q2, uvq.mean[1]/q2, 1./q2)  );
      line_list2.push_back(make_pair(p1,p2));

    }
    else{
      assert(false);
    }
  }

  IdStructObsMap globalid_structobs_map;
  {
    vector<SE3<> > ba_poses;

    ba_poses.push_back(true_pose_vec[0]);
    ba_poses.push_back(true_pose_vec[1]);
    vector<IdObs<2> > ba_obs;
    vector<Vector<3> > ba_points;
    vector<int> id_vec;
    SE3<>  kf1;
    SE3<>  kf2;
    kf1 = ba_poses[0];
    kf2 = ba_poses[1];


    int point_id = 0;
    for (map<int,UvqGauss>::const_iterator it=uvq_map.begin();
    it!=uvq_map.end();
    ++it)
    {
      if (it->second.num_updates<1)
        continue;

      int true_id = it->first;

      ba_obs.push_back(IdObs<2>(point_id,0,it->second.obs1));
      ba_obs.push_back(IdObs<2>(point_id,1,it->second.obs2));

      Vector<3>  xyz
          = transform(true_pose_vec[0].inverse(),
                      makeVector(it->second.mean[0],
                                 it->second.mean[1],1)/it->second.mean[2]);

      IdObsMap id_obs_map;
      id_obs_map.insert(make_pair(0,it->second.obs1));
      id_obs_map.insert(make_pair(1,it->second.obs2));

      int global_id = global_num;

      true_id_global_id_map[true_id] = global_id;
      global_id_true_id_map[global_id] = true_id;
      globalid_structobs_map.insert(make_pair(global_id,
                                              StructObs(xyz,id_obs_map)));

      id_vec.push_back(global_id);

      ba_points.push_back(xyz);
      ++point_id;
      ++global_num;
    }

    BA_SE3_XYZ::_TrackMap track_map;
    for (uint i_obs=0; i_obs<ba_obs.size();++i_obs)
    {
      IdObs<2> & id_obs = ba_obs[i_obs];

      BA_SE3_XYZ::_TrackMap::iterator it
          = track_map.find(id_obs.point_id);
      if (it==track_map.end())
      {
        list<IdObs<2> > obs_list;
        obs_list.push_back(id_obs);
        track_map.insert(make_pair(id_obs.point_id,obs_list));
      }
      else
      {
        it->second.push_back(id_obs);
      }
    }
    ba.calcFast(ba_poses,ba_points,se3xyz,track_map,1,ba_pars_k);

    for (uint i=0; i<id_vec.size(); ++i)
    {
      IdStructObsMap::iterator it = globalid_structobs_map.find(id_vec[i]);

      if (it!=globalid_structobs_map.end())
      {

        it->second.p3d = ba_points[i];
      }
      else
        assert(false);
    }
    cur_pose = ba_poses[1];
    keyframe_vec.push_back(ba_poses[0]);
    keyframe_vec.push_back(ba_poses[1]);
  }

  trackpoint_id_set.clear();
  uvq_map.clear();

  for(uint i=0; i<newobs_vec.size(); ++i)
  {
    const IdObs<2> & id_obs = newobs_vec[i];

    Matrix<3> Lambda = Zeros(3);
    Lambda(0,0) = Po2(cam_pars.F()[0]);
    Lambda(1,1) = Po2(cam_pars.F()[1]);
    Lambda(2,2) = 0;
    Vector<3> p;
    p.slice<0,2>() = cam_pars.unmap(id_obs.obs);
    p[2] = 0.9;

    uvq_map[id_obs.point_id] = UvqGauss(p,Lambda,id_obs.obs);

    trackpoint_id_set.insert(id_obs.point_id);
  }

  uint truepose_id = 2;
  DrawItems::LineList l_list;

  vector<pair<int,int> > loop_list;

  vector<double> median_vec;
  median_vec.push_back(1);
  median_vec.push_back(1);

  vector<Sim3<> > trans7_list;

  vector<SE3<> > trans6_list;
  vector<SE3<> > cor7_pose_vec;

  list<Constraint<SE3<>,6> > se3_list;
  list<Constraint<Sim3<>,7> > sim3_list;

  while (true){

    if (truepose_id<true_pose_vec.size())
    {
      set<double> depth_set;

      SE3<> old_pose = cur_pose;
      line_list2.clear();
      newobs_vec.clear();
      mappoint_id_set.clear();

      for (IdStructObsMap::iterator it=globalid_structobs_map.begin();
      it!=globalid_structobs_map.end(); ++it)
      {

        map<int,int>::iterator id_it = global_id_true_id_map.find(it->first);
        if (id_it!=global_id_true_id_map.end())
        {

          mappoint_id_set.insert(id_it->second);
        }
      }


      matchedmap_obs_vec.clear();
      matchedtrack_obs_vec.clear();

      getObsWith2IDs(cam_pars,
                     true_pose_vec[truepose_id],
                     true_point_vec,
                     mappoint_id_set,
                     trackpoint_id_set,
                     noise,

                     matchedmap_obs_vec,
                     matchedtrack_obs_vec,
                     newobs_vec);

      vector<Vector<3> > ba_point_vec;
      vector<IdObs<2> > ba_obs_vec;
      vector<SE3<> > ba_pose_vec;

      map<int,int> new_true_id_global_id_map;
      map<int,int> new_global_id_true_id_map;

      ba_pose_vec.push_back(cur_pose);
      //track features in map
      for (uint i=0; i<matchedmap_obs_vec.size();++i)
      {
        IdObs<2> id_obs = matchedmap_obs_vec[i];

        map<int,int>::iterator it = true_id_global_id_map.find(id_obs.point_id);

        if (it!=true_id_global_id_map.end())
        {
          int global_id = it->second;
          IdStructObsMap::iterator xyz_it
              = globalid_structobs_map.find(global_id);
          if (xyz_it!=globalid_structobs_map.end())
          {

            new_true_id_global_id_map[id_obs.point_id] = global_id;
            new_global_id_true_id_map[global_id] = id_obs.point_id;

            Vector<3> xyz = xyz_it->second.p3d;

            id_obs.point_id = ba_point_vec.size();
            id_obs.frame_id = 0;

            Vector<2> obs_d = se3xyz.map(cur_pose,xyz);

            double z = transform(cur_pose,xyz)[2];

            depth_set.insert(z);

            l_list.push_back(make_pair(obs_d,id_obs.obs+makeVector(2,2)));

            ba_point_vec.push_back(xyz);

            ba_obs_vec.push_back(id_obs);
            xyz_it->second.frameid_obs_map.insert(make_pair(keyframe_vec.size(),
                                                            id_obs.obs));
          }
          else
          {
            assert(false);
          }
        }
        else
        {
          assert(false);
        }
      }

      //optimise over cur_pose
      //      BA_SE3_XYZ::_TrackMap track_map;
      //       for (uint i_obs=0; i_obs<ba_obs_vec.size();++i_obs)
      //       {
      //         IdObs<2> & id_obs = ba_obs_vec[i_obs];

      //         BA_SE3_XYZ::_TrackMap::iterator it
      //             = track_map.find(id_obs.point_id);
      //         if (it==track_map.end())
      //         {
      //           list<IdObs<2> > obs_list;
      //           obs_list.push_back(id_obs);
      //           track_map.insert(make_pair(id_obs.point_id,obs_list));
      //         }
      //         else
      //         {
      //           it->second.push_back(id_obs);
      //         }
      //       }
      //ba.calcFastMotionOnly(ba_pose_vec,ba_point_vec,se3xyz,track_map,ba_pars_k);

      std::list<IdObs<2> > obs_list = std::list<IdObs<2> >(ba_obs_vec.begin(),
                                                           ba_obs_vec.end());
      ba.calcFastMotionOnly(cur_pose,ba_point_vec,se3xyz,obs_list,ba_pars_k);



      //ba.calcMotionOnly(ba_pose_vec,ba_point_vec,se3xyz,ba_obs_vec,0,ba_pars_k);
      // cur_pose = ba_pose_vec[0];

      ba_frame = cur_pose*old_pose.inverse();


      for(uint i=0; i<matchedtrack_obs_vec.size(); ++i)
      {
        const IdObs<2> & id_obs = matchedtrack_obs_vec[i];

        map<int, UvqGauss >::iterator uvq_it = uvq_map.find(id_obs.point_id);

        if(uvq_it!=uvq_map.end())
        {

          UvqGauss & uvq = uvq_it->second;

          ba.filterSingleFeatureOnly(ba_frame,
                                     uvq.mean,
                                     uvq.inv_cov,
                                     se3uvq,
                                     id_obs.obs,
                                     ba_pars_f);

          uvq.obs2 = id_obs.obs;
          ++(uvq.num_updates);

          Cholesky<3> Ch(uvq.inv_cov);

          Matrix<3> Sigma = Ch.get_inverse();


          double q1 = uvq.mean[2]+sqrt(Sigma(2,2))*3;
          double q2 = max(0.0000000001, uvq.mean[2]-sqrt(Sigma(2,2))*3);
          Vector<3> p1
              = transform(old_pose.inverse(),
                          makeVector(uvq.mean[0]/q1, uvq.mean[1]/q1, 1./q1)  );
          Vector<3> p2
              = transform(old_pose.inverse(),
                          makeVector(uvq.mean[0]/q2, uvq.mean[1]/q2, 1./q2)  );
          line_list2.push_back(make_pair(p1,p2));

          IdObsMap id_obs_map;

          Vector<3> xyz
              = transform(old_pose.inverse(),
                          makeVector(uvq.mean[0], uvq.mean[1],1)/uvq.mean[2]);

          id_obs_map.insert(make_pair(keyframe_vec.size()-1,
                                      uvq_it->second.obs1));
          id_obs_map.insert(make_pair(keyframe_vec.size(),
                                      uvq_it->second.obs2));

          int global_id = global_num;

          new_true_id_global_id_map[id_obs.point_id] = global_id;
          new_global_id_true_id_map[global_id] = id_obs.point_id;

          globalid_structobs_map.insert(make_pair(global_id,
                                                  StructObs(xyz,id_obs_map)));

          ++global_num;
        }
        else{
          assert(false);
        }
      }

      map<int,int>::iterator id_it = new_global_id_true_id_map.find(33);

      true_id_global_id_map = map<int,int>(new_true_id_global_id_map);
      global_id_true_id_map = map<int,int>(new_global_id_true_id_map);

      id_it = global_id_true_id_map.find(33);

      keyframe_vec.push_back(cur_pose);

      {
        vector<int> frame_id_vec;

        int id = truepose_id;

        uint num_frames_to_optimise = min(id+1,20);

        int num_fix_frames
            = std::max(2.,std::min(num_frames_to_optimise*0.25,5.));
        bool gauge_free = false;

        if (id<4)
        {
          num_fix_frames = 1;
          gauge_free = true;
        }

        for (uint i=0; i<num_frames_to_optimise; ++i)
          frame_id_vec.push_back(keyframe_vec.size()-num_frames_to_optimise+i);

        vector<Vector<3> > ba_point_vec;
        vector<SE3<> > ba_frame_vec;
        vector<IdObs<2> > ba_obs_vec;
        vector<int> globalid_vec;

        getBaParams(keyframe_vec,
                    globalid_structobs_map,
                    frame_id_vec,
                    ba_frame_vec,
                    ba_point_vec,
                    ba_obs_vec,
                    globalid_vec);

        BA_SE3_XYZ::_TrackMap track_map;
        for (uint i_obs=0; i_obs<ba_obs_vec.size();++i_obs)
        {
          IdObs<2> & id_obs = ba_obs_vec[i_obs];

          BA_SE3_XYZ::_TrackMap::iterator it
              = track_map.find(id_obs.point_id);
          if (it==track_map.end())
          {
            list<IdObs<2> > obs_list;
            obs_list.push_back(id_obs);
            track_map.insert(make_pair(id_obs.point_id,obs_list));
          }
          else
          {
            it->second.push_back(id_obs);
          }
        }
        ba.calcFast(ba_frame_vec,
                    ba_point_vec,
                    se3xyz,
                    track_map,
                    num_fix_frames,

                    ba_pars_k);

        //        ba.calcFull(ba_frame_vec,
        //                    ba_point_vec,
        //                    se3xyz,
        //                    ba_obs_vec,
        //                    num_fix_frames,
        //                    0,
        //                    ba_pars_k);

        for (uint i=0; i<ba_frame_vec.size(); ++i)
        {
          keyframe_vec.at(frame_id_vec.at(i)) = ba_frame_vec.at(i);
        }

        for (uint i=0; i<globalid_vec.size(); ++i)
        {
          globalid_structobs_map.find(globalid_vec.at(i))->second.p3d
              = ba_point_vec.at(i);
        }
      }

      cur_pose = keyframe_vec[keyframe_vec.size()-1];

      trackpoint_id_set.clear();
      uvq_map.clear();

      for(uint i=0; i<newobs_vec.size(); ++i)
      {
        const IdObs<2> & id_obs = newobs_vec[i];

        Matrix<3> Lambda = Zeros(3);
        Lambda(0,0) = Po2(cam_pars.F()[0]);
        Lambda(1,1) = Po2(cam_pars.F()[1]);
        Lambda(2,2) = 0;
        Vector<3> p;
        p.slice<0,2>() = cam_pars.unmap(id_obs.obs);
        p[2] = 0.9;

        uvq_map[id_obs.point_id] = UvqGauss(p,Lambda,id_obs.obs);

        trackpoint_id_set.insert(id_obs.point_id);
      }

      static bool match_last_time = false;
      static pair<int,int> cand_pair = make_pair(0,0);

      static double cand_mindist = 0.5;
      double min_dist=0.5;
      int arg_min_dist = 0;
      for (int i=0; i<max(0,int(keyframe_vec.size()-20)); ++i)
      {
        double n
            = norm((true_pose_vec[truepose_id]
                    *true_pose_vec[i].inverse()).get_translation());

        if (n<min_dist)
        {
          min_dist = n;
          arg_min_dist = i;
        }
      }
      if (min_dist<0.5)
      {

        if (!match_last_time || min_dist<cand_mindist)
        {
          cand_mindist = min_dist;
          cand_pair = make_pair(truepose_id,arg_min_dist);
          match_last_time = true;
        }

      }
      else
      {
        if (match_last_time)
        {
          loop_list.push_back(cand_pair);
          match_last_time = false;
        }
      }
      double m = median(depth_set);
      //cout << "MEDIAN" << m << endl;

      median_vec.push_back(m);

      if(truepose_id+1>=true_pose_vec.size())
      {
        loop_list.push_back(make_pair(0,true_pose_vec.size()-1));

        vector<SE3<> >  se3_vec;
        vector<Sim3<> >  sim3_vec;
        int trans_id = 0;
        for (uint i=0; i<keyframe_vec.size(); ++i)
        {
          SE3<> & pose = keyframe_vec[i];
          Sim3<>  sim(pose.get_rotation(),pose.get_translation(),1.);

          se3_vec.push_back(pose);
          sim3_vec.push_back(sim);

        }

        trans7_list.push_back(sim3_vec[0]);
        trans6_list.push_back(se3_vec[0]);

        for (uint i=1; i<keyframe_vec.size(); ++i)
        {
          trans7_list.push_back(sim3_vec[i]);
          trans6_list.push_back(se3_vec[i]);

          Matrix<6,6> inf6 = TooN::Identity;
          Matrix<7,7> inf7 = TooN::Identity;


          se3_list.push_back(
              Constraint<SE3<>,6 >(trans_id,
                                   trans_id+1,
                                   se3_vec[i]*se3_vec[i-1].inverse(),
                                   inf6));

          sim3_list.push_back(
              Constraint<Sim3<>,7 >(trans_id,
                                    trans_id+1,
                                    sim3_vec[i]*sim3_vec[i-1].inverse(),
                                    inf7));

          ++trans_id;
        }
        Matrix<6,6> inf6 = Identity(6);
        Matrix<7,7> inf7 = Identity(7);

        cout << "Number of loop constraints: " << loop_list.size() << endl;

        for (uint i=0; i<loop_list.size(); ++i)
        {
          pair<int,int> & int_pair = loop_list.at(i);

          SE3<> rel_pose
              = true_pose_vec.at(int_pair.first)
                *true_pose_vec.at(int_pair.second).inverse();

          double median_scale
              = median_vec.at(int_pair.first)/median_vec.at(int_pair.second);

          Sim3<> loop_constraint(rel_pose.get_rotation(),
                                 rel_pose.get_translation(),
                                 median_scale);

          se3_list.push_back(
              Constraint<SE3<>,6 >(int_pair.second,
                                   int_pair.first,
                                   rel_pose,
                                   inf6));

          sim3_list.push_back(
              Constraint<Sim3<>,7 >(int_pair.second,
                                    int_pair.first,
                                    loop_constraint,
                                    inf7));
        }

        SE3ConFun se3_confun;
        Sim3ConFun sim3_confun;

        GraphOptimizer<Sim3<>,7> opt7;
        StopWatch sw7;
        sw7.start();
        opt7.optimize(trans7_list,
                      sim3_list,
                      sim3_confun,
                      1,
                      2,
                      0.0000000000001);
        sw7.stop();

        for (uint i_f=0; i_f<trans7_list.size(); ++i_f)
        {
          cor7_pose_vec.push_back(
              SE3<>(trans7_list[i_f].get_rotation(),
                    trans7_list[i_f].get_translation()
                    /trans7_list[i_f].get_scale()));
        }

        SE3CompareModScale t;
        double s=1;
        double v1 = t.optimize(true_pose_vec,cor7_pose_vec,s,10);

        for (uint i_f=0; i_f<cor7_pose_vec.size(); ++i_f)
        {
          cor7_pose_vec[i_f].get_translation()
              = s*cor7_pose_vec[i_f].get_translation();
        }
        cout << "Sim3 rmse: " << sqrt(v1/true_pose_vec.size()) << endl;
        cout << "time in s: " << sw7.getStoppedTime() << endl;

        GraphOptimizer<SE3<>,6> opt6;
        StopWatch sw6;
        sw6.start();
        opt6.optimize(trans6_list,se3_list,se3_confun,1,2,0.0000000000001);
        sw6.stop();

        s=1;
        double v2 = t.optimize(true_pose_vec,trans6_list,s,10);

        for (uint i_f=0; i_f<trans6_list.size(); ++i_f)
        {
          trans6_list[i_f].get_translation()
              = s*trans6_list[i_f].get_translation();
        }
        cout << "Se3 rmse: " << sqrt(v2/true_pose_vec.size()) << endl;
        cout << "time in s: " << sw6.getStoppedTime() << endl;
        rmse_sim3 = sqrt(v1/true_pose_vec.size());
      }

      ++truepose_id;
    }

    glClearColor(1,1,1,1);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



    bool end = false;
    if(truepose_id>=true_pose_vec.size())
      end = true;

    for (uint i_view=0; i_view<4; ++i_view)
    {
      double ball_size = 0.02;
      if(i_view==0)
      {
#ifdef RV_OPENCV_SUPPORT
        view3d[0].activate2D();

        ColorImage rgb_img(ImageRef(640,480));
        rgb_img.fill(Rgb<byte>(255,255,255));
        cv::String str = "";

        double font_face = 2;

        if(end)
        {
          str = "7 DoF graph optimisation ";
          cv::putText(rgb_img.cv_mat,
                      str,
                      cv::Point(30,60),
                      font_face,
                      1.1,
                      cv::Scalar(0,255,0,255),
                      1.8,
                      CV_AA);

          str = "2 iterations";
          cv::putText(rgb_img.cv_mat,
                      str,
                      cv::Point(30,110),
                      font_face,
                      1.1,
                      cv::Scalar(0,255,0,255),
                      1.8,
                      CV_AA);
        }
        else{
          str = "Aircraft flying over a sphere";
          cv::putText(rgb_img.cv_mat,
                      str,
                      cv::Point(30,60),
                      font_face,
                      1.1,
                      cv::Scalar(0,0,255,255),
                      1.8,
                      CV_AA);

          str = "with downward looking camera";
          cv::putText(rgb_img.cv_mat,
                      str,
                      cv::Point(30,110),
                      font_face,
                      1.1,
                      cv::Scalar(0,0,255,255),
                      1.8,
                      CV_AA);
        }

        str = "Ground truth trajectory";
        cv::putText(rgb_img.cv_mat,
                    str,
                    cv::Point(30,200),
                    font_face,
                    0.8,
                    cv::Scalar(0,0,0,255),
                    1.8,
                    CV_AA);

        str = "Monocular visual odometry";
        cv::putText(rgb_img.cv_mat,
                    str,
                    cv::Point(30,260),
                    font_face,
                    0.8,
                    cv::Scalar(255,0,0,255),
                    1.8,
                    CV_AA);


        if(end)
        {
          char cstr[70];
          sprintf(cstr, "Root mean square error: %4.3f",rmse_sim3);
          str = cstr;
          cv::putText(rgb_img.cv_mat,
                      str,
                      cv::Point(30,290),
                      font_face,
                      0.8,
                      cv::Scalar(255,0,0,255),
                      1.8,
                      CV_AA);
        }

        glColor3f(1,1,1);

        view3d[0].drawTexture2D(rgb_img);
#endif
        continue;
      }
      else if (i_view==1)
      {

        view3d[i_view].activate3D(
            SE3<>(Identity,
                  makeVector(-0.1,-0.1,4))
            *true_pose_vec[keyframe_vec.size()-1]);

        ball_size = 0.015;
      }
      else
        view3d[i_view].activate3D();

      if((i_view==1 || i_view==3)    && !end)
      {

        glColor3f(1,0,0);
        for (IdStructObsMap::iterator it=globalid_structobs_map.begin();
        it!=globalid_structobs_map.end(); ++it)
        {
          view3d[i_view].drawBall3D(it->second.p3d,ball_size);
        }
      }

      if(i_view==2)
      {

        for (uint i=0; i<true_point_vec.size(); ++i)
        {
          glColor3f(0.5,0.5,0.5);
          view3d[i_view].drawBall3D(true_point_vec[i],ball_size);
        }
      }

      if(i_view==1 || i_view==2|| i_view==3)
      {
        glColor3f(0.75,0.75,0.75);
        for (uint i=0; i<keyframe_vec.size(); ++i)
        {
          view3d[i_view].drawPose3D(true_pose_vec[i]);
        }


      }

      if(!end && (i_view==1 || i_view==3))
      {
        glColor3f(0.75,0,0);
        for (uint i=0; i<keyframe_vec.size(); ++i)
        {
          view3d[i_view].drawPose3D(keyframe_vec[i]);
        }

      }


      if (end)
      {
        glColor3f(0,0,1);
        for (uint i=0; i<loop_list.size(); ++i)
        {
          pair<int,int> & int_pair = loop_list[i];
          if(i_view==1 || i_view==3)
          {
            view3d[i_view]
                .drawLine3D(trans2center(cor7_pose_vec[int_pair.first]),
                            trans2center(cor7_pose_vec[int_pair.second]));
          }
        }

        glColor3f(0.75,0,0);
        for (uint i=0; i<cor7_pose_vec.size(); ++i)
        {
          view3d[i_view].drawPose3D(cor7_pose_vec[i]);
        }

      }
    }

    glFlush();
    glwin.swap_buffers();

    usleep(10);
  }
}


int main(int argc, char *argv[])
{
  cout << endl;
  cout << "This program creates the figures of the paper:\n\n";
  cout << "> H. Strasdat, J.M.M. Montiel, A.J. Davison:\n"
      << "  'Scale Drift-Aware Large Scale Monocular SLAM',\n"
      << "  Proc. of Robotics: Science and Systems (RSS),\n"
      << "  Zaragoza, Spain, 2010.\n"
      << "  http://www.roboticsproceedings.org/rss06/p10.htm <\n\n\n";
  if (argc>1 && atoi(argv[1])==2)
    figure2();
  if (argc>1 && atoi(argv[1])==3)
    figure3();
  if (argc>1 && atoi(argv[1])==4)
    figure4();
  else
  {
    cout << "Please type: ./rss2010_demo [2|3|4]" << endl;
    cout << endl;
    cout << "Example: './rss2010_demo 2'" << endl << endl;
  }
  return 0;

}


