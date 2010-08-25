#include <vector>
#include <set>

#include <TooN/TooN.h>
#include <TooN/SVD.h>
#include <TooN/SymEigen.h>

#include <cvd/random.h>
#include <cvd/gl_helpers.h>

#include "gui_view.h"
#include "transformations.h"
#include "stopwatch.h"
#include "bundle_adjuster.h"
#include "graph_optimizer.h"



using namespace std;
using namespace CVD;
using namespace TooN;
using namespace RobotVision;



int main(int argc, char *argv[])
{


  const int WIDTH = 640;
  const int HEIGHT = 480;

  double image_noise = 0.5;
  double spurious_matches_prob = 0.05;
  bool robust_kernel = true;

  if (argc>3)
  {
    image_noise = atof(argv[1]);
    spurious_matches_prob = atof(argv[2]);
    robust_kernel = atof(argv[3]);
  }
  else
  {
    cout << endl;
    cout << "Type: ./ba_demo X Y Z" << endl;
    cout << endl;
    cout << "X := noise in the image e.g., X=0.5" << endl;
    cout << "Y :=  probability of a spurious match, e.g. Y=0.05 or Y=0.0"
        << endl;
    cout << "Z := use robust kernel? Z=0 (off) or Z=1 (on)"
        << endl;
    cout << endl;
    exit(0);
  }

  cout << endl;
  cout << "BUNDLE ADJUSTMENT EXAMPLE" << endl;
  cout << "Image noise: " << image_noise <<endl;
  cout << "Probability of a spurious match: "
      << spurious_matches_prob <<endl;
  cout << "use robust kernel: " << robust_kernel<<endl<<endl<<endl;

  ImageRef win_size(WIDTH,HEIGHT);
  GuiWindow glwin(win_size);

  LinearCamera cam(makeVector(484.69093,485.57926),
                   makeVector(312.22988,247.83456),
                   ImageRef(640,480));



  LinearCamera guiview_cam(makeVector(1000,1000),
                           makeVector(WIDTH/2,HEIGHT/2),
                           ImageRef(WIDTH,HEIGHT));

  SE3XYZ se3xyz(cam);

  GuiView view3d(ImageRef(WIDTH,HEIGHT),
                 guiview_cam,
                 makeVector(0,0,0),
                 makeVector(-0.5,0,3));

  view3d.init(&glwin,Rectangle(0,0,1,1));
  view3d.set_handler(new DoomViewHandler(&view3d));

  vector<Vector<3> > true_point_list;
  vector<Vector<3> > true_vis_point_list;
  vector<SE3<double> > true_frame_list;

  vector<Vector<3> > noisy_point_list;


  vector<SE3<double> > noisy_frame_list;


  // create list of true points
  for (int i=0; i<200; ++i)
  {
    true_point_list.push_back(makeVector(rand_u()*5,
                                         (rand_u()-0.5)*0.9,rand_u()+0.03));
  }

  // create list of noisy points
  for (uint i=0; i<true_point_list.size(); ++i)
  {
    noisy_point_list.push_back(true_point_list[i]
    + makeVector(0,0,rand_g()) * 0.03);
  }



  int num_frames = 50;
  double factor = 0.1;
  double initial_lambda = 0.05;

  //create list of true and noisy poses/frames
  for (int i=0; i<num_frames; ++i)
  {
    SO3<double> rot;
    true_frame_list.push_back(SE3<double>(rot,makeVector(-i*factor,0,0)));
    if (i<2)
    {
      noisy_frame_list.push_back(SE3<double>(rot,makeVector(-i*factor,0,0)));
    }
    else
    {
      noisy_frame_list.push_back(
          SE3<double>(rot,makeVector(-i*factor,0,0)
                      + makeVector(rand_g(),rand_g(),rand_g()) * 0.01 ));
    }
  }


  srand(5);


    vector<IdObs<2> > obs_vec;
    vector<Vector<3> > visible_point_list;
    vector<Vector<3> > cor_point_list;
    vector<Vector<4> > cor_hompoint_list;
    vector<SE3<double> > cor_frame_list;

   BundleAdjuster<SE3<> ,6,3,3,IdObs<2>,2>::_TrackMap track_list;

    //create list of observations and visible points
    int j_p = 0;
    for (uint j=0; j<true_point_list.size(); ++j)
    {
      bool matched = false;
      vector<IdObs<2> > tmp_vec;

      BundleAdjuster<SE3<> ,6,3,3,IdObs<2>,2>::_Track track;



      Vector<3> & p3d = true_point_list[j];
      for (uint i=0; i<true_frame_list.size(); ++i)
      {
        Vector<2> obs;

        obs = se3xyz.map(true_frame_list[i],p3d)
              + image_noise*makeVector(rand_g(),rand_g());

        if(obs[0]>=0 && obs[1]>=0 && obs[0]<=639 && obs[1] <= 479)
        {
          if (rand_u()<spurious_matches_prob)
          {
            obs = makeVector(rand_u()*640,rand_u()*480);
          }
          matched = true;
          tmp_vec.push_back(IdObs<2>(j_p,i,obs));


          track.push_back(IdObs<2>(j_p,i,obs));

        }
      }
      if (matched){

        obs_vec.insert(obs_vec.end(),tmp_vec.begin(),tmp_vec.end());
        visible_point_list.push_back(noisy_point_list[j]);
        cor_hompoint_list.push_back(unproject(noisy_point_list[j]));
        true_vis_point_list.push_back(noisy_point_list[j]);
        track_list.insert(make_pair(j_p,track));
          ++j_p;

      }
    }

    BundleAdjuster<SE3<> ,6,3,3,IdObs<2>,2> ba;

    ba.verbose = true;

    StopWatch sw;
    cor_frame_list = noisy_frame_list;
    cor_point_list = visible_point_list;


    uint num_iters = 5;
   BundleAdjusterParams opt_params(robust_kernel,1,num_iters,initial_lambda);




  cout << "Do " << num_iters << " iterations of bundle adjustment!" << endl;

   ba.verbose = true;
    sw.start();
       cout << ba.calcFast(cor_frame_list,
                   cor_point_list,
                   se3xyz,
                   track_list,
                   1,
                   opt_params,3) << endl;
       sw.stop();
       cout << "time in s: " << sw.getStoppedTime() << endl <<endl;


  



  




  while (true)
  {
    glClearColor(1,1,1,1);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    view3d.activate3D();

    glColor3f(0.8,0.8,0.8);
    for (uint i=0; i<noisy_point_list.size(); ++i)
    {
      view3d.drawBall3D(noisy_point_list[i],0.005);
    }
    for (uint i=0; i<noisy_frame_list.size(); ++i)
    {
      view3d.drawPose3D(noisy_frame_list[i]);
    }

    glColor3f(0,0,0.75);
    for (uint i=0; i<true_frame_list.size(); ++i)
    {
      view3d.drawPose3D(true_frame_list[i]);
    }
    for (uint i=0; i<true_vis_point_list.size(); ++i)
    {
      view3d.drawBall3D(true_vis_point_list[i],0.005);
    }


    glColor3f(0.75,0,0);
    view3d.drawBall3D(true_frame_list[num_frames-1].inverse().get_translation(),0.001);





        glColor3f(0.75,0,0);
        for (uint i=0; i<cor_frame_list.size(); ++i)
        {
          view3d.drawPose3D(cor_frame_list[i]);
        }
        for (uint i=0; i<cor_hompoint_list.size(); ++i)
        {
          view3d.drawLine3D(true_vis_point_list[i],project(cor_hompoint_list[i]));
          view3d.drawBall3D(project(cor_hompoint_list[i]),0.005);
        }

    glColor3f(0.3,0.3,0.3);
    view3d.drawLine3D(makeVector(0,0,0),makeVector(10000,0,0));
    view3d.drawLine3D(makeVector(0,0,0),makeVector(0,10000,0));
    view3d.drawLine3D(makeVector(0,0,0),makeVector(0,0,10000));

    glFlush();
    glwin.swap_buffers();

    glwin.handle_events_default();

    usleep(1000);
  }
}




