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
    cout << "Z := use robust kernel?. Z=0 (false) or Z=1 (true)"
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

  LinearCamera cam(makeVector(484.69093,-485.57926),
                   makeVector(312.22988,247.83456),
                   ImageRef(640,480));

  LinearCamera guiview_cam(makeVector(1000,-1000),
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
  vector<SE3<double> > true_frame_list;

  vector<Vector<3> > noisy_point_list;

  vector<Vector<3> > visible_point_list;
  vector<SE3<double> > noisy_frame_list;

  vector<Vector<3> > cor_point_list;
  vector<SE3<double> > cor_frame_list;

  // create list of true points
  for (int i=0; i<500; ++i)
  {
    true_point_list.push_back(makeVector(rand_u()*0.9,
                                         (rand_u()-0.5)*0.9,rand_u()));
  }

  // create list of noisy points
  for (uint i=0; i<true_point_list.size(); ++i)
  {
    noisy_point_list.push_back(true_point_list[i]
                               + makeVector(rand_g(),rand_g(),rand_g()) * 0.01);
  }

  vector<IdObs<2> > obs_vec;

  //create list of true and noisy poses/frames
  for (int i=0; i<10; ++i)
  {
    SO3<double> rot;
    true_frame_list.push_back(SE3<double>(rot,makeVector(-i*0.1,0,1)));
    if (i<1)
    {
      noisy_frame_list.push_back(SE3<double>(rot,makeVector(-i*0.1,0,1)));
    }
    else
    {
      noisy_frame_list.push_back(
          SE3<double>(rot,makeVector(-i*0.1,0,1)
                      + makeVector(rand_g(),rand_g(),rand_g()) * 0.1 ));
    }
  }

  //create list of observations and visible points
  int j_p = 0;
  for (uint j=0; j<true_point_list.size(); ++j)
  {
    bool matched = false;
    vector<IdObs<2> > tmp_vec;

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
      }
    }
    if (matched){
      ++j_p;
      obs_vec.insert(obs_vec.end(),tmp_vec.begin(),tmp_vec.end());
      visible_point_list.push_back(noisy_point_list[j]);

    }
  }

  BundleAdjuster<SE3<> ,6,3,3,IdObs<2>,2> ba;

  StopWatch sw;
  cor_frame_list = noisy_frame_list;
  cor_point_list = visible_point_list;

  uint num_iters = 50;
  BundleAdjusterParams opt_params(robust_kernel,1,num_iters);

  cout << "Do " << num_iters << " iterations of bundle adjustment!" << endl;

  sw.start();
  ba.calcFull(cor_frame_list,
              cor_point_list,
              se3xyz,
              obs_vec,
              1,
              0,
              opt_params,
              false);
  sw.stop();
  cout << "time in s: " << sw.getStoppedTime() << endl <<endl;

  SE3CompareModScale t;
  double s=1;
  double rmse = t.optimize(true_frame_list,cor_frame_list,s,10);
  cout << "RMS error: " << rmse << endl<<endl;
  cout << "True data is shown in blue." << endl;
  cout << "Noisy data before optimisation is shown in grey." << endl;
  cout << "Optimised data is shown in red." << endl;
  cout << endl <<endl;

  while (true)
  {
    glClearColor(1,1,1,1);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    view3d.activate3D();

    glColor3f(0.8,0.8,0.8);
    for (uint i=0; i<true_point_list.size(); ++i)
    {
      view3d.drawBall3D(noisy_point_list[i],0.005);
    }
        for (uint i=0; i<noisy_frame_list.size(); ++i)
    {
      view3d.drawPose3D(noisy_frame_list[i]);
    }

    glColor3f(0,0,0.75);
    for (uint i=0; i<true_frame_list.size(); ++i){
      view3d.drawPose3D(true_frame_list[i]);
    }
    for (uint i=0; i<true_point_list.size(); ++i)
    {
      view3d.drawBall3D(true_point_list[i],0.005);
    }

    glColor3f(0.75,0,0);
    for (uint i=0; i<cor_frame_list.size(); ++i){
      view3d.drawPose3D(cor_frame_list[i]);
    }
    for (uint i=0; i<cor_point_list.size(); ++i)
    {
      view3d.drawLine3D(true_point_list[i],cor_point_list[i]);
      view3d.drawBall3D(cor_point_list[i],0.005);
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




