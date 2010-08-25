/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * bundle_adjuster.h is part of RobotVision.
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

#ifndef RV_BUNDLE_ADJUSTER_H
#define RV_BUNDLE_ADJUSTER_H


#include <vector>
#include <map>

#include <TooN/TooN.h>
#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <TooN/se2.h>
#include <TooN/Cholesky.h>
#include <cvd/image_io.h>

#include "maths_utils.h"
#include "sparse_matrix.h"
#include "sparse_cholesky.h"
#include "transformations.h"
#include "stopwatch.h"



#ifdef RV_BUILD_RELEASE
//std::vector bounds check only for debug build

#define AT(x,y)  x[y]
//MARCRO definitions in headers are in general not a good idea,
//but here its is undefined later, so it is ok.

#else

#define AT(x,y)  x.at(y)
//MARCRO definitions in headers are in general not a good idea,
//but here its is undefined later, so it is ok.
#endif





namespace RobotVision
{
  /**
     * Parameter class for BundleAdjuster (see below)
     */
  class BundleAdjusterParams
  {
  public:
    BundleAdjusterParams(bool robust_kernel=true,
                         double kernel_param=1,
                         int num_iter=50,
                         double initial_mu=-1,
                         double final_mu_factor=999999.,
                         double tau = 0.00001):
    robust_kernel(robust_kernel),
    kernel_param(kernel_param),
    num_iter(num_iter),
    initial_mu(initial_mu),
    final_mu_factor(final_mu_factor),
    tau(tau)
    {
    }
    bool robust_kernel;
    double kernel_param;
    int num_iter;
    double initial_mu;
    double final_mu_factor;
    double tau;
  };



  /**
     * Helper class for BundleAdjuster (see below)
     */
  template <int FrameDoF, int PointDoF, int ObsDim> class JacData
  {
  public:
    JacData()
    {
      double nan = NAN;
      J_frame = TooN::Ones(ObsDim,FrameDoF)*nan;
      J_point = TooN::Ones(ObsDim,PointDoF)*nan;
      J_anchor = TooN::Ones(ObsDim,FrameDoF)*nan;
    }

  public:
    TooN::Matrix<ObsDim,FrameDoF> J_frame;
    TooN::Matrix<ObsDim,PointDoF> J_point;
    TooN::Matrix<ObsDim,FrameDoF> J_anchor;
  };





  /**
   *
   * This class performs Bundle adjustment (BA) using Levenberg-Marquardt (LM)
   * as described in
   *
   * > H. Strasdat, J.M.M. Montiel, A.J. Davison:
   *   "Scale Drift-Aware Large Scale Monocular SLAM",
   *   Proc. of Robotics: Science and Systems (RSS),
   *   Zaragoza, Spain, 2010.
   *   http://www.roboticsproceedings.org/rss06/p10.html <
   *
   * The the update strategy of mu follows
   * >M.I. A. Lourakis and A.A. Argyros, "The Design and Implementation
   *  of a Generic Sparse Bundle Adjustment Software Package Based on
   *  the Levenberg-Marquardt Algorithm", Technical Report, 2004.<
   * and the implementation of the Schur complement is inspired by
   * >C. Engels, H. Stewnius, D. Nister, "Bundle Adjustment Rules",
   * Photogrammetric Computer Vision (PCV), September 2006.
   * and
   * >K. Konolige, "Sparse Sparse Bundle Adjustment", BMVC 2010.
   *
   * Frame: How is the frame/pose represented? (e.g. SE3)
   * FrameDoF: How many DoF has the pose/frame? (e.g. 6 DoF, that is
   *           3 DoF translation, 3 DoF rotation)
   * PointParNum: number of parameters to represent a point
   *              (4 for a 3D homogenious point)
   * PointDoF: DoF of a point (3 DoF for a 3D homogenious point)
   * ObsDim: dimensions of observation (2 dim for (u,res) image
   *         measurement)
   */
  template <typename Frame,
  int FrameDoF,
  int PointParNum,
  int PointDoF,
  typename Obs,
  int ObsDim>
  class BundleAdjuster
  {
  public:

    typedef AbstractPrediction<Frame,FrameDoF,PointParNum,PointDoF,ObsDim>
        _AbsJac;
    typedef std::vector<JacData<FrameDoF, PointDoF,ObsDim> >
        _JacVec;
    typedef std::vector<TooN::Vector<PointParNum> >
        _PointVec;

    typedef std::list<Obs>
        _Track;
    typedef std::map<int,_Track >
        _TrackMap;



    BundleAdjuster(){
      verbose = 0;
    }

    int verbose;


    inline void SchurComp(const _Track & track,
                          int num_fix_frames,
                          int point_id,
                          const std::vector<TooN::Matrix<PointDoF,FrameDoF> > & H_p,
                          const TooN::Matrix<PointDoF,PointDoF> & inv_H_pp,
                          const TooN::Vector<PointDoF> & t_p,
                          std::vector<std::map<int,TooN::Matrix<FrameDoF,PointDoF> > > & T_,
                          TooN::Vector<> & B,
                          RowBlockMapVec<FrameDoF> & A_)
    {
      for (typename _Track::const_iterator it_cam = track.begin();
      it_cam!=track.end();
      ++it_cam)
      {
        const Obs & id_obs = *it_cam;
        int frame_id = id_obs.frame_id-num_fix_frames;
        if (frame_id>=0)
        {

          TooN::Matrix<FrameDoF,PointDoF> T_pc
              = AT(H_p,frame_id).T()*inv_H_pp;

          AT(T_,point_id).insert(std::make_pair(frame_id, T_pc));

          B.slice(frame_id*FrameDoF,FrameDoF)
              -= AT(H_p,frame_id).T()*t_p;

          for (typename _Track::const_iterator it_cam = track.begin();
          it_cam!=track.end();
          ++it_cam)
          {
            const Obs & id_obs2 = *it_cam;
            int frame_id2 = id_obs2.frame_id-num_fix_frames;
            if (frame_id2>=frame_id)
            {
              TooN::Matrix<FrameDoF,FrameDoF> tmp
                  = -(T_pc*AT(H_p,frame_id2));

              A_.add(tmp,frame_id,frame_id2);
              if (frame_id!= frame_id2)
              {
                A_.add(tmp.T(),frame_id2,frame_id);
              }
            }
          }
        }
      }
    }

    inline void SchurCompUpperTriag
        (const _Track  & track,
         int num_fix_frames,
         int point_id,
         const std::vector<TooN::Matrix<PointDoF,FrameDoF> > & H_p,
         const TooN::Matrix<PointDoF,PointDoF> & inv_H_pp,
         const TooN::Vector<PointDoF> & t_p,
         std::vector<std::map<int,TooN::Matrix<FrameDoF,PointDoF> > > & T_,
         TooN::Vector<> & B,
         RowBlockMapVec<FrameDoF> & A_)
    {
      for (typename _Track::const_iterator it_cam = track.begin();
      it_cam!=track.end();
      ++it_cam)
      {
        const Obs & id_obs = *it_cam;
        int frame_id = id_obs.frame_id-num_fix_frames;
        if (frame_id>=0)
        {

          TooN::Matrix<FrameDoF,PointDoF> T_pc
              = AT(H_p,frame_id).T()*inv_H_pp;

          AT(T_,point_id).insert(std::make_pair(frame_id, T_pc));

          B.slice(frame_id*FrameDoF,FrameDoF)
              -= AT(H_p,frame_id).T()*t_p;

          for (typename _Track::const_iterator it_cam = track.begin();
          it_cam!=track.end();
          ++it_cam)
          {
            const Obs & id_obs2 = *it_cam;
            int frame_id2 = id_obs2.frame_id-num_fix_frames;
            if (frame_id2>=frame_id)
            {
              TooN::Matrix<FrameDoF,FrameDoF> tmp
                  = -(T_pc*AT(H_p,frame_id2));

              A_.add(tmp,frame_id,frame_id2);
            }
          }
        }
      }
    }


    double calcFast(std::vector<Frame > & frame_list,
                    _PointVec & point_list,
                    _AbsJac & prediction,
                    const _TrackMap & track_list,
                    const int num_fix_frames,
                    const BundleAdjusterParams & ba_params,
                    int embedded_point_iter = 0,
                    bool use_qr= false)
    {
      assert(track_list.size()>0);
      assert(point_list.size()==track_list.size());

#ifndef RV_SUITESPARSE_SUPPORT
      use_qr = false; //QR is not implemented using CSparse (yet)
#endif


      int num_points = point_list.size();
      int num_frames = frame_list.size();



      double chi2 = 0;

      std::vector<Frame > new_frame_list(frame_list);
      _PointVec new_point_list(point_list);

      double norm_max_A = 0.;
      for (typename _TrackMap::const_iterator it_tl = track_list.begin();
      it_tl != track_list.end();
      ++it_tl)
      {
        int point_id = it_tl->first;
        const _Track & track = it_tl->second;
        TooN::Vector<PointParNum> & point = AT(point_list,point_id);

        for (typename _Track::const_iterator it_cam = track.begin();
        it_cam!=track.end();
        ++it_cam)
        {
          const Obs & id_obs = *it_cam;
          Frame & frame = AT(frame_list,id_obs.frame_id);

          TooN::Matrix<ObsDim,PointDoF> J_p;
          TooN::Matrix<ObsDim,FrameDoF> J_c;
          TooN::Vector<ObsDim> f = id_obs.obs
                                   - prediction.map_n_bothJac(frame,
                                                              point,
                                                              J_c,
                                                              J_p);

          norm_max_A
              = max(norm_max_A,
                    max(norm_max(mulW(J_c,J_c,id_obs).diagonal_slice()),
                        norm_max(mulW(J_p,J_p,id_obs).diagonal_slice())));

          chi2 += sqrW(f,id_obs);
        }
      }
      double mu = ba_params.initial_mu;


      if (ba_params.initial_mu==-1)
      {

        mu = ba_params.tau*norm_max_A;

      }
      double final_mu = ba_params.final_mu_factor*mu;

      if (verbose>0)
      {
        std::cout << "chi2: "<< chi2 << std::endl;
      }


      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;
      for (int i_g=0; i_g<ba_params.num_iter; i_g++)
      {
        if (verbose>0)
        {
          std::cout << "iteration: "<< i_g << std::endl;
        }


        double rho = 0; //Assign value, so the compiler is not complaining...
        do{

          if (verbose>0)
          {
            std::cout << "mu: " <<mu<< std::endl;
          }


          TooN::Matrix<FrameDoF,FrameDoF> muI_f = mu*TooN::Identity;
          TooN::Matrix<PointDoF,PointDoF> muI_p = mu*TooN::Identity;

          RowBlockMapVec<FrameDoF> A_(num_frames-num_fix_frames);
          TooN::Vector<> B = TooN::Zeros(FrameDoF*(num_frames-num_fix_frames));

          std::vector<TooN::Matrix<PointDoF,FrameDoF> >
              H_p(num_frames-num_fix_frames);

          std::vector<TooN::Vector<PointDoF> > t_(num_points);
          std::vector<std::map<int,TooN::Matrix<FrameDoF,PointDoF> > >
              T_(num_points);



          for (typename _TrackMap::const_iterator it_tl = track_list.begin();
          it_tl != track_list.end();
          ++it_tl)
          {

            int point_id = it_tl->first;
            const _Track & track = it_tl->second;
            TooN::Vector<PointParNum> & point = AT(point_list,point_id);
            TooN::Matrix<PointDoF,PointDoF> H_pp =  muI_p;
            TooN::Vector<PointDoF> b_p = TooN::Zeros;


            for (typename _Track::const_iterator it_cam = track.begin();
            it_cam!=track.end();
            ++it_cam)
            {



              const Obs & id_obs = *it_cam;
              Frame & frame = AT(frame_list,id_obs.frame_id);
              int frame_id = id_obs.frame_id-num_fix_frames;


              assert(id_obs.point_id==point_id);

              TooN::Matrix<ObsDim,PointDoF> J_p;
              TooN::Matrix<ObsDim,FrameDoF> J_c;

              TooN::Vector<ObsDim> f = id_obs.obs
                                       - prediction.map_n_bothJac(frame,
                                                                  point,
                                                                  J_c,
                                                                  J_p);
              J_c *= -1;
              J_p *= -1;


              if(ba_params.robust_kernel)
              {
                double nrm = std::max(0.00000000001,
                                      sqrt(sqrW(f,id_obs)));
                double w = sqrt(kernel(nrm,ba_params.kernel_param))/nrm;
                f *= w;
              }

              H_pp += mulW(J_p,J_p,id_obs);
              b_p -= mulW(J_p,f,id_obs);
              if (frame_id>=0)
              {
                A_.add(mulW(J_c,J_c,id_obs), muI_f, frame_id, frame_id);
                AT(H_p,frame_id) = mulW(J_p,J_c,id_obs);
                B.slice(frame_id*FrameDoF,FrameDoF) -= mulW(J_c,f,id_obs);
              }

            }

            TooN::Cholesky<PointDoF> Ch_H_pp(H_pp);
            TooN::Matrix<PointDoF,PointDoF> inv_H_pp = Ch_H_pp.get_inverse();

            TooN::Vector<PointDoF> t_p = inv_H_pp*b_p;
            AT(t_,point_id) = t_p;

            if (use_qr)
            {
              SchurComp(track,
                        num_fix_frames,
                        point_id,
                        H_p,
                        inv_H_pp,
                        t_p,
                        T_,
                        B,
                        A_);
            }
            else
            {
              SchurCompUpperTriag(track,
                                  num_fix_frames,
                                  point_id,
                                  H_p,
                                  inv_H_pp,
                                  t_p,
                                  T_,
                                  B,
                                  A_);
            }


          }

          int matrix_type = 1;//symmetric, upper_triangle
          if (use_qr)
          {
            matrix_type = 0;//unsymmetric
          }
          SparseMatrix<> sS(A_,matrix_type);


          TooN::Vector<> delta((num_frames-num_fix_frames)*FrameDoF);

          try{
            SparseSolver<> Ch(sS);

            delta = Ch.backsub(B);
          }catch (NotPosSemiDefException & e) {
            // not positive definite so increase mu and try again
            std::cout << "Not pose Def" << std::endl;

            mu *= nu;
            nu *= 2.;
            stop = (mu>final_mu);
            continue;
          }


          for (int i=0; i<num_frames-num_fix_frames; ++i)
          {
            AT(new_frame_list,i+num_fix_frames)
                = prediction.add(AT(frame_list,i+num_fix_frames),
                                 delta.slice(i*FrameDoF,FrameDoF));
          }

          for (typename _TrackMap::const_iterator it_tl = track_list.begin();
          it_tl != track_list.end();
          ++it_tl)
          {
            int point_id = it_tl->first;
            const _Track & track = it_tl->second;

            TooN::Vector<PointDoF> dp = AT(t_,point_id);

            for (typename _Track::const_iterator it_cam = track.begin();
            it_cam!=track.end();
            ++it_cam)
            {
              const Obs & id_obs = *it_cam;

              int frame_id = id_obs.frame_id-num_fix_frames;
              if (frame_id>=0)
              {
                dp -= AT(T_,point_id).find(frame_id)->second.T()
                      *delta.slice(frame_id*FrameDoF,FrameDoF);
              }

            }

            AT(new_point_list,point_id)
                = prediction.add(AT(point_list,point_id), dp);
          }
          double new_chi2 = 0;
          if (embedded_point_iter>0)
          {
            new_chi2 = calcFastStructureOnly(new_frame_list,
                                             new_point_list,
                                             prediction,
                                             track_list,
                                             BundleAdjusterParams
                                             (ba_params.robust_kernel,
                                              ba_params.kernel_param,
                                              embedded_point_iter,
                                              ba_params.initial_mu,
                                              ba_params.final_mu_factor));
          }
          else
          {


            for (typename _TrackMap::const_iterator it_tl = track_list.begin();
            it_tl != track_list.end();
            ++it_tl)
            {
              int point_id = it_tl->first;
              const _Track & track = it_tl->second;
              TooN::Vector<PointParNum> & point = AT(new_point_list,point_id);

              for (typename _Track::const_iterator it_cam = track.begin();
              it_cam!=track.end();
              ++it_cam)
              {
                const Obs & id_obs = *it_cam;
                Frame & frame = AT(new_frame_list,id_obs.frame_id);

                TooN::Vector<ObsDim> f = id_obs.obs - prediction.map(frame,
                                                                     point);
                new_chi2 +=  sqrW(f,id_obs);
              }
            }
          }


          if (isnan(new_chi2))
          {
            std::cerr << "res_new is NAN\n";
            exit(-1);
          }

          if (verbose>0)
          {
            std::cout << "new chi2: "<< new_chi2 << std::endl;
          }
          rho = (chi2-new_chi2)/(delta*(mu*delta+B));

          if(rho>0)
          {
            //std::cerr << mu << std::endl;

            frame_list = std::vector<Frame> (new_frame_list);
            point_list = _PointVec(new_point_list);

            chi2 = new_chi2;

            stop = norm_max(B)<=eps;
            mu *= std::max(1./3.,1-Po3(2*rho-1));
            nu = 2.;


          }
          else
          {
            if (verbose>0)
              std::cout << "no update: chi vs. new_chi2 "
                  << chi2
                  << " vs. "
                  << new_chi2
                  << std::endl;
            mu *= nu;
            nu *= 2.;

            stop = (mu>final_mu);
          }

        }while(!(rho>0 ||  stop));

        if (stop)
          break;

      }
      return chi2;
    }




    double calcFastStructureOnly(std::vector<Frame> & frame_list,
                                 _PointVec & point_list,
                                 _AbsJac & prediction,
                                 const _TrackMap & track_list,
                                 const BundleAdjusterParams & ba_params)
    {
      assert(point_list.size()==track_list.size());
      assert(track_list.size()>0);
      //int num_points = point_list.size();
      double chi2_sum=0;
      double old_chi2_sum=0;

      for (typename _TrackMap
           ::const_iterator it_tl = track_list.begin();
      it_tl != track_list.end();
      ++it_tl)
      {

        double nu = 2;
        double eps =0.000000000000001;
        bool stop = false;

        int point_id = it_tl->first;


        TooN::Vector<PointParNum> & point = AT(point_list,point_id);
        const _Track & track = it_tl->second;

        double chi2 = 0; 
        TooN::Vector<PointParNum> new_point = point;



        double norm_max_A = 0.;
        for (typename _Track::const_iterator it_cam = track.begin();
        it_cam!=track.end();
        ++it_cam)
        {
          const Obs & id_obs = *it_cam;
          Frame & frame = AT(frame_list,id_obs.frame_id);

          TooN::Matrix<ObsDim,PointDoF> J_p;
          TooN::Vector<ObsDim> f = id_obs.obs
                                   - prediction.map_n_pointJac(frame,
                                                               point,
                                                               J_p);
          norm_max_A
              = max(norm_max_A,
                    norm_max(mulW(J_p,J_p,id_obs).diagonal_slice()));
          chi2 += sqrW(f,id_obs);
        }
        old_chi2_sum += chi2;
        double mu = ba_params.initial_mu;


        if (ba_params.initial_mu==-1)
        {

          mu = ba_params.tau*norm_max_A;

        }
        double final_mu = ba_params.final_mu_factor*mu;

        for (int i_g=0; i_g<ba_params.num_iter; i_g++)
        {
          if (verbose>1)
          {
            std::cout << "iteration: "<< i_g << std::endl;
          }


          double rho = 0; //Assign value, so the compiler is not complaining...
          do{

            if (verbose>1)
            {
              std::cout << "mu: " <<mu<< std::endl;
            }

            TooN::Matrix<PointDoF,PointDoF> muI_p = mu*TooN::Identity;

            TooN::Matrix<PointDoF,PointDoF> H_pp =  muI_p;
            TooN::Vector<PointDoF> b_p = TooN::Zeros;


            for (typename _Track::const_iterator it_cam = track.begin();
            it_cam!=track.end();
            ++it_cam)
            {
              const Obs & id_obs = *it_cam;
              Frame & frame = AT(frame_list,id_obs.frame_id);
              TooN::Matrix<ObsDim,PointDoF> J_p;

              TooN::Vector<ObsDim> f = id_obs.obs
                                       - prediction.map_n_pointJac(frame,
                                                                   point,
                                                                   J_p);
              J_p *= -1;


              if(ba_params.robust_kernel)
              {
                double nrm = std::max(0.00000000001,
                                      sqrt(sqrW(f,id_obs)));
                double w = sqrt(kernel(nrm,ba_params.kernel_param))/nrm;
                f *= w;
              }

              H_pp += mulW(J_p,J_p,id_obs);
              b_p -= mulW(J_p,f,id_obs);

            }


            TooN::Vector<PointDoF> delta_p;



            TooN::Cholesky<PointDoF> Ch_H_pp(H_pp);

            delta_p = Ch_H_pp.backsub(b_p);





            new_point
                = prediction.add(point, delta_p);
            
            double new_chi2 = 0;


            for (typename _Track::const_iterator it_cam = track.begin();
            it_cam!=track.end();
            ++it_cam)
            {
              const Obs & id_obs = *it_cam;
              Frame & frame = AT(frame_list,id_obs.frame_id);

              TooN::Vector<ObsDim> f = id_obs.obs - prediction.map(frame,
                                                                   new_point);
              new_chi2 += sqrW(f,id_obs);
            }


            if (isnan(new_chi2))
            {
              std::cerr << "res_new is NAN\n";
              exit(-1);
            }
            rho = (chi2-new_chi2);//(delta*(mu*delta+g));

            if(rho>0)
            {


              point = new_point;

              chi2 = new_chi2;

              stop = norm_max(b_p)<=eps;
              mu *= std::max(1./3.,1-Po3(2*rho-1));
              nu = 2.;
            }
            else
            {
              if (verbose>1)
                std::cout << "no update: chi vs. new_chi2 "
                    << chi2
                    << " vs. "
                    << new_chi2
                    << std::endl;
              mu *= nu;
              nu *= 2.;

              stop = (mu>final_mu);
            }

          }while(!(rho>0 ||  stop));

          if (stop)
            break;

        }
        chi2_sum += chi2;

      }
      if (verbose>0)
        std::cout << " chi vs. new_chi2 "
            << old_chi2_sum
            << " vs. "
            << chi2_sum
            << std::endl;
      return chi2_sum;

    }





    double calcFastMotionOnly(Frame  & frame,
                              _PointVec & point_list,
                              _AbsJac & prediction,
                              const std::list<Obs> & obs_list,
                              const BundleAdjusterParams & ba_params)
    {
      // assert(point_list.size()==obs_list.size());
      assert(obs_list.size()>0);



      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;


      double chi2 = 0;

      Frame  new_frame = frame;
      double norm_max_A = 0.;
      for (typename std::list<Obs>
           ::const_iterator it_cam = obs_list.begin();
      it_cam != obs_list.end();
      ++it_cam)
      {
        int point_id = it_cam->point_id;

        TooN::Vector<PointParNum> & point = AT(point_list,point_id);


        const Obs & id_obs = *it_cam;

        TooN::Matrix<ObsDim,FrameDoF> J_c;
        TooN::Vector<ObsDim> f = id_obs.obs
                                 - prediction.map_n_frameJac(frame,
                                                             point,
                                                             J_c);
        norm_max_A = max(norm_max_A,
                         norm_max(mulW(J_c,J_c,id_obs).diagonal_slice()));
        chi2 += sqrW(f,id_obs);

      }
      double mu = ba_params.initial_mu;

      if (verbose>0)
      {
        std::cout << "init. chi2: "<< chi2 << std::endl;
      }
      if (ba_params.initial_mu==-1)
      {

        mu = ba_params.tau*norm_max_A;

      }
      double final_mu = ba_params.final_mu_factor*mu;

      for (int i_g=0; i_g<ba_params.num_iter; i_g++)
      {
        if (verbose>0)
        {
          std::cout << "iteration: "<< i_g << std::endl;
        }


        double rho = 0; //Assign value, so the compiler is not complaining...
        do{

          if (verbose>0)
          {
            std::cout << "mu: " <<mu<< std::endl;
          }




          TooN::Vector<FrameDoF> B = TooN::Zeros;


          TooN::Matrix<FrameDoF,FrameDoF> A = mu*TooN::Identity;

          for (typename std::list<Obs>
               ::const_iterator it_cam = obs_list.begin();
          it_cam != obs_list.end();
          ++it_cam)
          {
            int point_id = it_cam->point_id;

            TooN::Vector<PointParNum> & point = AT(point_list,point_id);




            const Obs & id_obs = *it_cam;



            TooN::Matrix<ObsDim,FrameDoF> J_c;

            TooN::Vector<ObsDim> f = id_obs.obs
                                     - prediction.map_n_frameJac(frame,
                                                                 point,
                                                                 J_c);
            J_c *= -1;

            if(ba_params.robust_kernel)
            {
              double nrm = std::max(0.00000000001,
                                    sqrt(sqrW(f,id_obs)));
              double w = sqrt(kernel(nrm,ba_params.kernel_param))/nrm;
              f *= w;
            }



            A += mulW(J_c,J_c,id_obs);
            B.slice(0,FrameDoF) -= mulW(J_c,f,id_obs);


          }



          TooN::Vector<FrameDoF> delta;

          TooN::Cholesky<FrameDoF> Ch(A);

          delta = Ch.backsub(B);

          //          std::cout<< A << std::endl<< std::endl;
          //          std::cout<< B << std::endl<< std::endl;



          //std::cout<< delta << std::endl<< std::endl;



          new_frame
              = prediction.add(frame,delta);


          double new_chi2 = 0;
          for (typename std::list<Obs>
               ::const_iterator it_cam = obs_list.begin();
          it_cam != obs_list.end();
          ++it_cam)
          {
            int point_id = it_cam->point_id;

            TooN::Vector<PointParNum> & point = AT(point_list,point_id);


            const Obs & id_obs = *it_cam;


            TooN::Vector<ObsDim> f = id_obs.obs - prediction.map(new_frame,
                                                                 point);
            new_chi2 += sqrW(f,id_obs);

          }

          if (isnan(new_chi2))
          {
            std::cerr << "res_new is NAN\n";
            exit(-1);
          }
          rho = (chi2-new_chi2)/(delta*(mu*delta+B));

          if(rho>0)
          {
            //std::cerr << mu << std::endl;


            frame = new_frame;

            chi2 = new_chi2;

            stop = norm_max(B)<=eps;
            mu *= std::max(1./3.,1-Po3(2*rho-1));
            nu = 2.;
          }
          else
          {
            if (verbose>0)
              std::cout << "no update: chi vs. new_chi2 "
                  << chi2
                  << " vs. "
                  << new_chi2
                  << std::endl;
            mu *= nu;
            nu *= 2.;

            stop = (mu>final_mu);
          }

        }while(!(rho>0 ||  stop));

        if (stop)
          break;


      }
      return chi2;
    }





    /**
     * This function implements an information filter for a single landmark
     * given known frame/pose. It can be used for a landmark initialisation
     * (without any depth prior) within a key-frame BA apporach as described
     * in:
     *
     * > H. Strasdat, J.M.M. Montiel, A.J. Davison:
     *   "Scale Drift-Aware Large Scale Monocular SLAM",
     *   Proc. of Robotics: Science and Systems (RSS),
     *   Zaragoza, Spain, 2010.
     *   http://www.roboticsproceedings.org/rss06/p10.html <
     *
     * frame: pose/frame
     * point: 3D point/landmark (normally in inverse depth coordinates)
     * Lambda: 3D point/landmark inverse covariance
     * prediction: prediction class
     * obs_vec: set of observations
     * ba_params: BA parameters
     */
    double filterSingleFeatureOnly(Frame & frame,
                                   TooN::Vector<PointParNum>  & point,
                                   TooN::Matrix<PointDoF,PointDoF> & Lambda,
                                   _AbsJac& prediction,
                                   const TooN::Vector<ObsDim> & obs,
                                   const BundleAdjusterParams & ba_params)
    {

      TooN::Vector<PointParNum> point_mean = point;

      TooN::Vector<ObsDim> residuals;

      TooN::Vector<PointDoF> residuals_dist(PointDoF);

      //this->ba_params = ba_params;

      TooN::Vector<PointDoF>  new_point;

      TooN::Matrix<ObsDim,PointDoF> J_point = prediction.pointJac(frame, point);

      TooN::Vector<ObsDim> delta_obs = obs - prediction.map(frame, point);

      double res =  delta_obs*delta_obs;

      residuals = delta_obs;

      TooN::Vector<> diff = (point_mean-point);
      TooN::Vector<> LambdaDiff = Lambda*diff;

      residuals_dist = LambdaDiff;

      res += diff*LambdaDiff;

      if (isnan(res))
      {
        std::cerr << "res is NAN\n";
        exit(-1);
      }
      if(verbose>0)
        std::cout << "res: " << res << std::endl;

      TooN::Matrix<PointDoF,PointDoF> V = J_point.T() * J_point;
      TooN::Vector<PointDoF> g= J_point.T() * residuals;
      g += residuals_dist;

      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;
      double mu = ba_params.initial_mu;


      if (ba_params.initial_mu==-1)
      {
        double norm_max_A = norm_max(V.diagonal_slice());


        mu = ba_params.tau*norm_max_A;

      }
      double final_mu = ba_params.final_mu_factor*mu;

      for (int i_g=0; i_g<ba_params.num_iter; i_g++)
      {
        if (verbose>0)
        {
          std::cout << "iteration: "<< i_g << std::endl;
        }
        double rho = 0; //Assign value, so the compiler is not complaining...
        do{
          if (verbose>0)
          {
            std::cout << "mu: " <<mu<< std::endl;
          }

          TooN::Matrix<PointDoF,PointDoF> H = Lambda+ V;
          H.diagonal_slice() += TooN::Ones(H.num_cols())*mu ;

          TooN::Cholesky<PointDoF> Ch_h(H);

          TooN::Vector<> delta = Ch_h.backsub(g);

          new_point =  point + delta;

          TooN::Vector<ObsDim> delta_obs
              = obs - prediction.map(frame, new_point);

          double res_new =  delta_obs*delta_obs;

          residuals = delta_obs;

          TooN::Vector<PointDoF> diff = (point_mean-new_point);
          TooN::Vector<PointDoF> LambdaDiff
              = Lambda*diff;

          residuals_dist  = LambdaDiff;

          res_new += diff*LambdaDiff;

          if (isnan(res_new))
          {
            std::cerr << "res_new is NAN\n";
            exit(-1);
          }

          rho = (res-res_new)/(delta*(mu*delta+g));

          if(rho>0)
          {
            if(verbose>0)
              std::cout << "res_new: " << res_new << std::endl;
            point= new_point;

            res = res_new;

            J_point = prediction.pointJac(frame, point);
            V = J_point.T() * J_point;
            g= J_point.T() * residuals;
            g += residuals_dist;

            stop = norm_max(g)<=eps;
            mu *= std::max(1./3.,1-Po3(2*rho-1));
            nu = 2.;
          }
          else
          {
            if (verbose>0)
              std::cout << "no update: res vs.res_new "
                  << res << " vs. " << res_new << std::endl;
            mu *= nu;
            nu *= 2.;

            if (verbose)
              std::cout << "mu" << mu << std::endl;

            stop = (mu>final_mu);
          }

        }while(!(rho>0 ||  stop));

        if (stop)
          break;

      }

      Lambda += V;
      return res;
    }

  protected:

    /** pseudo-huber cost function */
    double kernel(double delta, double kernel_param)
    {
      double b = kernel_param;
      return fabs(2*Po2(b)*(sqrt(1+Po2(delta/b))-1));
    }




    template <int Trans1DoF, int Trans2DoF>
        inline TooN::Matrix<Trans1DoF,Trans2DoF>
        mulW(const TooN::Matrix<ObsDim,Trans1DoF> & J_trans1,
             const TooN::Matrix<ObsDim,Trans2DoF> &  J_trans2,
             const IdObs<ObsDim> & obs )
    {
      return J_trans1.T() * J_trans2;
    }


    template <int Trans1DoF, int Trans2DoF>
        inline TooN::Matrix<Trans1DoF,Trans2DoF>
        mulW(const TooN::Matrix<ObsDim,Trans1DoF> & J_trans1,
             const TooN::Matrix<ObsDim,Trans2DoF> & J_trans2,
             const IdObsLambda<ObsDim> & obs )
    {
      return J_trans1.T() * obs.lambda *J_trans2;
    }


    template <int TransDoF>
        inline TooN::Vector<TransDoF>
        mulW(const TooN::Matrix<ObsDim,TransDoF> & J_trans,
             const TooN::Vector<ObsDim> & f,
             const IdObs<ObsDim> & obs )
    {
      return J_trans.T() * f;
    }

    template <int TransDoF>
        inline TooN::Vector<TransDoF>
        mulW(const TooN::Matrix<ObsDim,TransDoF> & J_trans,
             const TooN::Vector<ObsDim> & f,
             const IdObsLambda<ObsDim> & obs )
    {
      return J_trans.T() * obs.lambda * f;
    }

    inline double
        sqrW(const TooN::Vector<ObsDim> & f,
             const IdObs<ObsDim> & obs )
    {
      return f * f;
    }

    inline double
        sqrW(const TooN::Vector<ObsDim> & f,
             const IdObsLambda<ObsDim> & obs )
    {
      return f * obs.lambda * f;
    }









  };

  typedef BundleAdjuster<TooN::SE3<>,6,3,3,IdObs<2>,2> BA_SE3_XYZ;
  typedef BundleAdjuster<TooN::SE3<>,6,4,3,IdObs<2>,2> BA_SE3_XYZW;
  typedef BundleAdjuster<TooN::SE3<>,6,3,3,IdObs<2>,2> BA_SE3_XYZ_Lambda;
  typedef BundleAdjuster<TooN::SE3<>,6,4,3,IdObs<2>,2> BA_SE3_XYZW_Lambda;
  typedef BundleAdjuster<TooN::SE3<>,6,3,3,IdObs<3>,3> BA_SE3_XYZ_STEREO;


}


#undef AT

#endif // RV_BUNDLE_ADJUSTER_H
