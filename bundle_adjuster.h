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
#include <TooN/SVD.h>
#include <TooN/LU.h>
#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <TooN/se2.h>
#include <TooN/SymEigen.h>
#include <TooN/Cholesky.h>

#include "maths_utils.h"
#include "sparse_matrix.h"
#include "sparse_cholesky.h"
#include "transformations.h"

#define RV_USE_SPARSE_CHOLESKY 1



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
  class BundleAdjusterParams{
  public:
    BundleAdjusterParams(bool robust_kernel=true,
                         double kernel_param=1,
                         int num_iter=50,
                         double mu=-1):
    robust_kernel(robust_kernel),
    kernel_param(kernel_param),
    num_iter(num_iter),
    mu(mu)
    {
    }
    bool robust_kernel;
    double kernel_param;
    int num_iter;
    double mu;
  };

  /**
     * Helper class for BundleAdjuster (see below)
     */
  template <int FrameDof, int PointDof, int ObsDim> class JacData
  {
  public:
    JacData()
    {
      double nan = NAN;
      J_frame = TooN::Ones(ObsDim,FrameDof)*nan;
      J_point = TooN::Ones(ObsDim,PointDof)*nan;
    }

  public:
    TooN::Matrix<ObsDim,FrameDof>   J_frame;
    TooN::Matrix<ObsDim,PointDof>  J_point;
    int frame_id;
    int point_id;
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
   * The concrete implementation/notation of LM is partially inspired by
   * >M.I. A. Lourakis and A.A. Argyros, "The Design and Implementation
   *  of a Generic Sparse Bundle Adjustment Software Package Based on
   *  the Levenberg-Marquardt Algorithm", Technical Report, 2004.<
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
  int FrameDof,
  int PointParNum,
  int PointDof,
  typename Obs,
  int ObsDim>
  class BundleAdjuster
  {

    typedef AbstractPrediction<Frame,FrameDof,PointParNum,PointDof,ObsDim>
        _AbsJac;
    typedef std::vector<JacData<FrameDof, PointDof,ObsDim> >
        _JacVec;
    typedef std::vector<TooN::Vector<PointParNum> >
        _PointVec;
  public:
    BundleAdjuster(){
      verbose = 0;
    }

    int verbose;

    /**
     * perform full BA over points and frames
     *
     * frame_list: set of poses/frames
     * point_list: set of 3D points/landmarks
     * prediction: prediction class
     * obs_vec: set of observations
     * num_fix_frames: number of frames fixed during observations
     *                 (normally 1, the origin)
     * num_fix_points: number of points fixed (normally 0)
     * ba_params: BA parameters
     * alternation: use approximative, but very fast optimisation by using
     *              alternation between points and frames
     *              (normally not good, because of very slow convergence rate)
     */
    void calcFull(std::vector<Frame > & frame_list,
                  _PointVec & point_list,
                  _AbsJac & prediction,
                  const std::vector<Obs >  & obs_vec,
                  const int num_fix_frames,
                  const int num_fix_points,
                  const BundleAdjusterParams & ba_params,
                  bool alternation=false)
    {
      this->num_points = point_list.size();
      this->num_frames = frame_list.size();
      this->num_obs = obs_vec.size();

      this->num_fix_frames = num_fix_frames;
      this->num_fix_points = num_fix_points;

      this->ba_params = ba_params;

      std::vector<Frame > new_frame_list(num_frames);
      _PointVec new_point_list(num_points);

      residual_vec = std::vector<TooN::Vector<ObsDim> >(num_obs);
      jac_data_vec = _JacVec(num_obs);

      typename _JacVec::iterator it_jac
          = jac_data_vec.begin();
      for (typename std::vector<Obs >::const_iterator it_obs=obs_vec.begin();
      it_obs!=obs_vec.end();
      ++it_obs, ++it_jac)
      {
        it_jac->frame_id = it_obs->frame_id;
        it_jac->point_id = it_obs->point_id;
      }

      getJac(frame_list, point_list, prediction, jac_data_vec);


      double res;
      if(ba_params.robust_kernel)
        res = getRobustResidual(frame_list,
                                point_list,
                                prediction,
                                obs_vec,
                                residual_vec);
      else
        res = getResidual(frame_list,
                          point_list,
                          prediction,
                          obs_vec,
                          residual_vec);

      if (isnan(res))
      {
        std::cerr << "res is NAN\n";
        exit(-1);
      }
      if(verbose>0)
        std::cout << "res: " << res << std::endl;

      std::vector<TooN::Matrix<FrameDof,FrameDof> >
          U(num_frames-num_fix_frames,TooN::Zeros(FrameDof,FrameDof));
      std::vector<TooN::Vector<FrameDof> >
          eps_frame(num_frames-num_fix_frames,TooN::Zeros(FrameDof));

      std::vector<TooN::Matrix<PointDof,PointDof> >
          V(num_points-num_fix_points,TooN::Zeros(PointDof,PointDof));
      TooN::Vector<>  eps_point
          = TooN::Zeros(PointDof*(num_points-num_fix_points));

      calcUVeps(obs_vec,residual_vec,jac_data_vec, U,eps_frame,V,eps_point);

      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;
      double mu = ba_params.mu;

      if (ba_params.mu==-1)
      {
        double norm_max_A = norm_max(AT(U,0).diagonal_slice());
        for (int j=1; j<num_frames-num_fix_frames; ++j)
        {
          norm_max_A = max(norm_max_A,norm_max(AT(U,j).diagonal_slice()));
        }
        for (int i=0; i<num_points-num_fix_points; ++i)
        {
          norm_max_A = max(norm_max_A,norm_max(AT(V,i).diagonal_slice()));
        }
        double tau = 0.001;
        mu = tau*norm_max_A;
      }

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

          TooN::Matrix<> I_muFrame = TooN::Identity(FrameDof)*mu;
          TooN::Matrix<> I_muPoint = TooN::Identity(PointDof)*mu;
          std::vector<TooN::Matrix<FrameDof,FrameDof> >
              U_star(num_frames-num_fix_frames,I_muFrame);

          std::vector<TooN::Matrix<PointDof,PointDof> >
              V_inv(num_points-num_fix_points);

          for (uint i=0; i<U_star.size(); ++i)
          {
            U_star[i] += AT(U,i);
          }
          for (uint i=0; i<V.size(); ++i)
          {
            TooN::Cholesky<PointDof> Ch_V(AT(V,i)+I_muPoint);
            V_inv[i] = Ch_V.get_inverse();
          }

          std::map<std::pair<int,int>, TooN::Matrix<FrameDof,PointDof> > W;
          std::map<std::pair<int,int>, TooN::Matrix<FrameDof,PointDof> > Y;

          typename _JacVec::const_iterator
              jac_iter = jac_data_vec.begin();
          for (typename std::vector<Obs >::const_iterator obs_iter
               = obs_vec.begin();
          obs_iter != obs_vec.end();
          ++jac_iter, ++obs_iter)
          {
            const TooN::Matrix<ObsDim,PointDof> & J_point = jac_iter->J_point;
            if(!isnan(J_point(0,0)))
            {

              const TooN::Matrix<ObsDim,FrameDof> & J_frame = jac_iter->J_frame;
              if(!(isnan(J_frame(0,0))))
              {
                TooN::Matrix<FrameDof,PointDof> tmp = J_frame.T()*J_point;
                std::pair<int,int> id_pair = std::make_pair(obs_iter->frame_id,
                                                            obs_iter->point_id);
                W.insert(make_pair(id_pair,tmp));
                Y.insert(make_pair(id_pair,tmp
                                   *AT(V_inv,
                                       obs_iter->point_id-num_fix_points)));
              }
            }
          }
          TooN::Vector<> delta((num_frames-num_fix_frames)*FrameDof
                               + (num_points-num_fix_points)*PointDof);

          TooN::Vector<> e(FrameDof*(num_frames-num_fix_frames));
          if (alternation )
          {
            for (uint i=0; i<U_star.size(); ++i)
            {
              e.slice(i*FrameDof,FrameDof) = AT(eps_frame,i);
              TooN::LU<FrameDof> LU_U( AT(U,i)+I_muFrame);
              delta.slice(i*FrameDof,FrameDof) = LU_U.backsub(AT(eps_frame,i));

            }
          }
          else
          {
            std::map<std::pair<int,int>, TooN::Matrix<FrameDof,FrameDof> >
                YW_map;
            for (uint i=0; i<U.size(); ++i)
            {
              TooN::Matrix<FrameDof,FrameDof> & Us = AT(U_star,i);
              YW_map.insert(std::make_pair(std::make_pair(i,i),Us));
            }
            for (int j=0; j<num_frames-num_fix_frames; ++j)
            {
              e.slice(j*FrameDof,FrameDof) = AT(eps_frame,j);
            }
            for (int j=0; j<num_frames-num_fix_frames; ++j)
            {

              for (int i=0; i<num_points-num_fix_points; ++i)
              {
                typename std::map<std::pair<int,int>,
                TooN::Matrix<FrameDof,PointDof> >::iterator Y_ij
                = Y.find(std::make_pair(j+num_fix_frames,i+num_fix_points));
                if (Y_ij!=Y.end())
                {
                  for (int k=0; k<num_frames-num_fix_frames; ++k)
                  {
                    typename std::map<std::pair<int,int>,
                    TooN::Matrix<FrameDof,PointDof> >::iterator W_ik
                    = W.find(std::make_pair(k+num_fix_frames,i+num_fix_points));
                    if (W_ik!=W.end())
                    {
                      TooN::Matrix<FrameDof,FrameDof> YW
                          = -(Y_ij->second) * (W_ik->second).T();
                      typename std::map<std::pair<int,int>,
                      TooN::Matrix<FrameDof,FrameDof> >::iterator
                      it = YW_map.find(std::make_pair(j,k));
                      if(it!=YW_map.end())
                      {
                        it->second += YW;
                      }
                      else
                        YW_map.insert(std::make_pair(std::make_pair(j,k),YW));
                    }
                  }

                  e.slice(j*FrameDof,FrameDof)
                      -= (Y_ij->second) * eps_point.slice(i*PointDof, PointDof);
                }
              }
            }

#ifdef RV_USE_SPARSE_CHOLESKY
            TooN::TripletMatrix sparseS(FrameDof*(num_frames-num_fix_frames),
                                        FrameDof*(num_frames-num_fix_frames));
#else
            TooN::Matrix<> S
                = TooN::Zeros(FrameDof*(num_frames-num_fix_frames),
                              FrameDof*(num_frames-num_fix_frames));
#endif

            for (typename std::map<std::pair<int,int>,
                 TooN::Matrix<FrameDof,FrameDof> >::const_iterator it
                 = YW_map.begin(); it!=YW_map.end(); ++it)
            {
              const std::pair<int,int> & ids = it->first;
              const TooN::Matrix<FrameDof,FrameDof> & YW = it->second;

#ifdef RV_USE_SPARSE_CHOLESKY
              sparseS.add(ids.first*FrameDof,ids.second*FrameDof,YW);
#else
              S.slice(ids.first*FrameDof,ids.second*FrameDof,FrameDof,FrameDof)
                  = YW;
#endif
            }

#ifdef RV_USE_SPARSE_CHOLESKY

            TooN::SparseMatrix<> sS(sparseS);
            try{
              TooN::SparseCholesky<> Ch(sS);
              delta.slice(0,FrameDof*(num_frames-num_fix_frames))
                  = Ch.backsub(e);

            }catch (TooN::NotPosSemiDefException & e) {
              // not positive definite so increase mu and try again
              mu *= nu;
              nu *= 2.;
              stop = (mu>999999999.f);
              continue;
            }
#else
            TooN::Cholesky<> Ch(S);
            delta.slice(0,FrameDof*(num_frames-num_fix_frames))= Ch.backsub(e);
#endif
          }

          if (verbose>1)
            std::cout <<"delta: " << delta << std::endl;

          TooN::Vector<> g(FrameDof*(num_frames-num_fix_frames)
                           + PointDof*(num_points-num_fix_points));
          g.slice(0,FrameDof*(num_frames-num_fix_frames)) = e;

          for (int i=0; i<num_points-num_fix_points; ++i)
          {
            TooN::Vector<PointDof> tmp = eps_point.slice(PointDof*i,PointDof);
            for (int j=0; j<num_frames-num_fix_frames; ++j)
            {
              typename std::map<std::pair<int,int>,
              TooN::Matrix<FrameDof,PointDof> >::iterator W_ij
              = W.find(std::make_pair(j+num_fix_frames,i+num_fix_points));
              if (W_ij!=W.end())
              {
                tmp -= W_ij->second.T() * delta.slice(j*FrameDof,FrameDof);
              }
            }
            delta.slice((num_frames-num_fix_frames)
                        *FrameDof + i*PointDof, PointDof)
                = AT(V_inv,i) * tmp;
            g.slice((num_frames-num_fix_frames)*FrameDof + i*PointDof, PointDof)
                = eps_point.slice(PointDof*i,PointDof);
          }

          addToFrames(frame_list,
                      delta.slice(0,FrameDof*(num_frames-num_fix_frames)),
                      prediction,
                      new_frame_list);


          addToPoints(point_list,
                      delta.slice(FrameDof*(num_frames-num_fix_frames),
                                  PointDof*(num_points-num_fix_points)),
                      prediction,
                      new_point_list);

          double res_new;
          if(ba_params.robust_kernel)
            res_new = getRobustResidual(new_frame_list,
                                        new_point_list,
                                        prediction,
                                        obs_vec,
                                        residual_vec);
          else
            res_new = getResidual(new_frame_list,
                                  new_point_list,
                                  prediction,
                                  obs_vec,
                                  residual_vec);

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
            frame_list = std::vector<Frame> (new_frame_list);
            point_list = _PointVec(new_point_list);

            res = res_new;

            getJac(frame_list, point_list, prediction,jac_data_vec);

            U = std::vector<TooN::Matrix<FrameDof,FrameDof> >
                (num_frames-num_fix_frames,
                 TooN::Zeros(FrameDof,FrameDof) );
            eps_frame =  std::vector<TooN::Vector<FrameDof> >
                         (num_frames-num_fix_frames,
                          TooN::Zeros(FrameDof));

            V =  std::vector<TooN::Matrix<PointDof,PointDof> >
                 (num_points-num_fix_points,
                  TooN::Zeros(PointDof,PointDof) );

            eps_point = TooN::Zeros(PointDof*(num_points-num_fix_points));

            calcUVeps(obs_vec,
                      residual_vec,
                      jac_data_vec,
                      U,eps_frame,
                      V,
                      eps_point);

            stop = norm_max(g)<=eps;
            mu *= std::max(1./3.,1-Po3(2*rho-1));
            nu = 2.;
          }
          else
          {
            if (verbose>0)
              std::cout << "no update: res vs.res_new "
                  << res
                  << " vs. "
                  << res_new
                  << std::endl;
            mu *= nu;
            nu *= 2.;

            stop = (mu>999999999.f);
          }

        }while(!(rho>0 ||  stop));

        if (stop)
          break;

      }

    }

    /**
     * Struture-only BA
     *- optimise only wrt. to points and keep frames/poses fixed
     *
     * frame_list: set of poses/frames
     * point_list: set of 3D points/landmarks
     * prediction: prediction class
     * obs_vec: set of observations
     * num_fix_points: number of points fixed (normally 0)
     * ba_params: BA parameters
     */
    void calcStructOnly(std::vector<Frame > & frame_list,
                        _PointVec & point_list,
                        _AbsJac & prediction,
                        const std::vector<Obs >  & obs_vec,
                        const int num_fix_points,
                        const BundleAdjusterParams & ba_params)
    {
      this->num_points = point_list.size();
      this->num_frames = frame_list.size();
      this->num_obs = obs_vec.size();

      this->num_fix_frames = num_frames;
      this->num_fix_points = num_fix_points;

      this->ba_params = ba_params;

       _PointVec new_point_list(num_points);
      residual_vec = std::vector<TooN::Vector<ObsDim> >(num_obs);
      jac_data_vec = _JacVec(num_obs);

      typename _JacVec::iterator it_jac = jac_data_vec.begin();
      for (typename std::vector<Obs >::const_iterator it_obs=obs_vec.begin();
      it_obs!=obs_vec.end();
      ++it_obs, ++it_jac)
      {
        it_jac->frame_id = it_obs->frame_id;
        it_jac->point_id = it_obs->point_id;
      }

      getJac(frame_list, point_list,prediction, jac_data_vec);

      double res;
      if(ba_params.robust_kernel)
        res = getRobustResidual(frame_list,
                                point_list,
                                prediction,
                                obs_vec,
                                residual_vec);
      else
        res = getResidual(frame_list,
                          point_list,
                          prediction,
                          obs_vec,
                          residual_vec);

      if (isnan(res))
      {
        std::cout << "res is NAN\n";
        exit(-1);
      }

      if(verbose>0)
        std::cout << "res: " << res << std::endl;

      std::vector<TooN::Matrix<PointDof,PointDof> >
          V(num_points-num_fix_points,TooN::Zeros(PointDof,PointDof));
      TooN::Vector<>  eps_point
          = TooN::Zeros(PointDof*(num_points-num_fix_points));

      calcVeps(obs_vec,residual_vec,jac_data_vec, V,eps_point);

      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;

      double mu = ba_params.mu;

      if (ba_params.mu==-1)
      {

        double norm_max_A = norm_max(AT(V,0).diagonal_slice());

        for (int i=1; i<num_points-num_fix_points; ++i)
        {
          norm_max_A = max(norm_max_A,norm_max(AT(V,i).diagonal_slice()));
        }

        double tau = 0.001;
        mu = tau*norm_max_A;
      }


      for (int i_g=0; i_g<ba_params.num_iter; i_g++){

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

          TooN::Matrix<> I_muPoint = TooN::Identity(PointDof)*mu;
          TooN::Vector<> g(PointDof*(num_points-num_fix_points));

          TooN::Vector<> delta((num_points-num_fix_points)*PointDof);

          for (uint i=0; i<V.size(); ++i)
          {

            TooN::Cholesky<PointDof> Ch_V(AT(V,i)+I_muPoint);
            delta.slice(i*PointDof, PointDof)
                = Ch_V.backsub(eps_point.slice(PointDof*i,PointDof));
            g.slice(i*PointDof, PointDof)
                = eps_point.slice(PointDof*i,PointDof);

          }

          addToPoints(point_list,
                      delta,
                      prediction,
                      new_point_list);

          double res_new;
          if(ba_params.robust_kernel)
            res_new = getRobustResidual(frame_list,
                                        new_point_list,
                                        prediction,
                                        obs_vec,
                                        residual_vec);
          else
            res_new = getResidual(frame_list,
                                  new_point_list,
                                  prediction,
                                  obs_vec,
                                  residual_vec);

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

            point_list = _PointVec(new_point_list);

            res = res_new;

            getJac(frame_list, point_list,prediction,jac_data_vec);

            V =  std::vector<TooN::Matrix<PointDof,PointDof> >
                 (num_points-num_fix_points,
                  TooN::Zeros(PointDof,PointDof) );

            eps_point = TooN::Zeros(PointDof*(num_points-num_fix_points));

            calcVeps(obs_vec,residual_vec,jac_data_vec, V,eps_point);

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

            stop = (mu>999999999.f);
          }

        }while(!(rho>0 ||  stop));

        if (stop)
          break;
      }
    }


    /**
     * Motion-only BA
     * - optimise only wrt. to frames/poses and keep points fixed
     *
     * frame_list: set of poses/frames
     * point_list: set of 3D points/landmarks
     * prediction: prediction class
     * obs_vec: set of observations
     * num_fix_frames: number of frames fixed during observations
     *                 (normally 0)
     * ba_params: BA parameters
     */
    void calcMotionOnly(std::vector<Frame > & frame_list,
                        _PointVec & point_list,
                        _AbsJac & prediction,
                        const std::vector<Obs>  & obs_vec,
                        const int num_fix_frames,
                        const BundleAdjusterParams & ba_params)
    {
      this->num_points = point_list.size();
      this->num_frames = frame_list.size();
      this->num_obs = obs_vec.size();

      this->num_fix_frames = num_fix_frames;
      this->num_fix_points = num_points;

      this->ba_params = ba_params;

      std::vector<Frame > new_frame_list(num_frames);

      residual_vec = std::vector<TooN::Vector<ObsDim> >(num_obs);
      jac_data_vec = _JacVec(num_obs);

      typename _JacVec::iterator it_jac = jac_data_vec.begin();
      for (typename std::vector<Obs >::const_iterator it_obs=obs_vec.begin();
      it_obs!=obs_vec.end();
      ++it_obs, ++it_jac)
      {
        it_jac->frame_id = it_obs->frame_id;
        it_jac->point_id = it_obs->point_id;
      }

      getJac(frame_list, point_list,prediction, jac_data_vec);

      double res;
      if(ba_params.robust_kernel)
        res = getRobustResidual(frame_list,
                                point_list,
                                prediction,
                                obs_vec,
                                residual_vec);
      else
        res = getResidual(frame_list,
                          point_list,
                          prediction,
                          obs_vec,
                          residual_vec);

      if (isnan(res))
      {
        std::cerr << "res is NAN\n";
        exit(-1);
      }

      if(verbose>0)
        std::cout << "res: " << res << std::endl;

      std::vector<TooN::Matrix<FrameDof,FrameDof> >
          U(num_frames-num_fix_frames,TooN::Zeros(FrameDof,FrameDof));
      std::vector<TooN::Vector<FrameDof> >
          eps_frame(num_frames-num_fix_frames,TooN::Zeros(FrameDof));

      calcUeps(obs_vec,residual_vec,jac_data_vec, U,eps_frame);

      double norm_max_A = norm_max(AT(U,0).diagonal_slice());

      for (int j=1; j<num_frames-num_fix_frames; ++j)
      {
        norm_max_A = max(norm_max_A,norm_max(AT(U,j).diagonal_slice()));
      }

      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;
      double tau = 0.001;
      double mu = tau*norm_max_A;

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

          TooN::Matrix<> I_muFrame = TooN::Identity(FrameDof)*mu;

          std::vector<TooN::Matrix<FrameDof,FrameDof> >
              U_star(num_frames-num_fix_frames,I_muFrame);

          for (uint i=0; i<U_star.size(); ++i)
          {
            U_star[i] += AT(U,i);
          }

          TooN::Matrix<> S = TooN::Zeros(FrameDof*(num_frames-num_fix_frames),
                                         FrameDof*(num_frames-num_fix_frames));
          TooN::Vector<> g(FrameDof*(num_frames-num_fix_frames));
          for (uint i=0; i<U.size(); ++i)
          {
            S.slice(i*FrameDof,i*FrameDof,FrameDof,FrameDof) = AT(U_star,i);
          }

          TooN::Vector<> delta((num_frames-num_fix_frames)*FrameDof);

          for (uint i=0; i<U_star.size(); ++i)
          {
            g.slice(i*FrameDof,FrameDof) = AT(eps_frame,i);

            TooN::LU<FrameDof> LU_U( AT(U,i)+I_muFrame);
            delta.slice(i*FrameDof,FrameDof) = LU_U.backsub(AT(eps_frame,i));

            g.slice(i*FrameDof,FrameDof) = AT(eps_frame,i);

          }

          addToFrames(frame_list,
                      delta,
                      prediction,
                      new_frame_list);

          double res_new;
          if(ba_params.robust_kernel)
            res_new = getRobustResidual(new_frame_list,
                                        point_list,
                                        prediction,
                                        obs_vec,
                                        residual_vec);
          else
            res_new = getResidual(new_frame_list,
                                  point_list,
                                  prediction,
                                  obs_vec,
                                  residual_vec);

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

            frame_list = std::vector<Frame> (new_frame_list);

            res = res_new;

            getJac(frame_list, point_list,prediction,jac_data_vec);


            U = std::vector<TooN::Matrix<FrameDof,FrameDof> >
                (num_frames-num_fix_frames,
                 TooN::Zeros(FrameDof,FrameDof) );
            eps_frame =  std::vector<TooN::Vector<FrameDof> >
                         (num_frames-num_fix_frames,
                          TooN::Zeros(FrameDof));
            calcUeps(obs_vec,residual_vec,jac_data_vec, U,eps_frame);

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

            stop = (mu>999999999.f);
          }

        }while(!(rho>0 ||  stop));

        if (stop)
          break;

      }
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
    void filterSingleFeatureOnly(Frame & frame,
                                 TooN::Vector<PointParNum>  & point,
                                 TooN::Matrix<PointDof,PointDof> & Lambda,
                                 _AbsJac& prediction,
                                 const TooN::Vector<ObsDim> & obs,
                                 const BundleAdjusterParams & ba_params)
    {

      TooN::Vector<PointParNum> point_mean = point;

      TooN::Vector<ObsDim> residuals;

      TooN::Vector<PointDof> residuals_dist(PointDof);

      this->ba_params = ba_params;

      TooN::Vector<PointDof>  new_point;

      TooN::Matrix<ObsDim,PointDof> J_point = prediction.pointJac(frame, point);

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

      TooN::Matrix<PointDof,PointDof> V = J_point.T() * J_point;
      TooN::Vector<PointDof> g= J_point.T() * residuals;
      g += residuals_dist;

      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;
      double mu = ba_params.mu;

      if (ba_params.mu==-1)
      {
        double norm_max_A = norm_max(V.diagonal_slice());

        double tau = 0.001;
        mu = tau*norm_max_A;
      }

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

          TooN::Matrix<PointDof,PointDof> H = Lambda+ V;
          H.diagonal_slice() += TooN::Ones(H.num_cols())*mu ;

          TooN::Cholesky<PointDof> Ch_h(H);

          TooN::Vector<> delta = Ch_h.backsub(g);

          new_point =  point + delta;

          TooN::Vector<ObsDim> delta_obs
              = obs - prediction.map(frame, new_point);

          double res_new =  delta_obs*delta_obs;

          residuals = delta_obs;

          TooN::Vector<PointDof> diff = (point_mean-new_point);
          TooN::Vector<PointDof> LambdaDiff
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

            stop = (mu>999999999.f);
          }

        }while(!(rho>0 ||  stop));

        if (stop)
          break;

      }

      Lambda += V;
    }

  protected:
    void getJac(const std::vector<Frame > & frame_list,
                const _PointVec & point_list,
                const _AbsJac & prediction,
                std::vector<JacData<FrameDof,PointDof,ObsDim> > & jac_data_vec)
    {
      TooN::Matrix<ObsDim,PointDof> J_point;
      TooN::Matrix<ObsDim,FrameDof> J_frame;

      for (typename _JacVec::iterator jac_iter = jac_data_vec.begin();
      jac_iter!=jac_data_vec.end(); ++jac_iter)
      {
        int i_p = jac_iter->point_id;
        int i_f = jac_iter->frame_id;

        const TooN::Vector<PointParNum> & point = AT(point_list,i_p);
        const Frame & frame = AT(frame_list,i_f);

        if ((int)i_f>=num_fix_frames)
        {
          jac_iter->J_frame = prediction.frameJac(frame, point);
        }

        if ((int)i_p>=num_fix_points)
        {
          jac_iter->J_point = prediction.pointJac(frame, point);
        }
      }
    }


    /** standard residual function */
    double getResidual(const std::vector<Frame > & frame_list,
                       const _PointVec & point_list,
                       const _AbsJac & prediction,
                       const std::vector<IdObs<ObsDim> >  & obs_vec,
                       std::vector<TooN::Vector<ObsDim> > & residual_vec)
    {
      double sum = 0;
      typename std::vector<TooN::Vector<ObsDim> >::iterator residual_iter
          = residual_vec.begin();
      for (typename std::vector<Obs>::const_iterator obs_iter=obs_vec.begin();
      obs_iter!=obs_vec.end();
      ++obs_iter, ++residual_iter)
      {
        int i_p = obs_iter->point_id;
        int i_f = obs_iter->frame_id;
        const TooN::Vector<PointParNum> & point= AT(point_list,i_p);
        const Frame & frame = AT(frame_list,i_f);
        TooN::Vector<ObsDim> delta
            = obs_iter->obs - prediction.map(frame, point);
        sum +=  delta*delta;
        *residual_iter = delta;
      }
      return sum;
    }


    /** residual function using inverse observation uncertainty Lambda*/
    double getResidual(const std::vector<Frame > & frame_list,
                       const _PointVec & point_list,
                       const _AbsJac & prediction,
                       const std::vector<IdObsLambda<ObsDim> >  & obs_vec,
                       std::vector<TooN::Vector<ObsDim> > & residual_vec)
    {
      double sum = 0;
      typename std::vector<TooN::Vector<ObsDim> >::iterator residual_iter
          = residual_vec.begin();
      for (typename std::vector<Obs>::const_iterator obs_iter=obs_vec.begin();
      obs_iter!=obs_vec.end();
      ++obs_iter, ++residual_iter)
      {
        int i_p = obs_iter->point_id;
        int i_f = obs_iter->frame_id;
        const TooN::Vector<PointParNum> & point= AT(point_list,i_p);
        const Frame & frame = AT(frame_list,i_f);
        TooN::Vector<ObsDim> delta
            = obs_iter->obs - prediction.map(frame, point);
        sum +=  delta*obs_iter->lambda*delta;
        *residual_iter = delta;
      }
      return sum;
    }


    /** residual function using robust kernel*/
    double getRobustResidual(const std::vector<Frame > & frame_list,
                             const _PointVec & point_list,
                             const _AbsJac & prediction,
                             const std::vector<IdObs<ObsDim> >  & obs_vec,
                             std::vector<TooN::Vector<ObsDim> > & residual_vec)
    {
      double sum = 0;

      typename std::vector<TooN::Vector<ObsDim> >::iterator residual_iter
          = residual_vec.begin();

      for (typename std::vector<Obs>::const_iterator obs_iter=obs_vec.begin();
      obs_iter!=obs_vec.end();
      ++obs_iter, ++residual_iter)
      {
        int i_p = obs_iter->point_id;
        int i_f = obs_iter->frame_id;
        const TooN::Vector<PointParNum> & point= AT(point_list,i_p);
        const Frame & frame = AT(frame_list,i_f);
        TooN::Vector<ObsDim> delta
            = obs_iter->obs - prediction.map(frame, point);
        double sum_2=0;

        for (uint i=0; i<ObsDim;++i)
        {
          sum_2 += Po2(delta[i]);
        }
        double nrm = std::max(0.00000000001,sqrt(sum_2));
        double w = sqrt(kernel(nrm))/nrm;

        TooN::Vector<ObsDim> delta_r;
        for (uint i=0; i<ObsDim;++i)
        {
          delta_r[i] = w*delta[i];
        }
        sum +=  delta_r*delta_r;
        *residual_iter = delta_r;
      }
      return sum;
    }


    /** residual function using robust kernel and inverse covariance Lambda*/
    double getRobustResidual(const std::vector<Frame > & frame_list,
                             const _PointVec & point_list,
                             const _AbsJac & prediction,
                             const std::vector<IdObsLambda<ObsDim> >  & obs_vec,
                             std::vector<TooN::Vector<ObsDim> > & residual_vec)
    {
      double sum = 0;
      typename std::vector<TooN::Vector<ObsDim> >::iterator residual_iter
          = residual_vec.begin();

      for (typename std::vector<Obs>::const_iterator obs_iter=obs_vec.begin();
      obs_iter!=obs_vec.end();
      ++obs_iter, ++residual_iter)
      {
        int i_p = obs_iter->point_id;
        int i_f = obs_iter->frame_id;
        const TooN::Vector<PointParNum> & point= AT(point_list,i_p);
        const Frame & frame = AT(frame_list,i_f);

        TooN::Vector<ObsDim> delta
            = obs_iter->obs - prediction.map(frame, point);

        double sum_2=0;

        for (uint i=0; i<ObsDim;++i)
        {
          sum_2 += Po2(delta[i]);
        }

        double nrm = std::max(0.00000000001,sqrt(sum_2));
        double w = sqrt(kernel(nrm))/nrm;

        TooN::Vector<2> delta_r;
        for (uint i=0; i<ObsDim;++i)
        {
          delta_r[i] = w*delta[i];
        }
        sum +=  delta_r*obs_iter->lambda*delta_r;
        *residual_iter = delta_r;
      }
      return sum;
    }


    /** pseudo-huber cost function */
    double kernel(double delta)
    {
      double b = ba_params.kernel_param;
      return fabs(2*Po2(b)*(sqrt(1+Po2(delta/b))-1));
    }


    /** incremental updates of points */
    void addToPoints(const _PointVec & point_list,
                     const TooN::Vector<> & delta,
                     const _AbsJac & prediction,
                     _PointVec & new_point_list)
    {
      int i_par = 0;


      for (int i_p=0; i_p<num_fix_points; ++i_p){
        AT(new_point_list,i_p) = AT(point_list,i_p);
      }
      for (int i_p=num_fix_points; i_p<num_points; ++i_p){

        AT(new_point_list,i_p)
            =  prediction.add(AT(point_list,i_p),delta.slice(i_par,PointDof));
        i_par+=PointDof;
      }
      assert(i_par==delta.size());
    }


    /** incremental updates of frames */
    void addToFrames(const std::vector<Frame > & frame_list,
                     const TooN::Vector<> & delta,
                     const _AbsJac & prediction,
                     std::vector<Frame > & new_frame_list)
    {
      int i_par = 0;

      for (int i_f=0;i_f<num_fix_frames;++i_f)
      {
        AT(new_frame_list,i_f) = AT(frame_list,i_f);
      }
      for (int i_f=num_fix_frames;i_f<num_frames;++i_f){
        AT(new_frame_list,i_f)
            = prediction.add((AT(frame_list,i_f)), delta.slice(i_par,FrameDof));

        i_par+=FrameDof;
      }
      assert(i_par==delta.size());
    }


    void calcVeps(const std::vector<Obs >  & obs_vec,
                  const std::vector<TooN::Vector<ObsDim> >  & residual_vec,
                  const _JacVec & jac_data_vec,
                  std::vector<TooN::Matrix<PointDof,PointDof> > & V,
                  TooN::Vector<>  & eps_point)
    {
      typename  std::vector<Obs >::const_iterator obs_iter = obs_vec.begin();
      typename _JacVec::const_iterator
          jac_iter = jac_data_vec.begin();

      for (typename std::vector<TooN::Vector<ObsDim> >::const_iterator err_iter
           = residual_vec.begin();
      err_iter != residual_vec.end();
      ++err_iter , ++jac_iter, ++obs_iter)
      {
        int point_id = obs_iter->point_id - num_fix_points;
        if (point_id>=0)
        {
          const TooN::Matrix<ObsDim,PointDof> & J_point = jac_iter->J_point;
          AT(V,point_id) += J_point.T() * J_point;
          eps_point.slice(point_id*PointDof,PointDof) += J_point.T() * (*err_iter);
        }
      }

    }


    void calcUeps(const std::vector<Obs >  & obs_vec,
                  const std::vector<TooN::Vector<ObsDim> >  & residual_vec,
                  const _JacVec & jac_data_vec,
                  std::vector<TooN::Matrix<FrameDof,FrameDof> >  & U,
                  std::vector<TooN::Vector<FrameDof> >  & eps_frame)
    {
      typename std::vector<Obs >::const_iterator obs_iter = obs_vec.begin();
      typename _JacVec::const_iterator jac_iter = jac_data_vec.begin();

      for (typename std::vector<TooN::Vector<ObsDim> >::const_iterator err_iter
           = residual_vec.begin();
      err_iter != residual_vec.end();
      ++err_iter , ++jac_iter, ++obs_iter)
      {
        int frame_id = obs_iter->frame_id - num_fix_frames;

        if(frame_id>=0)
        {
          const TooN::Matrix<ObsDim,FrameDof> & J_frame = jac_iter->J_frame;
          AT(U,frame_id) += J_frame.T() * J_frame;
          AT(eps_frame,frame_id) += J_frame.T() * (*err_iter);
        }
      }
    }


    void calcUVeps(const std::vector<IdObs<ObsDim>  >  & obs_vec,
                   const std::vector<TooN::Vector<ObsDim> >  & residual_vec,
                   const _JacVec & jac_data_vec,
                   std::vector<TooN::Matrix<FrameDof,FrameDof> >  & U,
                   std::vector<TooN::Vector<FrameDof> >  & eps_frame,
                   std::vector<TooN::Matrix<PointDof,PointDof> > & V,
                   TooN::Vector<>  & eps_point)
    {
      typename std::vector<Obs >::const_iterator obs_iter = obs_vec.begin();
      typename _JacVec::const_iterator jac_iter = jac_data_vec.begin();

      for (typename std::vector<TooN::Vector<ObsDim> >::const_iterator err_iter
           = residual_vec.begin();
      err_iter != residual_vec.end();
      ++err_iter , ++jac_iter, ++obs_iter)
      {
        int frame_id = obs_iter->frame_id - num_fix_frames;
        int point_id = obs_iter->point_id - num_fix_points;
        if(frame_id>=0)
        {
          const TooN::Matrix<ObsDim,FrameDof> & J_frame = jac_iter->J_frame;
          AT(U,frame_id) += J_frame.T() * J_frame;
          AT(eps_frame,frame_id) += J_frame.T() * (*err_iter);
        }
        if (point_id>=0)
        {
          const TooN::Matrix<ObsDim,PointDof> & J_point = jac_iter->J_point;
          AT(V,point_id) += J_point.T() * J_point;
          eps_point.slice(point_id*PointDof,PointDof) += J_point.T() * (*err_iter);
        }
      }

    }




    void calcUVeps(const std::vector<IdObsLambda<ObsDim> >  & obs_vec,
                   const std::vector<TooN::Vector<ObsDim> >  & residual_vec,
                   const _JacVec & jac_data_vec,
                   std::vector<TooN::Matrix<FrameDof,FrameDof> >  & U,
                   std::vector<TooN::Vector<FrameDof> >  & eps_frame,
                   std::vector<TooN::Matrix<PointDof,PointDof> > & V,
                   TooN::Vector<>  & eps_point)
    {
      typename std::vector<IdObsLambda<ObsDim> >::const_iterator obs_iter
          = obs_vec.begin();
      typename _JacVec::const_iterator jac_iter = jac_data_vec.begin();

      for (typename std::vector<TooN::Vector<ObsDim> >::const_iterator err_iter
           = residual_vec.begin();
      err_iter != residual_vec.end();
      ++err_iter , ++jac_iter, ++obs_iter)
      {
        int frame_id = obs_iter->frame_id - num_fix_frames;
        int point_id = obs_iter->point_id - num_fix_points;
        if(frame_id>=0)
        {
          const TooN::Matrix<ObsDim,FrameDof> & J_frame = jac_iter->J_frame;
          AT(U,frame_id) += J_frame.T() * obs_iter->lambda * J_frame;
          AT(eps_frame,frame_id)
              += J_frame.T() * obs_iter->lambda * (*err_iter);


        }
        if (point_id>=0)
        {
          const TooN::Matrix<ObsDim,PointDof> & J_point = jac_iter->J_point;
          AT(V,point_id) += J_point.T() * obs_iter->lambda * J_point;
          eps_point.slice(point_id*PointDof,PointDof)
              += J_point.T() * obs_iter->lambda * (*err_iter);
        }
      }
    }


    _JacVec jac_data_vec;
    std::vector<TooN::Vector<ObsDim> > residual_vec;
    BundleAdjusterParams  ba_params;

    int num_points;
    int num_frames;
    int num_obs;
    int num_fix_frames;
    int num_fix_points;
  };

}

#undef AT

#endif // RV_BUNDLE_ADJUSTER_H
