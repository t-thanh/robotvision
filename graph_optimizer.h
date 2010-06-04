/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * grap_optimizer.h is part of RobotVision.
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


#ifndef RV_GRAPH_OPTIMIZER_H
#define RV_GRAPH_OPTIMIZER_H


#include <vector>
#include <list>
#include <set>

#include <TooN/se3.h>
#include <TooN/Cholesky.h>
#include <TooN/helpers.h>

#include "sim3.h"
#include "transformations.h"
#include "sparse_cholesky.h"

namespace RobotVision


{
  /** Compare two monocular SLAM trajctories (SE3) by minimising their
      *  difference wrt. the ambigious scale factor s. */
  class SE3CompareModScale
  {
  public:
    SE3CompareModScale()
    {
      verbose = 0;
    }

    double optimize(const std::vector<TooN::SE3<> > & trans_vec1,
                    const std::vector<TooN::SE3<> > & trans_vec2,
                    double & s,
                    int num_iter)
    {
      int num_trans = trans_vec1.size();


      double new_s;
      TooN::Vector<>residual_vec (num_trans*3);


      TooN::Matrix<> J = getJac(trans_vec1, trans_vec2, s);


      double res = getResidual(trans_vec1, trans_vec2, s, residual_vec);
      if (isnan(res))
      {
        std::cerr << "res is NAN\n";
        exit(-1);
      }

      TooN::Matrix<> A = J.T()*J;
      TooN::Vector<> g = -(J.T()*residual_vec);


      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;
      double mu = 0.0000000000000000001;

      if(verbose>0)
        std::cout << "res: " << res << std::endl;

      for (int i_g=0; i_g<num_iter; ++i_g){
        double rho = 0;
        if (verbose>0)
          std::cout << "iteration: "
              << i_g << " of "<< num_iter << std::endl;
        do
        {
          TooN::Matrix<> A_mu = A + TooN::Identity(A.num_cols())*mu;


          TooN::Cholesky<> Ch(A_mu);

          TooN::Vector<> delta = Ch.backsub(g);

          if (verbose>0)
            std::cout << "mu: " <<mu<< std::endl;

          if (verbose>1)
            std::cout << "delta: " << delta << std::endl;

          new_s = s+delta[0];


          double res_new = getResidual(trans_vec1,
                                       trans_vec2,
                                       new_s,
                                       residual_vec);

          rho = (res-res_new )/(delta*(mu*delta+g));

          if(rho>0)
          {
            if(verbose>0)
              std::cerr << "res_new: " << res_new << std::endl;
            s = new_s;

            res = res_new ;

            TooN::Matrix<> J = getJac(trans_vec1, trans_vec2, s);

            A = J.T()*J;

            g = -(J.T()*residual_vec);
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
      return res;
    }

    int verbose;

  private:
    TooN::Matrix<3,1> singleJac(const TooN::SE3<> & pose1,

                                const TooN::SE3<> & pose2,
                                double s)
    {
      double h = 0.000000000001;
      TooN::Matrix<3,1> J
          = TooN::Zeros;
      TooN::Vector<3> fun = diff(pose1,pose2,s);
      J.T()[0] = (diff(pose1,pose2,s+h) -fun)/h ;
      return J;
    }


    TooN::Matrix<> getJac(const std::vector<TooN::SE3<> > & trans_vec1,
                          const std::vector<TooN::SE3<> > & trans_vec2,
                          const double & s)
    {
      uint num_trans = trans_vec1.size();

      TooN::Matrix<> J(num_trans*3, 1);

      for (uint i=0 ; i<num_trans; ++i)
      {
        J.slice(i*3,0,3,1) = singleJac(trans_vec1[i], trans_vec2[i], s);
      }

      return J;
    }


    TooN::Vector<3> diff(const TooN::SE3<> & pose1,
                         const TooN::SE3<> & pose2,
                         const double & s)
    {
      return pose1.inverse().get_translation()
          - s*pose2.inverse().get_translation();
    }


    double getResidual(const std::vector<TooN::SE3<> > & trans_vec1,
                       const std::vector<TooN::SE3<> > & trans_vec2,
                       const double & s,
                       TooN::Vector<> & residual_vec)
    {
      double sum = 0;
      int trans_id=0;
      for (uint i=0 ; i<trans_vec1.size(); ++i)
      {

        TooN::Vector<3> delta = diff(trans_vec1[i], trans_vec2[i], s);
        residual_vec.slice(trans_id*3,3) = delta;
        sum += delta*delta;
        ++trans_id;
      }
      return sum;
    }
  };

  /** Helper class for pose-graph optimisation  (see below)*/
  template <typename Trans, int TransDoF> class Constraint
  {
  public:
    Constraint(int trans_id1,
               int trans_id2,
               const Trans & mean,
               const TooN::Matrix<TransDoF,TransDoF> & fisher_information)
                 : trans_id1(trans_id1),
                 trans_id2(trans_id2),
                 mean(mean),
                 fisher_information(fisher_information)
    {}

    int trans_id1;
    int trans_id2;
    Trans mean;
    TooN::Matrix<TransDoF,TransDoF> fisher_information;
  };



  /** This class perfoms pose-graph optimisation as described in
      * >H. Strasdat, J.M.M. Montiel, A.J. Davison: "Scale Drift-Aware
      *  Large Scale Monocular SLAM", Proc. of Robotics: Science and
      *  Systems (RSS), 2010.<
      *
      * The concrete implementation of Levenberg-Marquardt (LM) is
      * partially inspired by
      * >M.I. A. Lourakis and A.A. Argyros, "The Design and Implementation
      *  of a Generic Sparse Bundle Adjustment Software Package Based on
      *  the Levenberg-Marquardt Algorithm", Technical Report, 2004.<
      *
      * Trans: pose transformation (e.g. SE2, SE3, Sim3...)
      * TrandDoF: DoF of the corresponding transformation
      */
  template <typename Trans, int TransDoF> class GraphOptimizer
  {
    typedef std::list<Constraint<Trans,TransDoF> >
        _ConsList;
  public:
    GraphOptimizer()
    {
      verbose=0;
    }

    /** performs the pose-graph optimisation using LM
          * trans_vec: set of absolute pose transformations
          * contraint_list: set of relative pose-pose contraints
          * num_fix_trans: number of fixed transformations
          * num_iter: maximal number of iterations
          * mu: initial LM parameter which interpolates between
          *     Gauss-Newton (0) and gradient descent (infinity)
          */
    void optimize(std::vector<Trans> & trans_vec,
                  const std::list<Constraint<Trans, TransDoF> >
                  & constraint_list,
                  const AbstractConFun<Trans,TransDoF> & con_fun,
                  const int num_fixed_trans=0,
                  int num_iter=20,
                  double mu=0.0000000001)
    {
      int num_trans = trans_vec.size();
      int num_constraints = constraint_list.size();

      std::vector<Trans > new_trans_vec(num_trans);
      TooN::Vector<>residual_vec (num_constraints*TransDoF);

      std::pair<TooN::SparseMatrix<>,TooN::SparseMatrix<> >
          pair_J_JTLambda
          = getJac(trans_vec,
                   constraint_list,
                   con_fun,
                   num_fixed_trans);

      TooN::SparseMatrix<> & J = pair_J_JTLambda.first;
      TooN::SparseMatrix<> & J_T_Lambda = pair_J_JTLambda.second;


      double res = getResidual(trans_vec,
                               constraint_list,
                               con_fun,
                               residual_vec);
      if (isnan(res))
      {
        std::cerr << "res is NAN\n";
        exit(-1);
      }

      if(verbose>0)
        std::cout << "res: " << res << std::endl;

      TooN::SparseMatrix<> A = J_T_Lambda*J;
      TooN::Vector<> g = -(J_T_Lambda*residual_vec);
      TooN::SparseCholesky<> spCh (A);

      double nu = 2;
      double eps =0.000000000000001;
      bool stop = false;
      int i_g = 0;
      while(i_g<num_iter)
      {
        double rho = 0;
        do
        {
          if (verbose>0)
            std::cout << "iteration: "<< i_g << std::endl;
          if (verbose>0)
            std::cout << "mu: " <<mu<< std::endl;
          ++i_g;
          TooN::Matrix<> III = TooN::Identity(A.num_cols())*mu;
          TooN::SparseMatrix<> sIII(III);
          TooN::SparseMatrix<> A_mu = A + sIII;

          spCh.update(A_mu);
          TooN::Vector<> delta = spCh.backsub(g);
          addTrans(trans_vec,
                   delta,
                   con_fun,
                   new_trans_vec,
                   num_fixed_trans);
          double res_new = getResidual(new_trans_vec,
                                       constraint_list,
                                       con_fun,
                                       residual_vec);
          rho = (res-res_new)/(delta*(mu*delta+g));

          if(rho>0)
          {
            if(verbose>0)
              std::cerr << "res_new: " << res_new << std::endl;
            trans_vec = std::vector<Trans > (new_trans_vec);
            res = res_new;
            std::pair<TooN::SparseMatrix<>,TooN::SparseMatrix<> >
                pair_J_JTLambda
                = getJac(trans_vec,
                         constraint_list,
                         con_fun,
                         num_fixed_trans);

            TooN::SparseMatrix<> & J
                = pair_J_JTLambda.first;
            TooN::SparseMatrix<> & J_T_Lambda
                = pair_J_JTLambda.second;

            A = J_T_Lambda*J;
            g = -(J_T_Lambda*residual_vec);
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
        }while(!(rho>0 ||  stop || i_g>=num_iter));

        if (stop)
          break;
      }
    }

    int verbose;

  private:
    std::pair<TooN::SparseMatrix<>,TooN::SparseMatrix<> >
        getJac(const std::vector<Trans > & trans_vec,
               const _ConsList & constraint_list,
               const AbstractConFun<Trans,TransDoF> & con_fun,
               const int  num_fixed_trans )
    {
      int num_constraints = constraint_list.size();
      int num_trans = trans_vec.size();

      TooN::TripletMatrix
          J(TransDoF*num_constraints,
            TransDoF*(num_trans-num_fixed_trans ));
      TooN::TripletMatrix
          J_T_Lambda(TransDoF*(num_trans-num_fixed_trans ),
                     TransDoF*num_constraints);

      int constraint_id=0;
      for (typename _ConsList::const_iterator it
           = constraint_list.begin();
      it!=constraint_list.end();
      ++it)
      {
        int id1 = it->trans_id1;
        int id2 = it->trans_id2;
        const Trans & mean = it->mean;
        const TooN::Matrix<TransDoF,TransDoF> & inf
            = it->fisher_information;

        TooN::Matrix<TransDoF,TransDoF>  blockJac1
            = con_fun.d_diff_dT1(trans_vec[id1],
                                 mean,
                                 trans_vec[id2]);

        TooN::Matrix<TransDoF,TransDoF>  blockJac2
            = con_fun.d_diff_dT2(trans_vec[id1],
                                 mean,
                                 trans_vec[id2]);

        TooN::Matrix<TransDoF,TransDoF>  blockJac1_T_Lambda
            = blockJac1.T()*inf;
        TooN::Matrix<TransDoF,TransDoF>  blockJac2_T_Lambda
            = blockJac2.T()*inf;

        int col_id =  id1-num_fixed_trans;
        if (col_id>=0)
        {
          J.add(constraint_id*TransDoF,col_id*TransDoF,blockJac1);

          J_T_Lambda.add(col_id*TransDoF,
                         constraint_id*TransDoF,
                         blockJac1_T_Lambda);
        }
        col_id =  id2-num_fixed_trans;
        if (col_id>=0)
        {
          J.add(constraint_id*TransDoF,col_id*TransDoF, blockJac2);
          J_T_Lambda.add(col_id*TransDoF,
                         constraint_id*TransDoF,
                         blockJac2_T_Lambda);
        }
        ++constraint_id;
      }

      return std::make_pair(TooN::SparseMatrix<>(J),
                            TooN::SparseMatrix<>(J_T_Lambda));
    }


    double getResidual(const std::vector<Trans > & trans_vec,
                       const _ConsList & constraint_list,
                       const AbstractConFun<Trans,TransDoF> & con_fun,
                       TooN::Vector<> & residual_vec)
    {
      double sum = 0;
      int constraint_id=0;
      for (typename _ConsList::const_iterator it
           = constraint_list.begin();
      it!=constraint_list.end();
      ++it)

      {
        int id1 = it->trans_id1;
        int id2 = it->trans_id2;
        const Trans & mean = it->mean;
        const TooN::Matrix<TransDoF,TransDoF> & inf
            = it->fisher_information;

        TooN::Vector<TransDoF> delta
            = con_fun.diff(trans_vec[id1],
                           mean,
                           trans_vec[id2]);
        residual_vec.slice(constraint_id*TransDoF,TransDoF) = delta;
        sum += delta*inf*delta;
        ++constraint_id;
      }
      return sum;
    }

    void addTrans(const std::vector<Trans > & trans_vec,
                  const TooN::Vector<> &delta,
                  const AbstractConFun<Trans,TransDoF> & con_fun,
                  std::vector<Trans > & new_trans_vec,
                  int num_fix_trans)
    {
      for (int i=0; i<num_fix_trans; ++i)
        new_trans_vec[i] = trans_vec[i];

      for (uint i=num_fix_trans; i< trans_vec.size(); ++i)
      {
        new_trans_vec[i]
            = con_fun.add(trans_vec[i],
                          delta.slice((i-num_fix_trans)*TransDoF, TransDoF));

      }
    }
  };

}


#endif // RV_GRAPH_OPTIMIZER_H
