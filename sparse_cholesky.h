/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * sparse_cholesky.h is part of RobotVision.
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


#ifndef RV_SPARSE_CHOLESKY_H
#define RV_SPARSE_CHOLESKY_H

#ifdef RV_SUITESPARSE_SUPPORT
#include "SuiteSparseQR.hpp"
#endif



#include <list>
#include "math.h"

#include "sparse_matrix.h"

#include <stdexcept>



namespace RobotVision
{

  class NotPosSemiDefException : public std::runtime_error
  {
  public:
    NotPosSemiDefException()
      : std::runtime_error("Not positive semi-definite") { }
  };


  template<int Size=TooN::Dynamic>
  class SparseSolver
  {
  public:

    SparseSolver(const SparseMatrix<Size> & s_M) : s_M(s_M)
    {

#ifdef RV_SUITESPARSE_SUPPORT
      assert(s_M.data->stype==0||s_M.data->stype==1);
      if (s_M.data->stype==0)
      {
        use_qr = true;
      }
      else
      {
        use_qr = false;
      }


      if (!use_qr)
      {
        factor = cholmod_l_analyze(s_M.data,
                                   &(CholmodSingleton::getInstance().data));

        cholmod_l_factorize(s_M.data,
                            factor,
                            &(CholmodSingleton::getInstance().data));
        int status= CholmodSingleton::getInstance().data.status;

        if (status == CHOLMOD_NOT_POSDEF)
        {
          throw NotPosSemiDefException();
        }
      }

#else

      symbolic_structure = NULL;
      numeric_structure = NULL;

      symbolicDecomposition();
      decomposition();
      if (numeric_structure==NULL)
      {
        throw NotPosSemiDefException();
      }
#endif
    }

    TooN::Matrix<Size,Size> get_L()
    {
#ifdef RV_SUITESPARSE_SUPPORT
      TooN::Matrix<Size,Size> L
          = TooN::Zeros(factor->n,factor->n);

      int col_idx = -1;

      for (int r=0;r<(int)(s_M.data->nzmax); ++r)
      {
        if (reinterpret_cast<int *>(factor->p)[col_idx+1]==r)
          ++col_idx;
        L(reinterpret_cast<int *>(factor->i)[r],col_idx)
            = reinterpret_cast<double *>(factor->x)[r];
      }

      return L;
#else
      TooN::Matrix<Size,Size> L
          = TooN::Zeros(numeric_structure->L->m,numeric_structure->L->n);

      int col_idx = -1;

      for (int r=0;r<s_M.sparse_matrix.nzmax; ++r)
      {
        if (numeric_structure->L->p[col_idx+1]==r)
          ++col_idx;
        L(numeric_structure->L->i[r],col_idx) = numeric_structure->L->x[r];
      }

      return L;
#endif
    }


    TooN::Vector<Size> backsub (const TooN::Vector<Size>& b) const
    {
#ifdef RV_SUITESPARSE_SUPPORT
      DenseVector dense_b(b);

      cholmod_dense * dres;

      if (use_qr)
      {
        dres = SuiteSparseQR<double>(s_M.data,
                                     dense_b.data,
                                     &(CholmodSingleton::getInstance().data));
      }
      else
      {
        dres = cholmod_l_solve(CHOLMOD_A,
                               factor,
                               dense_b.data,
                               &(CholmodSingleton::getInstance().data));
      }
      DenseVector res(dres);



      return res.vec();
#else
      TooN::Vector<Size> tmp = b;
      TooN::Vector<Size> sol = b;

      cs_ipvec(symbolic_structure->pinv,&b[0],&tmp[0],b.size());
      //permute con. pivoting
      cs_lsolve(numeric_structure->L,&tmp[0]);
      cs_ltsolve(numeric_structure->L,&tmp[0]);
      cs_pvec(symbolic_structure->pinv,&tmp[0],&sol[0],b.size());
      //unpermute con. pivoting
      return sol;
#endif
    }


    TooN::Matrix<Size,Size> backsub (const TooN::Matrix<Size,Size>& M) const
    {

      TooN::Matrix<Size,Size> res(s_M.num_rows(),s_M.num_rows());

      for (int i=0;i<M.num_cols(); ++i)
      {
        res.T()[i] = backsub(M.T()[i]);
      }

      return res;
    }


    TooN::Matrix<Size,Size> get_inverse()const {

      TooN::Matrix<Size,Size> I = TooN::Identity(s_M.num_rows());


      return backsub(I);
    }


    ~SparseSolver()
    {
#ifdef RV_SUITESPARSE_SUPPORT
      if (!use_qr)
        cholmod_l_free_factor(&factor,&(CholmodSingleton::getInstance().data));

#else
      cs_nfree(numeric_structure);
      cs_sfree(symbolic_structure);
#endif
    }

  private:
    SparseMatrix<Size> s_M;
#ifdef RV_SUITESPARSE_SUPPORT
    cholmod_factor * factor;
    bool use_qr;
#else
    css * symbolic_structure;
    csn * numeric_structure;


    void symbolicDecomposition()
    {
      symbolic_structure = cs_schol (1, &(s_M.sparse_matrix)) ;
    }

    void decomposition()
    {
      numeric_structure = cs_chol (&(s_M.sparse_matrix),symbolic_structure) ;
    }
#endif
  };
}

#endif // RV_SPARSE_CHOLESKY_H
