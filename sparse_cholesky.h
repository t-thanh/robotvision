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

#include <list>
#include "math.h"

#include "sparse_matrix.h"

#include <stdexcept>

extern "C"{
#include "CSparse/cs.h"
}

namespace TooN
{

  class NotPosSemiDefException : public std::runtime_error
  {
  public:
    NotPosSemiDefException()
      : std::runtime_error("Not positive semi-definite") { }
  };


  template<int Size=TooN::Dynamic>
  class SparseCholesky
  {
  public:

    SparseCholesky(const SparseMatrix<Size> & s_M) : s_M(s_M)
    {
      symbolic_structure = NULL;
      numeric_structure = NULL;

      symbolicDecomposition();
      decomposition();
      if (numeric_structure==NULL)
      {
        throw NotPosSemiDefException();
      }
    }


    TooN::Matrix<Size,Size> get_L()
    {
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
    }


    TooN::Vector<Size> backsub (const TooN::Vector<Size>& x) const
    {
      TooN::Vector<Size> tmp = x;
      TooN::Vector<Size> sol = x;

      cs_ipvec(symbolic_structure->pinv,&x[0],&tmp[0],x.size());
      //permute con. pivoting
      cs_lsolve(numeric_structure->L,&tmp[0]);
      cs_ltsolve(numeric_structure->L,&tmp[0]);
      cs_pvec(symbolic_structure->pinv,&tmp[0],&sol[0],x.size());
      //unpermute con. pivoting
      return sol;
    }


    TooN::Matrix<Size,Size> backsub (const TooN::Matrix<Size,Size>& M) const {

      TooN::Matrix<Size,Size> res(s_M.sparse_matrix.m,s_M.sparse_matrix.m);

      for (int i=0;i<M.num_cols(); ++i)
      {
        res.T()[i] = backsub(M.T()[i]);
      }

      return res;
    }


    TooN::Matrix<Size,Size> get_inverse()const {

      TooN::Matrix<Size,Size> I = TooN::Identity(s_M.sparse_matrix.m);


      return backsub(I);
    }


    void update(const SparseMatrix<Size,Size> & other)
        //ATTENTION: other has to have the same sparseness structure!!
    {
      assert(s_M.sparse_matrix.nzmax == other.sparse_matrix.nzmax);
      assert(s_M.sparse_matrix.n == other.sparse_matrix.n);

      s_M.copy(&(other.sparse_matrix));

      cs_nfree(numeric_structure);
      numeric_structure = NULL;

      decomposition();
      if (numeric_structure==NULL)
      {
        throw NotPosSemiDefException();
      }
    }


    ~SparseCholesky()
    {
      cs_nfree(numeric_structure);

      cs_sfree(symbolic_structure);
    }

  private:

    css * symbolic_structure;
    csn * numeric_structure;
    SparseMatrix<Size> s_M;

    void symbolicDecomposition()
    {
      symbolic_structure = cs_schol (1, &(s_M.sparse_matrix)) ;
    }

    void decomposition()
    {
      numeric_structure = cs_chol (&(s_M.sparse_matrix),symbolic_structure) ;
    }
  };
}

#endif // RV_SPARSE_CHOLESKY_H
