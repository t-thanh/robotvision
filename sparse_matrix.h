/**
 * @author  Hauke Strasdat
 *
 * Copyright (C) 2010  Hauke Strasdat
 *                     Imperial College London
 *
 * sparse_matrix.h is part of RobotVision.
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


#ifndef RV_SPARSE_MATRIX_H
#define RV_SPARSE_MATRIX_H

#include <list>
#include <stdio.h>
#include <string.h>
#include "maths_utils.h"

extern "C"{
#include "CSparse/cs.h"
}


namespace TooN
{

  template<int Size>
  class SparseCholesky;

  template<int Size,int Size2>
  class SparseMatrix;

  /** TripletMatrix is the best way to fill SparseMatrix:
   *  - Create TripletMatrix(num_row,num_cols)
   *  - add non-zero data with add(i,j,v)
   *  - create SparseMatrix from TripletMatrix
   *  This is more memory- and time-efficient than
   * create SparseMatrix from a dense matrix!
   *  (Note: You cannot change data in SparseMatrix, once it is created.)
   */
  class TripletMatrix
  {
    template<int Size, int Size2>
    friend class SparseMatrix;

  public:
    /** create a sparse matrix in the triplet-form */
    TripletMatrix(int num_rows, int num_cols)
    {
      sparse_matrix.nzmax = 1;
      sparse_matrix.m = num_rows;
      sparse_matrix.n = num_cols;
      sparse_matrix.i = (int*)malloc(sizeof(int)*sparse_matrix.nzmax);
      sparse_matrix.p = (int*)malloc(sizeof(int)*(sparse_matrix.n+1));
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*sparse_matrix.nzmax);
      sparse_matrix.nz = 0;
    }

    /** add data to matrix in triplet form*/
    void add(int i, int j, double val)
    {
      int error = cs_entry(&sparse_matrix,i,j,val);
      if (error ==0)
      {
        std::cerr <<"CSparse error code: " << error << std::endl;
        exit(-1);
      }

    }

    /** add data to matrix in triplet form */
    //
    //ToDo:implement this also for Matrix<S1,S2,double,ColMajor>
    template <int S1, int S2>
        void add(int i, int j, const Matrix<S1,S2> & sub_m)
    {
      for (int ii=0; ii<sub_m.num_rows(); ++ii)
        for (int jj=0; jj<sub_m.num_cols(); ++jj)
        {
        int error = cs_entry(&sparse_matrix,i+ii, j+jj,sub_m(ii,jj));
        if (error ==0)
        {
          std::cerr <<"CSparse error code: " << error << std::endl;
          exit(-1);
        }
      }
    }

    ~TripletMatrix()
    {
      cs_free(sparse_matrix.i);
      cs_free(sparse_matrix.p);
      cs_free(sparse_matrix.x);
    }

  private:
    cs sparse_matrix;

  };

  template<int Size=TooN::Dynamic,int Size2=Dynamic>
  class SparseMatrix
  {
    friend class SparseCholesky<Size>;

  public:

    /** preferred constructure , (see above) */
    SparseMatrix(const TripletMatrix & t)
    {
      assert(t.sparse_matrix.nz>0);
      //triplet matrix must have at least one entry!

      cs * sm = cs_compress(&(t.sparse_matrix));
      sparse_matrix.i = (int*)malloc(sizeof(int)*sm->nzmax);
      sparse_matrix.p = (int*)malloc(sizeof(int)*(sm->n+1));
      sparse_matrix.x = (double*)malloc(sizeof(double)*sm->nzmax);
      copy(sm);
      cs_spfree(sm);
    }


    SparseMatrix(const TooN::Matrix<Size,Size2> & C,
                 const TooN::Matrix<Size,Size2,char> & mask)
    {
      std::list<int> row_list;
      std::list<int> col_list;
      std::list<double> content_list;

      for (int c=0; c<C.num_cols(); ++c)
      {
        col_list.push_back(row_list.size());
        for (int r=0; r<C.num_rows(); ++r)
        {
          if (mask(r,c)!=0)
          {
            row_list.push_back(r);
            content_list.push_back(C(r,c));
          }
        }
      }
      col_list.push_back(row_list.size());

      sparse_matrix.m = C.num_rows();
      sparse_matrix.n = C.num_cols();
      sparse_matrix.nzmax = content_list.size();
      sparse_matrix.i = (int*)malloc(sizeof(int)*row_list.size());
      sparse_matrix.p = (int*)malloc(sizeof(int)*col_list.size());
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*content_list.size());

      int r=0;
      for (std::list<int>::iterator it=row_list.begin();
      it!=row_list.end();
      ++it)
      {
        sparse_matrix.i[r] = *it;
        ++r;
      }

      int c=0;
      for (std::list<int>::iterator it=col_list.begin();
      it!=col_list.end();
      ++it)
      {
        sparse_matrix.p[c] = *it;
        ++c;
      }

      int e=0;
      for (std::list<double>::iterator it=content_list.begin();
      it!=content_list.end();
      ++it)
      {
        sparse_matrix.x[e] = *it;
        ++e;
      }

      sparse_matrix.nz = -1;
    }

    SparseMatrix(const TooN::Matrix<Size,Size2> & C)
        //ATTENTION: Make sure that all fill-in cells of the matrix
        //are non-zero!
        //For Example: C(i,j) might be sin(x). Make sure that C(i,j)!=0
        //even for x=0. E.g. set C(i,j)=0.
    {
      std::list<int> row_list;
      std::list<int> col_list;
      std::list<double> content_list;

      for (int c=0; c<C.num_cols(); ++c)
      {
        col_list.push_back(row_list.size());
        for (int r=0; r<C.num_rows(); ++r)
        {
          if (C(r,c)!=0)
          {
            row_list.push_back(r);
            content_list.push_back(C(r,c));
          }
        }
      }
      col_list.push_back(row_list.size());

      sparse_matrix.m = C.num_rows();
      sparse_matrix.n = C.num_cols();
      sparse_matrix.nzmax = content_list.size();
      sparse_matrix.i = (int*)malloc(sizeof(int)*row_list.size());
      sparse_matrix.p = (int*)malloc(sizeof(int)*col_list.size());
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*content_list.size());


      int r=0;
      for (std::list<int>::iterator it=row_list.begin();
      it!=row_list.end();
      ++it)
      {
        sparse_matrix.i[r] = *it;
        ++r;
      }

      int c=0;
      for (std::list<int>::iterator it=col_list.begin();
      it!=col_list.end();
      ++it)
      {
        sparse_matrix.p[c] = *it;
        ++c;
      }

      int e=0;
      for (std::list<double>::iterator it=content_list.begin();
      it!=content_list.end();
      ++it)
      {
        sparse_matrix.x[e] = *it;
        ++e;
      }
      sparse_matrix.nz = -1;
    }


    SparseMatrix(const SparseMatrix & other)
    {
      sparse_matrix.i
          = (int*)malloc(sizeof(int)*other.sparse_matrix.nzmax);
      sparse_matrix.p
          = (int*)malloc(sizeof(int)*(other.sparse_matrix.n+1));
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*other.sparse_matrix.nzmax);
      copy(&other.sparse_matrix);
    }


    SparseMatrix(const cs  * const sm)
    {
      assert(CS_CSC(sm)==true);
      sparse_matrix.i = (int*)malloc(sizeof(int)*sm->nzmax);
      sparse_matrix.p = (int*)malloc(sizeof(int)*(sm->n+1));
      sparse_matrix.x = (double*)malloc(sizeof(double)*sm->nzmax);
      copy(sm);
    }


    void  copy(const cs  * const sm)
    {
      sparse_matrix.m = sm->m;
      sparse_matrix.n = sm->n;
      sparse_matrix.nz = sm->nz;
      sparse_matrix.nzmax = sm->nzmax;

      memcpy(sparse_matrix.i,sm->i,sizeof(int)*sm->nzmax);
      memcpy(sparse_matrix.p,sm->p,sizeof(int)*(sm->n+1));
      memcpy(sparse_matrix.x,sm->x,sizeof(double)*sm->nzmax);
    }


    ~SparseMatrix()
    {
      cs_free(sparse_matrix.i);
      cs_free(sparse_matrix.p);
      cs_free(sparse_matrix.x);
    }

    void operator = (const SparseMatrix & other)
                    {
      cs_free(sparse_matrix.i);
      cs_free(sparse_matrix.p);
      cs_free(sparse_matrix.x);

      sparse_matrix.i
          = (int*)malloc(sizeof(int)*other.sparse_matrix.nzmax);
      sparse_matrix.p
          = (int*)malloc(sizeof(int)*(other.sparse_matrix.n+1));
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*other.sparse_matrix.nzmax);
      copy(&other.sparse_matrix);

    }


    //There might be a more efficient implementation for +/*/+=/*=/T
    // - working memory is allocated in cs_add/multiply...
    // - memory is allocated in SparseMatrix constructor, and data is copied.
    // - working memory is freed afterwards using cs_free.
    // It might be possible to allocate memory only once and avoid copying!
    // At least, the current version should have no memory leak
    // (tested with valgrind!).
    // Same holds for SparseMatrix(const TripletMatrix...).
    SparseMatrix operator + (const SparseMatrix & other) const
    {
      cs * sm = cs_add(&(this->sparse_matrix), &(other.sparse_matrix),1,1);
      SparseMatrix SM(sm);
      cs_spfree(sm);
      return SM;
    }

    SparseMatrix operator * (const SparseMatrix & other) const
    {
      cs * sm = cs_multiply(&(this->sparse_matrix), &(other.sparse_matrix));
      SparseMatrix SM(sm);
      cs_spfree(sm);
      return SM;
    }


    Vector<> operator * (const Vector<> & other) const
    {
      assert(other.size() == num_cols());
      Vector<> res = TooN::Zeros(num_rows());


      const double * y = &(other[0]);
      double * x = &(res[0]);
      cs_gaxpy(&sparse_matrix,y,x);
      return res;
    }


    void operator += (const SparseMatrix & other)
                     {
      cs * sm = cs_add(&(this->sparse_matrix), &(other.sparse_matrix),1,1);
      copy(sm);
      cs_spfree(sm);
    }

    void operator *= (const SparseMatrix & other)
                     {
      cs * sm = cs_multiply(&(this->sparse_matrix), &(other.sparse_matrix));
      copy(sm);
      cs_spfree(sm);
    }


    SparseMatrix T() const
    {
      cs * sm = cs_transpose(&sparse_matrix,1);
      SparseMatrix SM(sm);
      cs_spfree(sm);
      return SM;
    }


    Matrix<Size,Size2> get_dense() const
    {
      Matrix<Size,Size2> d_M = Zeros(sparse_matrix.m,sparse_matrix.n);
      int col_idx = -1;

      for (int r=0;r<sparse_matrix.nzmax; ++r)
      {
        if (sparse_matrix.p[col_idx+1]==r)
          ++col_idx;
        d_M(sparse_matrix.i[r],col_idx) = sparse_matrix.x[r];
      }

      return d_M;
    }


    int num_rows() const
    {
      return sparse_matrix.m;
    }


    int num_cols() const
    {
      return sparse_matrix.n;
    }


  private:
    cs sparse_matrix;

  };
}


#endif // RV_SPARSE_MATRIX_H
