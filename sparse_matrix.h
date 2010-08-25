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
#include <map>
#include <stdio.h>
#include <string.h>
#include "maths_utils.h"

extern "C"{
#ifdef RV_SUITESPARSE_SUPPORT
#include "cholmod.h"

#else
#include "CSparse/cs.h"
#endif
}



namespace RobotVision
{




  extern int global_int;


  template<int Size>
  class SparseSolver;

  template<int Size>
  class SparseMatrix;



#ifdef RV_SUITESPARSE_SUPPORT
  class CholmodSingleton
  {
  private:
    static CholmodSingleton _instance;


    CholmodSingleton() {
      cholmod_l_start(&data);
      data.print_function = printf;
      data.final_ll = true;
    }


    ~CholmodSingleton() {
      cholmod_l_finish(&data);
    }


    CholmodSingleton(const CholmodSingleton &);
    CholmodSingleton & operator=(const CholmodSingleton &);

  public:
    static CholmodSingleton &getInstance();
    cholmod_common data;
  };

#endif



#ifdef RV_SUITESPARSE_SUPPORT

  class DenseVector
  {

  public:
    DenseVector(const TooN::Vector<> & d)
    {
      data = cholmod_l_allocate_dense(d.size(),
                                      1,
                                      d.size(),
                                      CHOLMOD_REAL,
                                      &(CholmodSingleton::getInstance().data));
      
      memcpy(data->x,&(d[0]),sizeof(double)*d.size());
      
      own_data = true;
    }

    DenseVector(cholmod_dense * d)
    {
      assert(d->dtype == CHOLMOD_DOUBLE);
      assert(d->ncol == 1);
      assert(d->d = 1);
      assert(d->xtype==CHOLMOD_REAL);
      data = d;
      own_data = true;
    }

    TooN::Vector<> vec()
    {

      TooN::Vector<> res(data->nrow);

      for (uint i=0; i<data->nrow; ++i)
      {
        res[i] = reinterpret_cast<double *>(data->x)[i];
      }
      return res;
    }

    ~DenseVector()
    {
      if (own_data)
      {
        cholmod_l_free_dense(&data,&(CholmodSingleton::getInstance().data));
      }
      else
      {
        delete(data);
      }
    }

    cholmod_dense * data;
  private:
    bool own_data;
  };
#endif

  template<int BlockDim>
  class RowBlockMapVec
    : public std::vector<std::map<int,TooN::Matrix<BlockDim,BlockDim> > >
  {
  public:
    RowBlockMapVec(int num_cols)
      : std::vector<std::map<int,TooN::Matrix<BlockDim,BlockDim> > > (num_cols),
      num_elems(0)
    {
    }

    

    inline void add(const TooN::Matrix<BlockDim,BlockDim> & m,
                    uint row_idx,
                    uint col_idx)
    {
      assert(col_idx>=0);
      assert(col_idx<this->size());
      std::map<int,TooN::Matrix<BlockDim,BlockDim> > & col_ref
          = (*this)[col_idx];
      typename std::map<int,TooN::Matrix<BlockDim,BlockDim> >
          ::iterator map_it = col_ref.find(row_idx);
      if (map_it==col_ref.end())
      {
        col_ref.insert(std::make_pair(row_idx,m));
        num_elems += RobotVision::Po2(BlockDim);
      }
      else
      {
        map_it->second += m;
      }
    }

    inline void add(const TooN::Matrix<BlockDim,BlockDim> & m,
                    const TooN::Matrix<BlockDim,BlockDim> & add_once,
                    uint row_idx,
                    uint col_idx)
    {
      assert(col_idx>=0);
      assert(col_idx<this->size());
      std::map<int,TooN::Matrix<BlockDim,BlockDim> > & col_ref
          = (*this)[col_idx];
      typename std::map<int,TooN::Matrix<BlockDim,BlockDim> >
          ::iterator map_it = col_ref.find(row_idx);
      if (map_it==col_ref.end())
      {
        col_ref.insert(std::make_pair(row_idx,m+add_once));
        num_elems += RobotVision::Po2(BlockDim);
      }
      else
      {
        map_it->second += m;
      }
    }

    int num_elems;

  };


  class RowMapVec
    : public std::vector<std::map<int,double > >
  {
  public:
    RowMapVec(int num_cols)
      : std::vector<std::map<int,double> > (num_cols),
      num_elems(0)
    {
    }



    inline void add(double v,
                    uint row_idx,
                    uint col_idx)
    {
      assert(col_idx>=0);
      assert(col_idx<this->size());
      std::map<int,double > & col_ref
          = (*this)[col_idx];
      std::map<int,double>
          ::iterator map_it = col_ref.find(row_idx);
      if (map_it==col_ref.end())
      {
        col_ref.insert(std::make_pair(row_idx,v));
        ++num_elems;
      }
      else
      {
        map_it->second += v;
      }
    }

    inline void add(double v,
                    double add_to_v,
                    uint row_idx,
                    uint col_idx)
    {
      assert(col_idx>=0);
      assert(col_idx<this->size());
      std::map<int,double > & col_ref
          = (*this)[col_idx];
      std::map<int,double>
          ::iterator map_it = col_ref.find(row_idx);
      if (map_it==col_ref.end())
      {
        col_ref.insert(std::make_pair(row_idx,v+add_to_v));
        ++num_elems;
      }
      else
      {
        map_it->second += v;
      }
    }

    template<int S1,int S2, typename L>
    inline void add(const TooN::Matrix<S1,S2,double,L> & m,
                    uint row_idx,
                    uint col_idx)
    {
      int rows = m.num_rows();
      int cols = m.num_cols();
      assert(col_idx>=0);
      //  assert(col_idx+cols<this->size());


      for (int i_c=0;i_c<cols;++i_c)
      {
        for (int i_r=0;i_r<rows;++i_r)
        {

          add(m(i_r,i_c), row_idx+i_r, col_idx+i_c);
        }
      }
    }

    template<int S1,int S2, typename L>
    inline void add(const TooN::Matrix<S1,S2,double,L> & m,
                    const TooN::Matrix<S1,S2,double,L> & add_to_m,
                    uint row_idx,
                    uint col_idx)
    {
      int rows = m.num_rows();
      int cols = m.num_cols();
      assert(col_idx>=0);
      //  assert(col_idx+cols<this->size());


      for (int i_c=0;i_c<cols;++i_c)
      {
        for (int i_r=0;i_r<rows;++i_r)
        {
          add(m(i_r,i_c),
              add_to_m(i_r,i_c),
              row_idx+i_r,
              col_idx+i_c);
        }
      }
    }


    int num_elems;

  };

  template<int Size=TooN::Dynamic>
  class SparseMatrix
  {
    friend class SparseSolver<Size>;

  public:

    
    template<int Dim>
    SparseMatrix(const RowBlockMapVec<Dim> & C,
                 int matrix_type = 1)
    {
#ifdef RV_SUITESPARSE_SUPPORT

      data = cholmod_l_allocate_sparse(C.size()*Dim,
                                       C.size()*Dim,
                                       C.num_elems,
                                       false,
                                       true,
                                       matrix_type,
                                       CHOLMOD_REAL,
                                       &(CholmodSingleton::getInstance().data));

      int num_entries = 0;

      long long * p_p = reinterpret_cast<long long *>(data->p);
      long long * p_i = reinterpret_cast<long long *>(data->i);
      double * p_x = reinterpret_cast<double *>(data->x);

      for (typename
           std::vector<std::map<int,TooN::Matrix<Dim,Dim> > >::const_iterator
           it_col = C.begin(); it_col!=C.end(); ++it_col)
      {

        for (uint i_c=0; i_c<Dim; ++i_c)
        {
          *p_p = num_entries;
          p_p++;
          for (typename
               std::map<int,TooN::Matrix<Dim,Dim> >::const_iterator
               it_row = it_col->begin(); it_row!=it_col->end(); ++it_row)
          {
            int row_idx = it_row->first*Dim;
            const TooN::Matrix<Dim,Dim> & M = it_row->second;

            for (uint i_r=0; i_r<Dim; ++i_r)
            {
              *p_i = row_idx+i_r;
              p_i++;
              *p_x = M(i_r,i_c);
              p_x++;
            }
            num_entries+=Dim;
          }
        }
      }
      *p_p = num_entries;

#else

      sparse_matrix.m = C.size()*Dim;
      sparse_matrix.n = C.size()*Dim;
      sparse_matrix.nzmax = C.num_elems;
      sparse_matrix.i = (int*)malloc(sizeof(int)*sparse_matrix.nzmax);
      sparse_matrix.p = (int*)malloc(sizeof(int)*(sparse_matrix.n+1));
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*sparse_matrix.nzmax);



      int num_entries = 0;

      int * p_p = sparse_matrix.p;
      int * p_i = sparse_matrix.i;
      double * p_x = sparse_matrix.x;

      for (typename
           std::vector<std::map<int,TooN::Matrix<Dim,Dim> > >::const_iterator
           it_col = C.begin(); it_col!=C.end(); ++it_col)
      {

        for (uint i_c=0; i_c<Dim; ++i_c)
        {
          *p_p = num_entries;
          p_p++;
          for (typename
               std::map<int,TooN::Matrix<Dim,Dim> >::const_iterator
               it_row = it_col->begin(); it_row!=it_col->end(); ++it_row)
          {
            int row_idx = it_row->first*Dim;
            const TooN::Matrix<Dim,Dim> & M = it_row->second;
            for (uint i_r=0; i_r<Dim; ++i_r)
            {
              *p_i = row_idx+i_r;
              p_i++;
              *p_x = M(i_r,i_c);
              p_x++;
            }
            num_entries+=Dim;
          }
        }
      }
      *p_p = num_entries;

      sparse_matrix.nz = -1;
#endif
    }


    SparseMatrix(const RowMapVec & C,
                 int matrix_type = 1)
    {
#ifdef RV_SUITESPARSE_SUPPORT

      data = cholmod_l_allocate_sparse(C.size(),
                                       C.size(),
                                       C.num_elems,
                                       false,
                                       true,
                                       matrix_type,
                                       CHOLMOD_REAL,
                                       &(CholmodSingleton::getInstance().data));

      int num_entries = 0;

      long long * p_p = reinterpret_cast<long long *>(data->p);
      long long * p_i = reinterpret_cast<long long *>(data->i);
      double * p_x = reinterpret_cast<double *>(data->x);

      for (typename
           std::vector<std::map<int,double > >::const_iterator
           it_col = C.begin(); it_col!=C.end(); ++it_col)
      {


        *p_p = num_entries;
        p_p++;
        for (typename
             std::map<int,double >::const_iterator
             it_row = it_col->begin(); it_row!=it_col->end(); ++it_row)
        {
          //int row_idx = it_row->first;
          //const TooN::Matrix<Dim,Dim> & M = it_row->second;

          *p_i = it_row->first;
          p_i++;
          *p_x = it_row->second;
          p_x++;

          ++num_entries;
        }

      }
      *p_p = num_entries;

#else

      sparse_matrix.m = C.size();
      sparse_matrix.n = C.size();
      sparse_matrix.nzmax = C.num_elems;
      sparse_matrix.i = (int*)malloc(sizeof(int)*sparse_matrix.nzmax);
      sparse_matrix.p = (int*)malloc(sizeof(int)*(sparse_matrix.n+1));
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*sparse_matrix.nzmax);



      int num_entries = 0;

      int * p_p = sparse_matrix.p;
      int * p_i = sparse_matrix.i;
      double * p_x = sparse_matrix.x;

      for (typename
           std::vector<std::map<int,double > >::const_iterator
           it_col = C.begin(); it_col!=C.end(); ++it_col)
      {


        *p_p = num_entries;
        p_p++;
        for (typename
             std::map<int,double>::const_iterator
             it_row = it_col->begin(); it_row!=it_col->end(); ++it_row)
        {
          *p_i = it_row->first;
          p_i++;
          *p_x = it_row->second;
          p_x++;

          ++num_entries;
        }

      }
      *p_p = num_entries;

      sparse_matrix.nz = -1;
#endif
    }

    SparseMatrix(const SparseMatrix & other)
    {

#ifdef RV_SUITESPARSE_SUPPORT

      this->data = cholmod_l_copy_sparse(other.data,
                                         &(CholmodSingleton::getInstance().data));
#else

      sparse_matrix.i
          = (int*)malloc(sizeof(int)*other.sparse_matrix.nzmax);
      sparse_matrix.p
          = (int*)malloc(sizeof(int)*(other.sparse_matrix.n+1));
      sparse_matrix.x
          = (double*)malloc(sizeof(double)*other.sparse_matrix.nzmax);
      copy(&other.sparse_matrix);
#endif
    }




#ifdef RV_SUITESPARSE_SUPPORT
    SparseMatrix(cholmod_sparse  * data)
    {
      this->data = data;
    }


#else
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
#endif


    ~SparseMatrix()
    {
#ifdef RV_SUITESPARSE_SUPPORT
      cholmod_l_free_sparse(&data,&(CholmodSingleton::getInstance().data));
#else
      cs_free(sparse_matrix.i);
      cs_free(sparse_matrix.p);
      cs_free(sparse_matrix.x);
#endif
    }

    void operator = (const SparseMatrix & other)
                    {
#ifdef RV_SUITESPARSE_SUPPORT
      cholmod_l_free_sparse(&data,&(CholmodSingleton::getInstance().data));
      this->data = cholmod_l_copy_sparse(other.data);
#else
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
#endif

    }



    SparseMatrix operator + (const SparseMatrix & other) const
    {
#ifdef RV_SUITESPARSE_SUPPORT

      double alpha[2] = {1.,1.};
      double beta[2] = {1.,1.};
      return SparseMatrix(cholmod_add(this->data,
                                      other->data,
                                      alpha,
                                      beta,
                                      true,
                                      true,
                                      &(CholmodSingleton::getInstance().data)));

#else
      cs * sm = cs_add(&(this->sparse_matrix), &(other.sparse_matrix),1,1);
      SparseMatrix SM(sm);
      cs_spfree(sm);
      return SM;
#endif
    }

    SparseMatrix operator * (const SparseMatrix & other) const
    {
#ifdef RV_SUITESPARSE_SUPPORT
      return SparseMatrix(
          cholmod_l_ssmult(this->data,
                           other->data,
                           1,
                           true,
                           true,
                           &(CholmodSingleton::getInstance().data)));
#else
      cs * sm = cs_multiply(&(this->sparse_matrix), &(other.sparse_matrix));
      SparseMatrix SM(sm);
      cs_spfree(sm);
      return SM;
#endif
    }




    void operator += (const SparseMatrix & other)
                     {
#ifdef RV_SUITESPARSE_SUPPORT
      //inefficient implementation!!
      SparseMatrix C = *this + other;
      *this = C;
#else
      cs * sm = cs_add(&(this->sparse_matrix), &(other.sparse_matrix),1,1);
      copy(sm);
      cs_spfree(sm);
#endif
    }

    void operator *= (const SparseMatrix & other)
                     {
#ifdef RV_SUITESPARSE_SUPPORT
      //inefficient implementation!!
      SparseMatrix C = *this * other;
      *this = C;
#else
      cs * sm = cs_multiply(&(this->sparse_matrix), &(other.sparse_matrix));
      copy(sm);
      cs_spfree(sm);
#endif
    }


    SparseMatrix T() const
    {
#ifdef RV_SUITESPARSE_SUPPORT
      return SparseMatrix(
          cholmod_l_transpose(data,
                              true,
                              &(CholmodSingleton::getInstance().data)));;
#else
      cs * sm = cs_transpose(&sparse_matrix,1);
      SparseMatrix SM(sm);
      cs_spfree(sm);
      return SM;
#endif
    }


    TooN::Matrix<Size,Size> get_dense() const
    {  
      int col_idx = -1;

#ifdef RV_SUITESPARSE_SUPPORT
      TooN::Matrix<Size,Size> d_M = TooN::Zeros(data->nrow,data->ncol);
      for (int r=0;r<(int)(data->nzmax); ++r)
      {
        if (reinterpret_cast<int *>(data->p)[col_idx+1]==r)
          ++col_idx;
        d_M( reinterpret_cast<int *>(data->i)[r],col_idx)
            = reinterpret_cast<double *>(data->x)[r];
      }
#else
      TooN::Matrix<Size,Size> d_M
          = TooN::Zeros(sparse_matrix.m,sparse_matrix.n);
      for (int r=0;r<sparse_matrix.nzmax; ++r)
      {
        if (sparse_matrix.p[col_idx+1]==r)
          ++col_idx;
        d_M(sparse_matrix.i[r],col_idx) = sparse_matrix.x[r];
      }
#endif

      return d_M;
    }


    int num_rows() const
    {
#ifdef RV_SUITESPARSE_SUPPORT
      return data->nrow;
#else
      return sparse_matrix.m;
#endif
    }


    int num_cols() const
    {
#ifdef RV_SUITESPARSE_SUPPORT
      return data->ncol;
#else
      return sparse_matrix.n;
#endif
    }


  private:
#ifdef RV_SUITESPARSE_SUPPORT
    cholmod_sparse * data;
#else
    cs sparse_matrix;
#endif

  };
}


#endif // RV_SPARSE_MATRIX_H
