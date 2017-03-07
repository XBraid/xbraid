#ifndef BRAID_PFASST_HPP
#define BRAID_PFASST_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

namespace pfasst
{
   class NotImplementedYet
   : public runtime_error
   {
   public:
      /**
       * @param[in] msg component or algorithm the throwing function is required for.
       */
      explicit NotImplementedYet(const string& msg);
      virtual const char* what() const throw();
   };
   NotImplementedYet::NotImplementedYet(const string& msg)
   : runtime_error(msg)
   {}
   
   const char* NotImplementedYet::what() const throw()
   {
      return (string("Not implemented/supported yet, required for: ") + string(runtime_error::what())).c_str();
   }
   
   namespace quadrature
   {
      class Matrix
      {
      protected:
         double **data;
         int      nrows;
         int      ncols;
      public:
         Matrix(int m, int n);
         Matrix()  = default;
         ~Matrix();
         int get_nrows() const;
         int get_ncols() const;
         double& operator() (int i, int j);
      };
      
      Matrix::Matrix (int m, int n)
      {
         assert(m > 0 && n > 0);
         nrows = m;
         ncols = n;
         data = new double* [nrows];
         for (int i = 0; i < nrows; ++i)
            data[i] = new double [ncols];
         for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j)
               data[i][j] = 0.0;
      }
      //Matrix::Matrix()
      //: nrows(0)
      //{
      //   ncols = 0;
      //}
      Matrix::~Matrix()
      {
         for (int i = 0; i < nrows; ++i)
            delete data[i];
         delete[] data;
      }
      int Matrix::get_nrows() const
      {
         return nrows;
      }
      int Matrix::get_ncols() const
      {
         return ncols;
      }
      double& Matrix::operator() (int i, int j)
      {
         assert (i >= 0 && i < nrows);
         assert (j >= 0 && j < ncols);
         return data[i][j];
      }
      
      class Polynomial
      {
      protected:
         /**
          * Coefficients of the polynomial.
          *
          * The coefficient for the highest degree of \\( x \\) has index `0`.
          * The last coefficient is for \\( x^0 \\).
          */
         vector<double> c;
         
      public:
         //! @{
         Polynomial(size_t n);
         //! @}
         
         //! @{
         /**
          * Order of this polynomial.
          *
          * The order of this polynomial is one less the number of coefficients defined.
          *
          * @returns order of the polynomial
          */
         size_t order() const;
         
         /**
          * Access coefficient @p i
          *
          * @param[in] i coefficient index (zero-based)
          * @returns \\( i+1 \\)-th coefficient
          * @throws std::out_of_range if index is out of bounds, i.e. @p i >= Polynomial::order()
          *
          * @note The coefficients are stored in reversed order to the degree of the indeterminate,
          *   i.e. `i=0` corresponds to \\( c_n \\) while `i=n` corresponds to \\( c_0 \\).
          *   See also Polynomial::c.
          */
         double& operator[](const size_t i);
         
         /**
          * Differentiate this polynomial.
          *
          * Computes standard differential of this polynomial.
          *
          * @returns differentiated polynomial
          */
         Polynomial differentiate() const;
         
         /**
          * Integrates this polynomial.
          *
          * Computes integral of this polynomial.
          *
          * @returns integrated polynomial
          */
         Polynomial integrate() const;
         //! @}
         
         //! @{
         /**
          * Evaluate polynomial for given value.
          *
          * @tparam xtype numerical type of the value
          * @param[in] x value to evaluate polynomial at
          * @returns value of polynomial at @p x
          */
         template<typename xtype>
         xtype evaluate(const xtype x) const
         {
            size_t n = this->order();
            xtype v = c[n];
            for (int j = n - 1; j >= 0; j--) {
               v = x * v + c[j];
            }
            return v;
         }
         
         /**
          * Normalizes this polynomial with respect to \\( c_0 \\).
          *
          * @returns normalized polynomial
          */
         Polynomial normalize() const;
         
         /**
          * Computes the roots of this polynomial.
          *
          * @returns roots sorted with respect to their value
          */
         vector<double> roots(size_t num_iterations=100, double ztol=1.0e-20) const;
         //! @}
         
         //! @{
         /**
          * Computes the Legendre polynomial of given order.
          *
          * @param[in] order desired order of the Legendre polynomial
          * @returns Legendre polynomial of order @p order
          */
         static Polynomial legendre(const size_t order);
         //! @}
      };
      
      Polynomial::Polynomial(size_t n)
      : c(n, 0.0)
      {}
      
      size_t Polynomial::order() const
      {
         return c.size() - 1;
      }
      
      double& Polynomial::operator[](const size_t i)
      {
         return c.at(i);
      }
      
      Polynomial Polynomial::differentiate() const
      {
         Polynomial p(c.size() - 1);
         for (size_t j = 1; j < c.size(); j++) {
            p[j - 1] = j * c[j];
         }
         return p;
      }
      
      Polynomial Polynomial::integrate() const
      {
         Polynomial p(c.size() + 1);
         for (size_t j = 0; j < c.size(); j++) {
            p[j + 1] = c[j] / (j + 1);
         }
         return p;
      }
      
      Polynomial Polynomial::normalize() const
      {
         Polynomial p(c.size());
         for (size_t j = 0; j < c.size(); j++) {
            p[j] = c[j] / c.back();
         }
         return p;
      }
      
      vector<double> Polynomial::roots(size_t num_iterations, double ztol) const
      {
         assert(c.size() >= 1);
         size_t n = c.size() - 1;
         
         // initial guess
         vector<complex<double> > z0(n);
         for (size_t j = 0; j < n; j++) {
            z0[j] = pow(complex<double>(0.4, 0.9), j);
         }
         
         // durand-kerner-weierstrass iterations
         Polynomial p = this->normalize();
         for (size_t k = 0; k < num_iterations; k++) {
            complex<double> num, den;
            for (size_t i = 0; i < n; i++) {
               num = p.evaluate(z0[i]);
               den = 1.0;
               for (size_t j = 0; j < n; j++) {
                  if (j == i) { continue; }
                  den = den * (z0[i] - z0[j]);
               }
               z0[i] = z0[i] - num / den;
            }
         }
         
         vector<double> roots(n);
         for (size_t j = 0; j < n; j++) {
            roots[j] = abs(z0[j]) < ztol ? 0.0 : real(z0[j]);
         }
         
         sort(roots.begin(), roots.end());
         return roots;
      }
      
      Polynomial Polynomial::legendre(const size_t order)
      {
         if (order == 0) {
            Polynomial p(1);
            p[0] = 1.0;
            return p;
         }
         
         if (order == 1) {
            Polynomial p(2);
            p[0] = 0.0;
            p[1] = 1.0;
            return p;
         }
         
         Polynomial p0(order + 1), p1(order + 1), p2(order + 1);
         p0[0] = 1.0; p1[1] = 1.0;
         
         // (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
         for (size_t m = 1; m < order; m++) {
            for (size_t j = 1; j < order + 1; j++) {
               p2[j] = ((2 * m + 1) * p1[j - 1] - m * p0[j]) / (m + 1);
            }
            p2[0] = - int(m) * p0[0] / (m + 1);
            
            for (size_t j = 0; j < order + 1; j++) {
               p0[j] = p1[j];
               p1[j] = p2[j];
            }
         }
         
         return p2;
      }
      
      static Polynomial build_polynomial(const size_t node, const vector<double>& nodes)
      {
         const size_t num_nodes = nodes.size();
         Polynomial p(num_nodes + 1), p1(num_nodes + 1);
         p[0] = 1.0;
         
         for (size_t m = 0; m < num_nodes; ++m) {
            if (m == node) { continue; }
            
            // p_{m+1}(x) = (x - x_j) * p_m(x)
            p1[0] = 0.0;
            for (size_t j = 0; j < num_nodes;     ++j) { p1[j + 1]  = p[j]; }
            for (size_t j = 0; j < num_nodes + 1; ++j) { p1[j]     -= p[j] * nodes[m]; }
            for (size_t j = 0; j < num_nodes + 1; ++j) { p[j]       = p1[j]; }
         }
         
         return p;
      }
      
      
      /**
       * Compute quadrature matrix \\( Q \\) between two sets of nodes.
       *
       * Computing the quadrature matrix \\( Q \\) for polynomial-based integration from one set of
       * quadrature nodes (@p from) to another set of quadrature nodes (@p to).
       *
       * @tparam scalar precision of quadrature (i.e. `double`)
       * @param[in] from first set of quadrature nodes
       * @param[in] to second set of quadrature nodes
       * @returns quadrature matrix \\( Q \\) with `to.size()` rows and `from.size()` colums
       *
       * @pre For correctness of the algorithm it is assumed, that both sets of nodes are in the range
       *   \\( [0, 1] \\).
       *
       * @since v0.3.0
       */
      static Matrix *compute_q_matrix(const vector<double>& from, const vector<double>& to)
      {
         const size_t to_size = to.size();
         const size_t from_size = from.size();
         assert(to_size >= 1 && from_size >= 1);
         
         //std::cout << "compute_q_matrix(): matrix size " << from_size << "x" << to_size << std::endl;
         Matrix *q_mat = new Matrix::Matrix(to_size, from_size);
         for (size_t m = 0; m < from_size; ++m) {
            Polynomial p = build_polynomial(m, from);
            // evaluate integrals
            auto den = p.evaluate(from[m]);
            auto P = p.integrate();
            for (size_t j = 0; j < to_size; ++j) {
               (*q_mat)(j, m) = (P.evaluate(to[j]) - P.evaluate(0.0)) / den;
            }
         }
         
//         // print computed matrix
//         for (size_t m = 0; m < from_size; ++m)
//         {
//            for (size_t j = 0; j < to_size; ++j)
//               std::cout << (*q_mat)(j,m) << " ";
//            std::cout << std::endl;
//         }
         
         return q_mat;
      }
      
      
      /**
       * Compute quadrature matrix \\( Q \\) for one set of nodes.
       *
       * @tparam scalar precision of quadrature (i.e. `double`)
       * @param[in] nodes quadrature nodes to compute \\( Q \\) matrix for
       *
       * @since v0.3.0
       *
       * @overload
       */
      static Matrix *compute_q_matrix(const vector<double>& nodes)
      {
         return compute_q_matrix(nodes, nodes);
      }

      
      // Quadrature
      class IQuadrature
      {
      protected:
         static const bool LEFT_IS_NODE = false;
         static const bool RIGHT_IS_NODE = false;
         
         size_t num_nodes;
         Matrix *q_mat;
         vector<double> nodes;
         
      public:
         //! @{
         /**
          * @throws invalid_argument if number of nodes is invalid for quadrature type
          */
         explicit IQuadrature(const size_t num_nodes);
         /**
          * @throws invalid_argument if number of nodes is invalid for quadrature type
          */
         IQuadrature();
         virtual ~IQuadrature() = default;
         //! @}
         
         //! @{
         virtual const Matrix& get_q_mat() const;
         virtual const vector<double>& get_nodes() const;
         virtual size_t get_num_nodes() const;
         
         /**
          * @throws pfasst::NotImplementedYet if not overwritten by implementation;
          *   required for quadrature of any kind
          */
         virtual bool left_is_node() const;
         /**
          * @throws pfasst::NotImplementedYet if not overwritten by implementation;
          *   required for quadrature of any kind
          */
         virtual bool right_is_node() const;
         //! @}
         
         /**
          * Compute a rough estimate of the numerical error... XXX
          */
         double expected_error() const;
         
      protected:
         //! @{
         /**
          * @throws pfasst::NotImplementedYet if not overwritten by implementation;
          *   required for quadrature of any kind
          */
         virtual void compute_nodes();
         virtual void compute_weights();
         //! @}
      };
      IQuadrature::IQuadrature(const size_t num_nodes)
      : num_nodes(num_nodes)
      {
         if (this->num_nodes == 0) {
            throw invalid_argument("Any quadrature requires at least one quadrature nodes.");
         }
      }
      
      IQuadrature::IQuadrature()
      : num_nodes(0)
      {}
      
      const Matrix& IQuadrature::get_q_mat() const
      {
         return *(this->q_mat);
      }
      
      const vector<double>& IQuadrature::get_nodes() const
      {
         return this->nodes;
      }
      
      size_t IQuadrature::get_num_nodes() const
      {
         return this->num_nodes;
      }
      
      bool IQuadrature::left_is_node() const
      {
         throw NotImplementedYet("Quadrature");
         return LEFT_IS_NODE;
      }
      
      bool IQuadrature::right_is_node() const
      {
         throw NotImplementedYet("Quadrature");
         return RIGHT_IS_NODE;
      }
      
      void IQuadrature::compute_nodes()
      {
         throw NotImplementedYet("Quadrature");
      }
      
      void IQuadrature::compute_weights()
      {
         this->q_mat = compute_q_matrix(this->nodes);
//         std::cout << "compute_weights(): Q matrix of size " << this->q_mat->get_nrows();
//         std::cout << "x" << this->q_mat->get_ncols() << " computed:" << std::endl;
//         for (int i = 0; i < this->q_mat->get_nrows(); ++i)
//         {
//            for (int j = 0; j < this->q_mat->get_ncols(); ++j)
//               std::cout << (*(this->q_mat))(j,i) << " ";
//            std::cout << std::endl;
//         }
      }
      
      
      class GaussRadau
      : public IQuadrature
      {
      protected:
         //! @{
         static const bool LEFT_IS_NODE = false;
         static const bool RIGHT_IS_NODE = true;
         //! @}
         
      public:
         //! @{
         /**
          * @throws invalid_argument if less than two nodes are requested
          */
         explicit GaussRadau(const size_t num_nodes);
         GaussRadau() = default;
         virtual ~GaussRadau() = default;
         //! @}
         
         //! @{
         virtual bool left_is_node() const override;
         virtual bool right_is_node() const override;
         //! @}
         
      protected:
         //! @{
         virtual void compute_nodes() override;
         //! @}
      };
      
      GaussRadau::GaussRadau(const size_t num_nodes)
      : IQuadrature(num_nodes)
      {
         if (this->num_nodes < 2) {
            throw invalid_argument("Gauss-Radau quadrature requires at least two quadrature nodes.");
         }
         this->compute_nodes();
         this->compute_weights();
      }
      
      bool GaussRadau::left_is_node() const
      {
         return LEFT_IS_NODE;
      }
      
      bool GaussRadau::right_is_node() const
      {
         return RIGHT_IS_NODE;
      }
      
      void GaussRadau::compute_nodes()
      {
         this->nodes = vector<double>(this->num_nodes, 0.0);
         auto l   = Polynomial::legendre(this->num_nodes);
         auto lm1 = Polynomial::legendre(this->num_nodes - 1);
         
         for (size_t i = 0; i < this->num_nodes; i++) {
            l[i] += lm1[i];
         }
         auto roots = l.roots();
         for (size_t j = 1; j < this->num_nodes; j++) {
            this->nodes[j - 1] = 0.5 * (1.0 - roots[this->num_nodes - j]);
         }
         this->nodes.back() = 1.0;
      }
   }// ::pfasst::quadrature
}// ::pfasst

#endif
