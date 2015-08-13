#ifndef mfem_arnoldi_HEADER
#define mfem_arnoldi_HEADER

#include "mfem.hpp"

using mfem;

class DenseMatrixTimeDependentOperator : public TimeDependentOperator
{
    public:
        DenseMatrixTimeDependentOperator(DenseMatrix & A_) : A(A_) {}
        
        virtual void Mult(const Vector &x, Vector &y) const
        {
            A.Mult(x,y);
        }

    private:
        DenseMatrix & A;
};

class Arnoldi
{
    public:
        Arnoldi(int k_max_, MPI_Comm comm_)
        : H(k_max_,k_max_), Hop(H), comm(comm_), k_max(k_max_) 
        { 
            V = new Vector[k_max_];
            A = NULL;
        }

        ~Arnoldi()
        {
            delete[] V;
        }

        void SetOperator(Operator &A_)
        {
            A = &A_;
            int n = A.Width();

            for(int i = 0; i < k_max; i++)
            {
                V[i].SetSize(n);
            }
        }

        void ApplyV(const Vector & ubar, Vector & u)
        {
            u = 0.0;
            int k = H.Width();
            
            for(int i = 0; i < k; i++)
            {
                u.Add(ubar(i), V[i]);
            }

        }

        void ApplyVT(const Vector & u, Vector & ubar)
        {
            int k = H.Width();
            
            ubar.SetSize(k);
            for(int i = 0; i < k; i++)
            {
                ubar(i) = Dot(V[i], u);
            }
        }

        TimeDependentOperator & GetH( )
        {
            return Hop;
        }

        GenKrylovSpace(const Vector & u)
        {
            int k;
            double max_norm = 0.0;
            double tol = 1e-12;
            double norm_u = sqrt(Dot(u,u)); 
            DenseMatrix F(k_max, k_max);
            
            if(norm_u > 0.0)
            {
                V[0].Set( 1.0/norm_u, u );
            }
            else
            {
                MFEM_ABORT("Does Not Support All Zero Seed Vector\n");
            }

            for(k = 0; k < k_max-1; k++)
            {
                A.Mult(V[k], V[k+1]);

                for(int m = 0; m <=k; m++)
                {
                    F(m,k) = Dot(V[m], V[k+1]);
                    V[k+1].Add( -F(m,k), V[m] );
                }

                F(k+1,k) = sqrt( Dot(V[k+1], V[k+1]) );

                max_norm = max(F(k+1,k), max_norm);
                if (F(k+1,k) <= tol*max_norm)
                {
                    break;
                }

                V[k+1] /= F(k+1,k);
            }
            
            H.SetSize(k+1, k+1);
            H.CopyMN(F, k+1, k+1, 0, 0);
        }

    protected:
        Vector *V;
        DenseMatrix H;
        DenseMatrixTimeDependentOperator Hop;
        Operator *A;
        MPI_Comm comm;
        int k_max;

        double Dot(const Vector &x, const Vector &y) const
        {
            double local_dot = (x * y);
            double global_dot;

            MPI_Allreduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, comm);

            return global_dot;
        }


};

#endif // mfem_arnoldi_HEADER

