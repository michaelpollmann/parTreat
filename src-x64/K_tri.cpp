#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>

inline double intPow(double x,int n)
{
    double r = 1;
    for (int i = 0; i<n; i++)
    {
        r *= x; 
    }
    return r; 
}

inline double K_d0(double u)
{
    return ((35.0/32)*(intPow(1-intPow(u,2),3)));
}

inline double K_d1(double u)
{
    return ((3*35.0/32)*(intPow(1-intPow(u,2),2)*(-2*u)));
}

inline double K_d2(double u)
{
    return ((2*3*35.0/32)*((1-intPow(u,2))*(5*intPow(u,2)-1)));
}



// [[Rcpp::export]]
double K_tri(double x0, NumericVector x, NumericVector h, NumericVector h_pow_d, int d) {
    // link kernel function depending on derivative
    double (*K_func)(double);
    if (d == 0) {
        K_func = &K_d0;
    } else if (d==1) {
        K_func = &K_d1;
    } else if (d==2) {
        K_func = &K_d2;
    } else {
        return NA_REAL;
    }

    // output
    double K = 0.0;

    // if the same bandwidth for all observations
    if (h.length() == 1) {
        // should we stop at the next zero-kernel?
        bool stopNext0;
        stopNext0 = false;
        // normalized point to put into kernel function
        double u;
        for (int i=0; i<x.length(); i++) {
            u = (x0-x(i))/h(0);
            if (std::abs(u) <= 1) {
                K += K_func(u);
                stopNext0 = true;
            } else if ((std::abs(u) > 1) && (stopNext0)) {
                // we're past the non-zero kernels
                break;
            }
        }
        // divide by sample size 
        K = K/(x.length()*h_pow_d(0));
    }
    else
    {
        // normalized point to put into kernel function
        double u;
        for (int i=0; i<x.length(); i++)
        {
            u = (x0-x(i))/h(i);
            if (std::abs(u) <= 1)
            {
                K += K_func(u)/h_pow_d(i);
            }
        }
        // divide by sample size
        K = K/(x.length());
    }
    return K;
}
