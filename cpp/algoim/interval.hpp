#ifndef ALGOIM_INTERVAL_HPP
#define ALGOIM_INTERVAL_HPP

// algoim::Interval<N> methods which compute first-order Taylor series with bounded remainder.
// These methods are based on those described in the paper
//    R. I. Saye, High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes
//    in Hyperrectangles,SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015),
//    http://dx.doi.org/10.1137/140966290

#include "real.hpp"
#include "uvector.hpp"
#include "utility.hpp"

namespace algoim
{
    /* Interval arithmetic using a first-order Taylor series with remainder. A function's range of attainable values is
       evaluated as
         f(x_c + y) = alpha + beta.y + [-eps,eps] 
       where alpha is the value of the function at the centre x_c of an interval [x_c - delta, x_c + delta],
       beta is the n-dimensional gradient of the function evaluated at the centre, y is a placeholder for first-order
       variations, and eps bounds the remainder term.
   
       delta is a static (thread-local) variable and should be appropriately initialised before a sequence of calculations
       involving any Interval<N> objects. Thus, all intervals share the same delta bounds. This is thread-safe, however
       note this technique is not appropriate when using a complicated composition of algorithms (e.g., if a method to
       calculate a function's value internally depends on another method that applies separate interval arithmetic
       calculations).
   
       For a general C^2 function f, a first order Taylor series with remainder can be computed via
         f(alpha + beta.y + h(y)) = f(alpha) + f'(alpha)*(beta.y + h(y)) + 0.5 f''(xi) (beta.y + h(y))^2
       where xi is somewhere in the interval [alpha, alpha + beta.y + h(y)]. Thus, if C bounds |f''| on the same interval,
       then one can deduce
         f(alpha + beta.y + h(y)) = f(alpha) + f'(alpha)*beta.y + hh(y)
       where 
         hh(y) = f'(alpha) h(y) + 0.5 f''(xi) (beta.y + h(y))^2
       and
         |hh(y)| <= |f'(alpha)| eps + 0.5 C (|beta|.delta + eps)^2.
    */
    template<int N>
    struct Interval
    {
        real alpha;
        uvector<real,N> beta;
        real eps;

        // Default constructor corresponds to the zero function
        Interval() : alpha(0.0), beta(0.0), eps(0.0) {}

        // Constructor for a constant
        explicit Interval(real alpha) : alpha(alpha), beta(0.0), eps(0.0) {}

        // Constructor for the linear function f(y) = alpha + beta.y
        explicit Interval(real alpha, const uvector<real,N>& beta) : alpha(alpha), beta(beta), eps(0.0) {}

        // Constructor with all parameters given
        explicit Interval(real alpha, const uvector<real,N>& beta, real eps) : alpha(alpha), beta(beta), eps(eps) {}

        // Assignment from a constant
        Interval& operator=(real alpha)
        {
            this->alpha = alpha;
            beta = 0.0;
            eps = 0.0;
            return *this;
        }

        // Assigment from another Interval, i.e., copy constructor
        Interval& operator=(const Interval& other)
        {
            // Shouldn't need this, but if not done manually (i.e., cannot use "= default;"), Intel compiler
            // icpc 15 has uninitialised memory problem
            alpha = other.alpha;
            beta = other.beta;
            eps = other.eps;
            return *this;
        }

        // Maximum deviation of the interval's values from its value at the centre of the interval, i.e.,
        // place bounds on the linear plus remainder term.
        real maxDeviation() const
        {
            using std::abs;
            real b = eps;
            for (int dim = 0; dim < N; ++dim)
                b += abs(beta(dim)) * delta(dim);
            return b;
        }

        // Unary negation operator
        Interval operator-() const
        {
            return Interval(-alpha, -beta, eps);
        }

        // Addition by a constant
        Interval& operator+=(real c)
        {
            alpha += c;
            return *this;
        }

        // Addition by another Interval
        Interval& operator+=(const Interval& rhs)
        {
            alpha += rhs.alpha;
            beta += rhs.beta;
            eps += rhs.eps;
            return *this;
        }

        // Negation by a constant
        Interval& operator-=(real c)
        {
            alpha -= c;
            return *this;
        }

        // Negation by another Interval
        Interval& operator-=(const Interval& rhs)
        {
            // Note: if rhs == *this, the result could be considered the zero function, but this code may double eps
            alpha -= rhs.alpha;
            beta -= rhs.beta;
            eps += rhs.eps;
            return *this;
        }

        // Multiplication by a constant
        Interval& operator*=(real c)
        {
            using std::abs;
            alpha *= c;
            beta *= c;
            eps *= abs(c);
            return *this;
        }

        // Multiplication by another Interval
        Interval& operator*=(const Interval& rhs)
        {
            using std::abs;
            // This interval's function is:     (*this)(y) = alpha + beta . y + [-eps, eps]
            // The rhs's interval's function is:    rhs(y) = rhs.alpha + rhs.beta . y + [-rhs.eps, rhs.eps]
            // The product of the two functions takes the form
            //   (*this * rhs)(y) = alpha*rhs.alpha + (alpha*rhs.beta + rhs.alpha*beta).y + (beta . y)(rhs.beta . y) + [-eps, eps]*(rhs.alpha + rhs.beta . y) + [-rhs.eps, rhs.eps]*(alpha + beta . y) + [-eps,eps] * [-rhs.eps,rhs.eps]
            real ell_this = maxDeviation();
            real ell_rhs = rhs.maxDeviation();
            // "Reverse ordering" so old values are not prematurely overwritten; this works even if &rhs == this.
            eps = ell_this * ell_rhs + abs(alpha) * rhs.eps + abs(rhs.alpha) * eps;
            for (int dim = 0; dim < N; ++dim)
                beta(dim) = alpha * rhs.beta(dim) + rhs.alpha * beta(dim);
            alpha *= rhs.alpha;
            return *this;
        }

        // Division by a constant, presumed != 0
        Interval& operator/=(real c)
        {
            using std::abs;
            alpha /= c;
            beta /= c;
            eps /= abs(c);
            return *this;
        }

        // Division by another interval
        Interval& operator/=(const Interval& rhs)
        {
            using std::abs;
            real rhsxinv = 1.0 / rhs.alpha;
            real tau = rhs.maxDeviation()*abs(rhsxinv);
            if (tau >= 1.0)
                throw std::domain_error("Unable to perform operator/= with supplied argument");
            real r = util::sqr(tau)/util::cube(1.0 - tau);
            real bnum = maxDeviation();
            eps = (abs(alpha*rhsxinv)*rhs.eps + (abs(alpha) + bnum)*r + bnum*tau + eps) * abs(rhsxinv);
            for (int dim = 0; dim < N; ++dim)
                beta(dim) = beta(dim)*rhsxinv - alpha*rhs.beta(dim)*util::sqr(rhsxinv);
            alpha *= rhsxinv;
            return *this;
        }

        // If the interval's values are uniformly positive or negative, returns +1 or -1, respectively. If
        // the values cannot be guaranteed to have uniform sign, returns 0. Internally, this function uses
        // a small tolerance based on machine precision, so that, e.g., x -> x has uniform sign on the
        // interval [0,1]. The function x -> 0 is declared to have sign 0.
        int sign() const
        {
            real x = (real(1.0) - real(10.0)*real(std::numeric_limits<real>::epsilon()))*maxDeviation();
            if (alpha > x)
                return 1;
            else if (alpha < -x)
                return -1;
            else
                return 0;
        }

        // Returns whether the interval's value are guaranteed to have uniform sign
        bool uniformSign() const
        {
            return sign() != 0;
        }

        // Half-width of interval calculations, across all dimensions
        static uvector<real,N>& delta()
        {
            // thread_local to allow different threads to perform separate chains of Interval arithmetic.
            // Originally, delta was a static member variable of Interval<N> but GCC has compiler bugs 
            // concerning thread_local variables, possibly related to 
            // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81880. Changed to a static member function 
            // instead.
            static thread_local uvector<real,N> d{};
            return d;
        }

        // Half-width of interval calculations, for a specific dimension
        static real& delta(int dim)
        {
            return delta()(dim);
        }
    };

    // Interval<N> binary operators

    template<int N>
    Interval<N> operator+(Interval<N> alpha, const Interval<N>& y)
    {
        return alpha += y;
    }

    template<int N>
    Interval<N> operator+(Interval<N> alpha, real y)
    {
        return alpha += y;
    }

    template<int N>
    Interval<N> operator+(real alpha, Interval<N> y)
    {
        return y += alpha;
    }

    template<int N>
    Interval<N> operator-(Interval<N> alpha, const Interval<N>& y)
    {
        return alpha -= y;
    }

    template<int N>
    Interval<N> operator-(Interval<N> alpha, real y)
    {
        return alpha -= y;
    }

    template<int N>
    Interval<N> operator-(real alpha, Interval<N> y)
    {
        y *= -1.0;
        return y += alpha;
    }

    template<int N>
    Interval<N> operator*(Interval<N> alpha, const Interval<N>& y)
    {
        return alpha *= y;
    }

    template<int N>
    Interval<N> operator*(Interval<N> alpha, real y)
    {
        return alpha *= y;
    }

    template<int N>
    Interval<N> operator*(real alpha, Interval<N> y)
    {
        return y *= alpha;
    }

    template<int N>
    Interval<N> operator/(Interval<N> alpha, const Interval<N>& y)
    {
        return alpha /= y;
    }

    template<int N>
    Interval<N> operator/(Interval<N> alpha, real y)
    {
        return alpha /= y;
    }

    // Streaming operator for inspection
    template<int N>
    std::ostream& operator<< (std::ostream& o, const Interval<N>& alpha)
    {
        return o << "Interval<" << N << ">(" << alpha.alpha << ", " << alpha.beta << ", " << alpha.eps << ')';
    }

    // Various unary function evaluations of Interval objects

    template<int N>
    Interval<N> sin(const Interval<N>& i)
    {
        using std::sin;
        using std::cos;
        using std::abs;
        real c = cos(i.alpha);
        return Interval<N>(sin(i.alpha), c * i.beta, abs(c)*i.eps + 0.5*util::sqr(i.maxDeviation()));
    }

    template<int N>
    Interval<N> cos(const Interval<N>& i)
    {
        using std::sin;
        using std::cos;
        using std::abs;
        real s = -sin(i.alpha);
        return Interval<N>(cos(i.alpha), s * i.beta, abs(s)*i.eps + 0.5*util::sqr(i.maxDeviation()));
    }

    template<int N>
    Interval<N> sqrt(const Interval<N>& i)
    {
        using std::sqrt;
        real b = i.maxDeviation();
        if (b > i.alpha)
            throw std::domain_error("Unable to compute sqrt() with supplied argument");
        real sqrtx = sqrt(i.alpha);
        real sqrtxinv = 1.0 / sqrtx;
        real C = 0.25/((i.alpha - b)*sqrt(i.alpha - b));
        return Interval<N>(sqrtx, (0.5*sqrtxinv)*i.beta, 0.5 * C * util::sqr(b));
    }

    template<int N>
    Interval<N> exp(const Interval<N>& i)
    {
        using std::exp;
        real ex = exp(i.alpha);
        real b = i.maxDeviation();
        // Bounding (exp(alpha))'' using montonocity and thus argument [i.alpha + b]
        return Interval<N>(ex, ex * i.beta, ex*i.eps + 0.5 * exp(i.alpha + b) * util::sqr(b));
    }

    template<int N>
    Interval<N> cosh(const Interval<N>& z)
    {
        using std::sinh;
        using std::cosh;
        using std::abs;
        real coshz = cosh(z.alpha);
        real sinhz = sinh(z.alpha);
        real b = z.maxDeviation();
        // Second derivative of cosh(alpha) is cosh(alpha), and is even and monotone increasing for alpha >= 0. Therefore
        // for z.alpha - z.maxDeviation <= alpha <= z.alpha + z.maxDeviation, the maximum value for cosh(alpha) is cosh(abs(z.alpha) + z.maxDeviation)
        return Interval<N>(coshz, sinhz * z.beta, abs(sinhz)*z.eps + 0.5 * cosh(abs(z.alpha) + b) * util::sqr(b));
    }

    template<int N>
    Interval<N> tanh(const Interval<N>& z)
    {
        using std::tanh;
        using std::cosh;
        using std::abs;
        real tanhz = tanh(z.alpha);
        real sechzsqr = util::sqr(1.0/cosh(z.alpha));
        // The second derivative of tanh(alpha) is -2*tanh(alpha)*sech(alpha)^2, and has maximum absolute value 0.76980035892...
        return Interval<N>(tanhz, sechzsqr * z.beta, abs(sechzsqr)*z.eps + 0.5 * 0.77 * util::sqr(z.maxDeviation()));
    }

    inline real sech(real alpha)
    {
        using std::cosh;
        return real(1)/cosh(alpha);
    }

    template<int N>
    Interval<N> sech(const Interval<N>& z)
    {
        using std::tanh;
        using std::abs;
        real sechz = sech(z.alpha);
        real fprime = -sechz*tanh(z.alpha);
        uvector<real,N> beta = fprime * z.beta;
        // The second derivative of sech(alpha) is {sech(alpha)*tanh(alpha)^2 - sech(alpha)^3} and has maximum absolute value 1.0
        return Interval<N>(sechz, beta, abs(fprime)*z.eps + 0.5 * 1.0 * util::sqr(z.maxDeviation()));
    }
} // namespace algoim

#endif
