/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <vector>
#include <chrono>
#include <sstream>
#include <iostream>
#include <stack>
#include <valarray>
#include <map>
#include <cstring>
#include <limits>
#include "../common/RNM.hpp"
#include "../common/global.hpp"
#include <sys/stat.h>
#include <list>
#include <time.h>
#include <ctime>

typedef double R;

typedef unsigned int Uint;
typedef unsigned long int Ulint;
typedef unsigned short Usint;
typedef signed char byte;
typedef unsigned char Ubyte;

typedef KN<R> Rn;
typedef KNM<R> Rnm;
typedef double R;
typedef KN_<R> RN_;
typedef KN_<R> Rn_;
typedef KN<R> RN;
typedef KNM_<R> RNM_;
typedef KNMK_<R> RNMK_;

inline std::tuple<size_t, size_t> size(const std::map<std::pair<int, int>, double> &mat) {
    return (mat.size() == 0 ? std::make_pair(0, 0)
                            : std::make_pair(mat.rbegin()->first.first + 1, mat.rbegin()->first.second + 1));
}

namespace util {
inline byte sign(double d) { return d > 0. ? 1 : (d < 0. ? -1 : 0); }

inline byte fsign(double d) { return d > 0. ? 1 : (d < 0. ? -1 : 0); }

// inline bool changeSign(const Rn &v) {
//    R a = v(0);
//    for (int i = 1; i < v.size(); ++i) {
//       if (sign(a * v(i)) < 0)
//          return true;
//    }
//    return false;
// }

template <class Container> auto &at(Container &&c, std::size_t pos) {
    if (pos >= c.size())
        throw std::out_of_range("out of bounds");
    return c[pos];
}

inline bool changeSign(const R *v, int size) {
    R a = v[0];
    for (int i = 1; i < size; ++i) {
        if (sign(a * v[i]) < 0)
            return true;
    }
    return false;
}

///\brief Write the sign of the levelset function src to the sequence dst.
inline void copy_levelset_sign(const KN<double> &src, KN<byte> &dst) {
    dst.resize(src.size());
    std::transform(&src[0], &src[0] + src.size(), &dst[0], sign);
}

static bool contain(const KN<int> &v, int x) {
    for (int i = 0; i < v.size(); ++i)
        if (v(i) == x)
            return true;
    return false;
}
static bool contain(const std::list<int> &v, int x) {
    for (auto it = v.begin(); it != v.end(); ++it)
        if (*it == x)
            return true;
    return false;
}

} // namespace util

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
inline const std::string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d-%X", &tstruct);

    return buf;
}
inline bool IsPathExist(const std::string &s) {
    struct stat buffer;
    return (stat(s.c_str(), &buffer) == 0);
}

/// Get the address of the first element in a valarray
template <typename T> inline const T *Addr(const std::valarray<T> &x) {
    return &(const_cast<std::valarray<T> &>(x)[0]);
}

/// Get the address of the first element in a valarray
template <typename T> inline T *Addr(std::valarray<T> &x) { return &(x[0]); }

template <typename T> inline const T *Addr(const std::vector<T> &x) { return &(const_cast<std::vector<T> &>(x)[0]); }

/// Get the address of the first element in a vector
template <typename T> inline T *Addr(std::vector<T> &x) { return &(x[0]); }

static bool equal(double a, double b, double eps = globalVariable::Epsilon) { return (std::fabs(a - b) < eps); }
inline bool isEqual(double a, double b, double eps = globalVariable::Epsilon) { return (std::fabs(a - b) < eps); }
inline bool isNotEqual(double a, double b, double eps = globalVariable::Epsilon) { return (std::fabs(a - b) > eps); }
static bool nequal(double a, double b, double eps = globalVariable::Epsilon) { return (std::fabs(a - b) > eps); }

// stream buffer without output on screen
class NullStreamBuf : public std::streambuf {
  public:
    NullStreamBuf() {}
};

class MuteStdOstream {
  private:
    std::streambuf *bout_, *berr_, *blog_;
    NullStreamBuf devnull_;

  public:
    MuteStdOstream() : bout_(std::cout.rdbuf()), berr_(std::cerr.rdbuf()), blog_(std::clog.rdbuf()) {}

    void Mute(std::ostream &os) { os.rdbuf(&devnull_); }
    void Mute() {
        Mute(std::cout);
        Mute(std::clog);
    }
    void Recover() const {
        std::cout.rdbuf(bout_);
        std::cerr.rdbuf(berr_);
        std::clog.rdbuf(blog_);
    }
};

////========================================================////
////////////===== Barre de progression ======///////////////////

class progress {

    using time_t = std::chrono::time_point<std::chrono::high_resolution_clock>;

  private:
    const char *title;
    const int length;
    int it;
    int prg;
    time_t t_start;
    // clock_t t0;
    int verbose;

  public:
    progress(const char *aff, const int &l, int v) : title(aff), length(l) {
        t_start = std::chrono::high_resolution_clock::now();
        // t0 = clock();
        it      = 0;
        prg     = 0;
        verbose = v;

        if (verbose > 1) {
#pragma omp single nowait
            {
                std::cout << "\r";
                std::cout << title << ": \t";
                std::cout << 0 << "%";
                std::cout.flush();
            }
        }
    }
    // progress(const char *aff, const int &l, int v, int np) : progress(aff, l / np, v) {}

    void operator++(int n) {
#pragma omp atomic
        it++;

#pragma omp single nowait
        {
            if (int(it * 100. / length) > prg & verbose > 1) {
                prg = int(it * 100. / length);
                std::cout << "\r";
                std::cout << title << ": \t";
                std::cout << prg << "%";
                std::cout.flush();
            }
        }
    }

    void operator+=(int n) {
#pragma omp atomic
        it += n;
#pragma omp single nowait
        {
            if (int(it * 100. / length) > prg & verbose > 1) {
                prg = int(it * 100. / length);
                std::cout << "\r";
                std::cout << title << ": \t";
                std::cout << prg << "%";
                std::cout.flush();
            }
        }
    }

    void end() {
#pragma omp single nowait
        {
            auto t_end  = std::chrono::high_resolution_clock::now();
            double time = std::chrono::duration<double>(t_end - t_start).count();
            // t0         = clock() - t0;
            // time_t now;
            // time(&now);

            if (verbose > 1) {
                std::cout << "\r";
                std::cout << title << ": \t";
                std::cout << time << " sec." << std::endl;
            }
        }
    }
};

static void gather(std::vector<std::map<std::pair<int, int>, double>> &l) {
    for (int i = 1; i < l.size(); ++i) {
        auto &A(l[i]);
        for (auto &aij : A) {
            l[0][aij.first] += aij.second;
        }
        A.clear();
    }
}

inline decltype(auto) now() { return std::chrono::high_resolution_clock::now(); }

inline decltype(auto) seconds(const decltype(now()) &start, const decltype(now()) &end) {
    return std::chrono::duration<double>(end - start).count();
}

inline decltype(auto) milliseconds(const decltype(now()) &start, const decltype(now()) &end) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

////========================================================////
////////////=========     Timer     =========///////////////////

inline double CPUtime() {
#ifdef SYSTIMES
    struct tms buf;
    if (times(&buf) != -1)
        return ((double)buf.tms_utime + (double)buf.tms_stime) / (long)sysconf(_SC_CLK_TCK);
    else
#endif
        return ((double)clock()) / CLOCKS_PER_SEC;
}

inline double getTime() {
#ifdef USE_MPI
    // return MPIcf::Wtime();
    return CPUtime();
#else
    return CPUtime();
#endif
}

#endif
