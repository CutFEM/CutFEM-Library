#include <time.h>
////========================================================////
////////////=========     Timer     =========///////////////////

inline double CPUtime() {
#ifdef SYSTIMES
   struct tms buf;
   if (times(&buf) != -1)
      return ((double)buf.tms_utime + (double)buf.tms_stime) /
             (long)sysconf(_SC_CLK_TCK);
   else
#endif
      return ((double)clock()) / CLOCKS_PER_SEC;
}

double getTime() {
#ifdef USE_MPI
   return MPIcf::Wtime();
#else
   return CPUtime();
#endif