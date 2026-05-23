// Memory usage reporting
// Cross-platform: Linux (/proc/self/stat) and macOS (mach API)

#include "getMem.hpp"

#ifdef __APPLE__
#include <mach/mach.h>
#endif

void process_mem_usage(double& vm_usage, double& resident_set)
{
   vm_usage     = 0.0;
   resident_set = 0.0;

#ifdef __APPLE__
   // macOS: use mach kernel API
   struct mach_task_basic_info info;
   mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
   if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                 (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
       vm_usage = (double)info.virtual_size;
       resident_set = (double)info.resident_size;
   }
#else
   // Linux: read /proc/self/stat
   using std::ios_base;
   using std::ifstream;
   using std::string;

   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss;

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE);
   vm_usage     = vsize;
   resident_set = rss * page_size_kb;
#endif
}
