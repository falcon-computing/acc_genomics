#ifndef ACCLIB_TIMER_H
#define ACCLIB_TIMER_H

#include <map>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#ifndef TIMER_VERBOSE
#define TIMER_VERBOSE 1
#endif

#ifndef TIMER_REPORT
#define TIMER_REPORT true
#endif

extern std::map<std::string, uint64_t> g_total_time;
extern std::map<std::string, uint64_t> g_last_time;

inline uint64_t getUs() {
  struct timespec tr;
  clock_gettime(CLOCK_REALTIME, &tr);
  return (uint64_t)tr.tv_sec*1e6 + tr.tv_nsec/1e3;
}

static void print_global_timers() {
  for (auto p : g_total_time) {
    fprintf(stderr, "[Timer::Total]: %s takes %d us\n",
        p.first.c_str(),
        p.second);
  }
}

static uint64_t get_global_etime(std::string s) {
  if (g_total_time.count(s)) return g_total_time[s];
  else return 0;
}

// unit: us
static void add_global_etime(std::string s, uint64_t t) {
  if (g_total_time.count(s)) g_total_time[s] += t;
  else g_total_time[s] = t;
}

// unit: us
static uint64_t get_last_etime(std::string s) {
  if (g_last_time.count(s)) return g_last_time[s];
  else 0;
}

static void add_last_etime(std::string s, uint64_t t) {
  if (g_last_time.count(s)) g_last_time[s] += t;
  else g_last_time[s] = t;
}

class Timer {
 public: 
  Timer(std::string func = "-", 
      bool flag_report = TIMER_REPORT,
      int verbose = TIMER_VERBOSE): 
    func_(func), verbose_(verbose), flag_report_(flag_report)
  {
    if (func.empty()) {
      throw std::runtime_error("Timer::Timer(): timer name cannot be empty");
    }
    start_ts_ = getUs(); 
  }

  ~Timer() {
    uint64_t e_time = getUs()-start_ts_;
    if (verbose_ > 0) {
      fprintf(stderr, "[Timer]: %s takes %ld us\n", func_.c_str(), e_time);
    }
    if (flag_report_) {
      add_global_etime(func_, e_time);
      add_last_etime(func_, e_time);
    }
  }
  
 private:
  int verbose_;
  bool flag_report_;
  std::string func_;
  uint64_t start_ts_;
};

// only suppose to call once
#define DEFINE_GLOBAL_TIMER std::map<std::string, uint64_t> g_total_time; \
  std::map<std::string, uint64_t> g_last_time; 

#define PLACE_TIMER Timer __timer_obj(__func__);
#define CONCAT_FNAME(A, B) (std::string(A) + "::" + std::string(B))
#define PLACE_TIMER1(s) Timer __timer_obj(CONCAT_FNAME(__func__, s));

#endif
