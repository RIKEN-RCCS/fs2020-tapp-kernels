#include "ktrace_mt.h"

#ifndef UNDEF_KTRACE

void ktrace_enter_mt() {
#pragma omp parallel
  {
    ktrace_enter();
  }
}

void ktrace_exit_mt() {
#pragma omp parallel
  {
    ktrace_exit();
  }
}

void ktrace_start_mt(int i) {
#pragma omp parallel
  {
    ktrace_start(i);
  }
}

void ktrace_stop_mt(int i) {
#pragma omp parallel
  {
    ktrace_stop(i);
  }
}

void ktrace_enter_mt_() {
  ktrace_enter_mt();
}

void ktrace_exit_mt_() {
  ktrace_exit_mt();
}

void ktrace_start_mt_(int *i) {
  ktrace_start_mt(*i);
}

void ktrace_stop_mt_(int *i) {
  ktrace_stop_mt(*i);
}

#endif
