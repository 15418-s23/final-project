// From the 15-418 Assignment

#ifndef TIMING_H
#define TIMING_H

#include <chrono>

class Timer {
public:
  Timer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(clock_::now() - beg_).count();
  }

private:
  typedef std::chrono::high_resolution_clock clock_;
  std::chrono::time_point<clock_> beg_;
};

#endif