#include "stopwatch.h"

StopWatch::StopWatch()
{
}

void StopWatch::start() {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_begin);
    isStarted = true;
}


void StopWatch::stop() {
    if (!isStarted) {
        return;
    }

    // log current running time
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);
    double time = (ts_end.tv_sec - ts_begin.tv_sec) +
                  (ts_end.tv_nsec - ts_begin.tv_nsec) / 1e9;
    timeRunning += time;
}

void StopWatch::reset() {
    isStarted = false;
    timeRunning = 0.0;
}

double StopWatch::getTime() {
    return timeRunning;
}
