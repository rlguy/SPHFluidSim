#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <time.h>


// used to keep track of simulation times
class StopWatch
{
public:
    StopWatch();
    void start();
    void stop();
    void reset();
    double getTime();    // in seconds

private:
    bool isStarted = false;
    timespec ts_begin, ts_end;
    double timeRunning = 0.0;
};

#endif // STOPWATCH_H
