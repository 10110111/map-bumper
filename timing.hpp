#pragma once
#include <cmath>
#include <chrono>
#include <atomic>
#include <sstream>

template<typename T>
std::string formatDeltaTime(const std::chrono::time_point<T> timeBegin, const std::chrono::time_point<T> timeEnd)
{
    const auto microsecTaken=std::chrono::duration_cast<std::chrono::microseconds>(timeEnd-timeBegin).count();
    const auto secondsTaken=1e-6*microsecTaken;
    std::ostringstream ss;
    if(secondsTaken < 1e-3)
    {
        ss << microsecTaken << u8" \u03bcs";
    }
    else if(secondsTaken < 1)
    {
        ss << secondsTaken*1000 << " ms";
    }
    else if(secondsTaken < 60)
    {
        ss << secondsTaken << " s";
    }
    else
    {
        auto remainder=secondsTaken;
        const auto d = int(remainder/(24*3600));
        remainder -= d*(24*3600);
        const auto h = int(remainder/3600);
        remainder -= h*3600;
        const auto m = int(remainder/60);
        remainder -= m*60;
        const auto s = std::lround(remainder);
        if(d)
            ss << d << "d";
        if(d || h)
            ss << h << "h";
        if(d || h || m)
            ss << m << "m";
        ss << s << "s";
    }
    return ss.str();
}

void handleProgressReporting(const size_t totalItemCount,
                             const std::chrono::time_point<std::chrono::steady_clock> startTime,
                             std::chrono::time_point<std::chrono::steady_clock>& time0,
                             std::atomic_int& numThreadsReportedFirstProgress,
                             size_t& itemsDoneInThisThreadAfterLastUpdate,
                             std::atomic<unsigned>& itemsDone);
