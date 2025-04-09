#include "timing.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>

void handleProgressReporting(const size_t totalItemCount,
                             const std::chrono::time_point<std::chrono::steady_clock> startTime,
                             std::chrono::time_point<std::chrono::steady_clock>& time0,
                             std::atomic_int& numThreadsReportedFirstProgress,
                             size_t& itemsDoneInThisThreadAfterLastUpdate,
                             std::atomic<unsigned>& itemsDone)
{
    ++itemsDoneInThisThreadAfterLastUpdate;
    using namespace std::chrono;
    const auto time1 = steady_clock::now();
    if(time1 - time0 > seconds(5))
    {
        // ETA is reported only by the first thread in the iteration.
        // This seems to provide the most accurate estimate. Once the
        // reporting thread is determined, it will report after each
        // iteration.
        thread_local bool isTimeReportingThread = false;
        thread_local bool reportedProgressFirstTime = false;
        if(!reportedProgressFirstTime)
        {
            if(numThreadsReportedFirstProgress.fetch_add(1) == 0)
                isTimeReportingThread = true;
        }

        itemsDone += itemsDoneInThisThreadAfterLastUpdate;
        itemsDoneInThisThreadAfterLastUpdate = 0;
        const auto progress = double(itemsDone) / totalItemCount;
        const auto usecSincePrevReport = duration_cast<microseconds>(time1 - time0).count();
        const auto usecElapsed = duration_cast<microseconds>(time1 - startTime).count();
        const long secToEnd = std::lround((1. - progress) / progress * usecElapsed * 1e-6);
        std::ostringstream ss;
        const double progressPerReport = progress / usecElapsed * usecSincePrevReport;
        const int precDigits = std::clamp(int(std::ceil(-std::log10(progressPerReport))) - 1, 0, 6);
        ss << std::fixed << std::setprecision(precDigits) << progress*100 << "% done";
        // First ETA estimate is too far from reality, likely due
        // to overhead of thread startup, so skip first measurement.
        if(isTimeReportingThread && reportedProgressFirstTime)
        {
            const auto days = secToEnd/(3600*24);
            const auto hr = secToEnd/3600%24;
            const auto min = secToEnd/60%60;
            const auto sec = secToEnd%60;
            ss << ", ETA: ";
            if(days) ss << days << 'd';
            if(days || hr) ss << hr << 'h';
            if(days || hr || min) ss << min << 'm';
            ss << sec << "s\n";
        }
        else
        {
            ss << '\n';
        }
        reportedProgressFirstTime = true;
        std::cerr << ss.str();
        time0 = time1;
    }
}
