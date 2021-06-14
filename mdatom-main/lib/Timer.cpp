#include "Timer.h"
#include <iomanip>

using namespace std;

void Timer::programStart() {
    tProgramStart = chrono::system_clock::now(); // instant execution begins
}

void Timer::mdStart() {
    tMDStart = chrono::system_clock::now();
}

void Timer::mdEnd() {
    tMDEnd = chrono::system_clock::now();
}

void Timer::programEnd() {
    tProgramEnd = chrono::system_clock::now();
}

void Timer::output(std::ostream &out) const {
    /*
     * Timing.
     * The c++11 time_point objects are converted to time_t variables to allow for
     * formatted time output.
     */
    out << " Timing:" << endl << endl;
    auto time_0 = chrono::system_clock::to_time_t(tProgramStart);
    out << "  At start of program         : " << ctime(&time_0);
    auto time_start = chrono::system_clock::to_time_t(tMDStart);
    out << "  Before main iteration loop  : " << ctime(&time_start);
    auto time_end = chrono::system_clock::to_time_t(tMDEnd);
    out << "  After main iteration loop   : " << ctime(&time_end);
    auto now = chrono::system_clock::to_time_t(tProgramEnd);
    out << "  Now                         : " << ctime(&now);

    auto runmda_duration_ms = double(chrono::duration_cast<chrono::milliseconds>(tMDEnd - tMDStart).count());
    auto total_duration_ms = double(chrono::duration_cast<chrono::milliseconds>(tProgramEnd - tProgramStart).count());
    out << "  Time spent for MD iterations: " << setprecision(3) << runmda_duration_ms / 1000.0 << " sec" << endl;
    out << "  Time spent in total         : " << setprecision(3) << total_duration_ms / 1000.0 << " sec" << endl;
}

