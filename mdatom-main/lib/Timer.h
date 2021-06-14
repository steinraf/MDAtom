#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>

/*!
 * Class to time the MD program.
 */
class Timer {
public:
    void programStart();

    void mdStart();

    void mdEnd();

    void programEnd();

    void output(std::ostream &out) const;

private:
    std::chrono::system_clock::time_point tProgramStart,
            tProgramEnd,
            tMDStart,
            tMDEnd;
};

#endif // TIMER_H
