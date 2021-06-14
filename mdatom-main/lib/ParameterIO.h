#ifndef PARAMETERIO_H
#define PARAMETERIO_H

#include "MDParameters.h"
#include <string>
#include <iostream>

/*
 * Class to read the parameters from a file
 */
class ParameterIO {
public:
    static MDParameters readParameters(const std::string &fileName);

    static void outputParameters(std::ostream &out, const MDParameters &par);
};


#endif // PARAMETERIO_H