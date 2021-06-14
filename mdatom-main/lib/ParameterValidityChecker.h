#ifndef PARAMETERVALIDITYCHECKER_H
#define PARAMETERVALIDITYCHECKER_H

#include "MDParameters.h"
#include <string>

/*!
 * This class verifies the validity of the struct containing the MD parameters.
 */
class ParameterValidityChecker {
public:
    explicit ParameterValidityChecker(const MDParameters &parameters);

    bool valid() const;

    std::string getErrorMessage() const;

private:
    bool checkValidity();

    const MDParameters &par;
    bool parametersValid;
    std::string errorMessage;
};

#endif // PARAMETERVALIDITYCHECKER_H