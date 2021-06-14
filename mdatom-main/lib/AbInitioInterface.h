//
// Created by jonas on 5/28/21.
//

#ifndef MDATOM_ABINITIOINTERFACE_H
#define MDATOM_ABINITIOINTERFACE_H


#include <utility>
#include <vector>
#include <string>
#include <unordered_map>

struct CalculationResult {
    std::vector<double> gradient;
    double energy;
    std::unordered_map<size_t, std::pair<std::string, double>> atomTypes;
};

class AbInitioInterface {
public:
    AbInitioInterface();

    AbInitioInterface(const std::string &grad_path, const std::string &atm_list_path,
                      const std::string &atm_info_path);

    CalculationResult simulate(std::vector<double> positions, std::vector<size_t> atoms);


    void cleanFiles();

    void finalCleanup();

private:
    std::string gradient_path;
    std::string atom_list_path;
    std::string atom_info_path;

    std::pair<double, std::vector<double>> readGradientEnergy();

    std::unordered_map<size_t, std::pair<std::string, double>> readAtomInformation();

    void writeAtomList(std::vector<double> positions, std::vector<size_t> atoms);

};


#endif //MDATOM_ABINITIOINTERFACE_H
