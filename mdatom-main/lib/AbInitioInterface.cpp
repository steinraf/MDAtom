//
// Created by jonas on 5/28/21.
//
#include "AbInitioInterface.h"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <thread>
#include <chrono>
#include <unordered_map>

using namespace std::chrono_literals;

AbInitioInterface::AbInitioInterface()
        : gradient_path("../tmp/gradient.txt"), atom_list_path("../tmp/positions.txt"),
          atom_info_path("../tmp/atom_info.txt") {}


AbInitioInterface::AbInitioInterface(const std::string &grad_path, const std::string &atm_list_path,
                                     const std::string &atm_info_path)
        : gradient_path(grad_path), atom_list_path(atm_list_path), atom_info_path(atm_info_path) {}

std::pair<double, std::vector<double>> AbInitioInterface::readGradientEnergy() {

    auto file = std::ifstream(gradient_path);

    std::vector<double> gradient;

    std::string line;
    std::getline(file, line);

    double energy;
    auto lineStream = std::istringstream(line);
    lineStream >> energy;


    for (; std::getline(file, line);) {
        lineStream = std::istringstream(line);

        for (int j = 0; j < 3; j++) {
            std::string num;
            std::getline(lineStream, num, ',');

            gradient.push_back(std::stod(num));
        }
    }

    std::filesystem::remove(gradient_path);

    return std::pair<double, std::vector<double>>(energy, gradient);
}

void AbInitioInterface::writeAtomList(std::vector<double> positions, std::vector<size_t> atoms) {

    const auto N = atoms.size();

    auto of = std::ofstream(atom_list_path);

    if (!of.is_open()) {
        std::cerr
                << "Could not write to atom list path.\nBe sure to create a folder called tmp at the root of the project and run the program from the mdatom-main folder"
                << std::endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            of << positions[3 * i + j] << ",";
        }

        of << atoms[i] << std::endl;
    }

}

CalculationResult AbInitioInterface::simulate(std::vector<double> positions, std::vector<size_t> atoms) {

    writeAtomList(positions, atoms);

    while (!std::filesystem::is_regular_file(gradient_path) || !std::filesystem::is_regular_file(atom_info_path)) {
        std::this_thread::sleep_for(1000ms);
    }

    const auto[energy, gradient] = readGradientEnergy();
    const auto atom_information = readAtomInformation();

    CalculationResult result{gradient, energy, atom_information};

    return result;
}

std::unordered_map<size_t, std::pair<std::string, double>> AbInitioInterface::readAtomInformation() {
    auto file = std::ifstream(atom_info_path);

    std::unordered_map<size_t, std::pair<std::string, double>> atom_information;

    for (std::string line; std::getline(file, line);) {
        std::istringstream line_stream(line);

        unsigned int atomic_number;
        double atomic_mass;

        std::string symbol;

        line_stream >> atomic_number >> symbol >> atomic_mass;

        atom_information[atomic_number] = {symbol, atomic_mass};
    }

    std::filesystem::remove(atom_info_path);

    return atom_information;
}

void AbInitioInterface::cleanFiles() {
    std::filesystem::remove(gradient_path);
    std::filesystem::remove(atom_info_path);
    std::filesystem::remove(atom_list_path);
}

void AbInitioInterface::finalCleanup() {

    std::string ab_initio_path("../ab_initio/");

    //Remove all .clean files that are created by psi4
    for (auto &file: std::filesystem::directory_iterator(ab_initio_path)) {

        const std::string file_name = file.path();

        if (file_name.find(".clean") != std::string::npos)
            std::filesystem::remove(file_name);
    }

    //Remove all temporary files created to allow communication between c++ and python
    cleanFiles();


}


