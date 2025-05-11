#pragma once
#include "rotor.hpp"
#include "particle.hpp"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <ctime>
#include <iostream>
#include <vector>
#include <sstream>

namespace fs = std::filesystem;

namespace BET {
    struct TimeStepData {
        double time;
        std::vector<std::tuple<double, double, double, double, double, double, double, double, double>> bladeStates;
        std::vector<std::tuple<double, double, double, double, double, double>> particleStates;
    };

    class VTKWriter {
    private:
        std::string folderPath;
        int timeStepCounter;

    public:
        VTKWriter() : folderPath("vtk_files_BE_1"), timeStepCounter(0) {
            if (!fs::exists(folderPath)) {
                fs::create_directories(folderPath);
            }
        }

        void saveTimeStep(const Rotor& rotor, double time, const ParticleSystem* particleSystem) {
            // Create filename with padded time step number
            std::stringstream ss;
            ss << folderPath << "/simulation_" << std::setfill('0') << std::setw(6) << timeStepCounter << ".vtk";
            std::string filePath = ss.str();
            
            std::ofstream file(filePath);
            if (!file.is_open()) {
                std::cerr << "Error: Could not open file for writing VTK output." << std::endl;
                return;
            }

            // Write VTK header
            file << "# vtk DataFile Version 3.0\n";
            file << "Rotor Simulation Data\n";
            file << "ASCII\n";
            file << "DATASET POLYDATA\n";

            // Get current time step data
            auto bladeStates = rotor.getBladeElementStates(time);
            std::vector<std::tuple<double, double, double, double, double, double>> particleStates;
            if (particleSystem) {
                particleStates = particleSystem->getActiveParticles();
            }

            // Calculate total points
            size_t totalPoints = bladeStates.size() + particleStates.size();

            // Write points
            file << "POINTS " << totalPoints << " float\n";
            
            // Write blade element positions
            for (const auto& [x, y, z, UR, UT, UP, Ublade_x, Ublade_y, Ublade_z] : bladeStates) {
                file << std::fixed << std::setprecision(6) << x << " " << y << " " << z << "\n";
            }
            // Write particle positions
            for (const auto& [x, y, z, vx, vy, vz] : particleStates) {
                file << std::fixed << std::setprecision(6) << x << " " << y << " " << z << "\n";
            }

            // Write blade connectivity
            size_t elementsPerBlade = bladeStates.size() / 2; // Assuming 2 blades
            file << "LINES " << 2 << " " << (2 * (elementsPerBlade + 1)) << "\n";
            
            // Write connectivity for each blade
            for (int blade = 0; blade < 2; ++blade) {
                file << elementsPerBlade << " ";
                for (int elem = 0; elem < elementsPerBlade; ++elem) {
                    file << (blade * elementsPerBlade + elem) << " ";
                }
                file << "\n";
            }

            // Write point data
            file << "POINT_DATA " << totalPoints << "\n";
            
            // Write time data
            file << "SCALARS time float 1\n";
            file << "LOOKUP_TABLE default\n";
            // Write time for blade elements
            for (size_t i = 0; i < bladeStates.size(); ++i) {
                file << std::fixed << std::setprecision(6) << time << "\n";
            }
            // Write time for particles
            for (size_t i = 0; i < particleStates.size(); ++i) {
                file << std::fixed << std::setprecision(6) << time << "\n";
            }

            // Write velocities
            file << "VECTORS velocities float\n";
            // Write blade element velocities
            for (const auto& [x, y, z, UR, UT, UP, Ublade_x, Ublade_y, Ublade_z] : bladeStates) {
                file << std::fixed << std::setprecision(6) << Ublade_x << " " << Ublade_y << " " << Ublade_z << "\n";
            }
            // Write particle velocities
            for (const auto& [x, y, z, vx, vy, vz] : particleStates) {
                file << std::fixed << std::setprecision(6) << vx << " " << vy << " " << vz << "\n";
            }

            file.close();
            std::cout << "VTK file written to: " << filePath << std::endl;
            
            // Create or update the time series file
            std::string timeSeriesPath = folderPath + "/simulation.pvd";
            std::ofstream timeSeriesFile(timeSeriesPath, timeStepCounter == 0 ? std::ios::out : std::ios::app);
            
            if (timeStepCounter == 0) {
                timeSeriesFile << "<?xml version=\"1.0\"?>\n";
                timeSeriesFile << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
                timeSeriesFile << "<Collection>\n";
            }
            
            timeSeriesFile << "  <DataSet timestep=\"" << time << "\" file=\"simulation_" 
                          << std::setfill('0') << std::setw(6) << timeStepCounter << ".vtk\"/>\n";
            
            if (timeStepCounter == 100) { // Assuming 100 time steps
                timeSeriesFile << "</Collection>\n";
                timeSeriesFile << "</VTKFile>\n";
            }
            
            timeSeriesFile.close();
            timeStepCounter++;
        }
    };
} 