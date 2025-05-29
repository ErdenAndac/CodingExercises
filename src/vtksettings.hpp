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
        float time;
        std::vector<std::tuple<float, float, float, float, float, float>> bladeStates;
    };

    class VTKWriter {
    private:
        std::string folderPath;
        int timeStepCounter;

    public:
        VTKWriter() : folderPath("vtk_files_BETVPM"), timeStepCounter(0) {
            if (!fs::exists(folderPath)) {
                fs::create_directories(folderPath);
            }
        }

        void saveTimeStep(const Rotor& rotor, float time, const ParticleSystem& particleSystem = ParticleSystem()) {
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

            // Calculate total number of points (all blade elements + particles)
            size_t totalPoints = 0;
            for (const auto& blade : rotor.blades) {
                totalPoints += blade.elements.size();
            }
            totalPoints += particleSystem.particles.size();

            // Write points
            file << "POINTS " << totalPoints << " float\n";
            
            // Write blade element positions
            for (const auto& blade : rotor.blades) {
                for (const auto& element : blade.elements) {
                    file << std::fixed << std::setprecision(6) 
                         << element.position[0] << " " 
                         << element.position[1] << " " 
                         << element.position[2] << "\n";
                }
            }

            // Write particle positions
            for (const auto& particle : particleSystem.particles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.position[0] << " " 
                     << particle.position[1] << " " 
                     << particle.position[2] << "\n";
            }

            // Write blade connectivity
            size_t elementsPerBlade = rotor.blades[0].elements.size();
            file << "LINES " << rotor.numBlades << " " << (rotor.numBlades * (elementsPerBlade + 1)) << "\n";
            
            // Write connectivity for each blade
            for (int blade = 0; blade < rotor.numBlades; ++blade) {
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
            for (size_t i = 0; i < totalPoints; ++i) {
                file << std::fixed << std::setprecision(6) << time << "\n";
            }

            // Write circulation magnitude
            file << "SCALARS circulation_magnitude float 1\n";
            file << "LOOKUP_TABLE default\n";
            // Write blade element circulation
            for (const auto& blade : rotor.blades) {
                for (const auto& element : blade.elements) {
                    file << std::fixed << std::setprecision(6) << element.circulation.norm() << "\n";
                }
            }
            // Write particle circulation
            for (const auto& particle : particleSystem.particles) {
                file << std::fixed << std::setprecision(6) << particle.alpha.norm() << "\n";
            }

            // Write circulation vectors
            file << "VECTORS circulation_vectors float\n";
            // Write blade element circulation vectors
            for (const auto& blade : rotor.blades) {
                for (const auto& element : blade.elements) {
                    file << std::fixed << std::setprecision(6) 
                         << element.circulation[0] << " " 
                         << element.circulation[1] << " " 
                         << element.circulation[2] << "\n";
                }
            }
            // Write particle circulation vectors
            for (const auto& particle : particleSystem.particles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.alpha[0] << " " 
                     << particle.alpha[1] << " " 
                     << particle.alpha[2] << "\n";
            }

            file.close();
            // std::cout << "VTK file written to: " << filePath << std::endl;
            
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
            
            // Close the time series file properly
            timeSeriesFile.close();
            timeStepCounter++;
        }

        // Add a destructor to properly close the PVD file
        ~VTKWriter() {
            std::string timeSeriesPath = folderPath + "/simulation.pvd";
            std::ofstream timeSeriesFile(timeSeriesPath, std::ios::app);
            timeSeriesFile << "</Collection>\n";
            timeSeriesFile << "</VTKFile>\n";
            timeSeriesFile.close();
        }
    };
}