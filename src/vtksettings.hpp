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

        void saveTimeStep(const Rotor& rotor, float time, const ParticleSystem& particleSystem) {
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

            // Calculate total number of points (bound + wake particles)
            size_t totalPoints = particleSystem.boundParticles.size() + particleSystem.wakeParticles.size();

            // Write points
            file << "POINTS " << totalPoints << " float\n";
            
            // Write bound particle positions
            for (const auto& particle : particleSystem.boundParticles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.position[0] << " " 
                     << particle.position[1] << " " 
                     << particle.position[2] << "\n";
            }
            
            // Write wake particle positions
            for (const auto& particle : particleSystem.wakeParticles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.position[0] << " " 
                     << particle.position[1] << " " 
                     << particle.position[2] << "\n";
            }

            // Write point data
            file << "POINT_DATA " << totalPoints << "\n";
            
            // Write time data
            file << "SCALARS time float 1\n";
            file << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < totalPoints; ++i) {
                file << std::fixed << std::setprecision(6) << time << "\n";
            }

            // Write particle type (0 for bound, 1 for wake)
            file << "SCALARS particle_type int 1\n";
            file << "LOOKUP_TABLE default\n";
            // Write bound particle types
            for (const auto& particle : particleSystem.boundParticles) {
                file << "0\n";
            }
            // Write wake particle types
            for (const auto& particle : particleSystem.wakeParticles) {
                file << "1\n";
            }

            // Write circulation magnitude
            file << "SCALARS circulation_magnitude float 1\n";
            file << "LOOKUP_TABLE default\n";
            // Write bound particle circulation magnitudes
            for (const auto& particle : particleSystem.boundParticles) {
                file << std::fixed << std::setprecision(6) << particle.alpha.norm() << "\n";
            }
            // Write wake particle circulation magnitudes
            for (const auto& particle : particleSystem.wakeParticles) {
                file << std::fixed << std::setprecision(6) << particle.alpha.norm() << "\n";
            }

            // Write circulation vectors
            file << "VECTORS circulation_vectors float\n";
            // Write bound particle circulation vectors
            for (const auto& particle : particleSystem.boundParticles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.alpha[0] << " " 
                     << particle.alpha[1] << " " 
                     << particle.alpha[2] << "\n";
            }
            // Write wake particle circulation vectors
            for (const auto& particle : particleSystem.wakeParticles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.alpha[0] << " " 
                     << particle.alpha[1] << " " 
                     << particle.alpha[2] << "\n";
            }

            // Write velocity vectors
            file << "VECTORS velocity_vectors float\n";
            // Write bound particle velocities
            for (const auto& particle : particleSystem.boundParticles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.velocity[0] << " " 
                     << particle.velocity[1] << " " 
                     << particle.velocity[2] << "\n";
            }
            // Write wake particle velocities
            for (const auto& particle : particleSystem.wakeParticles) {
                file << std::fixed << std::setprecision(6) 
                     << particle.velocity[0] << " " 
                     << particle.velocity[1] << " " 
                     << particle.velocity[2] << "\n";
            }

            file.close();
            
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
            
            timeSeriesFile.close();
            timeStepCounter++;
        }

        ~VTKWriter() {
            std::string timeSeriesPath = folderPath + "/simulation.pvd";
            std::ofstream timeSeriesFile(timeSeriesPath, std::ios::app);
            timeSeriesFile << "</Collection>\n";
            timeSeriesFile << "</VTKFile>\n";
            timeSeriesFile.close();
        }
    };
}