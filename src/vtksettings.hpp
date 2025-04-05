#pragma once
#include "rotor.hpp"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <ctime>
#include <iostream>

namespace fs = std::filesystem;

namespace BET {
    // Function to save blade positions as .vtk file
    void saveVTK(const Rotor &rotor, double time, int timeIndex)
    {
        std::string folderPath = "vtk_files_BE_1";
        if (!fs::exists(folderPath))
        {
            fs::create_directories(folderPath);
        }

        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm tm = *std::localtime(&now_time);
        char buffer[100];
        std::strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S", &tm);
        std::string dateTime(buffer);

        std::string filePath = folderPath + "/blade_" + dateTime + "_" + std::to_string(timeIndex) + ".vtk";
        std::ofstream file(filePath);
        if (!file.is_open())
        {
            std::cerr << "Error: Could not open file for writing VTK output." << std::endl;
            return;
        }

        file << "# vtk DataFile Version 3.0\n";
        file << "Blade Element Positions\n";
        file << "ASCII\n";
        file << "DATASET POLYDATA\n";
        file << "FIELD FieldData 1\n";
        file << "TIME 1 1 float\n";
        file << time << "\n";

        auto positions = rotor.getBladeElementPositions(time);
        file << "POINTS " << positions.size() << " float\n";

        for (const auto &[x, y, z] : positions)
        {
            file << x << " " << y << " " << z << "\n";
        }

        file.close();
    }
} 