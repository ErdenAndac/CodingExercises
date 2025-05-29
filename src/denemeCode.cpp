#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <filesystem> // C++17 for filesystem operations
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace Eigen;
using vector = Vector3f;

int main()
{
    float rpm = 450.0;
    float dr = (12.0 - 0.48) / 10.0;
    float vi = 25.4232;
    vector omega = vector(0.0, 0.0, rpm * 2 * M_PI / 60.0);

    // Set fixed precision for all floating point numbers
    std::cout << std::fixed << std::setprecision(3);

    float psi = 0.0;
    float theta = -8.0 * M_PI / 180.0;
    float beta = 0.0;

    Matrix3f T;
    T = AngleAxis<float>(psi, vector::UnitZ()) * AngleAxis<float>(theta, vector::UnitY()) * AngleAxis<float>(beta, vector::UnitX());

    // Print table header
    std::cout << std::left << std::setw(20) << "Variable"
              << std::setw(15) << "X"
              << std::setw(15) << "Y"
              << std::setw(15) << "Z" << "\n";
    std::cout << std::string(65, '-') << "\n";

    for (int i = 0; i < 10; ++i)
    {

        std::cout << "Blade element " << i + 1 << ":\n";
        std::cout << std::string(65, '-') << "\n";

        // Update radius for this blade element
        vector radius = vector(0.0, 0.48 + i * dr + dr / 2.0, 0.0);
        std::cout << std::left << std::setw(20) << "Radius"
                  << std::setw(15) << radius[0]
                  << std::setw(15) << radius[1]
                  << std::setw(15) << radius[2] << "\n";

        // Calculate velocities
        vector BladeVelocity = T * omega.cross(radius);
        std::cout << std::left << std::setw(20) << "Blade velocity"
                  << std::setw(15) << BladeVelocity[0]
                  << std::setw(15) << BladeVelocity[1]
                  << std::setw(15) << BladeVelocity[2] << "\n";

        vector InducedVelocity = T * vector(0.0, 0.0, -vi);
        std::cout << std::left << std::setw(20) << "Induced velocity"
                  << std::setw(15) << InducedVelocity[0]
                  << std::setw(15) << InducedVelocity[1]
                  << std::setw(15) << InducedVelocity[2] << "\n";

        vector TotalVelocity = -BladeVelocity + InducedVelocity;
        std::cout << std::left << std::setw(20) << "Total velocity"
                  << std::setw(15) << TotalVelocity[0]
                  << std::setw(15) << TotalVelocity[1]
                  << std::setw(15) << TotalVelocity[2] << "\n\n";
    }

    std::cout << T << "\n";
}