#define _USE_MATH_DEFINES // for C++ Math Constants
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <ctime>

namespace fs = std::filesystem;

const double rhoAir = 0.0023769; // [slug/ft^3] density (sea-level standard)

// Blade Element Struct
struct BladeElement
{
    double radius; // [ft] radial position
    double chord;  // [ft] chord length
    double twist;  // [deg] twist angle

    BladeElement(double r, double c, double t)
        : radius(r), chord(c), twist(t) {}
};

// Blade Struct
struct Blade
{
    std::vector<BladeElement> elements; // list of blade elements
    double initialAngle;                // [rad] initial angular position
    double R;
    double CL_alpha;
    double CD_0;
    double twistRate;

    // Constructor with Twist Distribution
    Blade(double tipRadius, double rootRadius, int numElements, double chord, double twistRateInput, double angle)
        : initialAngle(angle),
          R(tipRadius),
          CL_alpha((2 * M_PI) * (M_PI / 180)),
          CD_0(0.0),
          twistRate(twistRateInput)
    {
        double dr = (tipRadius - rootRadius) / numElements; // Blade element spacing
        for (int i = 0; i < numElements; ++i)
        {
            double r = rootRadius + i * dr + dr / 2.0; // Midpoint of blade element

            double twist = twistRate * (r / R);

            elements.emplace_back(r, chord, twist);
        }
    }
};

// Rotor Struct
struct Rotor
{
    std::vector<Blade> blades; // store multiple blades
    int numBlades;             // number of blades
    double angularVelocity;    // [rad/s] angular (rotation) speed
    double discArea;           // [ft^2] rotor disk area
    double bladeArea;          // [ft^2] total blade area
    double solidity;           // [-] solidity
    double thrust;             // [lb] @ hover = MTOW
    double CT_mt;              // [-] thrust coefficient (momentum theory)
    double v_i;                // [ft/s] induced velocity (momentum theory)

    // Create a rotor with identical blades
    Rotor(int nBlades, double tipRadius, double rootRadius, int numElements, double chord, double twistRateInput, double rpm, double MTOW)
        : numBlades(nBlades),
          angularVelocity(rpm * 2.0 * M_PI / 60.0),
          discArea(M_PI * tipRadius * tipRadius),
          bladeArea(numBlades * (tipRadius - rootRadius) * chord),
          solidity(bladeArea / discArea),
          thrust(MTOW),
          CT_mt(thrust / (rhoAir * discArea * angularVelocity * angularVelocity * tipRadius * tipRadius)),
          v_i(sqrt(thrust / (2.0 * rhoAir * discArea)))
    {
        double angleSpacing = 2.0 * M_PI / numBlades;
        for (int i = 0; i < numBlades; ++i)
        {
            blades.emplace_back(tipRadius, rootRadius, numElements, chord, twistRateInput, i * angleSpacing);
        }
    }

    // Get blade element positions for a given time step
    std::vector<std::tuple<double, double, double>> getBladeElementPositions(double time) const
    {
        std::vector<std::tuple<double, double, double>> positions;
        for (const auto &blade : blades)
        {
            double theta = fmod(blade.initialAngle + angularVelocity * time, 2.0 * M_PI);
            for (const auto &element : blade.elements)
            {
                double x = element.radius * cos(theta);
                double y = element.radius * sin(theta);
                double z = 0.0;
                positions.emplace_back(x, y, z);
            }
        }
        return positions;
    }
    // Calculate thrust at each blade element
    std::vector<double> calculateCirculation(double V_c) const
    {
        std::vector<double> thrustValues;
        for (const auto &blade : blades)
        {
            for (const auto &element : blade.elements)
            {
                // Velocity components
                // Use momentum theory to calculate induced velocity
                double U_P = V_c + v_i;                              // [ft/s] out-of-plane velocity component
                double U_T = angularVelocity * element.radius;       // [ft/s] in-plane velocity component
                double phirad = atan(U_P / U_T);                     // [rad] inflow angle
                double phideg = phirad * 180.0 / M_PI;               // [deg] inflow angle
                double lambda = phirad * (element.radius / blade.R); // [rad] inflow ratio
                double lambdadeg = lambda * 180.0 / M_PI;            // [deg] inflow ratio
                double pitchAngle;
                if (blade.twistRate == 0.0)
                {
                    pitchAngle = (6.0 * CT_mt) / (solidity * blade.CL_alpha) + (3.0 / 2.0) * sqrt(CT_mt / 2.0);
                }
                else
                {
                    pitchAngle = ((6.0 * CT_mt) / (solidity * blade.CL_alpha)) - (3.0 * blade.twistRate / 4.0) + (3.0 * lambdadeg * 2.0);
                }
                double theta = pitchAngle;
                double CL = blade.CL_alpha * (theta - phideg);
                double boundCirc = 0.5 * CL * angularVelocity * element.radius * element.chord;

                /* double dCP_0 = (solidity / 2.0) * (blade.CD_0 * (element.radius / blade.R) * (element.radius / blade.R) * (element.radius / blade.R));
                double dCT = 0.5 * solidity * CL * (element.radius / blade.R) * (element.radius / blade.R);
                double dCP_i = (solidity / 2.0) * (CL * lambda *  (element.radius / blade.R) * (element.radius / blade.R) * (element.radius / blade.R)); */
                thrustValues.emplace_back(boundCirc);

            }
        }
        return thrustValues;
    }
};

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
    // std::cout << "Saved VTK file: " << filePath << std::endl;
}

int main()
{
    int numBlades = 2;
    double tipRadius = 12.0;       // [ft] tip radius
    double rootRadius = 0.48;      // [ft] cut-out radius
    int numElements = 100;         // number of elements per blade
    double chord = 0.6;            // [ft] chord length
    double twistRate = -8.0;       // [deg/ft] twist rate
    double rpm = 450.0;            // [rpm] angular (rotational) speed
    double MTOW = 1390.0;          // [lbf] total thrust
    double verticalVelocity = 0.0; // [ft/s] vertical velocity

    // Create Rotor
    Rotor myRotor(numBlades, tipRadius, rootRadius, numElements, chord, twistRate, rpm, MTOW);

    /// Calculate thrust at each blade element
    std::vector<double> circulationValues = myRotor.calculateCirculation(verticalVelocity);
    // Print thrust values
    std::cout << "Circulation values at each blade element:\n";
    for (size_t i = 0; i < circulationValues.size(); ++i)
    {
        std::cout << "Blade Element " << i + 1 << ": " << circulationValues[i] << " \n";
    }

    // Simulate and save VTK files at different time steps
    double timeStep = 0.01;
    for (int timeIndex = 0; timeIndex <= 100; ++timeIndex)
    {
        saveVTK(myRotor, timeIndex * timeStep, timeIndex);
    }

    return 0;
}


