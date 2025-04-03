#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>

// Constants
const double PI = 3.14159265359;
const double RHO = 0.002377; // Air density at sea level [slug/ft^3]
const double A0 = 2*PI;      // Lift curve slope [1/rad] (2π for thin airfoil theory)

// Blade Element struct to store local properties
struct BladeElement
{
    double radius;       // Local radius [ft]
    double chord;        // Local chord length [ft]
    double twist;        // Local twist angle [deg]
    double solidity;     // Local solidity ratio
    double inducedAngle; // Induced angle [rad]
    double inflowAngle;  // Inflow angle [rad]
    double alpha;        // Angle of attack [rad]
    double cl;           // Lift coefficient
    double cd;           // Drag coefficient
    double dL;           // Differential lift [lbf]
    double dD;           // Differential drag [lbf]
    double dT;           // Differential thrust [lbf]
    double dQ;           // Differential torque [ft-lbf]
};

// Blade struct to store blade properties and elements
struct Blade
{
    int numElements;                    // Number of elements per blade
    double rootRadius;                  // Root radius [ft]
    double tipRadius;                   // Tip radius [ft]
    double chord;                       // Chord length [ft]
    double twistRate;                   // Twist rate [deg/ft]
    std::vector<BladeElement> elements; // Vector of blade elements
};

// Rotor struct to store rotor properties and blades
struct Rotor
{
    int numBlades;             // Number of blades
    double rpm;                // Rotational speed [rpm]
    double omega;              // Angular velocity [rad/s]
    double verticalVelocity;   // Vertical velocity [ft/s]
    std::vector<Blade> blades; // Vector of blades
};

// Function to convert degrees to radians
double degToRad(double deg)
{
    return deg * PI / 180.0;
}

// Function to convert radians to degrees
double radToDeg(double rad)
{
    return rad * 180.0 / PI;
}

// Function to initialize a blade element
BladeElement initializeElement(const Blade &blade, int elementIndex)
{
    BladeElement element;

    // Calculate local radius
    double dr = (blade.tipRadius - blade.rootRadius) / blade.numElements;
    element.radius = blade.rootRadius + (elementIndex + 0.5) * dr;

    // Calculate local twist
    element.twist = degToRad(blade.twistRate * (element.radius - blade.rootRadius));

    // Set chord length (constant for this case)
    element.chord = blade.chord;

    // Calculate local solidity ratio
    element.solidity = (blade.numElements * element.chord) / (2.0 * PI * element.radius);

    // Initialize other parameters
    element.inducedAngle = 0.0;
    element.inflowAngle = 0.0;
    element.alpha = 0.0;
    element.cl = 0.0;
    element.cd = 0.0;
    element.dL = 0.0;
    element.dD = 0.0;
    element.dT = 0.0;
    element.dQ = 0.0;

    return element;
}

// Function to calculate induced velocity using momentum theory
double calculateInducedVelocity(const Rotor &rotor, double radius)
{
    double lambda_i = 0.0;
    double lambda_c = rotor.verticalVelocity / (rotor.omega * rotor.blades[0].tipRadius);

    // Simple momentum theory for hover
    if (std::abs(lambda_c) < 1e-6)
    {
        lambda_i = 0.5 * (std::sqrt(lambda_c * lambda_c + 4.0) - lambda_c);
    }
    else
    {
        // For vertical flight
        lambda_i = 0.5 * (std::sqrt(lambda_c * lambda_c + 4.0) - lambda_c);
    }

    return lambda_i * rotor.omega * rotor.blades[0].tipRadius;
}

// Function to calculate blade element forces
void calculateElementForces(BladeElement &element, const Rotor &rotor)
{
    // Calculate local velocity components
    double Vr = rotor.omega * element.radius;
    double Vi = calculateInducedVelocity(rotor, element.radius);
    double V = std::sqrt(Vr * Vr + (rotor.verticalVelocity + Vi) * (rotor.verticalVelocity + Vi));

    // Calculate inflow angle
    element.inflowAngle = std::atan2(rotor.verticalVelocity + Vi, Vr);

    // Calculate angle of attack
    element.alpha = element.twist - element.inflowAngle;

    // Calculate lift coefficient (simplified linear model)
    element.cl = A0 * element.alpha;

    // Calculate drag coefficient (simplified model)
    element.cd = 0.01 + 0.1 * element.cl * element.cl;

    // Calculate differential forces
    double dL = 0.5 * RHO * V * V * element.chord * element.cl;
    double dD = 0.5 * RHO * V * V * element.chord * element.cd;

    // Calculate thrust and torque components
    element.dT = dL * std::cos(element.inflowAngle) - dD * std::sin(element.inflowAngle);
    element.dQ = (dL * std::sin(element.inflowAngle) + dD * std::cos(element.inflowAngle)) * element.radius;
}

// Function to initialize and solve the rotor
Rotor initializeRotor()
{
    Rotor rotor;

    // Set rotor properties from given values
    rotor.numBlades = 2;
    rotor.rpm = 450.0;
    rotor.omega = 2.0 * PI * rotor.rpm / 60.0;
    rotor.verticalVelocity = 0.0;

    // Initialize blades
    for (int i = 0; i < rotor.numBlades; i++)
    {
        Blade blade;
        blade.numElements = 100;
        blade.rootRadius = 0.48;
        blade.tipRadius = 12.0;
        blade.chord = 0.6;
        blade.twistRate = 0.0;

        // Initialize elements for each blade
        for (int j = 0; j < blade.numElements; j++)
        {
            blade.elements.push_back(initializeElement(blade, j));
        }

        rotor.blades.push_back(blade);
    }

    return rotor;
}

// Function to solve the rotor performance
void solveRotor(Rotor &rotor)
{
    double totalThrust = 0.0;
    double totalTorque = 0.0;

    // Calculate forces for each blade element
    for (auto &blade : rotor.blades)
    {
        for (auto &element : blade.elements)
        {
            calculateElementForces(element, rotor);

            // Integrate forces
            double dr = (blade.tipRadius - blade.rootRadius) / blade.numElements;
            totalThrust += element.dT * dr;
            totalTorque += element.dQ * dr;
        }
    }

    // Multiply by number of blades
    totalThrust *= rotor.numBlades;
    totalTorque *= rotor.numBlades;

    // Print results
    std::cout << "Rotor Performance Results:" << std::endl;
    std::cout << "Total Thrust: " << totalThrust << " lbf" << std::endl;
    std::cout << "Total Torque: " << totalTorque << " ft-lbf" << std::endl;
    std::cout << "Power Required: " << totalTorque * rotor.omega / 550.0 << " hp" << std::endl;
}

// Function to get current timestamp for filename
std::string getTimestamp()
{
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return oss.str();
}

// Function to create VTK directory if it doesn't exist
void createVTKDirectory()
{
    std::filesystem::path dir("vtk_files_BE_1");
    if (!std::filesystem::exists(dir))
    {
        std::filesystem::create_directory(dir);
    }
}

// Function to generate VTK file for blade elements
void generateVTKFile(const Rotor &rotor)
{
    createVTKDirectory();

    // Time step parameters
    const int numTimeSteps = 100;
    const double timeStep = 0.01; // seconds
    const double totalTime = numTimeSteps * timeStep;

    // For each time step
    for (int step = 0; step < numTimeSteps; ++step)
    {
        double currentTime = step * timeStep;
        double baseRotationAngle = rotor.omega * currentTime;

        // Create filename with timestamp and step number
        std::string filename = "vtk_files_BE_1/blade_elements_" + getTimestamp() +
                               "_step" + std::to_string(step) + ".vtk";
        std::ofstream vtkFile(filename);

        if (!vtkFile.is_open())
        {
            std::cerr << "Error: Could not create VTK file: " << filename << std::endl;
            continue;
        }

        // Count total number of points and cells
        int totalPoints = 0;
        int totalCells = 0;
        for (const auto &blade : rotor.blades)
        {
            totalPoints += blade.elements.size();
            totalCells += blade.elements.size() - 1;
        }

        // Write VTK header
        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "Blade Element Data - Time: " << currentTime << " seconds\n";
        vtkFile << "ASCII\n";
        vtkFile << "DATASET POLYDATA\n\n";

        // Write points
        vtkFile << "POINTS " << totalPoints << " double\n";
        for (size_t bladeIndex = 0; bladeIndex < rotor.blades.size(); ++bladeIndex)
        {
            const auto &blade = rotor.blades[bladeIndex];
            // Calculate rotation angle for this blade
            double bladeRotationAngle = baseRotationAngle + (2.0 * PI * bladeIndex / rotor.numBlades);

            for (const auto &element : blade.elements)
            {
                // Calculate x, y coordinates with rotation
                double x = element.radius * std::cos(bladeRotationAngle);
                double y = element.radius * std::sin(bladeRotationAngle);
                double z = 0.0;
                vtkFile << x << " " << y << " " << z << "\n";
            }
        }
        vtkFile << "\n";

        // Write lines
        vtkFile << "LINES " << totalCells << " " << totalCells * 3 << "\n";
        int pointIndex = 0;
        for (const auto &blade : rotor.blades)
        {
            for (size_t i = 0; i < blade.elements.size() - 1; ++i)
            {
                vtkFile << "2 " << pointIndex << " " << pointIndex + 1 << "\n";
                pointIndex++;
            }
            pointIndex++; // Skip the last point of each blade
        }
        vtkFile << "\n";

        // Write point data
        vtkFile << "POINT_DATA " << totalPoints << "\n";

        // Write velocity components as vectors
        vtkFile << "VECTORS velocity double\n";
        for (size_t bladeIndex = 0; bladeIndex < rotor.blades.size(); ++bladeIndex)
        {
            const auto &blade = rotor.blades[bladeIndex];
            double bladeRotationAngle = baseRotationAngle + (2.0 * PI * bladeIndex / rotor.numBlades);

            for (const auto &element : blade.elements)
            {
                // Calculate local velocity components
                double Vr = rotor.omega * element.radius;
                double Vi = calculateInducedVelocity(rotor, element.radius);

                // Transform velocity components to global coordinates
                double Vx = Vr * std::cos(bladeRotationAngle);
                double Vy = Vr * std::sin(bladeRotationAngle);
                double Vz = rotor.verticalVelocity + Vi;

                vtkFile << Vx << " " << Vy << " " << Vz << "\n";
            }
        }

        vtkFile.close();
       // std::cout << "VTK file generated: " << filename << " (Time: " << currentTime << " seconds)" << std::endl;
    }
}

int main()
{
    // Initialize and solve the rotor
    Rotor rotor = initializeRotor();
    solveRotor(rotor);

    // Generate VTK file
    generateVTKFile(rotor);

    return 0;
}
