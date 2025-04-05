#include "rotor.hpp"
#include "vtksettings.hpp"
#include <iostream>

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
    BET::Rotor myRotor(numBlades, tipRadius, rootRadius, numElements, chord, twistRate, rpm, MTOW);

    /// Calculate thrust at each blade element
    std::vector<double> thrustValues = myRotor.calculateThrust(verticalVelocity);

    // Print values
    std::cout << "Thrust values at each blade element:\n";
    for (size_t i = 0; i < thrustValues.size(); ++i)
    {
        std::cout << "Blade Element " << i + 1 << ": " << thrustValues[i] << " \n";
    }

    // Simulate and save VTK files at different time steps
    double timeStep = 0.01;
    for (int timeIndex = 0; timeIndex <= 100; ++timeIndex)
    {
        BET::saveVTK(myRotor, timeIndex * timeStep, timeIndex);
    }

    return 0;
} 