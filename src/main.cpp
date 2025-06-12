#include "rotor.hpp"
#include "vtksettings.hpp"
#include "particle.hpp"
#include <iostream>

int main()
{
    std::cout << std::fixed << std::setprecision(5);
    // Define rotor parameters
    /*
    int numBlades = 2;
    float tipRadius = 12.0;       // [ft] tip radius
    float rootRadius = 0.48;      // [ft] cut-out radius
    int numElements = 10;         // number of elements per blade
    float chord = 0.6;            // [ft] chord length
    float twistRate = 0.0;       // [deg/ft] twist rate
    float rpm = 450.0;            // [rpm] angular (rotational) speed
    */

    // Sikorsky S-76
    int numBlades = 4;
    float tipRadius = 6.706;               // [m] tip radius
    float rootRadius = 0.1667 * tipRadius; // [m] cut-out radius
    int numElements = 10;                  // number of elements per blade
    float chord = 0.3937;                  // [m] chord length
    float twistRate = -1.4912;             // [deg/ft] twist rate
    float rpm = 293.0;                     // [rpm] angular (rotational) speed

    float omega = rpm * 2.0 * M_PI / 60.0;  // [rad/s] angular velocity
    float azimuthStep = 3.0 * M_PI / 180.0; // [rad]
    float timeStep = azimuthStep / omega;   // [s]

    BET::Rotor myRotor(numBlades, tipRadius, rootRadius, numElements, chord, twistRate, rpm);

    BET::VTKWriter vtkWriter;

    BET::ParticleSystem particleSystem;

    for (int index = 0; index <= 120 * 3; ++index)
    {
        float currentTime = index * timeStep;     // [s]
        float currentAngle = index * azimuthStep; // [rad]

        myRotor.getCirculation(currentTime);
        
        particleSystem.update(myRotor, timeStep);
       
        vtkWriter.saveTimeStep(myRotor, currentTime, particleSystem);

        
        if ((index) % 30 == 0)
        {
            std::cout << "Azimuth Angle: " << currentAngle << " radians" << std::endl;
            std::cout << "Time: " << currentTime << " seconds" << std::endl;
        }
        
    }

    return 0;
}