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
    float tipRadius = 6.706;       // [ft] tip radius
    float rootRadius = 0.1667 * tipRadius;      // [ft] cut-out radius
    int numElements = 10;         // number of elements per blade
    float chord = 0.3937;            // [ft] chord length
    float twistRate = -1.4912;       // [deg/ft] twist rate
    float rpm = 293.0;            // [rpm] angular (rotational) speed
    
    // Create Rotor (moved outside the time loop)
    BET::Rotor myRotor(numBlades, tipRadius, rootRadius, numElements, chord, twistRate, rpm);
    
    // Create VTK writer
    BET::VTKWriter vtkWriter;
    
    // Create Particle System
    BET::ParticleSystem particleSystem;
    
    // Simulate and save VTK files at different time steps
    float timeStep = 0.01; // [s] time step

    for (int timeIndex = 0; timeIndex <= 300; ++timeIndex)
    {
        float currentTime = timeIndex * timeStep;      
    
        myRotor.getCirculation(currentTime);

        particleSystem.generateParticles(myRotor, currentTime);

        particleSystem.calculateInducedVelocity();

        // Save time step data to VTK file
        vtkWriter.saveTimeStep(myRotor, currentTime, particleSystem);
        if ((timeIndex) % 100 == 0)
        {
            std::cout << "Time: " << currentTime << std::endl;
            
        }
    }

    return 0;
}