#include "rotor.hpp"
#include "vtksettings.hpp"
#include "particle.hpp"
#include <iostream>

int main()
{
    int numBlades = 2;
    double tipRadius = 12.0;       // [ft] tip radius
    double rootRadius = 0.48;      // [ft] cut-out radius
    int numElements = 10.0;         // number of elements per blade
    double chord = 0.6;            // [ft] chord length
    double twistRate = -8.0;       // [deg/ft] twist rate
    double rpm = 450.0;            // [rpm] angular (rotational) speed
    double MTOW = 1390.0;          // [lbf] total thrust
    
    // Create VTK writer
    BET::VTKWriter vtkWriter;
    
    // Simulate and save VTK files at different time steps
    double timeStep = 0.01; // [s] time step
    BET::ParticleSystem particleSystem; // Create particle system

    for (int timeIndex = 0; timeIndex <= 50; ++timeIndex)
    {
        double currentTime = timeIndex * timeStep;
        
        // Create Rotor
        BET::Rotor myRotor(numBlades, tipRadius, rootRadius, numElements, chord, twistRate, rpm, MTOW);

        // Generate particles at time index 0 and every 5 steps
        if (timeIndex == 0 || timeIndex % 1 == 0) {
            particleSystem.generateParticles(myRotor, currentTime);
        }
        
        // Get blade element states
        auto states = myRotor.getBladeElementStates(currentTime);
        
        // Calculate thrust values
        auto circulationValues = myRotor.calculateCirculation();

        // Save time step data to VTK file
        vtkWriter.saveTimeStep(myRotor, currentTime, &particleSystem);

        // Update existing particles
        particleSystem.updateParticles(timeStep);
    }

    return 0;
}