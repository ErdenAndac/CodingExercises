#pragma once
#include "rotor.hpp"
#include <tuple>
#include <cmath>
#include <algorithm>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

using vector = Eigen::Vector3f;

namespace BET
{
    // Gaussian kernel function for velocity calculation
    float q_sigma_gaussian(float rho)
    {
        return (1.0 / (4.0 * M_PI)) * (std::erf(rho / std::sqrt(2.0)) - std::sqrt(2.0 / M_PI) * rho * std::exp(-rho * rho / 2.0));
    }

    struct Particle
    {
        vector position;
        vector velocity;
        vector alpha;
        int timeStep; // Add time step tracking

        Particle(const vector &pos, const vector &vel, const vector &alpha, int step)
            : position(pos), velocity(vel), alpha(alpha), timeStep(step) {}
    };

    struct ParticleSystem
    {
        std::vector<Particle> particles;
        int currentTimeStep = 0; // Track current time step
        float timeStep; // Time step for simulation
        std::vector<vector> previousCirculations; // Store previous circulation values

        ParticleSystem(float dt = 0.01f) : timeStep(dt) {} // Constructor with default time step

        void generateParticles(const Rotor &rotor, float time)
        {
            // Store current circulation values for next time step
            std::vector<vector> currentCirculations;
            for (const auto &blade : rotor.blades)
            {
                const auto &tipElement = blade.elements.back();
                currentCirculations.push_back(tipElement.circulation);
            }

            // Skip particle generation in the first time step
            if (currentTimeStep > 0)
            {
                // Generate particles using previous circulation values
                for (size_t i = 0; i < rotor.blades.size(); ++i)
                {
                    const auto &blade = rotor.blades[i];
                    const auto &tipElement = blade.elements.back();
                    
                    // Use previous circulation values
                    vector circulation = previousCirculations[i];

                    // Create particle at blade tip position with current time step
                    particles.emplace_back(
                        tipElement.position,                    // Position from blade tip
                        vector::Zero(),                        // Initial velocity
                        circulation,                           // Circulation from previous time step
                        currentTimeStep                        // Current time step
                    );
                }
            }

            // Update previous circulations for next time step
            previousCirculations = currentCirculations;
            currentTimeStep++; // Increment time step after generating particles
        }

        auto calculateInducedVelocity()
        {
            // Remove oldest particles if total count exceeds 1500
            if (particles.size() > 1500) {
                // Sort particles by timeStep (oldest first)
                std::sort(particles.begin(), particles.end(),
                    [](const Particle& a, const Particle& b) {
                        return a.timeStep < b.timeStep;
                    }
                );
                // Remove oldest particles to get back to 1500
                particles.erase(particles.begin(), particles.begin() + (particles.size() - 800));
            }

            // Calculate induced velocities
            for (int particleIndex1 = 0; particleIndex1 < particles.size(); ++particleIndex1)
            {
                particles[particleIndex1].velocity = vector::Zero();
                
                for (int particleIndex2 = 0; particleIndex2 < particles.size(); ++particleIndex2)
                {
                    if (particleIndex1 == particleIndex2) continue;

                    vector distanceVector = (particles[particleIndex1].position - particles[particleIndex2].position);
                    float distanceScalar = distanceVector.norm();
                    
                    // Skip if particles are too close to avoid singularity
                    if (distanceScalar < 1e-6 || distanceScalar > 15)
                        continue;

                    float radius = distanceScalar * 1.2;
                    float rho = distanceScalar / radius;
                    float q_sigma = q_sigma_gaussian(rho);
                    float factor = q_sigma / (distanceScalar * distanceScalar * distanceScalar);
                    vector inducedVelocity = factor * distanceVector.cross(particles[particleIndex2].alpha);

                    particles[particleIndex1].velocity += inducedVelocity;
                }

                // Update particle position using the calculated velocity
                particles[particleIndex1].position += particles[particleIndex1].velocity * timeStep;
            }
        }
    };
};