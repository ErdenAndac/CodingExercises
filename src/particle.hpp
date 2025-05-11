#pragma once
#include "rotor.hpp"
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

namespace BET {
    // Gaussian kernel function for velocity calculation
    double q_sigma_gaussian(double rho) {
        return (1.0 / (4.0 * M_PI)) * (std::erf(rho / std::sqrt(2.0)) - std::sqrt(2.0 / M_PI) * rho * std::exp(-rho * rho / 2.0));
    }

    struct Particle {
        double x, y, z;           // Position
        double vx, vy, vz;        // Velocity
        double radius;            // Distance from motor hub (0,0,0)
        double circulation;       // Bound circulation value from blade element
        double initial_vx, initial_vy, initial_vz;        // Store initial z-velocity

        Particle(double pX, double pY, double pZ, double pUx, double pUy, double pUz, double circulation_ = 0.0)
            : x(pX), y(pY), z(pZ), vx(pUx), vy(pUy), vz(pUz), circulation(circulation_), initial_vx(pUx), initial_vy(pUy), initial_vz(pUz) {
            // Calculate radius as distance from motor hub (0,0,0)
            // radius = std::sqrt(x * x + y * y + z * z);
             radius = 0.5;
        }
    };

    class ParticleSystem {
    private:
        std::vector<Particle> particles;

        // Calculate induced velocities using Gaussian kernel
        void calculateInducedVelocities() {
            const size_t numParticles = particles.size();

            for (size_t i = 0; i < numParticles; i++) {
                // Reset velocities for this particle
                particles[i].vx = particles[i].initial_vx;
                particles[i].vy = particles[i].initial_vy;
                particles[i].vz = particles[i].initial_vz;  // Start with initial z-velocity

                for (size_t j = 0; j < numParticles; j++) {
                    if (i == j) continue;

                    // std::cout << "Ciculation: " << particles[i].circulation << std::endl;

                    // Calculate distance vector
                    double dx = particles[i].x - particles[j].x;
                    double dy = particles[i].y - particles[j].y;
                    double dz = particles[i].z - particles[j].z;
                    
                    // Calculate distance scalar
                    double distanceScalar = std::sqrt(dx*dx + dy*dy + dz*dz);
                    
                    if (distanceScalar < 1e-10) continue;  // Avoid self-interaction

                    // Calculate induced velocity using Gaussian kernel
                    double rho = distanceScalar / particles[j].radius;
                    double q_sigma = q_sigma_gaussian(rho);
                    double factor = q_sigma / (distanceScalar * distanceScalar * distanceScalar);

                    // Cross product of distance vector with circulation
                    // Using circulation as the strength of the vortex

                    double cross_x = dy * particles[j].circulation;
                    double cross_y = -dx * particles[j].circulation;
                    double cross_z = 0.0;

                    

                    // Add induced velocity to existing velocity
                    particles[i].vx += factor * cross_x;
                    particles[i].vy += factor * cross_y;
                    particles[i].vz += factor * cross_z;
                }
            }
        }

    public:
        ParticleSystem() = default;

        // Generate particles at all blade element positions
        void generateParticles(const Rotor& rotor, double time) {
            auto states = rotor.getBladeElementStates(time);
            auto circulationValues = rotor.calculateCirculation();
            
            // Generate particles at each blade element position
            for (size_t i = 0; i < states.size(); ++i) {
                const auto& [x, y, z, UR, UT, UP, Ublade_x, Ublade_y, Ublade_z] = states[i];
                const auto& [radius, circulation] = circulationValues[i];
                
                // Create particle with position, initial z-velocity from blade element, and circulation
                // Radius will be calculated in the Particle constructor
                particles.emplace_back(x, y, z, 0.0, 0.0, 0.0, circulation);
            }
        }

        // Update particle positions based on their velocities
        void updateParticles(double dt) {
            // Calculate induced velocities for all particles
            calculateInducedVelocities();

            for (auto& particle : particles) {
                // Update position using calculated velocities
                particle.x += particle.vx * dt;
                particle.y += particle.vy * dt;
                particle.z += particle.vz * dt;
                
                // Update radius after position update
                // particle.radius = std::sqrt(particle.x * particle.x + particle.y * particle.y + particle.z * particle.z);

                // std::cout << "particle.radius: " << particle.radius << std::endl;

            }
        }

        // Get all particles for VTK output
        std::vector<std::tuple<double, double, double, double, double, double>> getActiveParticles() const {
            std::vector<std::tuple<double, double, double, double, double, double>> activeParticles;
            for (const auto& particle : particles) {
                activeParticles.emplace_back(
                    particle.x, particle.y, particle.z,
                    particle.vx, particle.vy, particle.vz
                );
            }
            return activeParticles;
        }
    };
} 