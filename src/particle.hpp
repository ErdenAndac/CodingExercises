#pragma once
#include "rotor.hpp"
#include <tuple>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
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

    // Helper struct for grid cell coordinates
    struct GridCell {
        int x, y, z;
        
        bool operator==(const GridCell& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    // Hash function for GridCell
    struct GridCellHash {
        std::size_t operator()(const GridCell& cell) const {
            return std::hash<int>()(cell.x) ^ 
                   (std::hash<int>()(cell.y) << 1) ^ 
                   (std::hash<int>()(cell.z) << 2);
        }
    };

    struct ParticleSystem
    {
        std::vector<Particle> particles;
        int currentTimeStep = 0; // Track current time step
        float cellSize = 2.0f; // Size of each grid cell

        // Convert position to grid cell coordinates
        GridCell positionToGridCell(const vector& position) const {
            return {
                static_cast<int>(std::floor(position[0] / cellSize)),
                static_cast<int>(std::floor(position[1] / cellSize)),
                static_cast<int>(std::floor(position[2] / cellSize))
            };
        }

        // Get neighboring cells for a given cell
        std::vector<GridCell> getNeighboringCells(const GridCell& cell) const {
            std::vector<GridCell> neighbors;
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        neighbors.push_back({cell.x + dx, cell.y + dy, cell.z + dz});
                    }
                }
            }
            return neighbors;
        }

        void generateParticles(const Rotor &rotor, float time)
        {
            // Generate particles only at blade tips
            for (const auto &blade : rotor.blades)
            {
                // Get the last element (tip) of each blade
                const auto &tipElement = blade.elements.back();

                // Create particle at blade tip position with current time step
                particles.emplace_back(
                    tipElement.position,                    // Position from blade tip
                    vector::Zero(), // Initial velocity
                    tipElement.circulation,                 // Circulation from blade tip
                    currentTimeStep                         // Current time step
                );
            }
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

            // Create spatial hash map
            std::unordered_map<GridCell, std::vector<int>, GridCellHash> spatialGrid;
            
            // Assign particles to grid cells
            for (int i = 0; i < particles.size(); ++i) {
                GridCell cell = positionToGridCell(particles[i].position);
                spatialGrid[cell].push_back(i);
            }

            // Calculate induced velocities using spatial partitioning
            for (int particleIndex1 = 0; particleIndex1 < particles.size(); ++particleIndex1)
            {
                GridCell cell1 = positionToGridCell(particles[particleIndex1].position);
                auto neighboringCells = getNeighboringCells(cell1);

                // Check particles in neighboring cells
                for (const auto& neighborCell : neighboringCells) {
                    auto it = spatialGrid.find(neighborCell);
                    if (it == spatialGrid.end()) continue;

                    for (int particleIndex2 : it->second) {
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
                }

                // Update particle position using the calculated velocity
                float timeStep = 0.01; // Using the same time step as in main.cpp
                particles[particleIndex1].position += particles[particleIndex1].velocity * timeStep;
            }
        }
    };
};