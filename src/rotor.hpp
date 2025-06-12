#pragma once
#include "blade.hpp"
#include "constants.hpp"
#include <vector>
#include <cmath>
#include <tuple>
#include <array>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "SC1095PolarData.h"

using namespace Eigen;
using vector = Vector3f;

namespace BET
{
    struct Rotor
    {
        int numBlades;             // number of blades
        std::vector<Blade> blades; // store multiple blades
        vector omega;              // [rad/s] angular (rotation) speed
        vector v_i;                // [ft/s] induced velocity
        vector Ubody;              // body velocity
        float solidity;

        // Create a rotor with identical blades
        Rotor(int nBlades, float tipRadius, float rootRadius, int numElements, float chord, float twistRate, float rpm)
            : numBlades(nBlades),
              omega(vector(0.0, 0.0, rpm * 2.0 * M_PI / 60.0)),
              v_i(vector(0.0, 0.0, -12.141)), // [m/s] 
              Ubody(vector(0.0, 0.0, 0.0)),
              solidity((numBlades * chord) / (tipRadius * M_PI))
        {
            float angleSpacing = 2.0 * M_PI / numBlades;
            for (int i = 0; i < numBlades; ++i)
            {
                float initialAngle = 0.0 * M_PI / 180.0 + i * angleSpacing;
                blades.emplace_back(tipRadius, rootRadius, numElements, chord, twistRate, initialAngle);
            }
        }

        // Helper function to interpolate SC1095 polar data
        float interpolateCL(float alphaDeg) const {
            // Clamp alpha to the range of available data
            alphaDeg = std::max(-15.0f, std::min(16.75f, alphaDeg));
            
            // Find the two closest data points
            int lowerIndex = 0;
            while (lowerIndex < SC1095PolarDataSize - 1 && SC1095PolarData[lowerIndex + 1][0] < alphaDeg) {
                lowerIndex++;
            }
            
            // If we're at the last point, return its CL value
            if (lowerIndex == SC1095PolarDataSize - 1) {
                return SC1095PolarData[lowerIndex][1];
            }
            
            // Linear interpolation
            float alpha1 = SC1095PolarData[lowerIndex][0];
            float alpha2 = SC1095PolarData[lowerIndex + 1][0];
            float CL1 = SC1095PolarData[lowerIndex][1];
            float CL2 = SC1095PolarData[lowerIndex + 1][1];
            
            return CL1 + (CL2 - CL1) * (alphaDeg - alpha1) / (alpha2 - alpha1);
        }

        // Get blade element positions and velocities for a given time step
        auto getCirculation(float time) const
        {
            std::vector<std::tuple<float, float, float, float, float, float>> states;

            for (int bladeIndex = 0; bladeIndex < blades.size(); ++bladeIndex)
            {
                auto &blade = blades[bladeIndex];
                float psi = ((blade.initialAngle + omega[2] * time));
                float beta = -(blade.flappingAngle);

                for (int elementIndex = 0; elementIndex < blade.elements.size(); ++elementIndex)
                {
                    auto &element = blade.elements[elementIndex];

                    float theta = -(blade.elements[elementIndex].twist);

                    // First rotation matrix for blade position (CCW rotation around z-axis)
                    Matrix3f Rz = AngleAxis<float>(psi, vector::UnitZ()).toRotationMatrix();

                    // Second rotation matrix for blade pitch and flapping
                    Matrix3f Ryx = (AngleAxis<float>(theta, vector::UnitY()) *
                                    AngleAxis<float>(beta, vector::UnitX()))
                                       .toRotationMatrix();

                    Matrix3f Rzx = (AngleAxis<float>(psi, vector::UnitZ()) *
                                    AngleAxis<float>(beta, vector::UnitX()))
                                       .toRotationMatrix();

                    // Combined transformation matrix
                    Matrix3f R = Rz * Ryx;

                    // Calculate element position in global frame
                    // For initial angle 0°, blade is along x-axis
                    vector elementPosition;

                    elementPosition = Rzx * vector(element.radius, 0.0, 0.0);

                    // Update element position
                    blade.elements[elementIndex].updatePosition(elementPosition);

                    // Velocities in blade frame
                    vector bladeVelocity = Ryx * omega.cross(vector(0.0, element.radius, 0.0));
                    vector inducedVelocity = Ryx * v_i;
                    vector totalVelocity = -bladeVelocity + inducedVelocity;
                    vector Ublade = totalVelocity + R * Ubody;

                    blade.elements[elementIndex].updateUblade(Ublade);

                    float alphaRad = atan2(Ublade[2], Ublade[0]);
                    float alphaDeg = alphaRad * 180.0 / M_PI;

                    // Use SC1095 polar data instead of linear approximation
                    float CL = interpolateCL(alphaDeg);
                    
                    float circulationMagnitude = 0.5 * CL * Ublade.norm() * element.chord[0];
                    vector circulation = circulationMagnitude * elementPosition.normalized();                               

                    /*
                    std::cout << "Blade number: " << bladeIndex + 1 << std::endl;
                    std::cout << "Element number: " << elementIndex + 1 << std::endl;                    
                    std::cout << "Circulation magnitude: " << circulationMagnitude << std::endl;
                    std::cout << "Circulation vector: " << circulation.transpose() << std::endl;                    
                    */
                   
                    blade.elements[elementIndex].updateCirculation(circulation);

                    states.emplace_back(elementPosition[0], elementPosition[1], elementPosition[2], Ublade[0], Ublade[1], Ublade[2]);
                }
            }
            return states;
        }
    };
}