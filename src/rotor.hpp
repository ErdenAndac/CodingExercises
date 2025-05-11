#pragma once
#include "blade.hpp"
#include "constants.hpp"
#include <vector>
#include <cmath>
#include <tuple>
#include <array>
#include <iostream>

namespace BET {
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
        Rotor(int nBlades, double tipRadius, double rootRadius, int numElements, double chord, double twistRate, double rpm, double MTOW)
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
                double initialAngle = 90.0 * M_PI / 180.0 + i * angleSpacing;
                blades.emplace_back(tipRadius, rootRadius, numElements, chord, twistRate, initialAngle);
            }
        }

        // Transformation matrix for rotation around z-axis
        std::array<std::array<double, 3>, 3> getT1(double psi) const {
            return {{
                {cos(psi), sin(psi), 0},
                {-sin(psi), cos(psi), 0},
                {0, 0, 1}
            }};
        }

        // Transformation matrix for rotation around y-axis
        std::array<std::array<double, 3>, 3> getT2(double theta) const {
            return {{
                {cos(theta), 0, -sin(theta)},
                {0, 1, 0},
                {sin(theta), 0, cos(theta)}
            }};
        }

        // Transformation matrix for rotation around x-axis
        std::array<std::array<double, 3>, 3> getT3(double beta) const {
            return {{
                {1, 0, 0},
                {0, cos(beta), sin(beta)},
                {0, -sin(beta), cos(beta)}
            }};
        }

        // Matrix multiplication helper
        std::array<double, 3> matrixMultiply(const std::array<std::array<double, 3>, 3>& mat, const std::array<double, 3>& vec) const {
            std::array<double, 3> result;
            for (int i = 0; i < 3; ++i) {
                result[i] = 0;
                for (int j = 0; j < 3; ++j) {
                    result[i] += mat[i][j] * vec[j];
                }
            }
            return result;
        }

        // Get blade element positions and velocities for a given time step
        std::vector<std::tuple<double, double, double, double, double, double, double, double, double>> getBladeElementStates(double time) const {
            std::vector<std::tuple<double, double, double, double, double, double, double, double, double>> states;
            
            for (const auto& blade : blades) {
                double psi = -((blade.initialAngle + angularVelocity * time) - blade.initialAngle);
                double theta = -(blade.elements[0].twist);
                double beta = -(blade.flappingAngle);              
                
                // Get transformation matrices
                auto T1 = getT1(psi);
                auto T2 = getT2(theta);
                auto T3 = getT3(beta);
                
                for (const auto& element : blade.elements) {
                    // Calculate position in body frame
                    // Azimuth rotation
                    double x = element.radius * cos(blade.initialAngle + angularVelocity * time);
                    double y = element.radius * sin(blade.initialAngle + angularVelocity * time) * cos(blade.flappingAngle);
                    double z = element.radius * sin(blade.flappingAngle);
                   
                    
                    // Calculate velocities in blade frame
                    std::array<double, 3> Ubody = {0.0, 0.0, 0.0};                      
                    std::array<double, 3> Ublade = matrixMultiply(T1, matrixMultiply(T2, matrixMultiply(T3, Ubody)));
                    
                    // Calculate local velocity components [UR, UT, UP]
                    double UR = 0.0;
                    double UT = element.radius * angularVelocity;
                    double UP = v_i;

                    // std::cout << "twist: " << element.twist << std::endl;
                    
                    Ublade[0] = Ublade[0] + UT * cos(element.twist) + UP * sin(element.twist);  // x
                    Ublade[1] = Ublade[1]; // y
                    Ublade[2] = Ublade[2] + UT * sin(element.twist) - UP * cos(element.twist); //z                                                      // z 
                    
                    
                    states.emplace_back(x, y, z, UR, UT, UP, Ublade[0], Ublade[1], Ublade[2]);
                }
            }
            return states;
        }

        // Calculate circulation at each blade element
        std::vector<std::tuple<double, double>> calculateCirculation() const {
            std::vector<std::tuple<double, double>> circulationValues;
            auto states = getBladeElementStates(0.0);  // Get states at current time
            
            for (const auto& [x, y, z, UR, UT, UP, Ublade_x, Ublade_y, Ublade_z] : states) {                
                double phirad = atan(UP / UT);
                double phideg = phirad * 180.0 / M_PI;

                double alpharad = (blades[0].elements[0].twist - phirad);
                double alphadeg = alpharad * 180.0 / M_PI;         

                // std::cout << "alphadeg: " << alphadeg << std::endl;

                double CL = blades[0].CL_alpha * alpharad;
                double Umag = sqrt(Ublade_x*Ublade_x + Ublade_y*Ublade_y + Ublade_z*Ublade_z);
                double boundCirc = 0.5 * CL * Umag * blades[0].elements[0].chord;

                // std::cout << "Circulation: " << boundCirc << std::endl;
                
                circulationValues.emplace_back(std::sqrt(x*x + y*y + z*z), boundCirc);
            }
            return circulationValues;
        }


    };
} 