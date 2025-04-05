#pragma once
#include "blade.hpp"
#include <vector>
#include <cmath>
#include <tuple>

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
        Rotor(int nBlades, double tipRadius, double rootRadius, int numElements, double chord, double twistRateInput, double rpm, double MTOW)
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
                blades.emplace_back(tipRadius, rootRadius, numElements, chord, twistRateInput, i * angleSpacing);
            }
        }

        // Get blade element positions for a given time step
        std::vector<std::tuple<double, double, double>> getBladeElementPositions(double time) const
        {
            std::vector<std::tuple<double, double, double>> positions;
            for (const auto &blade : blades)
            {
                double theta = fmod(blade.initialAngle + angularVelocity * time, 2.0 * M_PI);
                for (const auto &element : blade.elements)
                {
                    double x = element.radius * cos(theta);
                    double y = element.radius * sin(theta);
                    double z = 0.0;
                    positions.emplace_back(x, y, z);
                }
            }
            return positions;
        }

        // Calculate thrust at each blade element
        std::vector<double> calculateThrust(double V_c) const
        {
            std::vector<double> thrustValues;
            for (const auto &blade : blades)
            {
                for (const auto &element : blade.elements)
                {
                    // Velocity components
                    // Use momentum theory to calculate induced velocity
                    double U_P = V_c + v_i;                              // [ft/s] out-of-plane velocity component
                    double U_T = angularVelocity * element.radius;       // [ft/s] in-plane velocity component
                    double phirad = atan(U_P / U_T);                     // [rad] inflow angle
                    double phideg = phirad * 180.0 / M_PI;               // [deg] inflow angle
                    double lambda = phirad * (element.radius / blade.R); // [rad] inflow ratio
                    double lambdadeg = lambda * 180.0 / M_PI;            // [deg] inflow ratio
                    double pitchAngle;
                    if (blade.twistRate == 0.0)
                    {
                        pitchAngle = (6.0 * CT_mt) / (solidity * blade.CL_alpha) + (3.0 / 2.0) * sqrt(CT_mt / 2.0);
                    }
                    else
                    {
                        pitchAngle = ((6.0 * CT_mt) / (solidity * blade.CL_alpha)) - (3.0 * blade.twistRate / 4.0) + (3.0 * lambdadeg * 2.0);
                    }
                    double theta = pitchAngle;
                    double CL = blade.CL_alpha * (theta - phideg);
                    double boundCirc = 0.5 * CL * angularVelocity * element.radius * element.chord;
                    thrustValues.emplace_back(boundCirc);
                }
            }
            return thrustValues;
        }
    };
} 