#pragma once
#include "bladeelement.hpp"
#include "constants.hpp"
#include <vector>
#include <cmath>

namespace BET {
    struct Blade
    {
        std::vector<BladeElement> elements; // list of blade elements
        double initialAngle;                // [rad] initial angular position
        double R;
        double CL_alpha;
        double CD_0;
        double pitchAngle;
        double flappingAngle;
        double twistRate;       

        // Constructor with Twist Distribution
        Blade(double tipRadius, double rootRadius, int numElements, double chord, double twistRateInput, double angle)
            : initialAngle(angle),
              R(tipRadius),
              CL_alpha((2 * M_PI)),
              CD_0(0.0),
              pitchAngle(5.0 * M_PI / 180.0),
              flappingAngle(0.0 * M_PI / 180.0),
              twistRate(twistRateInput * M_PI / 180.0)
        {
            double dr = (tipRadius - rootRadius) / numElements; // Blade element spacing
            for (int i = 0; i < numElements; ++i)
            {
                double r = rootRadius + i * dr + dr / 2.0; // Midpoint of blade element
                double twist = pitchAngle + twistRate * (r / R);
                elements.emplace_back(r, chord, twist);
            }
        }
    };
} 