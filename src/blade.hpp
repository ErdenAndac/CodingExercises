#pragma once
#include "bladeelement.hpp"
#include "constants.hpp"
#include <vector>
#include <cmath>
#include <Eigen/Geometry>
#include <Eigen/Core>

using vector = Eigen::Vector3f;


namespace BET {
    struct Blade
    {
        std::vector<BladeElement> elements; // list of blade elements
        float initialAngle;                // [rad] initial angular position
        float R;
        float CL_alpha;
        float CD_0;
        float pitchAngle;
        float flappingAngle;
        float twistRate;        


        // Constructor with Twist Distribution
        Blade(float tipRadius, float rootRadius, int numElements, float chord, float twistRateInput, float angle)
            : initialAngle(angle),
              R(tipRadius),
              CL_alpha((2 * M_PI)),                     // [1/rad] lift slope
              CD_0(0.0),                                // [-] zero-lift drag coefficient    
              pitchAngle(6.0 * M_PI / 180.0),           // [rad] pitch angle
              flappingAngle(0.0 * M_PI / 180.0),        // [rad] flapping angle
              twistRate(twistRateInput * M_PI / 180.0)  // [rad/ft] twist rate              
        {
            float dr = (tipRadius - rootRadius) / numElements; // Blade element spacing
            for (int i = 0; i < numElements; ++i)
            {
                float r = rootRadius + i * dr + dr / 2.0;          // Midpoint of blade element
                Eigen::Vector3f c = Eigen::Vector3f(chord, 0.0, 0.0);
                float t = pitchAngle + twistRate * (r);    // [rad] twist angle
                elements.emplace_back(r, c, t);
            }            
        }
    };
} 