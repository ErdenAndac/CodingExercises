#pragma once
#include <cmath>
#include "constants.hpp"

namespace BET {
    struct BladeElement
    {
        double radius; // [ft] radial position
        double chord;  // [ft] chord length
        double twist;  // [deg] twist angle

        BladeElement(double r, double c, double t)
            : radius(r), chord(c), twist(t) {}
    };
} 