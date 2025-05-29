#pragma once
#include <cmath>
#include "constants.hpp"
#include <Eigen/Geometry>
#include <Eigen/Core>

using vector = Eigen::Vector3f;

namespace BET {
    struct BladeElement
    {
        float radius; // [ft] radial position
        vector chord;  // [ft] chord length
        float twist;  // [deg] twist angle
        mutable vector circulation; // [ft^2/s] circulation
        mutable vector position; // [ft] position
        mutable vector Ublade; // [ft/s] blade velocity

        // Constructor with optional circulation
        BladeElement(float r, vector c, float t, vector circ = vector::Zero(), vector U = vector::Zero())
            : radius(r), chord(c), twist(t), circulation(circ), position(r, 0.0, 0.0), Ublade(U) {}

        // Method to update circulation
        void updateCirculation(const vector& newCirculation) const {
            circulation = newCirculation;
        }

        void updatePosition(const vector& newPosition) const {
            position = newPosition;
        }

        void updateUblade(const vector& newUblade) const {
                Ublade = newUblade;
            }
    };
} 