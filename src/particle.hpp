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
    float q_sigma_gaussian(float rho)
    {
        return (1.0 / (4.0 * M_PI)) * (std::erf(rho / std::sqrt(2.0)) - std::sqrt(2.0 / M_PI) * rho * std::exp(-rho * rho / 2.0));
    }

    float q_sigma_loworder(const float rho)
    {
        return (std::pow(rho, 3)) / (4.0 * M_PI * std::pow((std::pow(rho, 2) + 1), 1.5));
    }

    float q_sigma_highorder(const float rho)
    {
        return (std::pow(rho, 3) * (std::pow(rho, 2) + 5.0 / 2.0)) / (4.0 * M_PI * std::pow((std::pow(rho, 2) + 1), 2.5));
    }

    enum class ParticleType
    {
        BOUND,
        WAKE
    };

    struct Particle
    {
        vector position;
        vector velocity;
        vector alpha; // circulation vector
        ParticleType type;

        Particle(const vector &pos, const vector &vel, const vector &alpha, ParticleType t)
            : position(pos), velocity(vel), alpha(alpha), type(t) {}
    };

    struct ParticleSystem
    {
        std::vector<Particle> boundParticles;
        std::vector<Particle> wakeParticles;

        void generateBoundParticles(const Rotor &rotor)
        {
            boundParticles.clear();
            for (const auto &blade : rotor.blades)
            {
                for (const auto &element : blade.elements)
                {
                    boundParticles.emplace_back(
                        element.position,    // Position from blade element
                        vector::Zero(),      // Initial velocity
                        element.circulation, // Circulation from blade element
                        ParticleType::BOUND);
                }
            }
        }

        void generateWakeParticles(const Rotor &rotor)
        {
            for (const auto &blade : rotor.blades)
            {
                const auto &tipElement = blade.elements.back();
                wakeParticles.emplace_back(
                    tipElement.position,    // Position from blade tip
                    vector::Zero(),         // Initial velocity
                    tipElement.circulation, // Circulation from blade tip
                    ParticleType::WAKE);
            }
        }

        vector calculateInducedVelocity(const vector &targetPos, const vector &sourcePos, const vector &sourceAlpha)
        {
            vector distanceVector = (targetPos - sourcePos);
            float distanceScalar = distanceVector.norm();

            // Skip if particles are too close to avoid singularity
            if (distanceScalar < 1e-6)
                return vector::Zero();

            float radius = 1.0f;
            float rho = distanceScalar / radius;

            // Gaussian kernel
            float q_sigma = q_sigma_gaussian(rho);
            float factor_gaussian = q_sigma / (distanceScalar * distanceScalar * distanceScalar);

            // Low-order kernel
            // float q_sigma_low = q_sigma_loworder(rho);
            // float factor_low = q_sigma_low / (distanceScalar * distanceScalar * distanceScalar);

            // High-order kernel
            // float q_sigma_high = q_sigma_highorder(rho);
            // float factor_high = q_sigma_high / (distanceScalar * distanceScalar * distanceScalar);

            return factor_gaussian * distanceVector.cross(sourceAlpha);
        }

        // Calculate bound-bound interactions
        void updateBoundParticles(const Rotor &rotor)
        {
            for (auto &particle : boundParticles)
            {
                particle.velocity = vector::Zero();

                for (const auto &otherParticle : boundParticles)
                {
                    if (&particle == &otherParticle)
                        continue;

                    particle.velocity += calculateInducedVelocity(
                        particle.position,
                        otherParticle.position,
                        otherParticle.alpha);
                }
            }
        }

        // Calculate wake-wake interactions
        void updateWakeParticles(const Rotor &rotor, float timeStep)
        {
            for (auto &particle : wakeParticles)
            {
                // particle.velocity = vector::Zero();

                for (const auto &otherWakeParticle : wakeParticles)
                {
                   
                    if (&particle == &otherWakeParticle)
                        continue;

                    particle.velocity += calculateInducedVelocity(
                        particle.position,
                        otherWakeParticle.position,
                        otherWakeParticle.alpha);
                }                
                particle.position += (particle.velocity) * timeStep;
            }
        }

        void update(const Rotor &rotor, float timeStep)
        {
            generateBoundParticles(rotor);
            generateWakeParticles(rotor);
            updateBoundParticles(rotor);
            updateWakeParticles(rotor, timeStep);
        }
    };
};