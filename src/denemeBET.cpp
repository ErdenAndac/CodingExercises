#define _USE_MATH_DEFINES // for C++ Math Constants
#include <iostream>
#include <vector>
#include <cmath>

// Blade Element Struct
struct BladeElement {
    double radius;  // Radial position (ft)
    double chord;   // Chord length (ft)
    double twist;   // Twist angle (deg)
    double pitch;   // Blade pitch angle (deg)

    // Constructor
    BladeElement(double r, double c, double t, double p)
        : radius(r), chord(c), twist(t), pitch(p) {}
};

// Blade Struct
struct Blade {
    std::vector<BladeElement> elements; // List of blade elements
    double initialAngle;                // Initial angular position (rad)

    // Constructor: Creates a blade with uniform properties
    Blade(double tipRadius, double rootRadius, int numElements, double chord, double twist, double pitch, double angle)
        : initialAngle(angle) {
        double dr = (tipRadius - rootRadius) / numElements; // Element spacing
        for (int i = 0; i < numElements; ++i) {
            double r = rootRadius + i * dr + dr / 2.0; // Midpoint of element
            elements.emplace_back(r, chord, twist, pitch);
        }
    }
};

// Rotor Struct (Manages Multiple Blades)
struct Rotor {
    std::vector<Blade> blades; // Stores multiple blades
    int numBlades;             // Number of blades
    double angularVelocity;     // Rotation speed (rad/s)

    // Constructor: Creates a rotor with N identical blades
    Rotor(int nBlades, double tipRadius, double rootRadius, int numElements, double chord, double twist, double pitch, double rpm)
        : numBlades(nBlades), angularVelocity(rpm * 2.0 * M_PI / 60.0) {
        double angleSpacing = 2.0 * M_PI / numBlades; // Evenly space blades around rotor
        for (int i = 0; i < numBlades; ++i) {
            blades.emplace_back(tipRadius, rootRadius, numElements, chord, twist, pitch, i * angleSpacing);
        }
    }

    // Compute and print blade positions at a given time
    void printBladePositions(double time) const {
        std::cout << "Blade positions at t = " << time << " sec:\n";
        for (size_t i = 0; i < blades.size(); ++i) {
            double theta = fmod(blades[i].initialAngle + angularVelocity * time, 2.0 * M_PI); // Keep within 0 to 2π
            std::cout << " Blade " << i + 1 << " angle: " << theta * 180.0 / M_PI << " degrees\n";
        }
    }
};

// Example usage
int main() {
    int numBlades = 3;      // 2-blade rotor
    double tipRadius = 12.0;  // ft
    double rootRadius = 0.48; // ft
    int numElements = 100;    // Elements per blade
    double chord = 2.0;       // Chord length (ft)
    double twist = 0.0;       // No twist for now
    double pitch = 0.0;       // No pitch for now
    double rpm = 450.0;       // Rotational speed in RPM

    // Create Rotor
    Rotor myRotor(numBlades, tipRadius, rootRadius, numElements, chord, twist, pitch, rpm);

    // Simulate rotation at different time steps
    double timeStep = 0.01; // Time step in seconds
    for (double t = 0; t <= 0.1; t += timeStep) {
        myRotor.printBladePositions(t);
    }

    return 0;
}
