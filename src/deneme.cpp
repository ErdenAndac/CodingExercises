#define _USE_MATH_DEFINES // for C++
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <filesystem> // C++17 for filesystem operations
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace Vortex_RT {
using Real   = float;
using Vector = Eigen::Vector3f;
} // namespace Vortex_RT

namespace Vortex_RT {
struct Particle {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vector position;
  Vector velocity;
  Vector alpha;
  float  radius{};

  // Ctor
  Particle() = delete;
  Particle(Vector pos, Vector alph, float rad)
      : position(pos), velocity(Vector::Zero()), alpha(alph), radius(rad) {}
};

struct BladeElement {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // Rotor Settings
  Vector BEposition;
  Vector BEvelocity;
  Vector BEradius;
  Vector theta0;
  Vector theta;
  double nb      = 2;     // [-] number of blades
  double rr      = 12;    // [ft] rotor radius
  double r0      = 0.48;  // [ft] blade root cut-out
  double twist   = 0.0;   // [degree] blade twist
  double chord   = 0.6;   // [ft] chord length
  double avrpm   = 450.0; // [rpm] angular velocity
  double CLalpha = 0.1320;
  double CD0     = 0.0056;

  // Other Constants
  double rho  = 2.3769e-3; // [lb/ft^3] density
  double g    = 32.174;    // [ft/s^2] gravitational acceleration
  double MTOW = 1390;      // [lb] maximum take-off weight
  double Vc   = 0;         // [ft/s] forward velocity

  // Calculations
  double avrad  = avrpm * 2 * (M_PI / 60.0);                      // [rad/s] angular velocity
  double A      = (M_PI * (rr * rr));                             // [ft^2] disc area
  double Ablade = nb * (rr - r0) * chord;                         // [ft^2] total blade area
  double sigma  = Ablade / A;                                     // [-] solidity
  double Vtip   = avrad * rr;                                     // [ft/s] rotor blade tip speed
  double Vi     = sqrt(MTOW / (2.0 * rho * A));                   // [ft/s] induced velocity
  double CTmt   = MTOW / (rho * A * (avrad * rr) * (avrad * rr)); // [-] thrust coeff. from momentum theory

  BladeElement()
      : BEposition(Vector::Zero()),
        BEvelocity(Vector::Zero()),
        BEradius(Vector::Zero()),
        theta0(Vector::Zero()),
        theta(Vector::Zero()) {}
};

} // namespace Vortex_RT

namespace Vortex_RT {
void BladeElementTheory(std::vector<Vortex_RT::BladeElement> &elements, const size_t numElement, const double dt, const double tf)
{
  Vortex_RT::BladeElement elem;
  double                  numStep     = (tf - 0.0) / dt;
  double                  dr          = (elem.rr - elem.r0) / numElement;
  double                  integralCT  = 0.0;
  double                  integralCPi = 0.0;

  for (size_t j = 0; j <= 0; j++) {
    elements.clear();

    double prevCT       = 0.0;
    double prevCPi      = 0.0;
    double prevBEradius = elem.BEradius.y() / elem.rr;

    for (size_t i = 0; i < numElement; i++) {

      elem.BEradius.x() = elem.r0 + (dr * (i + 1) - (dr / 2.0));
      elem.BEradius.y() = elem.BEradius.x() / elem.rr;

      // Update blade element position
      elem.BEposition.x() = elem.BEradius.x() * cos(0 + elem.avrad * dt * j);
      elem.BEposition.y() = elem.BEradius.x() * sin(0 + elem.avrad * dt * j);

      // Update blade element velocity
      elem.BEvelocity.x() = elem.avrad * elem.BEradius.x() * sin(0 + elem.avrad * dt * j);
      elem.BEvelocity.y() = elem.avrad * elem.BEradius.x() * cos(0 + elem.avrad * dt * j);
      elem.BEvelocity.z() = -(elem.Vc + elem.Vi);

      // theta(r) calculation
      elem.theta.x() = (6.0 * elem.CTmt) / (elem.sigma * elem.CLalpha) + ((3.0 / 2.0) * sqrt(elem.CTmt / 2.0));
      elem.theta.y() = elem.theta.x() + elem.BEradius.y() * elem.twist;

      // local inflow ratio
      double phirad = (elem.Vc + elem.Vi) / (elem.avrad * elem.BEradius.x()); // [1/rad] inflow angle
      double phideg = phirad * (180.0 / M_PI);                                // [deg] inflow angle
      double lambda = phirad * elem.BEradius.y();

      // local effective angle of attack
      double alpha = elem.theta.y() - phideg;

      // Cl and cDp0 calculation
      double Cl   = elem.CLalpha * alpha;
      double dCP0 = (elem.sigma / 2.0) * (elem.CD0 * elem.BEradius.y() * elem.BEradius.y() * elem.BEradius.y());

      // incremental thrust coefficient
      double dCT = 0.5 * elem.sigma * Cl * elem.BEradius.y() * elem.BEradius.y();

      // incremental induced power coefficient
      double dCPi = (elem.sigma / 2.0) * (Cl * lambda * elem.BEradius.y() * elem.BEradius.y());

      // thrust coefficient integration
      if (i > 0) {
        double dY = elem.BEradius.y() - prevBEradius;
        integralCT += 0.5 * (prevCT + dCT) * dY;
        integralCPi += 0.5 * (prevCPi + dCPi) * dY;
      }

      // update previous values
      prevCT       = dCT;
      prevCPi      = dCPi;
      prevBEradius = elem.BEradius.y();
      elements.push_back(elem);
    }
  }
  std::cout << "Integral of CT with respect to BEradius.y(): " << integralCT << "\n";
  std::cout << "Integral of CPi with respect to BEradius.y(): " << integralCPi << "\n";
}
} // namespace Vortex_RT

void initRandom(std::vector<Vortex_RT::Particle> &particles, const size_t numPart)
{
  // Define the range for random numbers
  double lower_bound = 0.0;
  double upper_bound = 1.0;

  // Create a random number generator and a distribution
  std::random_device                     rd;              // Seed for the random number generator
  std::mt19937                           generator(rd()); // Mersenne Twister RNG
  std::uniform_real_distribution<double> distribution(lower_bound, upper_bound);

  particles.reserve(numPart);
  Vortex_RT::Particle randomParticle(Vortex_RT::Vector::Zero(), Vortex_RT::Vector::Zero(), 0.5);
  for (size_t i = 0; i < numPart; i++) {
    randomParticle.position = Vortex_RT::Vector(distribution(generator), distribution(generator), distribution(generator));
    randomParticle.alpha    = Vortex_RT::Vector(distribution(generator), distribution(generator), distribution(generator));
    particles.push_back(randomParticle);
  }
}

namespace Vortex_RT {

struct Gaussian {
};
struct LowOrder {
};
struct HighOrder {
};

Real q_sigma_gaussian(const float rho)
{
  // return 1.0 / rho * std::erf(rho / std::sqrt(2.0));
  return (1.0 / (4.0 * M_PI)) * (std::erf(rho / std::sqrt(2.0)) - std::sqrt(2.0 / M_PI) * rho * std::exp(-rho * rho / 2.0));
}

Real q_sigma_loworder(const float rho)
{
  return (std::pow(rho, 3)) / (4.0 * M_PI * std::pow((std::pow(rho, 2) + 1), 1.5));
}

Real q_sigma_highorder(const float rho)
{
  return (std::pow(rho, 3) * (std::pow(rho, 2) + 5.0 / 2.0)) / (4.0 * M_PI * std::pow((std::pow(rho, 2) + 1), 2.5));
}

void resetInducedVelocities(std::vector<Vortex_RT::Particle> &particles)
{
  for (auto &particle : particles) {
    particle.velocity(Vortex_RT::Vector::Zero());
  }
}

template <typename T>
void calculateInducedVelocities(std::vector<Particle> &particles, T kernel)
{
  const auto numParticles = particles.size();

  for (size_t i = 0; i < numParticles; i++) {
    for (size_t j = 0; j < numParticles; j++) {
      if (i == j) {
        continue;
      }
      const auto distanceVector = (particles[i].position - particles[j].position);
      // std::cout << "Distance vector = " << distanceVector << "\n" ;

      const auto distanceScalar = distanceVector.norm();
      // std::cout << "Distance scalar = " << distanceScalar << "\n" ;
      
      // std::cout << "Radius = " << particles[j].position << "\n" ;
      
      if constexpr (std::is_same_v<T, Gaussian>) {
        particles[i].velocity += q_sigma_gaussian(distanceScalar / particles[j].radius) / (distanceScalar * distanceScalar * distanceScalar) * distanceVector.cross(particles[j].alpha);
      }

      
      if constexpr (std::is_same_v<T, LowOrder>) {
        particles[i].velocity += q_sigma_loworder(distanceScalar / particles[j].radius) / (distanceScalar * distanceScalar * distanceScalar) * distanceVector.cross(particles[j].alpha);
      }
      if constexpr (std::is_same_v<T, HighOrder>) {
        particles[i].velocity += q_sigma_highorder(distanceScalar / particles[j].radius) / (distanceScalar * distanceScalar * distanceScalar) * distanceVector.cross(particles[j].alpha);
      }
    }
  }
}

template <typename T>
void integrate(std::vector<Particle> &particles, const float dt, T kernel)
{
  calculateInducedVelocities(particles, kernel);
  for (auto &particle : particles) {
    particle.position += particle.velocity * dt;
  }
  resetInducedVelocities(particles);
}

} // namespace Vortex_RT

// Function to write particles to a VTK file
void writeParticlesToVTK(const std::vector<Vortex_RT::Particle> &particles, const std::string &filename)
{

  // Define the output directory inside the build folder
  std::string outputDir = "vtk_files_2";

  // Check if the directory exists, if not, create it
  if (!std::filesystem::exists(outputDir)) {
    std::filesystem::create_directories(outputDir);
  }

  // Construct the full file path
  std::string fullPath = outputDir + "/" + filename;

  std::ofstream vtkFile(fullPath);
  if (!vtkFile.is_open()) {
    throw std::runtime_error("Could not open file for writing: " + fullPath);
  }

  // Write VTK header
  vtkFile << "# vtk DataFile Version 3.0\n";
  vtkFile << "Particles data\n";
  vtkFile << "ASCII\n";
  vtkFile << "DATASET POLYDATA\n";

  // Write points (positions of particles)
  vtkFile << "POINTS " << particles.size() << " float\n";
  for (const auto &particle : particles) {
    vtkFile << particle.position.x() << " "
            << particle.position.y() << " "
            << particle.position.z() << "\n";
  }

  // Write velocities as point data
  vtkFile << "POINT_DATA " << particles.size() << "\n";
  vtkFile << "VECTORS velocities float\n";
  for (const auto &particle : particles) {
    vtkFile << particle.velocity.x() << " "
            << particle.velocity.y() << " "
            << particle.velocity.z() << "\n";
  }

  vtkFile.close();
  std::cout << "VTK file written to: " << fullPath << std::endl;
}

int main()
{
  std::vector<Vortex_RT::Particle> MyParticles;
  initRandom(MyParticles, 1000);

   for (size_t i = 0; i < 10; i++)
  {

      Vortex_RT::integrate(MyParticles, 1, Vortex_RT::Gaussian{});
      writeParticlesToVTK(MyParticles, "Gaussian_part" + std::to_string(i) + ".vtk");

      // Vortex_RT::integrate(MyParticles, 1, Vortex_RT::LowOrder{});
      // writeParticlesToVTK(MyParticles, "LowOrder_part" + std::to_string(i) + ".vtk");

      // Vortex_RT::integrate(MyParticles, 1, Vortex_RT::HighOrder{});
      // writeParticlesToVTK(MyParticles, "HighOrder_part" + std::to_string(i) + ".vtk");

  } 

  /* for (const auto &particle : MyParticles)
  {
      std::cout << particle.position.transpose() << " " << particle.alpha.transpose() << "\n";
  } */

  // std::vector<Vortex_RT::BladeElement> MyElements;
  // BladeElementTheory(MyElements, 100, 0.01, 0.1);

  /* for (const auto &element : MyElements) {
    std::cout << std::fixed << std::setprecision(4) << element.BEradius.transpose() << "   " << element.BEposition.transpose() << "   "
              << element.BEvelocity.transpose() << "   " << element.theta.transpose() << "\n";
  } */
}