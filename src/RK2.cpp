#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

// Function definition (dy/dt)
double dydt(double t, double y){  
  return 30 * cos(t - y) + 3 * sin(t);
}

int main() {
  // Define time and time step size
  double h = 0.1;
  double start = 0.0;
  double end = 1.0;

  // Dynamic array (vector)
  vector<double> t;
  for (double i = start; i <= end; i += h) {
    t.push_back(i);
  }

  for (double print : t) {
    cout << print << " ";
  }

  // Initial value
  vector<double> y(t.size());
  y[0] = 1;

  double F1 = 0;
  double F2 = 0;

  // Runge-Kutta 2
  for (double i = 0; i < 11; i++) {
    F1 = h * dydt(t[i], y[i]);
    F2 = h * dydt(t[i] + h, y[i] + F1);
    y[i + 1] = y[i] + 0.5 * (F1 + F2);
  }

  for (double print : y) {
    cout << "\n"
         << print << " ";
  }

  // Plot
  plt::plot(t,y);
  plt::xlabel("t");
  plt::ylabel("y");
  plt::grid(true);
  plt::show();
}