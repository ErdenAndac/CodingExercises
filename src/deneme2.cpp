#include <vector>

class Blade {
private:
    double blade_radius;
    double cutoutradius;
public:
    double getBladeRadius() const { return blade_radius; }
    double getCutoutRadius() const { return cutoutradius; }
};

class BladeElement {
private:
    double bladeelement_radius;  // Single double, not a vector
public:
    double getBladeElementRadius() const { return bladeelement_radius; }
    void setBladeElementRadius(double value) { bladeelement_radius = value; }
};

void calculateBladeElementRadius(std::vector<BladeElement>& elements, std::vector<Blade>& blades, int numElem) {
    double dr = (blades[0].getBladeRadius() - blades[0].getCutoutRadius()) / numElem;
    for (int i = 0; i < numElem; i++) {
        elements[i].setBladeElementRadius(blades[0].getCutoutRadius() + (dr * (i + 1) - (dr / 2.0)));
    }
}