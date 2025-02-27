#ifndef BOUND_H
#define BOUND_H

#include <iostream>

class Bound {
public:
    float lower; // Lower bound
    float upper; // Upper bound
    float mid;   // Midpoint

    // Default constructor
    Bound();

    // Parameterized constructor
    Bound(float lower, float upper);

    // Recalculates the midpoint
    void recalibrateMid();

    // Prints the bound values
    void print();
};

#endif // BOUND_H