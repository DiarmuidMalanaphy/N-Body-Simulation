
#include <stdio.h>
#include <iostream>


#include "Bound.h"

// Default constructor implementation
Bound::Bound() : lower(0.0f), upper(0.0f), mid(0.0f) {}

// Parameterized constructor implementation
Bound::Bound(float lower, float upper)
    : lower(lower), upper(upper) {
    recalibrateMid();
}

// Recalculates the midpoint
void Bound::recalibrateMid() {
    this->mid = (lower + upper) / 2;
}

// Prints the bound values
void Bound::print() {
    std::cout << "Lower : " << lower << ", Middle : " << mid << " , Upper : " << upper << std::endl;
}




