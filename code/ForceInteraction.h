#ifndef FORCEINTERACTION_H
#define FORCEINTERACTION_H

#include <glm/glm.hpp>

// Struct to represent force interactions between particles
typedef struct ForceInteraction {
    glm::vec2 forceVector;      // The force vector applied
    glm::vec2 directionalVector; // The directional vector between particles
    bool interacts;              // Whether the particles interact
} ForceInteraction;

#endif