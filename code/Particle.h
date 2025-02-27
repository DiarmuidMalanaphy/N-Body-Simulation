#ifndef PARTICLE_H
#define PARTICLE_H
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glad/glad.h>
#include "shader_s.h"



#include "Bound.h"
#include "SphereRender.h"
#include "ForceInteraction.h"

extern const float GLOBAL_DT;

class Particle {
public:
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 oldAcceleration;
    float mass;
    SphereRenderData sphereVAOData;
    unsigned int particleID;
    Bound xBounds;
    Bound yBounds;
    Shader* shader;

    // Constructors
    Particle(glm::vec3 position, Bound xBounds, Bound yBounds, float mass, 
             unsigned int particleID, SphereRenderData sphereVAOData, Shader* shader);
    Particle();

    // Member functions
    void printPosition();
    ForceInteraction getForceInteraction(const Particle& terminalParticle) const;
    void updateParticle(glm::vec2 force);
    void clampPosition(glm::vec3 position);
    double getScalarAngle(Particle terminalParticle);
    void drawParticle();
    void updateAndDrawParticle(glm::vec2 movement);
};

#endif