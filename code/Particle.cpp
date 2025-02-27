


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glad/glad.h>
#include "shader_s.h"



#include "Bound.h"
#include "SphereRender.h"
#include "ForceInteraction.h"
#include "SphereRenderData.h"
#include "Particle.h"


const float GLOBAL_DT = 0.05;

Particle::Particle(glm::vec3 position, Bound xBounds, Bound yBounds, float mass, unsigned int particleID, SphereRenderData sphereVAOData, Shader* shader)
    : position(position), xBounds(xBounds), yBounds(yBounds), mass(mass), particleID(particleID), sphereVAOData(sphereVAOData), shader(shader) {
    velocity = glm::vec3(0.0f, 0.0f, 0.0f);
    oldAcceleration = glm::vec3(0.0f, 0.0f, 0.0f);
}

// Default constructor
Particle::Particle()
    : position(glm::vec3(0.0f, 0.0f, 0.0f)), mass(1.0f), particleID(0), sphereVAOData(), shader(nullptr) {
    velocity = glm::vec3(0.0f, 0.0f, 0.0f);
    oldAcceleration = glm::vec3(0.0f, 0.0f, 0.0f);
}


void Particle::printPosition() {
    std::cout << "(" << position.x << ", " << position.y << ", " << position.z << ")" << std::endl;
}

ForceInteraction Particle::getForceInteraction(const Particle& terminalParticle) const {
    float dx = position.x - terminalParticle.position.x;
    float dy = position.y - terminalParticle.position.y;
    glm::vec2 directionalVector = glm::normalize(glm::vec2(dx, dy));
    float euclideanDistance = dx * dx + dy * dy;
    float dampening_distance = 0.1f;
    float distance = (std::max)(euclideanDistance / 50.0f, dampening_distance);
    float forceMagnitude = (mass * terminalParticle.mass) * 4 / distance;
    bool collides = (abs(dx) < 1.1 && abs(dy) < 1.1);
    return ForceInteraction{ -directionalVector * forceMagnitude, directionalVector, collides };
}

void Particle::updateParticle(glm::vec2 force) {
    float dt = GLOBAL_DT;

    glm::vec3 acceleration = glm::vec3(force.x, force.y, 0.0f) / mass;
    velocity = velocity + ((acceleration + oldAcceleration) * dt) / 2.0f;
    position = position + (velocity * dt);
    clampPosition(position);

    // Update model matrix
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, position);
    shader->setMat4("model", model);
}


void Particle::clampPosition(glm::vec3 position) {
    glm::vec3 newPosition = position;
    bool clamped = false;

    // Clamp x position
    if (newPosition.x < xBounds.lower) {
        newPosition.x = xBounds.upper;
    }
    else if (newPosition.x > xBounds.upper) {
        newPosition.x = xBounds.lower;
    }

    // Clamp y position
    if (newPosition.y < yBounds.lower) {
        newPosition.y = yBounds.upper;
    }
    else if (newPosition.y > yBounds.upper) {
        newPosition.y = yBounds.lower;
    }

    // Update position
    position = newPosition;
}

// Calculate scalar angle between two particles
double Particle::getScalarAngle(Particle terminalParticle) {
    double dotProduct = glm::dot(this->position, terminalParticle.position);
    double magnitudeVec1 = glm::length(this->position);
    double magnitudeVec2 = glm::length(terminalParticle.position);
    float angleInRadians = std::acos(dotProduct / (magnitudeVec1 * magnitudeVec2));
    double angleInDegrees = glm::degrees(angleInRadians);
    return angleInDegrees;
}

// Draw the particle
void Particle::drawParticle() {
    glDrawElements(GL_TRIANGLE_STRIP, sphereVAOData.indexCount, GL_UNSIGNED_INT, 0);
}

// Update and draw the particle
void Particle::updateAndDrawParticle(glm::vec2 movement) {
    updateParticle(movement);
    drawParticle();
}

// class Particle {
//     public:
//     glm::vec3 position;
//     glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
//     glm::vec3 oldAcceleration = glm::vec3(0.0f, 0.0f, 0.0f);
//     float mass;
//     SphereRenderData sphereVAOData;
//     unsigned int particleID;

//     Bound xBounds;
//     Bound yBounds;
//     Shader* shader;  // Change to pointer

//     Particle(glm::vec3 position, Bound xBounds, Bound yBounds, float mass, unsigned int particleID, SphereRenderData sphereVAOData, Shader* shader)
//         : position(position), xBounds(xBounds), yBounds(yBounds), mass(mass), particleID(particleID), sphereVAOData(sphereVAOData), shader(shader)
//     {}

//     Particle()
//         : position(glm::vec3(0.0f, 0.0f, 0.0f)),  // Default position at origin
//           mass(1.0f),
//           particleID(0),
//           sphereVAOData(),                        // Default initialization for sphereVAOData
//           shader(nullptr)
//     {}


//     void printPosition() 
//     {
//         std::cout << "(" << position.x << ", " << position.y << ", " << position.z << ")" << std::endl;
//     }

//     ForceInteraction getForceInteraction(const Particle& terminalParticle) const 
//     {
//         float dx = position.x - terminalParticle.position.x;
//         float dy = position.y - terminalParticle.position.y;
//         glm::vec2 directionalVector = glm::normalize(glm::vec2(dx, dy));
//         float euclideanDistance = dx * dx + dy * dy;
//         float dampening_distance = 0.1f;
//         float distance = (std::max)(euclideanDistance/50.0f, dampening_distance);
//         float forceMagnitude = (mass * terminalParticle.mass)*4 / distance;
//         bool collides = (abs(dx) < 1.1 && abs(dy) < 1.1);
//         return ForceInteraction{-directionalVector * forceMagnitude, directionalVector, collides}; 
//     }


//     void updateParticle(glm::vec2 force)
//     {
//         float dt = GLOBAL_DT;

//         glm::vec3 acceleration = glm::vec3(force.x, force.y, 0.0f) / mass;
//         velocity = velocity + ((acceleration + oldAcceleration) * dt) / 2.0f; 
//         position = position + (velocity * dt);
//         clampPosition(position);

        
//         // Update model matrix
//         glm::mat4 model = glm::mat4(1.0f);
//         model = glm::translate(model, position);
//         shader->setMat4("model", model);
//     }

//     void clampPosition(glm::vec3 position)
//     {
//         glm::vec3 newPosition = position;
//         bool clamped = false;
//         // Clamp x position
//         if (newPosition.x < xBounds.lower)
//         {
//             newPosition.x = xBounds.upper;
//         }
//         else if (newPosition.x > xBounds.upper)
//         {
//             newPosition.x = xBounds.lower;
//         }

//         if (newPosition.y < yBounds.lower)
//         {

            
//             newPosition.y = yBounds.upper;
//         }
//         else if (newPosition.y > yBounds.upper)
//         {

//             newPosition.y = yBounds.lower;
//         }

//         // Update position
//         position = newPosition;

//     }

//     double getScalarAngle(Particle terminalParticle){
//         double dotProduct = glm::dot(this ->position, terminalParticle.position);
//         double magnitudeVec1 = glm::length(this->position);
//         double magnitudeVec2 = glm::length(this ->position);
//         float angleInRadians = std::acos(dotProduct / (magnitudeVec1 * magnitudeVec2)) ;
//         double angleInDegrees = glm::degrees(angleInRadians);
//         return angleInDegrees;
//     }
    
    
    
//     void drawParticle() 
//     {
//         glDrawElements(GL_TRIANGLE_STRIP, sphereVAOData.indexCount, GL_UNSIGNED_INT, 0);
//     }



//     void updateAndDrawParticle(glm::vec2 movement) 
//     {
//         updateParticle(movement);
//         drawParticle();
//     }
// };

