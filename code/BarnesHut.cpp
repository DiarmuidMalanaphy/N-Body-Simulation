
#include <vector>
#include <memory>
#include "Bound.h"
#include <glm/glm.hpp>
#include <unordered_map>

#include "BarnesHut.h"
#include "Particle.h"
#include "Node.h"


BarnesHut::BarnesHut(std::unordered_map<unsigned int, std::shared_ptr<Particle>>& particles, Bound xBounds, Bound yBounds, float theta)
    : theta(theta), xBounds(xBounds), yBounds(yBounds) {
    // Transfer ownership of particles using std::move
    for (auto& particlePair : particles) {
        this->particles[particlePair.first] = std::move(particlePair.second);
    }
    calculateRootNode();
}

void BarnesHut::calculateRootNode() {
    this->rootNode = InternalNode(this->xBounds, this->yBounds);
    for (auto& entry : this->particles) {
        Particle* particle = entry.second.get();
        rootNode.addParticle(particle);
    }
}

std::vector<Node*> BarnesHut::getNodes(Particle* particle) {
    std::vector<Node*> nodes;
    getNodesRecursive(&rootNode, particle, nodes);
    return nodes;
}


void BarnesHut::updateAndDrawParticles() {
    std::unordered_map<int, glm::vec2> forces;

    // Calculate forces for each particle
    for (const auto& [id, particle] : particles) {
        forces[id] = glm::vec2(0.0f, 0.0f);
        std::vector<Node*> relevantNodes = getNodes(particle.get());

        for (const auto* node : relevantNodes) {
            if (node->isInternal) {
                const InternalNode* internalNode = static_cast<const InternalNode*>(node);
                glm::vec2 force = calculateForceFromNode(*particle, internalNode);
                forces[id] += force;
            } else {
                const ExternalNode* externalNode = static_cast<const ExternalNode*>(node);
                if (externalNode->getParticle() != particle.get() && externalNode->isFilled()) {
                    ForceInteraction force = particle->getForceInteraction(*externalNode->particle);
                    forces[id] += force.forceVector;
                }
            }
        }
    }

    // Update and draw particles
    for (auto& [id, particle] : particles) {
        glm::vec2 movement = forces[id];
        particle->updateAndDrawParticle(movement);
    }

    // Recalculate the root node after updating particles
    calculateRootNode();
}

// Recursively get nodes for a particle
void BarnesHut::getNodesRecursive(Node* node, Particle* particle, std::vector<Node*>& nodes) {
    if (node == nullptr) return;

    if (!node->isInternal) {
        ExternalNode* externalNode = static_cast<ExternalNode*>(node);
        if (externalNode->particle != particle && externalNode->filled) {
            nodes.push_back(node);
        }
        return;
    }

    InternalNode* internalNode = static_cast<InternalNode*>(node);
    float s = (std::max)(node->xBounds.upper - node->xBounds.lower,
                         node->yBounds.upper - node->yBounds.lower);
    float d = glm::distance(internalNode->centreOfMass, particle->position);

    if (s / d < theta) {
        // Use the internal node as an approximation
        nodes.push_back(node);
        return;
    }

    // Traverse deeper into the tree
    for (const auto& subNode : internalNode->subNodes) {
        getNodesRecursive(subNode.get(), particle, nodes);
    }
}

// Calculate force from an internal node
glm::vec2 BarnesHut::calculateForceFromNode(const Particle& particle, const InternalNode* node) const {
    glm::vec2 direction = node->centreOfMass - particle.position;
    float distance = glm::length(direction);
    float forceMagnitude = (particle.mass * node->mass) / (distance * distance);
    return glm::normalize(direction) * forceMagnitude;
}

// class BarnesHut 
// {
//     InternalNode rootNode;
//     float theta;
//     Bound xBounds;
//     Bound yBounds;
//     std::unordered_map<unsigned int, std::shared_ptr<Particle>> particles;
//     public:

//     BarnesHut(std::unordered_map<unsigned int, std::shared_ptr<Particle>>& particles, Bound xBounds, Bound yBounds, float theta = 0.7f)
//         : theta(theta), xBounds(xBounds), yBounds(yBounds), particles(particles)
//     {
//         for (auto& particlePair : particles) {
//             this->particles[particlePair.first] = std::move(particlePair.second); // Transfer ownership using std::move
//         }
//         calculateRootNode();
//     }

//     void calculateRootNode()

//     {

//         this->rootNode = InternalNode(this->xBounds, this->yBounds);
//         for (auto& entry : this->particles) 
//         {
//             Particle* particle = entry.second.get();
//             rootNode.addParticle(particle); 
//         }
//     }

//     std::vector<Node*> getNodes(Particle* particle)
//     {
//         std::vector<Node*> nodes;
//         getNodesRecursive(&rootNode, particle, nodes);
//         return nodes;
//     }

//     void updateAndDrawParticles()
//     {
//         std::unordered_map<int, glm::vec2> forces;
//         for (const auto& [id, particle] : particles)
//         {
//             forces[id] = glm::vec2(0.0f, 0.0f);
//             std::vector<Node*> relevantNodes = getNodes(particle.get());
//             for (const auto* node : relevantNodes)
//             {
//                 if (node->isInternal)
//                 {
//                     const InternalNode* internalNode = static_cast<const InternalNode*>(node);
//                     glm::vec2 force = calculateForceFromNode(*particle, internalNode);
//                     forces[id] += force;
//                 }
//                 else
//                 {
//                     const ExternalNode* externalNode = static_cast<const ExternalNode*>(node);
//                     if (externalNode->getParticle() != particle.get() && externalNode->isFilled())
//                     {
//                         ForceInteraction force = particle->getForceInteraction(*externalNode->particle);
//                         forces[id] += force.forceVector;
//                     }
//                 }
//             }
//         }

//         // Update and draw particles
//         for (auto& [id, particle] : particles)
//         {
//             glm::vec2 movement = forces[id];
//             particle->updateAndDrawParticle(movement);
//         }

//         calculateRootNode();
//     }

    

// private:
    
//     void getNodesRecursive(Node* node, Particle* particle, std::vector<Node*>& nodes)
//     {
//         if (node == nullptr) return;

//         if (!node->isInternal)
//         {
//             ExternalNode* externalNode = static_cast<ExternalNode*>(node);
//             if (externalNode->particle != particle && externalNode->filled)
//             {
//                 nodes.push_back(node);
//             }
//             return;
//         }

//         InternalNode* internalNode = static_cast<InternalNode*>(node);
//         float s = (std::max)(node->xBounds.upper - node->xBounds.lower, 
//                            node->yBounds.upper - node->yBounds.lower);
//         float d = glm::distance(internalNode->centreOfMass, particle->position);

//         if (s / d < theta)
//         {
//             // Use the internal node as an estimation as it is far enough away
//             nodes.push_back(node);
//             return;
//         }
//         // We need to go deeper into the tree
//         for (const auto& subNode : internalNode->subNodes)
//         {
//             getNodesRecursive(subNode.get(), particle, nodes);
//         }
//     }

//     glm::vec2 calculateForceFromNode(const Particle& particle, const InternalNode* node) const
//     {
//         glm::vec2 direction = node->centreOfMass - particle.position;
//         float distance = glm::length(direction);
//         float forceMagnitude = (particle.mass * node->mass) / (distance * distance);
//         return glm::normalize(direction) * forceMagnitude;
//     }

// };

    