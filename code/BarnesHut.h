#ifndef BARNESHUT_H
#define BARNESHUT_H

#include <vector>
#include <memory>
#include "Bound.h"
#include <glm/glm.hpp>
#include <unordered_map>


#include "Particle.h"
#include "Node.h"


class BarnesHut 
{
    InternalNode rootNode;
    float theta;
    Bound xBounds;
    Bound yBounds;
    std::unordered_map<unsigned int, std::shared_ptr<Particle>> particles;

public:
    BarnesHut(std::unordered_map<unsigned int, std::shared_ptr<Particle>>& particles, Bound xBounds, Bound yBounds, float theta = 0.7f);
    
    void calculateRootNode();
    std::vector<Node*> getNodes(Particle* particle);
    void updateAndDrawParticles();

private:
    void getNodesRecursive(Node* node, Particle* particle, std::vector<Node*>& nodes);
    glm::vec2 calculateForceFromNode(const Particle& particle, const InternalNode* node) const;
};

#endif // BARNESHUT_H
