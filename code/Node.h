#ifndef NODE_H
#define NODE_H

#include "Bound.h"
#include "Particle.h"
#include <glm/glm.hpp>


#include <array>
#include <memory>

class Node {
public:
    Bound xBounds;
    Bound yBounds;
    int depth;
    bool filled;
    bool isInternal;

    Node::Node(Bound xBounds, Bound yBounds, int depth);
    Node::Node();
    virtual ~Node() = default;

    virtual void addParticle(Particle* particle) = 0;
};

class ExternalNode : public Node {
public:
    Particle* particle;

    ExternalNode(Bound xBounds, Bound yBounds, int depth);
    bool isFilled() const;
    Particle* getParticle() const;

    void addParticle(Particle* particle) override;
};

class InternalNode : public Node {
public:
    float mass;
    std::array<std::unique_ptr<Node>, 4> subNodes;
    glm::vec3 centreOfMass;

    InternalNode();
    InternalNode(Bound xBounds, Bound yBounds);
    InternalNode(ExternalNode* externalNode);

    void addParticle(Particle* particle) override;
    int getBoundIndex(Particle* particle);
    void recalibrateCentreOfMass(Particle* particle);
    void initialiseSubNodes();
};

#endif // NODE_H