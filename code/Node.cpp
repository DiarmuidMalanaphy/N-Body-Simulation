#include "Bound.h"
#include "Particle.h"
#include <glm/glm.hpp>


#include <array>
#include <memory>


#include "Node.h"

// Node constructor
Node::Node(Bound xBounds, Bound yBounds, int depth)
    : xBounds(xBounds), yBounds(yBounds), depth(depth) {
    xBounds.recalibrateMid();
    yBounds.recalibrateMid();
}

// Default Node constructor
Node::Node() : xBounds(), yBounds() {}

// ExternalNode constructor
ExternalNode::ExternalNode(Bound xBounds, Bound yBounds, int depth)
    : Node(xBounds, yBounds, depth), particle(nullptr) {
    this->isInternal = false;
    this->filled = false;
}

// Check if the ExternalNode is filled
bool ExternalNode::isFilled() const {
    return filled;
}

// Get the particle in the ExternalNode
Particle* ExternalNode::getParticle() const {
    return particle;
}

// Add a particle to the ExternalNode
void ExternalNode::addParticle(Particle* particle) {
    this->filled = true;
    this->particle = particle;
}

// InternalNode default constructor
InternalNode::InternalNode()
    : Node(), centreOfMass(0.0f, 0.0f, 0.0f), mass(0) {}

// InternalNode constructor with bounds
InternalNode::InternalNode(Bound xBounds, Bound yBounds)
    : Node(xBounds, yBounds, 0), centreOfMass(0.0f, 0.0f, 0.0f), mass(0) {
    initialiseSubNodes();
}

// InternalNode constructor from an ExternalNode
InternalNode::InternalNode(ExternalNode* externalNode)
    : Node(externalNode->xBounds, externalNode->yBounds, (externalNode->depth) + 1),
      centreOfMass(0.0f, 0.0f, 0.0f), mass(0) {
    initialiseSubNodes();
    addParticle(externalNode->particle); // Add the particle from the external node
}

// Add a particle to the InternalNode
void InternalNode::addParticle(Particle* particle) {
    // Recalibrate the centre of mass
    recalibrateCentreOfMass(particle);

    // Determine which sub-node to add the particle to
    int boundIndex = getBoundIndex(particle);
    Node* targetNode = subNodes[boundIndex].get();

    // If the target node is an external node and already filled, convert it to an internal node
    if (!targetNode->isInternal && targetNode->filled) {
        ExternalNode* node = static_cast<ExternalNode*>(targetNode);
        subNodes[boundIndex] = std::make_unique<InternalNode>(node);
        targetNode = subNodes[boundIndex].get();
    }

    // Add the particle to the target node
    targetNode->addParticle(particle);
}

// Determine which sub-node a particle belongs to
int InternalNode::getBoundIndex(Particle* particle) {
    bool top = (particle->position.y > yBounds.mid);
    bool right = (particle->position.x > xBounds.mid);

    if (top) {
        return right ? 1 : 0; // Top-right (1) or top-left (0)
    } else {
        return right ? 3 : 2; // Bottom-right (3) or bottom-left (2)
    }
}

// Recalibrate the centre of mass when adding a particle
void InternalNode::recalibrateCentreOfMass(Particle* particle) {
    glm::vec3 weightedSum = centreOfMass * mass;
    weightedSum += particle->position * particle->mass;
    mass += particle->mass;
    centreOfMass = weightedSum / mass;
}

// Initialize the sub-nodes of the InternalNode
void InternalNode::initialiseSubNodes() {
    Bound low_x_bound = Bound{xBounds.lower, xBounds.mid};
    Bound high_x_bound = Bound{xBounds.mid, xBounds.upper};
    Bound low_y_bound = Bound{yBounds.lower, yBounds.mid};
    Bound high_y_bound = Bound{yBounds.mid, yBounds.upper};
    int new_depth = depth + 1;

    // Initialize the four sub-nodes
    subNodes[0] = std::make_unique<ExternalNode>(low_x_bound, high_y_bound, new_depth); // Top-left
    subNodes[1] = std::make_unique<ExternalNode>(high_x_bound, high_y_bound, new_depth); // Top-right
    subNodes[2] = std::make_unique<ExternalNode>(low_x_bound, low_y_bound, new_depth); // Bottom-left
    subNodes[3] = std::make_unique<ExternalNode>(high_x_bound, low_y_bound, new_depth); // Bottom-right
}

// class Node 
// {
// public:
//     Bound xBounds;
//     Bound yBounds;
//     int depth;
//     bool filled = false;
//     bool isInternal = true;
//     Node(Bound xBounds, Bound yBounds, int depth)
//         : xBounds(xBounds), yBounds(yBounds), depth(depth)
//     {
//         xBounds.recalibrateMid();
//         yBounds.recalibrateMid();
//     }
//     Node() 
//         : xBounds(), yBounds() // Initializes both bounds using their default constructor
//     {
//     }
//     virtual void addParticle(Particle* particle) = 0;
// };

// class ExternalNode : public Node 
// {

//     public:
//     bool isFilled() const { return filled; }
//     Particle* getParticle() const { return particle; }
//     Particle* particle;
//     ExternalNode(Bound xBounds, Bound yBounds, int depth) 
//         : Node{xBounds, yBounds, depth}, particle(nullptr)
//     {
//         this->isInternal = false;
//         this->filled = false;
//     }

//     void addParticle(Particle* particle)
//     {
//         this->filled = true;
//         this->particle = particle;
//     }
// };

// class InternalNode : public Node
// {
//     // 0 is top-left
//     // 1 is top-right
//     // 2 is bottom-left
//     // 3 is bottom-right
    
//     public:

//     float mass;
//     std::array<std::unique_ptr<Node>, 4> subNodes;
//     glm::vec3 centreOfMass; 

//     void addParticle(Particle* particle)
//     {
//         //Every time a particle passes through here we must recalibrate the centre of mass 
//         recalibrateCentreOfMass(particle);
//         int boundIndex = getBoundIndex(particle);
//         Node* targetNode = subNodes[boundIndex].get();
        
//         if (targetNode->isInternal == false && targetNode->filled == true) 
//         {
//             ExternalNode* node = static_cast<ExternalNode*>(targetNode); 
//             subNodes[boundIndex] = std::make_unique<InternalNode>(node);
//             targetNode = subNodes[boundIndex].get();
//         }
//         targetNode->addParticle(particle);
//     }

//     int getBoundIndex(Particle *particle)
//     {
//         bool top = (particle->position.y > yBounds.mid);
//         bool right = (particle->position.x > xBounds.mid);
//         if (top)
//         {
//             if (!right) 
//             {
//                 //Top left
//                 return 0;
//             }
//             //Top-right
//             return 1;
//         }
//         if (!right) 
//             {
//                 //Bottom left
//                 return 2;
//             }
//         //Bottom-right
//         return 3;
//     }

//     void recalibrateCentreOfMass(Particle *particle) 
//     {
//         glm::vec3 weightedSum = centreOfMass * mass;
//         weightedSum += particle->position * particle->mass; 
//         mass += particle->mass;
//         centreOfMass = weightedSum / mass;
//     }


//     InternalNode()
//         : Node(),
//           centreOfMass(0.0f, 0.0f, 0.0f),
//           mass(0) 

//     {}

//     InternalNode(Bound xBounds, Bound yBounds)
//         : Node{xBounds, yBounds, 0}, 
//           centreOfMass(0.0f, 0.0f, 0.0f), 
//           mass(0)
//     {
//         initialiseSubNodes();

//     }

//     InternalNode(ExternalNode* externalNode)
//         : Node{externalNode->xBounds, externalNode->yBounds, (externalNode->depth)+1}, 
//           centreOfMass(0.0f, 0.0f, 0.0f), 
//           mass(0)
//     {
//         initialiseSubNodes();
//         addParticle(externalNode->particle);  // Assuming `addParticle` is implemented correctly
//     }
    
//     void initialiseSubNodes()
//     {
//         Bound low_x_bound = Bound
//                         {xBounds.lower, xBounds.mid};

//         Bound high_x_bound = Bound
//                         {xBounds.mid, xBounds.upper};

//         Bound low_y_bound = Bound
//                         {yBounds.lower, yBounds.mid};

//         Bound high_y_bound = Bound
//                         {yBounds.mid, yBounds.upper};
//         int new_depth = depth + 1;
//         //printf("%d newDepth\n", new_depth);
//         // 0 is top-left
//         // 1 is top-right
//         // 2 is bottom-left
//         // 3 is bottom-right
//         subNodes[0] = std::make_unique<ExternalNode>(low_x_bound, high_y_bound, new_depth);
//         subNodes[1] = std::make_unique<ExternalNode>(high_x_bound, high_y_bound, new_depth);
//         subNodes[2] = std::make_unique<ExternalNode>(low_x_bound, low_y_bound, new_depth);
//         subNodes[3] = std::make_unique<ExternalNode>(high_x_bound, low_y_bound, new_depth);
//     }
// };

