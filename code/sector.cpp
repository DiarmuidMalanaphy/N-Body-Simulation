#include "shader_s.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <chrono>  // For high-precision timing
#include <iterator>
#include <random>
#include <limits>
#include <execution>
#include <array>         // For std::array
#include <memory>        // For std::unique_ptr
#include <stdio.h>
#include <windows.h>
#include <algorithm> 
#include <iostream>
#include <cstdlib>
#include <ctime>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <math.h>

typedef struct SphereRenderData{
    unsigned int sphereVAO;
    unsigned int indexCount;
}   SphereRenderData;


class Profiler {
    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> start;
        std::string name;

    public:
        // Constructor that optionally takes a name for the profiled section
        explicit Profiler(const std::string& profileName = "Unnamed") 
            : start(std::chrono::high_resolution_clock::now()), name(profileName) {}

        // Destructor automatically ends profiling
        ~Profiler() {
            //end();
        }

        // Manually end profiling and print results
        void end() {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            //std::cout << name << " -> Time taken: " << duration.count() << " nanoseconds" << std::endl;
            printf("%s -> Time taken: %lld milliseconds\n", name.c_str(), duration.count());
            printf("\n");
        }

        // Restart the profiler
        void restart() {
            start = std::chrono::high_resolution_clock::now();
        }

        // Get elapsed time without ending the profile
        long long elapsed() const {
            auto now = std::chrono::high_resolution_clock::now();
            return std::chrono::duration_cast<std::chrono::microseconds>(now - start).count();
        }
    };




class Particle {
    public:
    glm::vec3 position; 
    int mass;
    SphereRenderData sphereVAOData;
    unsigned int particleID;
    unsigned int sectorID;

    Shader* shader;  // Change to pointer

    Particle(float x, float y, float z, int mass, unsigned int particleID, unsigned int sectorID, SphereRenderData sphereVAOData, Shader* shader)
        : position(glm::vec3(x, y, z)), mass(mass), particleID(particleID), sectorID(sectorID), sphereVAOData(sphereVAOData), shader(shader)
    {}
    Particle()
        : position(glm::vec3(0.0f, 0.0f, 0.0f)),  // Default position at origin
          mass(1),
          sectorID(0),
          particleID(0),
          sphereVAOData(),                        // Default initialization for sphereVAOData
          shader(nullptr)
    {}


    void printPosition() 
    {
        std::cout << "(" << position.x << ", " << position.y << ", " << position.z << ")" << std::endl;
    }

    void updateSectorID(unsigned int newSectorID)
    {
        sectorID = newSectorID;
    }

    glm::vec2 getForceInteraction(const Particle& terminalParticle) const 
        {
            float dx = position.x - terminalParticle.position.x;
            float dy = position.y - terminalParticle.position.y;
            glm::vec2 directionalVector = glm::normalize(glm::vec2(dx, dy));
            float euclideanDistance = dx * dx + dy * dy;
            float forceMagnitude;
            if (euclideanDistance != 0.0f) {
                printf("%f eucl\n", euclideanDistance);
                forceMagnitude = (mass * terminalParticle.mass) / (euclideanDistance * euclideanDistance);
            } else
            {
                forceMagnitude = (mass * terminalParticle.mass) / 0.001;
            }
            
            
            return (directionalVector * forceMagnitude);

        }

    void updateParticle(glm::vec2 movement){

        glBindVertexArray(sphereVAOData.sphereVAO);
        glm::mat4 model = glm::mat4(1.0f);
        glm::vec3 Dmovement = glm::vec3(movement.x, movement.y, 0.0f);
        position += Dmovement;
        model = glm::translate(model, position);
        shader->setMat4("model", model);

        }
    
    void drawParticle() {
        glDrawElements(GL_TRIANGLE_STRIP, sphereVAOData.indexCount, GL_UNSIGNED_INT, 0);
    }

    void updateAndDrawParticle(glm::vec2 movement) {
        updateParticle(movement);
        drawParticle();
    }
};



class Bound
{
    public:

    float lower;
    float upper;
    float mid;
    Bound() : lower(0.0f), upper(0.0f), mid(0.0f) {}

    Bound(float lower, float upper)
        : lower(lower), upper(upper)
    {
        recalibrateMid();
    }

    void recalibrateMid()
    {
        this->mid = (lower+upper)/2;
    }

    void print()
    {

        std::cout << "Lower : " << lower << ", Middle : " << mid << " , Upper : " << upper << std::endl;
    }

};

class Node 
{
public:
    Bound xBounds;
    Bound yBounds;
    bool filled = false;
    bool internal = true;
    Node(Bound xBounds, Bound yBounds)
        : xBounds(xBounds), yBounds(yBounds)
    {
        xBounds.recalibrateMid();
        yBounds.recalibrateMid();
    }
    Node() 
        : xBounds(), yBounds() // Initializes both bounds using their default constructor
    {}
    virtual void addParticle(Particle* particle) = 0;
    virtual ~Node() = default;
};

class ExternalNode : public Node 
{

    public:
    Particle* particle;
    ExternalNode(Bound xBounds, Bound yBounds) 
        : Node{xBounds, yBounds}, particle(nullptr)
    {
        this->internal = false;
        this->filled = false;
    }

    void addParticle(Particle* particle)
    {
        this->filled = true;
        this->particle = particle;
    }
};




class InternalNode : public Node
{
    // 0 is top-left
    // 1 is top-right
    // 2 is bottom-left
    // 3 is bottom-right
    
    public:

    float mass;
    std::array<std::unique_ptr<Node>, 4> subNodes;
    glm::vec3 centreOfMass;
    void addParticle(Particle* particle)
    {
        //Every time a particle passes through here we must recalibrate the centre of mass 
        recalibrateCentreOfMass(particle);
        int boundIndex = getBoundIndex(particle);
        Node* targetNode = subNodes[boundIndex].get();
        

        if (targetNode->internal == false && targetNode->filled == true) 
        {
            Particle* existingParticle = static_cast<ExternalNode*>(subNodes[boundIndex].get())->particle;
            subNodes[boundIndex] = std::make_unique<InternalNode>(subNodes[boundIndex]->xBounds, subNodes[boundIndex]->yBounds);
            
            targetNode = subNodes[boundIndex].get();
            targetNode->addParticle(existingParticle);           
        }
        targetNode->addParticle(particle);
    }

    int getBoundIndex(Particle *particle)
    {
        bool top = (particle->position.y > yBounds.mid);
        bool right = (particle->position.x > xBounds.mid);

        if (top)
        {
            if (!right) 
            {
                //Top left
                return 0;
            }
            //Top-right
            return 1;
        }
        if (!right) 
            {
                //Bottom left
                return 2;
            }
        //Bottom-right
        return 3;
    }

    void recalibrateCentreOfMass(Particle *particle) 
    {
        glm::vec3 weightedSum = centreOfMass * mass;
        weightedSum += particle->position * static_cast<float>(particle->mass); 
        mass += static_cast<float>(particle->mass);
        centreOfMass = weightedSum / mass;
    }


    InternalNode()
        : Node(),
          centreOfMass(0.0f, 0.0f, 0.0f),
          mass(0) 

    {}

    InternalNode(Bound xBounds, Bound yBounds)
        : Node{xBounds, yBounds}, 
          centreOfMass(0.0f, 0.0f, 0.0f), 
          mass(0)
    {
        initialiseSubNodes();
    }
    
    void initialiseSubNodes()
    {
        Bound low_x_bound = Bound
                        {xBounds.lower, xBounds.mid};

        Bound high_x_bound = Bound
                        {xBounds.mid, xBounds.upper};

        Bound low_y_bound = Bound
                        {yBounds.lower, yBounds.mid};

        Bound high_y_bound = Bound
                        {yBounds.mid, yBounds.upper};

        // 0 is top-left
        // 1 is top-right
        // 2 is bottom-left
        // 3 is bottom-right
        subNodes[0] = std::make_unique<ExternalNode>(low_x_bound, high_y_bound);
        subNodes[1] = std::make_unique<ExternalNode>(high_x_bound, high_y_bound);
        subNodes[2] = std::make_unique<ExternalNode>(low_x_bound, low_y_bound);
        subNodes[3] = std::make_unique<ExternalNode>(high_x_bound, low_y_bound);

    }
};




class BarnesHut 
{
    InternalNode rootNode;
    float theta;
    Bound xBounds;
    Bound yBounds;
    std::vector<Particle>& particles;
    public:

    BarnesHut(std::vector<Particle>& particles, Bound xBounds, Bound yBounds, float theta = 0.5f)
        : theta(theta), xBounds(xBounds), yBounds(yBounds), particles(particles)
    {
        calculateRootNode();
    }

    void calculateRootNode()
    {
        this->rootNode = InternalNode(this->xBounds, this->yBounds);
        for (auto &particle : this->particles) 
            {
                rootNode.addParticle(&particle); 
            }
    }

    std::vector<Node*> getNodes(Particle* particle)
    {
        std::vector<Node*> nodes;
        getNodesRecursive(&rootNode, particle, nodes);
        return nodes;
    }

    void updateAndDrawParticles(float deltaTime)
    {
        std::vector<glm::vec2> forces(particles.size(), glm::vec2(0.0f, 0.0f));

        for (size_t i = 0; i < particles.size(); ++i)
        {
            auto& particle = particles[i];
            std::vector<Node*> relevantNodes = getNodes(&particle);

            for (const auto* node : relevantNodes)
            {
                if (node->internal)
                {
                    const InternalNode* internalNode = static_cast<const InternalNode*>(node);
                    glm::vec2 force = calculateForceFromNode(particle, internalNode);
                    forces[i] += force;
                }
                else
                {
                    const ExternalNode* externalNode = static_cast<const ExternalNode*>(node);
                    if (externalNode->particle != &particle && externalNode->filled)
                    {
                        glm::vec2 force = particle.getForceInteraction(*externalNode->particle);
                        forces[i] += force;
                    }
                }
            }
        }

        // Update and draw particles
        for (size_t i = 0; i < particles.size(); ++i)
        {
            glm::vec2 movement = forces[i] * deltaTime;
            particles[i].updateAndDrawParticle(movement);
        }

        // Recalculate the tree for the next frame
        calculateRootNode();
    }

    

private:
    
    void getNodesRecursive(Node* node, Particle* particle, std::vector<Node*>& nodes)
    {
        if (node == nullptr) return;

        if (!node->internal)
        {
            ExternalNode* externalNode = static_cast<ExternalNode*>(node);
            if (externalNode->particle != particle && externalNode->filled)
            {
                nodes.push_back(node);
            }
            return;
        }

        InternalNode* internalNode = static_cast<InternalNode*>(node);
        float s = (std::max)(node->xBounds.upper - node->xBounds.lower, 
                           node->yBounds.upper - node->yBounds.lower);
        float d = glm::distance(internalNode->centreOfMass, particle->position);

        if (s / d < theta)
        {
            // Use the internal node as an estimation as it is far enough away
            nodes.push_back(node);
            return;
        }
        // We need to go deeper into the tree
        for (const auto& subNode : internalNode->subNodes)
        {
            getNodesRecursive(subNode.get(), particle, nodes);
        }
    }
    glm::vec2 calculateForceFromNode(const Particle& particle, const InternalNode* node) const
    {
        glm::vec2 direction = node->centreOfMass - particle.position;
        float distance = glm::length(direction);
        float forceMagnitude = (particle.mass * node->mass) / (distance * distance);
        return glm::normalize(direction) * forceMagnitude;
    }

};

class Sector;

typedef struct DistanceRelationship{
    Sector *terminalSector;
    float distance;

} DistanceRelationship;


class Sector {
    std::vector<DistanceRelationship> distanceVectors;
    
    glm::vec2 x_bounds;
    glm::vec2 y_bounds;
    glm::vec2 z_bounds; 
    unsigned int depth;
    public:

    glm::vec2 centreOfMass; // (x,y)
    std::vector<Particle*> particles;
    unsigned int sectorID;
    // Default constructor
    Sector() : particles(), x_bounds(), y_bounds(), z_bounds(), sectorID(0), depth(1) {}

    Sector(unsigned int id, unsigned int depthLevel, glm::vec2 x_bounds, glm::vec2 y_bounds, glm::vec2 z_bounds) 
        : sectorID(id), depth(depthLevel), x_bounds(x_bounds), y_bounds(y_bounds), z_bounds(z_bounds){}
    
    //Returning boolean 0 1 if it adds
    unsigned int checkAndAddParticle(Particle* particle) 
    {
        if (particle->position.x > x_bounds.x && particle->position.x < x_bounds.y && 
            particle->position.y > y_bounds.x && particle->position.y < y_bounds.y) 
        {
            particles.push_back(particle);
            particle->sectorID = sectorID;
            return 1;
        }
        return 0;
    }

    void printCentreOfMass()
    {
        printf("(%f, %f)\n", centreOfMass.x, centreOfMass.y);
    }

    unsigned int appendID(unsigned int baseID, unsigned int appendNum) 
    {
        unsigned int factor = 1;
        while (appendNum >= factor) {
            factor *= 10;  // Increase factor to match the digit count of appendNum
        }
        return baseID * factor + appendNum;
    }


    inline void calculateCentreOfMass() 
    {
        glm::vec2 cumulativeVector(0.0f, 0.0f); 
        for (auto & particle : particles) {
            cumulativeVector.x += particle->position.x;
            cumulativeVector.y += particle->position.y;
        }
        
        if (!particles.empty()) {
            centreOfMass = cumulativeVector / static_cast<float>(particles.size());
        }
        
        
    }

    int calculateMass() 
    {
        return particles.size(); 
    }

    float calculateDistance(const Sector& terminalSector) 
    {
        float dx = centreOfMass.x - terminalSector.centreOfMass.x;
        float dy = centreOfMass.y - terminalSector.centreOfMass.y;
        
        return (dx * dx + dy * dy);
    }

    glm::vec2 calculateAttractiveForceFromSectorsParticles(const Particle& startingParticle) 
    {
        glm::vec2 modifiedVec = glm::vec2(0.0f, 0.0f);
        for (Particle* & terminalParticle : particles) 
        {

            modifiedVec += startingParticle.getForceInteraction(*terminalParticle);
            //modifiedVec += startingParticle.getForceInteraction(terminalParticle);
        }
        return modifiedVec;

    }

    std::vector<glm::vec2> getAttractiveParticleLevelForces()
    {
        const size_t endIDx = static_cast<size_t>(distanceVectors.size() * 0.05f);
        
        const size_t particleCount = particles.size();
        std::vector<glm::vec2> ourParticlesForces;
        ourParticlesForces.resize(particles.size());
        
        glm::vec2 modifiedVec = glm::vec2(0.0f, 0.0f);
        std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),[&](const auto& particle) 
        {
                    size_t i = &particle - &particles[0];
                    glm::vec2 modifiedVec(0.0f, 0.0f);
                    //Profiler r("Other ");


                    for (size_t j = 0; j < endIDx; ++j) 
                    {
                        const auto& distanceVector = distanceVectors[j];
                        modifiedVec += distanceVector.terminalSector->calculateAttractiveForceFromSectorsParticles(*particle);
                    }
                    //r.end();
                    //Profiler d("Inner ");

                    for (size_t j = i + 1; j < particleCount; ++j) 
                    {
                        glm::vec2 force = particle->getForceInteraction(*particles[j]);
                        modifiedVec += force;
                        ourParticlesForces[j] -= force;
                    }

                    //d.end();
                    ourParticlesForces[i] += modifiedVec;

        });
        return ourParticlesForces;
   
    }

    
    glm::vec2 getAttractiveSectorLevelForce()
    {
        //Top 5% we actually do particle stuff on -> Defer this to the particles themselves.
        // Top 50% we simulate with sector calculations -> This function handles this bit we essentially just want to calculate and return the force pull on this sector from every other sector..
        // Bottom 50% we either ignore or do hyper optimised calculations with.
        //
        //
        int sectorMass = calculateMass();
        glm::vec2 totalForce(0.0f, 0.0f);
        glm::vec2 direction;
        size_t startIDx = static_cast<size_t>(distanceVectors.size() * 0.05f);
        size_t endIDx = static_cast<size_t>(distanceVectors.size() * 0.950f);

        for (size_t i = startIDx; i < endIDx; ++i) {
            DistanceRelationship& distanceVector = distanceVectors[i];
            float forceMagnitude = (sectorMass * distanceVector.terminalSector->calculateMass()) / (distanceVector.distance);
            direction = glm::normalize(distanceVector.terminalSector->centreOfMass - centreOfMass);
            //printf("magnitude %f\n", forceMagnitude);
            totalForce += direction * forceMagnitude;
        }


        return totalForce;
    }

    void updateAndDrawParticles() 
    {

        glm::vec2 sectorLevelAttraction = getAttractiveSectorLevelForce();

        std::vector<glm::vec2> particleLevelForces = getAttractiveParticleLevelForces();
        
        glm::vec2 modifiedVec;
        for (size_t i = 0; i < particles.size(); ++i) 
        {
            modifiedVec = sectorLevelAttraction * 10.0f + particleLevelForces[i];
            particles[i]->updateAndDrawParticle(modifiedVec);

        }


    }

    void updateParticles() 
    {

        glm::vec2 sectorLevelAttraction = getAttractiveSectorLevelForce();

        //std::cout << "BABE" << std::endl;           
        std::vector<glm::vec2> particleLevelForces = getAttractiveParticleLevelForces();
        
        //std::cout << "BABE" << std::endl;           
        glm::vec2 modifiedVec;
        for (size_t i = 0; i < particles.size(); ++i) 
        {
            modifiedVec = sectorLevelAttraction + particleLevelForces[i];
            particles[i]->updateParticle(modifiedVec);

        }


    }

    void drawParticles() 
    {
        for (size_t i = 0; i < particles.size(); ++i) 
        {
            particles[i]->drawParticle();

        }
    }


    void addDistanceCached(DistanceRelationship distance) {
        distanceVectors.push_back(distance);    
    }

    void addDistanceFromCalculation(Sector* terminalSector) 
    {
        float distance = calculateDistance(*terminalSector);
        
        DistanceRelationship mine = DistanceRelationship{terminalSector, distance};
        
        //
        DistanceRelationship theirs = DistanceRelationship{this, distance};
        distanceVectors.push_back(mine);
        terminalSector->addDistanceCached(theirs);

    }

    inline void sortDistanceVectors() {
        std::sort(distanceVectors.begin(), distanceVectors.end(), 
            [](const DistanceRelationship& a, const DistanceRelationship& b) {
                return a.distance < b.distance; // Sort by distance (ascending)
        });
    }

    std::vector<Sector> getSplitBoundaries() 
    {
        float mid_x_point = (x_bounds.y+x_bounds.x)/2; 
        glm::vec2 low_x_bound = glm::vec2(x_bounds.x, mid_x_point);
        glm::vec2 high_x_bound = glm::vec2(mid_x_point, x_bounds.y);
        
        float mid_y_point = (y_bounds.y+y_bounds.x)/2;
        glm::vec2 low_y_bound = glm::vec2(y_bounds.x, mid_y_point);
        glm::vec2 high_y_bound = glm::vec2(mid_y_point, y_bounds.y);
        Sector topLeft = Sector(appendID(sectorID, 1), depth + 1, low_x_bound, high_y_bound, low_y_bound);

        Sector topRight = Sector(appendID(sectorID, 2), depth + 1, high_x_bound, high_y_bound, low_y_bound);

        Sector bottomLeft = Sector(appendID(sectorID, 3), depth + 1, low_x_bound, low_y_bound, low_y_bound);

        Sector bottomRight = Sector(appendID(sectorID, 4), depth + 1, high_x_bound, low_y_bound, low_y_bound);
        
        return {topLeft, topRight, bottomLeft, bottomRight};
    }

    std::vector<Sector> fillSectors(std::vector<Sector>& sectors, std::vector<Particle>& particles) {
        
        std::vector<bool> particleProcessed(particles.size(), false);
        for (Sector& sector : sectors) {
            for (size_t i = 0; i < particles.size(); ++i) {
                if (!particleProcessed[i] && sector.checkAndAddParticle(&particles[i])) {
                    particleProcessed[i] = true;
                }
            }
        }

        return sectors;
    }

    std::vector<Sector> split(int max_depth, int min_particles) {
        
        //printf("%d, \n", particles.size());

        if (depth >= max_depth) {
            return std::vector<Sector>();
        }
        if (particles.size() < min_particles) {
            return std::vector<Sector>();
        }
        std::vector<Sector> childSectors = getSplitBoundaries();
        for (auto& particle : particles) {
            for (Sector& sector : childSectors) {
                if (sector.checkAndAddParticle(particle)) {
                    break;
                }
            }
        }

        return childSectors;
    }

    void printBounds()
    {
        std::cout << "X Bound : (" << x_bounds.x << " , " << x_bounds.y << ")"  <<std::endl;
        std::cout << "Y Bound : (" << y_bounds.x << " , " << y_bounds.y << ")"  <<std::endl;
    }

    void print() 
    {
        std::cout << "Sector : " << sectorID << std::endl;
        std::cout << "Depth : " << depth << std::endl;
        std::cout << "Number of Particles : " << particles.size() << std::endl;
        printBounds();
        printf("\n\n");
    }

    


};


class SectorMap    
{

    private:

    std::vector<Particle> particles;
    std::vector<Sector> storedSectors;
    int max_depth;
    int bounds;
    int min_particles;
    public:
    
    SectorMap() {}
    
    SectorMap(std::vector<Particle> particles, int number_of_particles, int bounds, int max_depth, int min_particles) : max_depth(max_depth), min_particles(min_particles), particles(particles), bounds(bounds) 
     
    {
    }


    Sector createInitialSector(int bounds) 
        {
            float halfBound = bounds/2;
            Sector sector = Sector(1, 1, glm::vec2(-halfBound, halfBound), glm::vec2(-halfBound, halfBound), glm::vec2(0 - halfBound, halfBound));
            return sector; 
        }


    void calculateSectors()
    {
            //printf("Particles : %d\n", particles.size());
            Profiler p("Testing ");
            Sector initialSector = createInitialSector(bounds);
            std::vector<Sector> sectorVector = std::vector<Sector>();
            sectorVector.push_back(initialSector);
            sectorVector = initialSector.fillSectors(sectorVector, particles);

            std::vector<Sector> resultantVectors;

            storedSectors.clear();
            
            while (sectorVector.size() != 0) 
            {
                Sector currentSector = sectorVector.back();

                resultantVectors = currentSector.split(max_depth, min_particles);

                sectorVector.pop_back();

                if (resultantVectors.size() == 0 || resultantVectors.size() == 1) { 
                    
                    currentSector.calculateCentreOfMass();
                    storedSectors.push_back(currentSector);

                } else {
                    
                    sectorVector.insert(sectorVector.end(), resultantVectors.begin(), resultantVectors.end());
                }
            }
    
            //for (auto& sector : storedSectors)
            //{
            //    sector.printCentreOfMass();
            //}
            p.end();

            int count = 0;
            int n = storedSectors.size(); 
            for (int i = 0; i < n; ++i) 
            {
                // No need to check processedMask since we only process each element once in this loop
                for (int j = i + 1; j < n; ++j)  // Only compare with the elements ahead of i
                {
                    count += 1;
                    storedSectors[i].addDistanceFromCalculation(&storedSectors[j]);
                }
            }

            for (int i = 0; i < storedSectors.size(); ++i)
            {
                storedSectors[i].sortDistanceVectors();
            }
            
        }

    void updateAndDrawParticles()
    {
        for (auto& sector : storedSectors) {

            sector.updateAndDrawParticles();
        }

    }

    void updateParticles()
    {
        for (auto& sector : storedSectors) {
            sector.updateParticles();
        }

    }

    void drawParticles()
    {
        for (auto& sector : storedSectors) {
            sector.drawParticles();
        }

    }
        
        
    
   
    
};

    


SphereRenderData renderSphere(float sizeOfSpheres);
//std::vector<Sector> splitIntoSectors(Particle *particles, unsigned int numParticles, int bounds);
void processInput(GLFWwindow *window);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
glm::vec3 getDirectionalVector(Particle startingParticle, Particle terminalParticle);
double getScalarAngle(Particle startingParticle, Particle terminalParticle);
float getAttractiveForce(Particle startingParticle, Particle terminalParticle);

#define _USE_MATH_DEFINES 
const float cameraSpeed = 50.0f;
const unsigned int SCR_WIDTH = 1600;
const unsigned int SCR_HEIGHT = 1000;
const float SIZE_OF_SPHERES = 1.0f;
const int bounds = 2000;
const int numParticles = 2000;
const int min_particles = 50;
const int max_depth = 15;

glm::vec3 cameraPos   = glm::vec3(1.0f, 0.0f,  1.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* window = glfwCreateWindow(SCR_HEIGHT, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    std::cout<< "HERE" << std::endl;
    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    Shader shader("vertexshader.vert", "fragshader.frag"); 
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glEnable(GL_DEPTH_TEST);  

    SphereRenderData sphereData = renderSphere(SIZE_OF_SPHERES); 

    //
    std::vector<Particle> particles;
    particles.reserve(numParticles);

    
    // Seed the random number generator for different positions
    // Generate 100 particles with random X and Y, all with Z = -100
    //Sector(0, 0, glm::vec2(0, 1000), glm::vec2(0, 1000), glm::vec2(0, 1000));
    Bound xBound = Bound(-bounds/2, bounds/2);
    Bound yBound = Bound(-bounds/2, bounds/2);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(xBound.lower, xBound.upper);
    std::vector<float> randomXVals(numParticles);
    std::vector<float> randomYVals(numParticles);
    std::vector<float> randomZVals(numParticles);

    std::generate(randomXVals.begin(), randomXVals.end(), [&]() { return dis(gen); });
    std::generate(randomYVals.begin(), randomYVals.end(), [&]() { return dis(gen); });
    std::generate(randomZVals.begin(), randomZVals.end(), [&]() { return -1000.0f - (static_cast<float>(std::rand()) / RAND_MAX) * 40.0f; });
    
    for (int i = 0; i < numParticles; ++i) {


        particles.emplace_back(randomXVals[i], randomYVals[i], randomZVals[i], 1, i, 0, sphereData, &shader);
    }
    
    SectorMap secMap = SectorMap(particles, numParticles, bounds, max_depth, min_particles);
    secMap.calculateSectors();      
    BarnesHut bhut = BarnesHut(particles, xBound, yBound);
    

   
    glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);
    glm::vec3 cameraDirection = glm::normalize(cameraPos - cameraTarget);
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f); 
    glm::vec3 cameraRight = glm::normalize(glm::cross(up, cameraDirection));
    glm::vec3 cameraUp = glm::cross(cameraDirection, cameraRight);
    glm::mat4 view;



    glm::mat4 model;
    glm::mat4 projection;
    unsigned int modelLoc;
    unsigned int viewLoc;
    unsigned int timeLoc;
    float lastFrame = 0.0f;
    float deltaTime = 0.0f;
    while(!glfwWindowShouldClose(window))
        {
           
            processInput(window);
            glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

            float timeValue = glfwGetTime();
            timeLoc = glGetUniformLocation(shader.ID, "time");
            shader.use(); 
            model = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
            //const float radius = 40.0f;
            glm::mat4 view;
            view = glm::lookAt(cameraPos, cameraPos + cameraFront ,cameraUp);
            projection = glm::mat4(1.0f);
            
            projection = glm::perspective(glm::radians(60.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 10000.0f);
            shader.setMat4("projection", projection);
            shader.setFloat("time", timeValue);
            glUniform1f(timeLoc, timeValue);
            
            modelLoc = glGetUniformLocation(shader.ID, "model");
            viewLoc  = glGetUniformLocation(shader.ID, "view");
            glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
            glUniformMatrix4fv(viewLoc, 1, GL_FALSE, &view[0][0]);

            
            float currentFrame = glfwGetTime();
            deltaTime = currentFrame - lastFrame;
            lastFrame = currentFrame;
            Profiler p("Barnes");
            bhut.updateAndDrawParticles(deltaTime*1500);
            p.end();
            

            glfwSwapBuffers(window);
            glfwPollEvents();  
           
        }

    glfwTerminate();

    return 0;
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}  

void processInput(GLFWwindow *window)
    {
        if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            {
            glfwSetWindowShouldClose(window, true);
            }

        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            {
                cameraFront += cameraSpeed * glm::vec3(-0.05f, 0, 0);
                if (cameraFront.x < -2.5){
                    cameraFront.x = 2.5;
                }
            }
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            {
                cameraFront += cameraSpeed * glm::vec3(0, 0.05f, 0);
                if (cameraFront.y > 2.5) {
                    cameraFront.y = -2.5;
                }
            }
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            {
                cameraFront += cameraSpeed * glm::vec3(0, -0.05f, 0);
                if (cameraFront.y < -2.5) {
                    cameraFront.y = 2.5;
                }
            }

        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            {
                cameraFront += cameraSpeed * glm::vec3(0.05f, 0, 0);
                if (cameraFront.x > 2.5) {
                    cameraFront.x = -2.5;
                }
            }
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
            cameraPos += cameraSpeed * cameraFront;
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
            cameraPos -= cameraSpeed * cameraFront;
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
            cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
            cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;


    }

unsigned int sphereVAO = 0;
unsigned int indexCount;
//Size is what the size is multiplied (if it's below zero it will break 
SphereRenderData renderSphere(float sizeOfSpheres)
{
    if (sizeOfSpheres < 0) {
        abort();
    }
    //Would add a check to verify if number of spheres is below zero, but it's an impossibility.

    if (sphereVAO == 0)
    {
        glGenVertexArrays(1, &sphereVAO);
        unsigned int vbo, ebo;
        glGenBuffers(1, &vbo);
        glGenBuffers(1, &ebo);

        std::vector<glm::vec3> positions;
        std::vector<glm::vec3> colors;
        std::vector<glm::vec3> normals;
        std::vector<unsigned int> indices;
        const unsigned int X_SEGMENTS = 64;
        const unsigned int Y_SEGMENTS = 64;
        const float PI = 3.14159265359f;
        const glm::vec3 lightBlueColor(0.68f, 0.85f, 0.9f);
        for (unsigned int x = 0; x <= X_SEGMENTS; ++x)
        {
            for (unsigned int y = 0; y <= Y_SEGMENTS; ++y)
            {
                float xSegment = (float)x / (float)X_SEGMENTS;
                float ySegment = (float)y / (float)Y_SEGMENTS;
                float xPos = std::cos(xSegment * 2.0f * PI) * std::sin(ySegment * PI) * sizeOfSpheres ;
                float yPos = std::cos(ySegment * PI) * sizeOfSpheres ;
                float zPos = std::sin(xSegment * 2.0f * PI) * std::sin(ySegment * PI) * sizeOfSpheres;

                positions.push_back(glm::vec3(xPos, yPos, zPos));
                colors.push_back(lightBlueColor);
                normals.push_back(glm::vec3(xPos, yPos, zPos));
            }
        }

        bool oddRow = false;
        for (unsigned int y = 0; y < Y_SEGMENTS; ++y)
        {
            if (!oddRow) // even rows: y == 0, y == 2; and so on
            {
                for (unsigned int x = 0; x <= X_SEGMENTS; ++x)
                {
                    indices.push_back(y * (X_SEGMENTS + 1) + x);
                    indices.push_back((y + 1) * (X_SEGMENTS + 1) + x);
                }
            }
            else
            {
                for (int x = X_SEGMENTS; x >= 0; --x)
                {
                    indices.push_back((y + 1) * (X_SEGMENTS + 1) + x);
                    indices.push_back(y * (X_SEGMENTS + 1) + x);
                }
            }
            oddRow = !oddRow;
        }
        indexCount = static_cast<unsigned int>(indices.size());
        std::vector<float> data;
        for (unsigned int i = 1; i < positions.size(); ++i)
        {
            data.push_back(positions[i].x);
            data.push_back(positions[i].y);
            data.push_back(positions[i].z);
            // Push color (light blue)
            data.push_back(colors[i].r);
            data.push_back(colors[i].g);
            data.push_back(colors[i].b);

            if (normals.size() > 0)
            {
                data.push_back(normals[i].x);
                data.push_back(normals[i].y);
                data.push_back(normals[i].z);
            }
        }
        glBindVertexArray(sphereVAO);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(float), &data[0], GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
        unsigned int stride = (3 + 3 + 3) * sizeof(float);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, (void*)(6 * sizeof(float)));
    }
    SphereRenderData data = {sphereVAO, indexCount};
    return data;
}





double getScalarAngle(Particle startingParticle, Particle terminalParticle){
    double dotProduct = glm::dot(startingParticle.position, terminalParticle.position);
    double magnitudeVec1 = glm::length(startingParticle.position);
    double magnitudeVec2 = glm::length(terminalParticle.position);
    float angleInRadians = std::acos(dotProduct / (magnitudeVec1 * magnitudeVec2)) ;
    double angleInDegrees = glm::degrees(angleInRadians);
    return angleInDegrees;
}



