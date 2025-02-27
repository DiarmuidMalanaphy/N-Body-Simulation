#include "shader_s.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iterator>
#include <random>
#include <limits>
#include <array>         // For std::array
#include <memory>        // For std::unique_ptr
#include <stdio.h>
#include <windows.h>
#include <algorithm> 
#include <iostream>
#include <unordered_map>
#include <cstdlib>
#include <ctime>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <math.h>

#include "Particle.h"

#include "SphereRender.h"
#include "BarnesHut.h"
#include "Bound.h"

void processInput(GLFWwindow *window);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

#define _USE_MATH_DEFINES 
const float cameraSpeed = 50.0f;
const unsigned int SCR_WIDTH = 1600;
const unsigned int SCR_HEIGHT = 1000;
const float SIZE_OF_SPHERES = 6.0f;
const int BOUNDS_SIZE = 1000;
const float CIRCLE_RADIUS = BOUNDS_SIZE * 0.3f;
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

    GLFWmonitor* monitor = glfwGetPrimaryMonitor(); // Get the primary monitor
    const GLFWvidmode* mode = glfwGetVideoMode(monitor); // Get the monitor's resolution

    GLFWwindow* window = glfwCreateWindow(mode->width, mode->height, "Particle Simulator", monitor, NULL);

    // GLFWwindow* window = glfwCreateWindow(SCR_HEIGHT, SCR_HEIGHT, "Particle simulator", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
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

    std::unordered_map<unsigned int, std::shared_ptr<Particle>> particles;
    
    // Seed the random number generator for different positions
    // Generate 100 particles with random X and Y, all with Z = -100
    //
    Bound xBound = Bound :: Bound(-BOUNDS_SIZE, BOUNDS_SIZE);
    Bound yBound = Bound :: Bound(-BOUNDS_SIZE, BOUNDS_SIZE);
    std::random_device rd;
    std::mt19937 gen(rd());


    std::uniform_real_distribution<float> disRadius(0.0f, BOUNDS_SIZE); // Radius within the circle
    std::uniform_real_distribution<float> disAngle(0.0f, 2 * M_PI);    // Angle for polar coordinates
    std::uniform_real_distribution<float> disZ(-1040.0f, -1000.0f);    // Z-coordinate range

    std::vector<float> randomXVals(numParticles);
    std::vector<float> randomYVals(numParticles);
    std::vector<float> randomZVals(numParticles);

    // Generate random positions in a 2D circle
    for (int i = 0; i < numParticles; ++i) {
        float radius = std::sqrt(disRadius(gen)); // Square root for uniform distribution
        float angle = disAngle(gen);             // Random angle
        randomXVals[i] = radius * 20 * std::cos(angle); // Convert polar to Cartesian X
        randomYVals[i] = radius * 20 * std::sin(angle); // Convert polar to Cartesian Y
        randomZVals[i] = disZ(gen);              // Random Z-coordinate
    }



    // std::uniform_real_distribution<float> dis(xBound.lower/2, xBound.upper/2);
    // std::vector<float> randomXVals(numParticles);
    // std::vector<float> randomYVals(numParticles);
    // std::vector<float> randomZVals(numParticles);
    // std::generate(randomXVals.begin(), randomXVals.end(), [&]() { return dis(gen); });
    // std::generate(randomYVals.begin(), randomYVals.end(), [&]() { return dis(gen); });
    // std::generate(randomZVals.begin(), randomZVals.end(), [&]() { return -1000.0f - (static_cast<float>(std::rand()) / RAND_MAX) * 40.0f; });

    for (int particleID = 0; particleID < numParticles; ++particleID) {
        glm::vec3 position = glm::vec3(randomXVals[particleID], randomYVals[particleID], -500);
        particles[particleID] = std::make_unique<Particle>(position, xBound, yBound, 1.0f, particleID, sphereData, &shader);
    }
  
    
    BarnesHut bhut = BarnesHut :: BarnesHut(particles, xBound, yBound);
    

   
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
            
            projection = glm::perspective(glm::radians(60.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100000.0f);
            shader.setMat4("projection", projection);
            shader.setFloat("time", timeValue);
            glUniform1f(timeLoc, timeValue);
            
            modelLoc = glGetUniformLocation(shader.ID, "model");
            viewLoc  = glGetUniformLocation(shader.ID, "view");
            glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
            glUniformMatrix4fv(viewLoc, 1, GL_FALSE, &view[0][0]);

            
            //Profiler p("Barnes");
            bhut.updateAndDrawParticles();
            //p.end();
            

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








