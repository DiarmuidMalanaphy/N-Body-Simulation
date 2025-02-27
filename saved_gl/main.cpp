#include "shader_s.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <math.h>

void processInput(GLFWwindow *window);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

unsigned int SCR_WIDTH = 800;
unsigned int SCR_HEIGHT = 600;

int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* window = glfwCreateWindow(800, 600, "LearnOpenGL", NULL, NULL);
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

    Shader ourShader("vertexshader.vert", "fragshader.frag"); 
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glEnable(GL_DEPTH_TEST);  
    float vertices[] = {
        // positions          // texture coords  // colors (R, G, B)
        -0.5f, -0.5f, -0.5f,        1.0f, 0.0f, 0.0f, // red
         0.5f, -0.5f, -0.5f,          0.0f, 1.0f, 0.0f, // green
         0.5f,  0.5f, -0.5f,          0.0f, 0.0f, 1.0f, // blue
         0.5f,  0.5f, -0.5f,          0.0f, 0.0f, 1.0f, // blue
        -0.5f,  0.5f, -0.5f,          1.0f, 1.0f, 0.0f, // yellow
        -0.5f, -0.5f, -0.5f,          1.0f, 0.0f, 0.0f, // red

        -0.5f, -0.5f,  0.5f,        1.0f, 0.0f, 1.0f, // magenta
         0.5f, -0.5f,  0.5f,          0.0f, 1.0f, 1.0f, // cyan
         0.5f,  0.5f,  0.5f,        1.0f, 0.5f, 0.0f, // orange
         0.5f,  0.5f,  0.5f,          1.0f, 0.5f, 0.0f, // orange
        -0.5f,  0.5f,  0.5f,          0.5f, 0.0f, 1.0f, // violet
        -0.5f, -0.5f,  0.5f,        1.0f, 0.0f, 1.0f, // magenta

        -0.5f,  0.5f,  0.5f,         0.0f, 1.0f, 0.0f, // green
        -0.5f,  0.5f, -0.5f,          0.5f, 0.5f, 0.5f, // gray
        -0.5f, -0.5f, -0.5f,          1.0f, 0.0f, 0.0f, // red
        -0.5f, -0.5f, -0.5f,          1.0f, 0.0f, 0.0f, // red
        -0.5f, -0.5f,  0.5f,         1.0f, 0.0f, 1.0f, // magenta
        -0.5f,  0.5f,  0.5f,         0.0f, 1.0f, 0.0f, // green

         0.5f,  0.5f,  0.5f,         0.0f, 1.0f, 1.0f, // cyan
         0.5f,  0.5f, -0.5f,         0.5f, 0.5f, 1.0f, // light blue
         0.5f, -0.5f, -0.5f,         0.0f, 0.5f, 1.0f, // sky blue
         0.5f, -0.5f, -0.5f,          0.0f, 0.5f, 1.0f, // sky blue
         0.5f, -0.5f,  0.5f,          0.0f, 1.0f, 1.0f, // cyan
         0.5f,  0.5f,  0.5f,          0.0f, 1.0f, 1.0f, // cyan

        -0.5f, -0.5f, -0.5f,          1.0f, 0.0f, 0.0f, // red
         0.5f, -0.5f, -0.5f,          0.0f, 1.0f, 0.0f, // green
         0.5f, -0.5f,  0.5f,          0.0f, 0.0f, 1.0f, // blue
         0.5f, -0.5f,  0.5f,          0.0f, 0.0f, 1.0f, // blue
        -0.5f, -0.5f,  0.5f,          1.0f, 0.0f, 1.0f, // magenta
        -0.5f, -0.5f, -0.5f,          1.0f, 0.0f, 0.0f, // red

        -0.5f,  0.5f, -0.5f,         1.0f, 1.0f, 0.0f, // yellow
         0.5f,  0.5f, -0.5f,         1.0f, 0.5f, 0.0f, // orange
         0.5f,  0.5f,  0.5f,         0.0f, 1.0f, 0.0f, // green
         0.5f,  0.5f,  0.5f,         0.0f, 1.0f, 0.0f, // green
        -0.5f,  0.5f,  0.5f,          1.0f, 1.0f, 0.0f, // yellow
        -0.5f,  0.5f, -0.5f,          1.0f, 1.0f, 0.0f  // yellow
    };

    glm::vec3 cubePositions[] = {
        glm::vec3( 1.0f,  0.0f,  0.0f), 
        glm::vec3( 2.0f,  5.0f, -15.0f), 
        glm::vec3(-1.5f, -2.2f, -2.5f),  
        glm::vec3(-3.8f, -2.0f, -12.3f),  
        glm::vec3( 2.4f, -0.4f, -3.5f),  
        glm::vec3(-1.7f,  3.0f, -7.5f),  
        glm::vec3( 1.3f, -2.0f, -2.5f),  
        glm::vec3( 1.5f,  2.0f, -2.5f), 
        glm::vec3( 1.5f,  0.2f, -1.5f), 
        glm::vec3(-1.3f,  1.0f, -1.5f)  
    };



       
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    //glGenBuffers(1, &EBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // Color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
        // position attribute
    
    

    


    glm::mat4 model;
    glm::mat4 view;
    glm::mat4 projection;
    unsigned int modelLoc;
    unsigned int viewLoc;
    unsigned int timeLoc;
    GLUquadric* quad = gluNewQuadric();
    if (quad != NULL) {
    }
    while(!glfwWindowShouldClose(window))
        {
           

            processInput(window);


            glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

            float timeValue = glfwGetTime();
            timeLoc = glGetUniformLocation(ourShader.ID, "time");
            
            ourShader.use(); 
            model         = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
            view          = glm::mat4(1.0f);
            projection    = glm::mat4(1.0f);
            
            model = glm::rotate(model, glm::radians(-55.0f), glm::vec3(1.0f, 0.0f, 0.0f));
            view  = glm::translate(view, glm::vec3(0.0f, 0.0f, -3.0f));
            projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
            ourShader.setFloat("time", timeValue);
            glUniform1f(timeLoc, timeValue);
            
            modelLoc = glGetUniformLocation(ourShader.ID, "model");
            viewLoc  = glGetUniformLocation(ourShader.ID, "view");
            glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
            glUniformMatrix4fv(viewLoc, 1, GL_FALSE, &view[0][0]);
            ourShader.setMat4("projection", projection);
            


            glBindVertexArray(VAO);
            for(unsigned int i = 0; i < 10; i++)
            {
                // calculate the model matrix for each object and pass it to shader before drawing
                glm::mat4 model = glm::mat4(1.0f);

                float oscillation = sin(glfwGetTime());
                glm::vec3 oscillatingPosition = glm::vec3(
                    oscillation * cubePositions[i].x,
                    oscillation * cubePositions[i].y,
                    oscillation * cubePositions[i].z
                );
                model = glm::translate(model, oscillatingPosition);
                float angle = 20.0f * i; 
                model = glm::rotate(model, glm::radians(angle), glm::vec3(1.0f, 0.3f, 0.5f));
                ourShader.setMat4("model", model);
                
                glDrawArrays(GL_TRIANGLES, 0, 36);           
            }    

            gluSphere(quad, 1.0, 32, 32);  // Draw a sphere using the quadrics object
            glfwSwapBuffers(window);
            glfwPollEvents();  
           
            //glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        }

    gluDeleteQuadric(quad);        // Clean up when done
    glfwTerminate();

    return 0;
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}  

void processInput(GLFWwindow *window)
    {
        if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS){

            glfwSetWindowShouldClose(window, true);
        
            }

        if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
            {
            glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            }
        if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
            {
            glDrawArrays(GL_TRIANGLES, 0, 3);
            }
        if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
            {
            glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
            }


    
}
