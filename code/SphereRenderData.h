#ifndef SPHERERENDERDATA_H
#define SPHERERENDERDATA_H

#include <vector>
#include <cmath>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Structure to hold sphere rendering data
typedef struct SphereRenderData
{
    unsigned int sphereVAO;
    unsigned int indexCount;
} SphereRenderData;



#endif // SPHERE_RENDERER_H