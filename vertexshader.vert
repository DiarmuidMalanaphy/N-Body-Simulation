
#version 330 core
layout (location = 0) in vec3 aPos; // the position variable has attribute position 0
layout (location = 1) in vec3 aColor;

out vec4 vertexColor; // specify a color output to the fragment shader


uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 transform;

void main()
{

    gl_Position = projection * view * model * vec4(aPos, 1.0);
    
    vertexColor = vec4(aColor, 1.0);  // Pass the vertex color to the fragment shader
}

