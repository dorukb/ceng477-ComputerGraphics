#include "EclipseMap.h"

using namespace std;

/* Geometry variables */
glm::mat4 MVP, M_model, M_view, M_projection;

/* Uniform variable locations */
int MVP_location, heightFactor_location, cameraPos_location, heightmap_location,
texture_location, lightPos_location, textureOffset_location;
int activeProgram[9] = {0};


void EclipseMap::initMatrices() {

	/* Initialize Cam vectors first */
	// cameraLeft = glm::cross(cameraUp, cameraDirection);

	/* Now Set MVP */
	// M_model = glm::rotate(M_model, (float) glm::radians(-60.0), glm::vec3(1, 0, 0));

    // glMatrixMode(GL_MODELVIEW);
    // glLoadIdentity();
    // gluLookAt(cameraPosition,cameraPosition + cameraDirection,cameraUp);

    M_model = glm::mat4();
	M_view = glm::lookAt(cameraPosition, cameraPosition + cameraDirection, cameraUp);
	M_projection = glm::perspective(projectionAngle, aspectRatio, near, far);
	MVP = M_projection * M_view * M_model;
}
void EclipseMap::initBuffers() {
	/* Init VAO */

	// /* Init VBOs */

	/* glVertexAttribPointer(array_index, #_of_coods_per_Vertex, type, need_normalization?,
	 * Byte_offset_between_consecutive_vertices, offset_of_fields_in_vertex_structure) */
	// glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) nullptr);
	// glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) offsetof(vertex, normal));
	// glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) (offsetof(vertex, texture)));

}
void EclipseMap::createSphereVertices() {
	float x, y, z, xy, u, v, alpha, beta;
	for (int i = 0; i <= verticalSplitCount; i++) {
		beta = M_PI / 2 - (i * M_PI / verticalSplitCount); // pi/2 to -pi/2
		xy = radius * cosf(beta);
		z = radius * sinf(beta);
		for (int j = 0; j <= horizontalSplitCount; j++) {
			alpha = j * 2 * M_PI / horizontalSplitCount;  // 0 to 2pi
			x = xy * cosf(alpha);
			y = xy * sinf(alpha);
			u = (float)j / horizontalSplitCount;
			v = (float)i / verticalSplitCount;

			vertex vertex(glm::vec3(x,y,z), 
                          glm::normalize(glm::vec3(x/radius, y/radius, z/radius)), 
                          glm::vec2(u, v));
			vertices.push_back(vertex);
		}
	}
}
void EclipseMap::createMoonVertices() {
	float x, y, z, xy, u, v, alpha, beta;
	for (int i = 0; i <= verticalSplitCount; i++) {
		beta = M_PI / 2 - (i * M_PI / verticalSplitCount); // pi/2 to -pi/2
		xy = moonRadius * cosf(beta);
		z = moonRadius * sinf(beta);
		for (int j = 0; j <= horizontalSplitCount; j++) {
			alpha = j * 2 * M_PI / horizontalSplitCount;  // 0 to 2pi
			x = xy * cosf(alpha);
			y = xy * sinf(alpha);
			u = (float)j / horizontalSplitCount;
			v = (float)i / verticalSplitCount;

			vertex vertex(glm::vec3(x,y,z), 
                          glm::normalize(glm::vec3(x/moonRadius, y/moonRadius, z/moonRadius)), 
                          glm::vec2(u, v));
            
            glm::mat4 Ty = glm :: translate(glm::mat4(1.0f),glm::vec3(0.0f,2660.0f,0.0f));
            glm::vec4 temp = Ty *glm::vec4 (vertex.position.x,vertex.position.y,vertex.position.z,1.0f);
            vertex.position = glm::vec3(temp.x,temp.y,temp.z);

			moonVertices.push_back(vertex);
		}
	}
}
void EclipseMap::createMoonIndices() {
	/* Initialize indices per pixel, be careful about the winding order! */
	int k1, k2;
	for(int i = 0; i < verticalSplitCount; ++i) {
		k1 = i * (horizontalSplitCount + 1);
		k2 = k1 + horizontalSplitCount + 1;
		for(int j = 0; j < horizontalSplitCount; ++j, ++k1, ++k2) {
			// 2 triangles per sector excluding first and last stacks
			/* Provide indices for the first triangle with a correct winding order that suits RH-rule */
			if(i != 0) {
				moonIndices.push_back(k1);
				moonIndices.push_back(k2);
				moonIndices.push_back(k1 + 1);
			}
			/* Provide indices for the second triangle with a correct winding order that suits RH-rule */
			if(i != (horizontalSplitCount - 1)) {
				moonIndices.push_back(k1 + 1);
				moonIndices.push_back(k2);
				moonIndices.push_back(k2 + 1);
			}
		}
	}
}

void EclipseMap::createSphereIndices() {
	/* Initialize indices per pixel, be careful about the winding order! */
	int k1, k2;
	for(int i = 0; i < verticalSplitCount; ++i) {
		k1 = i * (horizontalSplitCount + 1);
		k2 = k1 + horizontalSplitCount + 1;
		for(int j = 0; j < horizontalSplitCount; ++j, ++k1, ++k2) {
			// 2 triangles per sector excluding first and last stacks
			/* Provide indices for the first triangle with a correct winding order that suits RH-rule */
			if(i != 0) {
				indices.push_back(k1);
				indices.push_back(k2);
				indices.push_back(k1 + 1);
			}
			/* Provide indices for the second triangle with a correct winding order that suits RH-rule */
			if(i != (horizontalSplitCount - 1)) {
				indices.push_back(k1 + 1);
				indices.push_back(k2);
				indices.push_back(k2 + 1);
			}
		}
	}
}

void EclipseMap::setUniforms(GLuint worldShaderId, GLuint moonShaderId) {
	/* Set the uniform variables so that our shaders can access them */
	MVP_location = glGetUniformLocation(worldShaderId, "MVP");
	glUniformMatrix4fv(MVP_location, 1, GL_FALSE, glm::value_ptr(MVP));

	heightFactor_location = glGetUniformLocation(worldShaderId, "heightFactor");
	glUniform1f(heightFactor_location, heightFactor);

	textureOffset_location = glGetUniformLocation(worldShaderId, "textureOffset");
	glUniform1i(textureOffset_location, textureOffset);

	cameraPos_location = glGetUniformLocation(worldShaderId, "cameraPosition");
	glUniform3fv(cameraPos_location, 1, glm::value_ptr(cameraPosition));

	lightPos_location = glGetUniformLocation(worldShaderId, "lightPosition");
	glUniform3fv(lightPos_location, 1, glm::value_ptr(lightPos));

	// heightmap_location = glGetUniformLocation(worldShaderId, "TexGrey");
	// glUniform1i(heightmap_location, 0);

	// texture_location = glGetUniformLocation(worldShaderId, "TexColor");
	// glUniform1i(texture_location, 1);
}


void EclipseMap::Render(const char *coloredTexturePath, const char *greyTexturePath, const char *moonTexturePath) {
    // Open window
    GLFWwindow *window = openWindow(windowName, screenWidth, screenHeight);

    // Moon commands
    // Load shaders
    GLuint moonShaderID = initShaders("moonShader.vert", "moonShader.frag");

    initMoonColoredTexture(moonTexturePath, moonShaderID);

    
    // TODO: Set moonVertices
    createMoonVertices();
    createMoonIndices();
    // TODO: Configure Buffers
    

    // World commands
    // Load shaders
    GLuint worldShaderID = initShaders("worldShader.vert", "worldShader.frag");

    initColoredTexture(coloredTexturePath, worldShaderID);
    initGreyTexture(greyTexturePath, worldShaderID);

    // TODO: Set worldVertices

    //edit sphere drawing
    setUniforms(worldShaderID, moonShaderID);

    createSphereVertices();
    createSphereIndices();
    initMatrices();
    // TODO: Configure Buffers
    // Earth buffers
    unsigned int VAO;
    glGenVertexArrays(1, &VAO); 
    glBindVertexArray(VAO);
    
    unsigned int VBO;
    glGenBuffers(1, &VBO);  
    glBindBuffer(GL_ARRAY_BUFFER, VBO);  
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertex), vertices.data(), GL_STATIC_DRAW);
    
    GLuint earthIndicesVBO;
    glGenBuffers(1, &earthIndicesVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, earthIndicesVBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), indices.data(), GL_STATIC_DRAW);

	/* glVertexAttribPointer(array_index, #_of_coods_per_Vertex, type, need_normalization?,
	 * Byte_offset_between_consecutive_vertices, offset_of_fields_in_vertex_structure) */
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) nullptr);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) offsetof(vertex, normal));
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) (offsetof(vertex, texture)));

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);

    // Moon buffers
    unsigned int moonVAO;
    glGenVertexArrays(1, &moonVAO); 
    glBindVertexArray(moonVAO);
    
    unsigned int moonVBO;
    glGenBuffers(1, &moonVBO);  
    glBindBuffer(GL_ARRAY_BUFFER, moonVBO);  
    glBufferData(GL_ARRAY_BUFFER, moonVertices.size() * sizeof(vertex), moonVertices.data(), GL_STATIC_DRAW);
    
    GLuint moonIndicesVBO;
    glGenBuffers(1, &moonIndicesVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, moonIndicesVBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, moonIndices.size() * sizeof(int), moonIndices.data(), GL_STATIC_DRAW);

	/* glVertexAttribPointer(array_index, #_of_coods_per_Vertex, type, need_normalization?,
	 * Byte_offset_between_consecutive_vertices, offset_of_fields_in_vertex_structure) */
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) nullptr);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) offsetof(vertex, normal));
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) (offsetof(vertex, texture)));

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
    glUseProgram(worldShaderID);


    M_projection = glm::perspective(projectionAngle, aspectRatio, near, far);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);

    // Main rendering loop
    do {
        glViewport(0, 0, screenWidth, screenHeight);

        glClearStencil(0);
        glClearDepth(1.0f);
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);



        // TODO: Handle key presses
        handleKeyPress(window);

        // TODO: Manipulate rotation variables
        
        // TODO: Bind textures
        // glUniform1i(glGetUniformLocation(worldShaderID, "TexColor"), 0);

        
        // TODO: Use moonShaderID program
        // glClientActiveTexture(GL_TEXTURE0 + 1);
        // glEnableClientState(GL_TEXTURE_COORD_ARRAY);

		// glUniform1i(glGetUniformLocation(worldShaderID, "TexColor"), 1);
        
        // TODO: Update camera at every frame

        updateCamera();
        // TODO: Update uniform variables at every frame
        
        // TODO: Bind moon vertex array        

        // TODO: Draw moon object
        
        /*************************/

        // TODO: Use worldShaderID program
        
        // TODO: Update camera at every frame

        // TODO: Update uniform variables at every frame
        
        // TODO: Bind world vertex array
        
        // TODO: Draw world object
        // glDrawArrays(GL_TRIANGLES, 0, 3);
        /* Now render the frame */
        float pitchDiff =  pitch-startPitch;
        float yawDiff = yaw-startYaw;

        if(abs(pitchDiff) < FLT_EPSILON) pitchDiff = 0.0f;
        if(abs(yawDiff) < FLT_EPSILON) yawDiff = 0.0f;
        cameraUp = glm::normalize(glm::rotate(cameraStartUp, pitchDiff, cameraLeft));
        cameraDirection = glm::normalize(glm::rotate(cameraStartDirection, pitchDiff, cameraLeft));

		cameraDirection = glm::normalize(glm::rotate(cameraDirection, yawDiff, cameraUp));
        cameraLeft = glm::normalize(glm::rotate(cameraStartLeft, yawDiff, cameraUp));

        cameraPosition += speed * cameraDirection;
        M_view = glm::lookAt(cameraPosition, cameraPosition + cameraDirection, cameraUp); // gluLookAt(eye, center, up)
        MVP = M_projection * M_view * M_model;

        // Do not forget the update the uniforms of geometry too
        glUniformMatrix4fv(MVP_location, 1, GL_FALSE, glm::value_ptr(MVP));
        glUniform3fv(cameraPos_location, 1, glm::value_ptr(cameraPosition));
        glUniform1f(heightFactor_location, heightFactor);

        // std::cout<<"indices size: " << indices.size() << std::endl;
        // draw earth
        glBindVertexArray(VAO);
	    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, nullptr);

        //draw moon
        glBindVertexArray(moonVAO);
	    glDrawElements(GL_TRIANGLES, moonIndices.size(), GL_UNSIGNED_INT, nullptr);

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (!glfwWindowShouldClose(window));

    // Delete buffers
    glDeleteBuffers(1, &moonVAO);
    glDeleteBuffers(1, &moonVBO);
    glDeleteBuffers(1, &moonEBO);

    
    // Delete buffers
    glDeleteBuffers(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
   
    glDeleteProgram(moonShaderID);
    glDeleteProgram(worldShaderID);

    // Close window
    glfwTerminate();
}


void EclipseMap::updateCamera()
{

        if (activeProgram[0] == 1) {
            pitch+= 0.05;
            if (pitch > 360.0)
                pitch = 360.0;

            // cameraUp = glm::rotate(cameraUp, 0.05f, cameraLeft);
		    // cameraDirection = glm::rotate(cameraDirection, 0.05f, cameraLeft);
        }
        if (activeProgram[1] == 1) {

            //pitch decrease
            pitch -= 0.05;
            if (pitch < 0.0)
                pitch = 0.0;
            
            // cameraUp = glm::rotate(cameraUp, -0.05f, cameraLeft);
		    // cameraDirection = glm::rotate(cameraDirection, -0.05f, cameraLeft);
            //heightFactor -= 0.5;
            //glUniform1f(heightFactor_location, heightFactor);
        }
        if (activeProgram[2] == 1) {
            //right
            yaw += 0.05;
            if (360.0<yaw)
                yaw -= 360.0;
            // cameraLeft = glm::rotate(cameraLeft, 0.05f, cameraUp);
		    // cameraDirection = glm::rotate(cameraDirection, 0.05f, cameraUp);
        }
        if (activeProgram[3] == 1){

            yaw -= 0.05;
            if (yaw < 0.0)
                yaw += 360.0;
            // cameraLeft = glm::rotate(cameraLeft, -0.05f, cameraUp);
		    // cameraDirection = glm::rotate(cameraDirection, -0.05f, cameraUp);
        }
        if (activeProgram[4] == 1) {

            speed += 0.01;

        }
        if (activeProgram[5] == 1) {
            speed -= 0.01;
            if(speed < 0){
                speed = 0;
            }
        }


        if (activeProgram[6] == 1) { 
            speed = 0;
        }
        if (activeProgram[7] == 1) {
            // reset to initial configs
            speed = startSpeed;
            pitch = startPitch;
            yaw = startYaw;
            cameraPosition = cameraStartPosition;
            cameraDirection = cameraStartDirection;
            cameraUp = cameraStartUp;
  
            /*
            pos = glm::vec3(textureWidth / 2.0, textureWidth / 10.0, -textureWidth / 4.0);
            textureOffset = 0;
            glUniform1i(textureOffset_location, textureOffset);
            heightFactor = 10.0;
            glUniform1f(heightFactor_location, heightFactor);
        
             */

        }
        // if (activeProgram[8]){

        //     const GLFWvidmode* mode = glfwGetVideoMode(monitor);

        //     glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);

        //     //make full screen
        //     //glfwSetWindowMonitor(window, nullptr, windowX, windowY, windowX, windowY, 0);
        //     //glViewport(0, 0, windowX, windowY);
        // }
}

void EclipseMap::handleKeyPress(GLFWwindow *window) {

    int key = glfwGetKey(window, GLFW_KEY_ESCAPE);

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    else if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
    {
        activeProgram[0] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
    {
        // cout << "pressed S" << endl;
        activeProgram[1] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
    {
        activeProgram[2] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    {
        activeProgram[3] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
    {
        activeProgram[4] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS)
    {
        activeProgram[5] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
    {
        activeProgram[6] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS)
    {
        // cout<<"I pressed" << endl;
        activeProgram[7] = 1;
    }
    else if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    {
        activeProgram[8]= 1;
    }
    ///// release

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_RELEASE)
    {
        activeProgram[0] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_RELEASE)
    {
        activeProgram[1] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_RELEASE)
    {
        activeProgram[2] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_RELEASE)
    {
        activeProgram[3] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_RELEASE)
    {
        activeProgram[4] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_RELEASE)
    {
        activeProgram[5] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_RELEASE)
    {
        activeProgram[6] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_RELEASE)
    {
        // cout << "I released" << endl;
        activeProgram[7] = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_RELEASE)
    {
        activeProgram[8]= 0;
    }

}

GLFWwindow *EclipseMap::openWindow(const char *windowName, int width, int height) {
    if (!glfwInit()) {
        getchar();
        return 0;
    }

    const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow *window = glfwCreateWindow(width, height, windowName, NULL, NULL);
    glfwSetWindowMonitor(window, NULL, 1, 31, screenWidth, screenHeight, mode->refreshRate);

    if (window == NULL) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwMakeContextCurrent(window);

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glClearColor(0, 0, 0, 0);

    return window;
}


void EclipseMap::initColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &textureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, textureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Colered Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexColor"), 0);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initGreyTexture(const char *filename, GLuint shader) {

    glGenTextures(1, &textureGrey);
    glBindTexture(GL_TEXTURE_2D, textureGrey);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    int width, height;

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Grey Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
  



    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexGrey"), 1);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initMoonColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &moonTextureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, moonTextureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Moon COlor Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "MoonTexColor"), 2);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}
