#include <stdio.h>
#include <string>

#include <SDL2/SDL.h>

#if defined(__APPLE__)
  #include <OpenGL/OpenGL.h>
  #include <OpenGL/gl3.h>
#else
  #define GL_GLEXT_PROTOTYPES
  #include <GL/glx.h>
#endif

/*
  you must define one of these values, which defines what values are supplied to the
  fragmant shader:
    - FRAG_COORDS_NORMALISED
      inputs are in the range of 0,0 to 1,1.
    - FRAG_COORDS_SCREEN
      inputs are in the range of -1,-1 to 1,1.
    - FRAG_COORDS_VIEWPORT
      inputs are in the range of 0,0 to {size x},{size y}, where size x and size y are
      the current dimensions of the viewport, i.e, the size of the window.
*/
//#define FRAG_COORDS_NORMALISED
#define FRAG_COORDS_SCREEN
//#define FRAG_COORDS_VIEWPORT

struct float2 {
  float x = 0.f;
  float y = 0.f;
};

namespace window {
  struct parameters {
    float2 size = (float2){
      .x = 0,
      .y = 0,
    };
    bool   fullscreen = false;
  };

  static SDL_Window *
  build
  (const parameters &parameters)
  {
    SDL_Window *result = SDL_CreateWindow(
      "Shade",
      SDL_WINDOWPOS_UNDEFINED,
      SDL_WINDOWPOS_UNDEFINED,
      parameters.size.x,
      parameters.size.y,
      SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI
    );

    if (NULL == result)
    {
      printf("Could not initialise SDL: %s\n", SDL_GetError());
      return NULL;
    }

    if (parameters.fullscreen)
      SDL_SetWindowFullscreen(result, SDL_WINDOW_FULLSCREEN_DESKTOP);

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

    SDL_GLContext glContext = SDL_GL_CreateContext(result);
    if (NULL == glContext)
    {
      printf("Could not initialise OpenGL: %s\n", SDL_GetError());
      return NULL;
    }

    if(SDL_GL_SetSwapInterval(1) < 0)
    {
      printf("Could not set VSync: %s\n", SDL_GetError());
      return NULL;
    }

    return result;
  }
}

namespace gl {
  static void
  init
  (void)
  {
    GLuint vertexArray = 0;
    glGenVertexArrays(1, &vertexArray);
    glBindVertexArray(vertexArray);

    glClearColor(0.f, 0.f, 0.f, 1.f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  namespace geometry {
    static GLuint
    build
    (void)
    {
      static const float vertexPositions[] = {
        -1, 1, -1,-1,  1, 1,
         1,-1,  1, 1, -1,-1,
      };

      GLuint result = 0;
      glGenBuffers(1, &result);
      glBindBuffer(GL_ARRAY_BUFFER, result);
      glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(vertexPositions),
        vertexPositions,
        GL_STATIC_DRAW
      );

      return result;
    }
  }

  namespace program {
    namespace shader {
      static GLuint
      build
      (const char *shaderSource, const GLenum shaderType)
      {
        GLuint result = glCreateShader(shaderType);
        glShaderSource(result, 1, &shaderSource, NULL);
        glCompileShader(result);

        GLint ret = GL_FALSE;
        glGetShaderiv(result, GL_COMPILE_STATUS, &ret);

        if (GL_FALSE == ret)
        {
          GLint len = 0; glGetShaderiv(result, GL_INFO_LOG_LENGTH, &len);
          char errorMessage[len+1];
          glGetShaderInfoLog(result, len, NULL, errorMessage);
          printf("Error: Couldn't build shader:\n %s", errorMessage);
          return 0;
        }

        return result;
      }
    }

    static GLuint
    link
    (GLuint vertexShader, GLuint fragmentShader)
    {
      if (0 == vertexShader || 0 == fragmentShader)
        return 0;

      GLuint result = glCreateProgram();
      glAttachShader(result, vertexShader);
      glAttachShader(result, fragmentShader);
      glLinkProgram(result);

      glDetachShader(result, vertexShader);
      glDetachShader(result, fragmentShader);
      glDeleteShader(vertexShader);
      glDeleteShader(fragmentShader);

      GLint ret = GL_FALSE;
      glGetProgramiv(result, GL_LINK_STATUS, &ret);

      if (GL_FALSE == ret)
      {
        GLint len = 0; glGetProgramiv(result, GL_INFO_LOG_LENGTH, &len);
        char errorMessage[len+1];
        glGetProgramInfoLog(result, len, NULL, errorMessage);
        printf("Error: Unable to link shaders:\n %s", errorMessage);
      }

      return result;
    }

    GLuint
    build
    (void)
    {
      static const std::string fragmentShaderSource(
        R"(
          #version 330 core
          in  vec2 fPosition;
          out vec4 color;

          uniform float uTime;
          uniform vec2 uSize;
        )"

        #include "build/fragment.glsl.h"
        
        R"(
          void main() { color.rgba = fragment(fPosition); }
        )"
      );

      static const std::string vertexShaderSource =
        R"(
          #version 330 core
          layout(location = 0) in vec2 vPosition;
          uniform vec2 uSize;

          out vec2 fPosition;

          vec2 getFragmentCoord(vec2 vc) { )"

        #if defined(FRAG_COORDS_SCREEN)
          "return vc;"
        #elif defined(FRAG_COORDS_VIEWPORT)
          "return ((vc/2.f) + vec2(.5f,.5f))*uSize;"
        #elif defined(FRAG_COORDS_NORMALISED)
          "return ((vc/2.f) + vec2(.5f,.5f));"
        #else
          #error Please define one of FRAG_COORDS_SCREEN, FRAG_COORDS_VIEWPORT, or FRAG_COORDS_NORMALISED
        #endif

        R"( }

            void main()
            {
              gl_Position.xy = vPosition;
              gl_Position.zw = vec2(0, 1);
              fPosition      = getFragmentCoord(vPosition);
            }
      )";

      return program::link(
        shader::build(vertexShaderSource.data(), GL_VERTEX_SHADER), 
        shader::build(fragmentShaderSource.data(), GL_FRAGMENT_SHADER)
      );
    }
  }
}

static float2
getDrawableSize
(SDL_Window *window)
{
  int sx = 0, sy = 0;
  SDL_GL_GetDrawableSize(window, &sx, &sy);

  return (float2) {
    .x = (float) sx,
    .y = (float) sy
  };
}

static void
setShaderProgram
(GLuint program, float time, const float2 size)
{
  glUseProgram(program);
  glUniform1f(glGetUniformLocation(program, "uTime"), time);
  glUniform2f(glGetUniformLocation(program, "uSize"), size.x, size.y);
}

int
main
(int argc, char **argv)
{
  const auto parameters = (window::parameters) {
    .size = (float2) {
      .x = 512,
      .y = 512,
    },
    .fullscreen = false,
  };

  SDL_Window *window = window::build(parameters);
  if (NULL == window)
    return 1;

  gl::init();

  GLuint geometry = gl::geometry::build();
  GLuint program = gl::program::build();

  //not real time. each frame increments this by 1/60. might need to adjust for fancy monitors..
  float time = 0.f;
  bool running = (program != 0);
  while(running) {
    const auto size = getDrawableSize(window);

    glViewport(0, 0, size.x, size.y);

    glClear(GL_COLOR_BUFFER_BIT);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, geometry);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);

    setShaderProgram(program, time, size);

    glDrawArrays(GL_TRIANGLES, 0, 12);
    SDL_GL_SwapWindow(window);

    //don't want to run into precision issues! :D
    time = (time < 3600.f) ? (time + 1.f/60.f) : 0.f;

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT)
        running = false;

      if (event.type == SDL_KEYUP)
      if (event.key.keysym.sym == SDLK_ESCAPE)
        running = false;
    }
  }
}
