#pragma once

#include <stdexcept>
#include <iostream>

// // --- legacy OpenGL headers ---
// #if defined(_WIN32)
//     #include <windows.h>
//     #include <GL/gl.h>
//     #include <GL/glu.h>
// #elif defined(__APPLE__)
//     #include <OpenGL/gl.h>
//     #include <OpenGL/glu.h>
// #else
//     #include <GL/gl.h>
//     #include <GL/glu.h>
// #endif


class MechanicalIllustrativeShader
{
public:
    MechanicalIllustrativeShader()
    {
        compile();
    }

    ~MechanicalIllustrativeShader()
    {
        if (program_ != 0)
            glDeleteProgram(program_);
    }

    void use() const
    {
        glUseProgram(program_);
    }

    // ------- Uniform setters (float only) --------

    void setBaseColor(float r,float g,float b) const
    {
        glUseProgram(program_);
        glUniform3f(glGetUniformLocation(program_,"uBaseColor"),r,g,b);
    }

    void setLineColor(float r,float g,float b) const
    {
        glUseProgram(program_);
        glUniform3f(glGetUniformLocation(program_,"uLineColor"),r,g,b);
    }

    void setLightDirView(float x,float y,float z) const
    {
        glUseProgram(program_);
        glUniform3f(glGetUniformLocation(program_,"uLightDirView"),x,y,z);
    }

    void setEdgeWidth(float w) const
    {
        glUseProgram(program_);
        glUniform1f(glGetUniformLocation(program_,"uEdgeWidth"),w);
    }

    void setRim(float strength,float power) const
    {
        glUseProgram(program_);
        glUniform1f(glGetUniformLocation(program_,"uRimStrength"),strength);
        glUniform1f(glGetUniformLocation(program_,"uRimPower"),power);
    }

    GLuint id() const { return program_; }


private:
    GLuint program_ = 0;

    static const char* vertexShader()
    {
        return
        " #version 120\n"
        " varying vec3 vNormalView;\n"
        " varying vec3 vViewPos;\n"
        " void main()\n"
        " {\n"
        "   vec4 viewPos = gl_ModelViewMatrix * gl_Vertex;\n"
        "   vViewPos = viewPos.xyz;\n"
        "   vNormalView = normalize(gl_NormalMatrix * gl_Normal);\n"
        "   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
        " }\n";
    }

    static const char* fragmentShader()
{
    return
    " #version 120\n"
    " varying vec3 vNormalView;\n"
    " varying vec3 vViewPos;\n"
    " uniform vec3 uBaseColor;\n"
    " uniform vec3 uLineColor;\n"
    " uniform vec3 uLightDirView;\n"
    " uniform float uEdgeWidth;\n"
    " uniform float uRimStrength;\n"
    " uniform float uRimPower;\n"
    " void main()\n"
    " {\n"
    "   vec3 N = normalize(vNormalView);\n"
    "   vec3 L = normalize(-uLightDirView);\n"
    "   vec3 V = normalize(-vViewPos);\n"
    "\n"
    "   // combinazione luce + vista (così il “dorso” verso la camera è più chiaro)\n"
    "   float dl = max(dot(N, L), 0.0);\n"
    "   float dv = max(dot(N, V), 0.0);\n"
    "   float t  = 0.55 * dl + 0.45 * dv; // [0,1]\n"
    "\n"
    "   // quantizzazione forte in 5 toni\n"
    "   float tone;\n"
    "   if      (t < 0.20) tone = 0.18;  // ombra\n"
    "   else if (t < 0.35) tone = 0.32;  // shadow-mid\n"
    "   else if (t < 0.55) tone = 0.50;  // mid\n"
    "   else if (t < 0.80) tone = 0.72;  // light\n"
    "   else               tone = 0.94;  // highlight\n"
    "\n"
    "   // base clay: mix tra scuro e chiaro, poi tinto da uBaseColor\n"
    "   vec3 clay = mix(vec3(0.3), vec3(1.0), tone);\n"
    "   vec3 color = clay * uBaseColor;\n"
    "\n"
    "   // spec “soft band” per un po’ di carattere automotive\n"
    "   vec3 H = normalize(L + V);\n"
    "   float spec = pow(max(dot(N, H), 0.0), 50.0);\n"
    "   float specBand = smoothstep(0.10, 0.16, spec);\n"
    "   color += vec3(0.35) * specBand; // leggera banda chiara\n"
    "\n"
    "   // rim chiaro sul bordo opposto alla luce\n"
    "   float rim = pow(1.0 - max(dot(N, V), 0.0), uRimPower) * uRimStrength;\n"
    "   color += vec3(rim);\n"
    "\n"
    "   // silhouette: linee nere spesse\n"
    "   float ndotv = abs(dot(N, V));\n"
    "   float edge = 1.0 - smoothstep(0.0, uEdgeWidth, ndotv);\n"
    "   color = mix(color, uLineColor, edge);\n"
    "\n"
    "   gl_FragColor = vec4(color, 1.0);\n"
    " }\n";
}




    GLuint compileShader(GLenum type,const char* src)
    {
        GLuint shader = glCreateShader(type);
        glShaderSource(shader,1,&src,nullptr);
        glCompileShader(shader);

        GLint ok = 0;
        glGetShaderiv(shader,GL_COMPILE_STATUS,&ok);
        if(!ok)
        {
            char log[1024];
            glGetShaderInfoLog(shader,1024,nullptr,log);
            std::cerr << "Shader compile error:\n" << log << std::endl;
            throw std::runtime_error("Shader compilation failed");
        }
        return shader;
    }

    void compile()
    {
        GLuint vs = compileShader(GL_VERTEX_SHADER,vertexShader());
        GLuint fs = compileShader(GL_FRAGMENT_SHADER,fragmentShader());

        program_ = glCreateProgram();
        glAttachShader(program_,vs);
        glAttachShader(program_,fs);
        glLinkProgram(program_);

        glDeleteShader(vs);
        glDeleteShader(fs);

        GLint ok = 0;
        glGetProgramiv(program_,GL_LINK_STATUS,&ok);
        if(!ok)
        {
            char log[1024];
            glGetProgramInfoLog(program_,1024,nullptr,log);
            std::cerr << "Link error:\n" << log << std::endl;
            throw std::runtime_error("Program linking failed");
        }
    }
};

