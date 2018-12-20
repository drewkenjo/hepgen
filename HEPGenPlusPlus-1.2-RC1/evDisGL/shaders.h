#ifndef SHADERS_HH
#define SHADERS_HH
#include <QString>

const QString vertexShaderSource= 
    "#version 330 core\n"
    "in highp vec4 vertex;\n"
    "in mediump vec4 color;\n"
    "uniform highp mat4 matrix;\n"
    "uniform highp mat4 worldMatrix;\n"
    "out mediump vec4 fragColor;\n"
    "void main(void)\n"
    "{\n"
    "   gl_Position = matrix * vertex;\n"
    "   fragColor = color;\n"
    "}";




const QString fragmentShaderSource=
   "#version 330 core\n"
   "in mediump vec4 fragColor;\n"
   "out mediump vec4 gl_FragColor;\n"
   "uniform lowp vec2 screenDim;\n"
   "void main(void)\n"
   "{\n"
   "   float x = gl_FragCoord.x / screenDim.x;\n"
   "   float y = gl_FragCoord.y / screenDim.y;\n"
   "   x -= 0.5; \n y-= 0.5;\n"
   "   float alpha = 1.2 * sqrt(x*x+y*y);\n"
   "   gl_FragColor = fragColor;\n"
   "   gl_FragColor.a -= alpha;\n"
   "}";
    //th


#endif

