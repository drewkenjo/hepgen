#ifndef vstructs
#define vstructs
#include <vector>


typedef struct {
    float posVertex1[3];
    float colorVertex[3];
    float posVertex2[3];
    int size;
    float width;
} line;

inline void setBlue(float* _in){
  _in[0]=0.0;
  _in[1]=0.0;
  _in[2]=1.0;
}

inline void setGreen(float* _in){
  _in[0]=0.0;
  _in[1]=1.0;
  _in[2]=0.0;
}

inline void setRed(float* _in){
  _in[0]=1.0;
  _in[1]=0.0;
  _in[2]=0.0;
}

inline void setYellow(float* _in){
  _in[0]=1.0;
  _in[1]=1.0;
  _in[2]=0.0;  
}
inline void setPurple(float* _in){
  _in[0]=139.0/255.0;
  _in[1]=0.0;
  _in[2]=139.0/255.0;  
}
inline void setTurk(float* _in){
  _in[0]=0.0f;
  _in[1]=0.0f;
  _in[2]=1.0f;  
}
inline void setOrange(float* _in){
  _in[0]=1.0f;
  _in[1]=0.5f;
  _in[2]=0.0f;  
}


#endif



