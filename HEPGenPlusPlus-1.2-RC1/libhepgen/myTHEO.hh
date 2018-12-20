#ifndef MYTHEO
#define MYTHEO

extern "C" {
#include <string.h>
#include <math.h>
#include <stdlib.h>
}
#include <iostream>

#ifndef SQERR
#define SQERR

using namespace std;


inline  double sq(double a) {
    return a*a;
}
inline void error(const char* c)
{
    cout<<c;
    exit(3);
}

#endif

class  Complexe {
    double re, im;
public:
    Complexe(double r=0.) {
        re=r;
        im=0.;
    }
    Complexe(double r, double i) {
        re=r;
        im=i;
    }
    ~Complexe() {
        ;
    }

    double imag() const {
        return im;
    }
    double real() const {
        return re;
    }
    double Norm() const {
        return sqrt(re*re+im*im);
    }
    double Norm2() const {
        return re*re+im*im;
    }
    double SignRe() const {
        return re/sqrt(re*re);
    }
    double SignIm() const {
        return im/sqrt(im*im);
    }
    double AbsRe() const {
        return re*re/sqrt(re*re);
    }
    double AbsIm() const {
        return im*im/sqrt(im*im);
    }

    Complexe& operator=(const Complexe&);
    int operator==(const Complexe& c) const {
        return (re==c.re && im==c.im);
    }

    Complexe operator+(const Complexe& b) const {
        return Complexe(re+b.re,im+b.im);
    }
    Complexe operator+(double b) const {
        return Complexe(re+b,im);
    }
    friend Complexe  operator+(const double&,const Complexe&);
    Complexe operator-(const Complexe& b) const {
        return Complexe(re-b.re,im-b.im);
    }
    Complexe operator-(double b) const {
        return Complexe(re-b,im);
    }
    friend Complexe  operator-(const double& ,const Complexe& );
    Complexe operator*(const Complexe& b)  const
    {
        return Complexe(re*b.re-im*b.im,im*b.re+re*b.im);
    }
    Complexe operator*(double b) const
    {
        return Complexe(re*b,im*b);
    }
    friend Complexe  operator*(const double&,const Complexe&);
    Complexe operator/(Complexe) const;
    Complexe operator/(double) const;
    friend Complexe  operator/(double,Complexe);
    Complexe sq() const  {
        return Complexe(re*re-im*im,im*re+re*im);
    }
    Complexe exp(int n) const  {
        Complexe a=1.;
        for (int j=0; j<n; j++) a=a*(*this);
        return a;
    }
    Complexe sqroot() const;
    Complexe Conj() const;
};
// fonctions associees
inline   ostream& operator<<(ostream&,const Complexe&);
inline  Complexe cos(const Complexe& a) {
    return Complexe(cos(a.real())*cosh(a.imag()) , -1.*sin(a.real())*sinh(a.imag()) );
}
inline  Complexe sin(const Complexe& a) {
    return Complexe(sin(a.real())*cosh(a.imag()) , cos(a.real())*sinh(a.imag()) );
}
inline Complexe exp(const Complexe& a) {
    return Complexe( exp(a.real())*cos(a.imag()) , exp(a.real())*sin(a.imag()) );
}
inline  Complexe lambda(const Complexe& a,const Complexe& b,const Complexe& c)
{
    return a.sq()+b.sq()+c.sq()-2.*a*b-2.*a*c-2.*b*c;
}

template<class T> class Matrix;

template<class T>
class Vector {
protected:
    T *v;
    int sz;
    int sh; // column 1, line 0
public:
    Vector(const T& a=T(0.),int s=4, int sc=1);
    Vector(const Vector&);
    Vector(const T&,const T&,const T&, int sc=1);
    Vector(const T&,const T&,const T&,const T&, int sc=1);
    ~Vector();

    T GetElem(int i) const {
        return v[i];
    }
    T GetSum() const {
        return v[0]+v[1]+v[3];
    }
    int GetSize() const {
        return sz;
    }
    int GetGeom() const {
        return sh;
    }
    void SetGeom(int i) {
        sh=i;
    }
    void SetElem(int i,const T& a) {
        v[i]=a;
    }

    Vector& operator=(const Vector&);
    int operator==(const Vector&) const;
    Vector operator+(const Vector&) const;
    Vector operator-(const Vector&) const;
    Vector operator*(const Matrix<T>&) const; // et son equiv. dans Matrix<T>
    Matrix<T> operator*(const Vector&) const; // vscal et vprod pour reste
    Vector operator*(const T&) const;  // et son equiv. T*V
    Vector operator/(const Complexe&) const;

    T vscal(const Vector&) const;
    Vector Conj() const;
    Vector Trans() const;
    Vector JRotate3(const Complexe&,const Complexe&) const;

};


template<class T>
class Matrix {
protected:
    T **m;
    int sz;
public:
    Matrix(const Matrix&);
    Matrix(const T& =T(0.),int size=4);
    Matrix(const char* const,int );
    ~Matrix();

    T GetElem(int i,int j) const {
        return m[i][j];
    }
    int GetSize() const {
        return sz;
    }
    void SetElem(int i,int j,const T& a) {
        m[i][j]=a;
    }

    Matrix& operator=(const Matrix&);
    int operator==(const Matrix&) const;

    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix sq() const {
        return this->operator*(*this);
    }
    Vector<T> operator*(const Vector<T>&) const;// et son equiv. dans Vector<T>
    Matrix operator/(const Complexe&) const;
    Matrix operator*(const T&) const; // et son equivalent inverse

    Matrix Conj() const;
    Matrix Trans() const;

    //    Complexe det();
    //    Matrix operator/(const Matrix&);
    //    friend Matrix operator/(const Complexe&,const Matrix&);
};


typedef Matrix<Complexe> CMatrix ;

class  Index
{
protected:
    char *ind;
    int co; // 1:covariant 0:contravariant
public:
    Index(const char* const c="ukwn",int ty=1);
    Index(const Index&);
    ~Index();

    Index& operator=(const Index&);
    int operator==(const Index&) const;

    const char* GetName() const {
        return ind;
    }
    int GetType() const {
        return co;
    }
    void SetName(const char* c) {
        strcpy(ind,c);
    }
    void SetType(int ty) {
        co=ty;
    }
    void SetIndex(const Index& a) {
        co=a.co;
        strcpy(ind,a.ind);
    }
};

template<class T> class Lorentz_Tensor;

template<class T>
class Lorentz_Vector:public Vector<T>,public Index
{
public:
    Lorentz_Vector();
    Lorentz_Vector(const Lorentz_Vector&);
    Lorentz_Vector(const Index&,const T&,const T&,const T&,const T&);
    Lorentz_Vector(const char*,int,const T&,const T&,const T&,const T&);
    Lorentz_Vector(const Index&,const Vector<T>&);
    Lorentz_Vector(const char*,int,const Vector<T>&);
    ~Lorentz_Vector() {}

    Index GetIndex() const {
        return Index(GetName(),GetType());
    }
    Vector<T> Get3Vector() const {
        return Vector<T>(Vector<T>::v[1],Vector<T>::v[2],Vector<T>::v[3]);
    }
    Vector<T> Get4Vector() const {
        return Vector<T>(Vector<T>::v[0],Vector<T>::v[1],Vector<T>::v[2],Vector<T>::v[3]);
    }

    T Mass2() const {
        return Vector<T>::v[0].sq()-Vector<T>::v[1].sq()-Vector<T>::v[2].sq()-Vector<T>::v[3].sq();
    }
    T GetE() const {
        return Vector<T>::v[0];
    }
    T GetX() const {
        return Vector<T>::v[1];
    }
    T GetY() const {
        return Vector<T>::v[2];
    }
    T GetZ() const {
        return Vector<T>::v[3];
    }
    T GetElem(int jj) const {
        return Vector<T>::v[jj];
    }
    void SetElem(int jj,const T& a) {
        Vector<T>::v[jj]=a;
    }

    Lorentz_Vector g(const char*,int,const char*,int) const;
    Lorentz_Vector Idxd(const char*,int) const;
    Lorentz_Vector Idxd(const Index&) const;


    Lorentz_Vector& operator=(const Lorentz_Vector&);
    Lorentz_Vector& operator=(const Vector<T>&);
    Lorentz_Vector& operator=(const Index&);
    int operator==(const Lorentz_Vector&) const;

    Lorentz_Vector operator+(const Lorentz_Vector&) const;
    Lorentz_Vector operator-(const Lorentz_Vector&) const;

    Lorentz_Vector operator*(const T&) const;// et son equiv. T*LV !! il faut un constructeur C->T
    Lorentz_Tensor<T> operator*(const Lorentz_Vector&) const;
    Lorentz_Vector operator/(const Complexe&) const;

    T C(const Lorentz_Vector&) const;
    Lorentz_Vector C(const Lorentz_Tensor<T>& b) const;

    Lorentz_Vector Conj() const
    {
        return Lorentz_Vector( ind,co,this->Vector<T>::Conj() );
    }
};

template<class T>
class Lorentz_Tensor:public Matrix<T>
{
protected:
    Index id1;
    Index id2;
public:
    Lorentz_Tensor();
    Lorentz_Tensor(const Lorentz_Tensor&);
    Lorentz_Tensor(const Index&,const Index&,const Matrix<T>&);
    Lorentz_Tensor(const Lorentz_Vector<T>&,const Lorentz_Vector<T>&);
    ~Lorentz_Tensor() {}

    Index GetIndex1() const {
        return Index(id1);
    }
    Index GetIndex2() const {
        return Index(id2);
    }
    void SetIndex1(const Index& a) {
        id1=a;
    }
    void SetIndex2(const Index& a) {
        id2=a;
    }
    T GetElem(int jj,int kk) const {
        return Matrix<T>::m[jj][kk];
    }
    void SetElem(int jj,int kk,const T& a) {
        Matrix<T>::m[jj][kk]=a;
    }

    Lorentz_Tensor& operator=(const Lorentz_Tensor&);
    Lorentz_Tensor& operator=(const Matrix<T>&);
    int operator==(const Lorentz_Tensor&) const;

    Lorentz_Tensor operator+(const Lorentz_Tensor&) const;
    Lorentz_Tensor operator-(const Lorentz_Tensor&) const;
    Lorentz_Tensor operator*(const T&) const;// et son equiv. T*LT ! il faut un constructeur C->T
    Lorentz_Tensor operator/(const Complexe&) const;

    Lorentz_Tensor C(const Lorentz_Tensor&) const;
    Lorentz_Vector<T> C(const Lorentz_Vector<T>&) const; // et son equiv. ds LV
    Lorentz_Tensor Trans() const {
        return Lorentz_Tensor(id2,id1,this->Matrix<T>::Trans());
    }
};


class  C_Lorentz_Vector;
class  CM_Lorentz_Vector;
class  C_Lorentz_Tensor;
class  CM_Lorentz_Tensor;
class  Spinor;

class  C_Lorentz_Vector:public Lorentz_Vector<Complexe>
{
public:
    C_Lorentz_Vector();
    C_Lorentz_Vector(const Lorentz_Vector<Complexe>&);
    C_Lorentz_Vector(const char*,int,const Complexe&,const Complexe&,const Complexe&,const Complexe&);
    C_Lorentz_Vector(const Index&,const Complexe&,const Complexe&,const Complexe&,const Complexe&);
    C_Lorentz_Vector(const Index&,const Vector<Complexe>&);
    C_Lorentz_Vector(const char*,int,const Vector<Complexe>&);
    C_Lorentz_Vector(const Complexe&,const char*,int,const Complexe&,const Complexe&,const Complexe&);
    C_Lorentz_Vector(const Complexe&,const char*,int,const Vector<Complexe>&);
    C_Lorentz_Vector(const Spinor&,const CM_Lorentz_Vector&, const Spinor&);
    ~C_Lorentz_Vector() {}

    C_Lorentz_Vector g(const char* n1,int t1,const char* n2,int t2) const
    {
        return C_Lorentz_Vector(Lorentz_Vector<Complexe>::g(n1,t1,n2,t2));
    }
    C_Lorentz_Vector Idxd(const char* n1,int t1) const {
        return C_Lorentz_Vector(Lorentz_Vector<Complexe>::Idxd(n1,t1));
    }
    C_Lorentz_Vector Idxd(const Index& i1) const {
        return C_Lorentz_Vector(Lorentz_Vector<Complexe>::Idxd(i1));
    }

    Complexe Mass() const {
        return (v[0].sq()-v[1].sq()-v[2].sq()-v[3].sq()).sqroot();
    }
    Complexe GetP() const {
        return (v[1].sq()+v[2].sq()+v[3].sq()).sqroot();
    }
    Complexe GetTheta() const;
    Complexe GetPhi() const;

    CM_Lorentz_Vector operator*(const CMatrix&) const; // et equiv.
    C_Lorentz_Vector operator*(const Complexe&) const; // et equiv.
    C_Lorentz_Vector operator*(double) const; // et equiv.
    C_Lorentz_Vector operator/(const Complexe& c) const
    {
        return C_Lorentz_Vector( this->Lorentz_Vector<Complexe>::operator/(c) );
    }
    C_Lorentz_Vector operator/(double c) const
    {
        return C_Lorentz_Vector( this->Lorentz_Vector<Complexe>::operator/(c) );
    }

    C_Lorentz_Vector operator+(const C_Lorentz_Vector& CLV) const
    {
        return C_Lorentz_Vector(this->Lorentz_Vector<Complexe>::operator+(CLV));
    }
    C_Lorentz_Vector operator-(const C_Lorentz_Vector& CLV) const
    {
        return C_Lorentz_Vector(this->Lorentz_Vector<Complexe>::operator-(CLV));
    }

    CM_Lorentz_Vector operator+(const CM_Lorentz_Vector&);
    CM_Lorentz_Vector operator-(const CM_Lorentz_Vector&);


// les contractions:
    Complexe C(const C_Lorentz_Vector&) const;
    CMatrix C(const CM_Lorentz_Vector&) const;
    C_Lorentz_Vector C(const C_Lorentz_Tensor&) const;
    CM_Lorentz_Vector C(const CM_Lorentz_Tensor&) const;

    CMatrix Slash() const;

    CM_Lorentz_Vector Id() const;
};

class  CM_Lorentz_Vector:public Lorentz_Vector<CMatrix>
{
public:
    CM_Lorentz_Vector();
    CM_Lorentz_Vector(const Lorentz_Vector<CMatrix>&);
    //CM_Lorentz_Vector(const C_Lorentz_Vector&);
    CM_Lorentz_Vector(const char*,int,const CMatrix&,const CMatrix&,const CMatrix&,const CMatrix&);
    CM_Lorentz_Vector(const Index&,const CMatrix&,const CMatrix&,const CMatrix&,const CMatrix&);
    CM_Lorentz_Vector(const Index&,const Vector<CMatrix>&);
    CM_Lorentz_Vector(const char*,int,const Vector<CMatrix>&);
    ~CM_Lorentz_Vector() {}

    CM_Lorentz_Vector g(const char* n1,int t1,const char* n2,int t2) const
    {
        return CM_Lorentz_Vector(Lorentz_Vector<CMatrix>::g(n1,t1,n2,t2));
    }
    CM_Lorentz_Vector Idxd(const char* n1,int t1) const
    {
        return CM_Lorentz_Vector(this->Lorentz_Vector<CMatrix>::Idxd(n1,t1));
    }
    CM_Lorentz_Vector Idxd(const Index& i1) const
    {
        return CM_Lorentz_Vector(this->Lorentz_Vector<CMatrix>::Idxd(i1));
    }

    CM_Lorentz_Vector operator*(const CMatrix&) const; // et equiv.
    CM_Lorentz_Vector operator*(const Complexe&) const; // et equiv.
    CM_Lorentz_Vector operator*(double) const; // et equiv.
    CM_Lorentz_Vector operator/(const Complexe& c) const
    {
        return CM_Lorentz_Vector( this->Lorentz_Vector<CMatrix>::operator/(c) );
    }
    CM_Lorentz_Vector operator/(double c) const
    {
        return CM_Lorentz_Vector( this->Lorentz_Vector<CMatrix>::operator/(c) );
    }

    CM_Lorentz_Vector operator+(const CM_Lorentz_Vector& CMLV) const
    {
        return CM_Lorentz_Vector(this->Lorentz_Vector<CMatrix>::operator+(CMLV));
    }
    CM_Lorentz_Vector operator-(const CM_Lorentz_Vector& CMLV) const
    {
        return CM_Lorentz_Vector(this->Lorentz_Vector<CMatrix>::operator-(CMLV));
    }

    CM_Lorentz_Vector operator+(const C_Lorentz_Vector&);
    CM_Lorentz_Vector operator-(const C_Lorentz_Vector&);

// les contractions:
    CMatrix C(const CM_Lorentz_Vector&) const;
    CMatrix C(const C_Lorentz_Vector&) const;
    CM_Lorentz_Vector C(const CM_Lorentz_Tensor&) const;
    CM_Lorentz_Vector C(const C_Lorentz_Tensor&) const;

};

typedef C_Lorentz_Vector Current ;

// fonctions associees aux Lorentz_Vector :

inline C_Lorentz_Vector  operator*(const Complexe& a, const C_Lorentz_Vector& CLV)
{
    return C_Lorentz_Vector(CLV*a);
}
inline C_Lorentz_Vector  operator*(double a, const C_Lorentz_Vector& CLV)
{
    return C_Lorentz_Vector(CLV*Complexe(a));
}
inline CM_Lorentz_Vector  operator*(const CMatrix& M, const C_Lorentz_Vector& CLV)
{
    return CM_Lorentz_Vector(CLV*M);
}
inline CM_Lorentz_Vector  operator*(double a, const CM_Lorentz_Vector& CMLV)
{
    return CM_Lorentz_Vector(CMLV*Complexe(a));
}
inline CM_Lorentz_Vector  operator*(const Complexe& a, const CM_Lorentz_Vector& CMLV)
{
    return CM_Lorentz_Vector(CMLV*a);
}
inline CM_Lorentz_Vector  operator*(const CMatrix& , const CM_Lorentz_Vector& );
inline CM_Lorentz_Vector  GammaDirac(const char*,int);
CMatrix  Gamma5();

class  C_Lorentz_Tensor:public Lorentz_Tensor<Complexe>
{
public:
    C_Lorentz_Tensor();
    C_Lorentz_Tensor(const Lorentz_Tensor<Complexe>&);
    C_Lorentz_Tensor(const C_Lorentz_Vector&,const C_Lorentz_Vector&);
    C_Lorentz_Tensor(const char*,int,const char*,int);
    C_Lorentz_Tensor(const char*,int,const char*,int,const Complexe&,const Complexe&);
    C_Lorentz_Tensor(const Spinor&,const CM_Lorentz_Tensor&, const Spinor&);
    ~C_Lorentz_Tensor() {}


    CM_Lorentz_Tensor operator*(const CMatrix&) const; // et equiv.
    C_Lorentz_Tensor operator*(const Complexe&) const; // et equiv.
    C_Lorentz_Tensor operator*(double) const; // et equiv.
    C_Lorentz_Tensor operator/(const Complexe& c) const
    {
        return C_Lorentz_Tensor( this->Lorentz_Tensor<Complexe>::operator/(c) );
    }
    C_Lorentz_Tensor operator/(double d) const
    {
        return C_Lorentz_Tensor( this->Lorentz_Tensor<Complexe>::operator/(Complexe(d)) );
    }

    C_Lorentz_Tensor operator+(const C_Lorentz_Tensor& CLT) const
    {
        return C_Lorentz_Tensor(this->Lorentz_Tensor<Complexe>::operator+(CLT));
    }
    C_Lorentz_Tensor operator-(const C_Lorentz_Tensor& CLT) const
    {
        return C_Lorentz_Tensor(this->Lorentz_Tensor<Complexe>::operator-(CLT));
    }

    CM_Lorentz_Tensor operator+(const CM_Lorentz_Tensor&);
    CM_Lorentz_Tensor operator-(const CM_Lorentz_Tensor&);

// les contractions:
    C_Lorentz_Tensor C(const C_Lorentz_Tensor&) const;
    CM_Lorentz_Tensor C(const CM_Lorentz_Tensor&) const;
    C_Lorentz_Vector C(const C_Lorentz_Vector&) const;
    CM_Lorentz_Vector C(const CM_Lorentz_Vector&) const;

    friend C_Lorentz_Tensor  gT(const char*,int,const char*,int);
};


void  epsil(double e[4][4][4][4]);
C_Lorentz_Tensor   LeviCita(const char*,int,const char*,int,const C_Lorentz_Tensor&);
C_Lorentz_Tensor   TLorentz(const C_Lorentz_Vector&,const char*,int,const char*,int,const char* dir="direct");
C_Lorentz_Tensor  gT(const char*,int,const char*,int);

typedef  C_Lorentz_Tensor RotTens;

class  CM_Lorentz_Tensor:public Lorentz_Tensor<CMatrix>
{
public:
    CM_Lorentz_Tensor();
    CM_Lorentz_Tensor(const Lorentz_Tensor<CMatrix>&);
    //CM_Lorentz_Tensor(const C_Lorentz_Tensor&);
    CM_Lorentz_Tensor(const CM_Lorentz_Vector&,const CM_Lorentz_Vector&);
    CM_Lorentz_Tensor(const char*,int,const char*,int);
    CM_Lorentz_Tensor(const char*,int,const char*,int,const char*);
    ~CM_Lorentz_Tensor() {}

    CM_Lorentz_Tensor operator*(const CMatrix&) const; // et equiv.
    CM_Lorentz_Tensor operator*(const Complexe&) const; // et equiv.
    CM_Lorentz_Tensor operator*(double) const; // et equiv.

    CM_Lorentz_Tensor operator+(const CM_Lorentz_Tensor& CMLV) const
    {
        return CM_Lorentz_Tensor(this->Lorentz_Tensor<CMatrix>::operator+(CMLV));
    }
    CM_Lorentz_Tensor operator-(const CM_Lorentz_Tensor& CMLV) const
    {
        return CM_Lorentz_Tensor(this->Lorentz_Tensor<CMatrix>::operator-(CMLV));
    }
    CM_Lorentz_Tensor operator+(const C_Lorentz_Tensor&);
    CM_Lorentz_Tensor operator-(const C_Lorentz_Tensor&);

// les contractions:
    CM_Lorentz_Tensor C(const CM_Lorentz_Tensor&) const;
    CM_Lorentz_Tensor C(const C_Lorentz_Tensor&) const;
    CM_Lorentz_Vector C(const CM_Lorentz_Vector&) const;
    CM_Lorentz_Vector C(const C_Lorentz_Vector&) const;

};
// fonctions associees aux Lorentz_Tensor
inline C_Lorentz_Tensor  operator*(const C_Lorentz_Vector& CLV1,const C_Lorentz_Vector& CLV2)
{
    return C_Lorentz_Tensor( CLV1.Lorentz_Vector<Complexe>::operator*(CLV2) );
}
inline CM_Lorentz_Tensor  operator*(const CM_Lorentz_Vector& CMLV1,const CM_Lorentz_Vector& CMLV2)
{
    return CM_Lorentz_Tensor(CMLV1.Lorentz_Vector<CMatrix>::operator*(CMLV2));
}

C_Lorentz_Tensor  operator*(double a, const C_Lorentz_Tensor& CLT);
C_Lorentz_Tensor  operator*(const Complexe& a, const C_Lorentz_Tensor& CLT);
CM_Lorentz_Tensor  operator*(const CMatrix& a, const C_Lorentz_Tensor& b);

CM_Lorentz_Tensor  operator*(double a, const CM_Lorentz_Tensor& b);
CM_Lorentz_Tensor  operator*(const Complexe& a, const CM_Lorentz_Tensor& b);
CM_Lorentz_Tensor  operator*(const CMatrix& a, const CM_Lorentz_Tensor& b);

CM_Lorentz_Tensor  TSigmaDirac(const char*,int,const char*,int);

enum  helicity {neg,pos};
/* enum helicity {zero,neg,pos,fin}; */
enum  particule {part,anti};
enum  invar {heli,spin};

class  Spinor:public Vector<Complexe>
{
private:
    helicity h;
    particule pa;
    invar inv;
public:
    Spinor();
    Spinor(helicity,particule,invar,const Vector<Complexe>& );
    Spinor(double);
    Spinor(particule,const C_Lorentz_Vector&,helicity,Complexe ma,invar inv=spin);
    Spinor(particule,const C_Lorentz_Vector&,double,Complexe ma,invar inv=spin);
    ~Spinor() {
        ;
    }

    helicity GetHelicity() const {
        return h;
    }
    particule GetParticule() const {
        return pa;
    }
    invar GetInvar() const {
        return inv;
    }
    Spinor bar() const;

    Spinor operator*(Complexe a) {
        return Spinor(h,pa,inv,Vector<Complexe>::operator*(a));
    }
    Spinor operator+(const Spinor& s) {
        return Spinor(h,pa,inv,Vector<Complexe>::operator+(s));
    }
    C_Lorentz_Vector mult(const CM_Lorentz_Vector&,const Spinor&);
    C_Lorentz_Tensor mult(const CM_Lorentz_Tensor&,const Spinor&);

    friend ostream& operator<<(ostream&,const Spinor&);
};

class  Dspinor:public Lorentz_Vector<Spinor>
{
public:
    Dspinor();
    Dspinor(const Lorentz_Vector<Spinor>&);
    Dspinor(const char*,int,const C_Lorentz_Vector&, double, Complexe);
    ~Dspinor() {
        ;
    }

    Dspinor bar();

    Spinor C(const C_Lorentz_Vector&);

};



//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                          TEMPLATE VECTOR
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


template<class T>
Vector<T>::Vector(const T& a,int s, int sc)
{
    sz=s;
    sh=sc;
    v = new T[sz];
    for ( int ii=0; ii<sz; ii++ ) v[ii]=a;
}

template<class T>
Vector<T>::Vector(const Vector<T>& b)
{
    sz=b.sz;
    sh=b.sh;
    v = new T[sz];
    for (int i=0; i<sz; i++) v[i]=b.v[i];
}

template<class T>
Vector<T>::Vector(const T& a,const T& b,const T& c, int sc)
{
    sz=3;
    sh=sc;
    v = new T[3];
    v[0]=a;
    v[1]=b;
    v[2]=c;
}

template<class T>
Vector<T>::Vector(const T& a,const T& b,const T& c,const T& d, int sc)
{
    sz=4;
    sh=sc;
    v = new T[4];
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
}

template<class T>
Vector<T>::~Vector()
{
    if (v!=NULL) delete[] v;
}

template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& b)
{
    sz=b.sz;
    sh=b.sh;
    if (v!=NULL) delete[] v;
    v = new T[sz];
    for (int i=0; i<sz; i++) v[i]=b.v[i];
    return *this;
}

template<class T>
int Vector<T>::operator==(const Vector<T>& b) const
{
    int comp=0;
    if ( sz==b.sz && sh==b.sh) {
        int j=0;
        while (j<sz && v[j]==b.v[j]) {
            j++;
        }
        if (j==sz) comp=1;
    }
    return comp;
}

template<class T>
Vector<T> Vector<T>::operator+(const Vector<T>& b) const
{
    Vector<T> c(T(0.),sz,sh);

    if (sz!=b.sz) error("Bad Vector<T> size for addition\n");
    if (sh!=b.sh) error("Bad Vector<T> geometry for addition\n");

    for (int i=0; i<sz; i++) c.v[i]=v[i]+b.v[i];
    return c;
}

template<class T>
Vector<T> Vector<T>::operator-(const Vector<T>& b) const
{
    Vector<T> c(T(0.),sz,sh);

    if (sz!=b.sz) error("Bad Vector<T> size for substraction\n");
    if (sh!=b.sh) error("Bad Vector<T> geometry for substraction\n");

    for (int i=0; i<sz; i++) c.v[i]=v[i]-b.v[i];
    return c;
}

template<class T>
Vector<T> Vector<T>::operator*(const T& b) const
{
    Vector<T> c(T(0.),sz,sh);

    for (int i=0; i<sz; i++) c.v[i]=b*v[i];
    return c;
}

template<class T>
Vector<T> operator*(const T& b,Vector<T> a)
{
    Vector<T> c(T(0.),a.GetSize(),a.GetGeom());
    for (int i=0; i<a.GetSize(); i++) c.SetElem(i,b*a.GetElem(i));
    return c;
}

template<class T>
Vector<T> operator*(double b,Vector<T> a)
{
    Vector<T> c(T(0.),a.GetSize(),a.GetGeom());
    for (int i=0; i<a.GetSize(); i++) c.SetElem(i,b*a.GetElem(i));
    return c;
}

template<class T>
Vector<T> Vector<T>::operator*(const Matrix<T>& ma) const
{
    if (sz!=ma.GetSize()) error("Bad Vector<T> size for multiplication befor Matrix\n");
    if (sh!=0) error("Bad Vector<T> Type for multiplication befor Matrix\n");
    Vector<T> mprod(T(0.),sz,0);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            mprod.v[i]=mprod.v[i]+v[j]*ma.GetElem(j,i);
    return mprod;
}

template<class T>
Matrix<T> Vector<T>::operator*(const Vector<T>& va) const
{
    if (sz!=va.GetSize()) error("Bad Vector<T> size for multiplication befor Vector\n");
    if (sh!=1 || va.GetGeom()!=0) error("Bad Vector<T> Type for multiplication befor Vector\n");
    Matrix<T> mprod(T(0.),sz);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            mprod.SetElem(i,j,v[i]*va.v[j]);
    return mprod;
}

template<class T>
Vector<T> Vector<T>::operator/(const Complexe& b) const
{
    Vector<T> c(T(0.),sz,sh);
    if ( b==Complexe(0.) ) error("Division of a Vector<T> by Complexe zero\n");

    for (int i=0; i<sz; i++) c.v[i]=v[i]/b;
    return c;
}

template<class T>
Vector<T> Vector<T>::Trans() const
{
    Vector<T> b=*this;
    b.SetGeom( ! this->GetGeom() );
    return b;
}

template<class T>
Vector<T> Vector<T>::Conj() const
{
    Vector<T> b(T(0.),this->GetSize(),this->GetGeom());
    for (int i=0; i<this->GetSize(); i++) b.SetElem( i , (this->v[i]).Conj() );
    return b;
}

template<class T>
T Vector<T>::vscal(const Vector<T>& va) const
{
    T sum=v[0];
    sum=sum-sum;
    for (int i=0; i<sz; i++) sum = sum + v[i]*va.v[i];
    return sum;
}

template<class T>
Vector<T> Vector<T>::JRotate3(const Complexe& t,const Complexe& p) const
{
//rotation of a vector in a fixed frame

    Vector<T> vv=*this;
    vv.v[0]=(cos(t.real())*cos(p.real())*v[0]-sin(p.real())*v[1]
             +sin(t.real())*cos(p.real())*v[2]);
    vv.v[1]=(sin(p.real())*cos(t.real())*v[0]+cos(p.real())*v[1]
             +sin(p.real())*sin(t.real())*v[2]);
    vv.v[2]=-sin(t.real())*v[0]+cos(t.real())*v[2];
    return vv;
}



/* template<class T> */
/* ostream& operator<<(ostream& f,Vector<T> vec) */
/* { */
/*  if (vec.GetGeom()==0) f<<"line"; */
/*  else f<<"column"<<"\n"; */
/*  for (int i=0;i<vec.GetSize();i++)  f<<"\n"<<vec.GetElem(i); */
/*  f<<"\n"; */
/*  f.clear(ios::goodbit); */
/*  return f; */
/* } */

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                           TEMPLATE  Matrix
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

template<class T>
Matrix<T>::Matrix(const T& c,int size)
{
    if (size<1.) error("Bad Matrix construction size");
    sz=size;
    m = new T*[sz];
    for (int i=0; i<sz; i++) {
        m[i] = new T[sz];
        for (int j=0; j<sz; j++)
            m[i][j]=T(0.);
    }
    for (int i=0; i<sz; i++) {
        m[i][i]=c;
    }
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& ma)
{
    sz=ma.sz;
    m = new T*[sz];
    for (int i=0; i<sz; i++) {
        m[i] = new T[sz];
        for (int j=0; j<sz; j++)
            m[i][j]=ma.m[i][j];
    }
}

template<class T>
Matrix<T>::Matrix(const char* const n,int t)
{
    sz=4;
    m = new T*[sz];
    for (int i=0; i<sz; i++) {
        m[i] = new T[sz];
        for (int j=0; j<sz; j++)
            m[i][j]=T(0.);
    }
    switch(t) {
    case 0:
        m[0][0]=T(1.);
        m[1][1]=T(1.);
        m[2][2]=T(-1.);
        m[3][3]=T(-1.);
        break;
    case 1:
        m[0][3]=T(1.);
        m[1][2]=T(1.);
        m[2][1]=T(-1.);
        m[3][0]=T(-1.);
        break;
    case 2:
        m[0][3]=Complexe(0.,-1.);
        m[1][2]=Complexe(0.,1.);
        m[2][1]=Complexe(0.,1.);
        m[3][0]=Complexe(0.,-1.);
        break;
    case 3:
        m[0][2]=T(1.);
        m[1][3]=T(-1.);
        m[2][0]=T(-1.);
        m[3][1]=T(1.);
        break;
    case 5:
        m[0][2]=T(1.);
        m[1][3]=T(1.);
        m[2][0]=T(1.);
        m[3][1]=T(1.);
        break;
    }
}

template<class T>
Matrix<T>::~Matrix()
{
    // cout<<"destruction Matrix";
    if (m!=NULL) {
        for(int i=0; i<sz; i++) delete[] m[i];
        delete[] m;
    }
// cout<<"\t finie\n";
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& ma)
{
    if (m!=NULL) {
        for(int i=0; i<sz; i++) delete[] m[i];
        delete[] m;
    }
    sz=ma.sz;
    m = new T*[sz];
    for (int i=0; i<sz; i++) {
        m[i] = new T[sz];
        for (int j=0; j<sz; j++)
            m[i][j]=ma.m[i][j];
    }
    return *this;
}

template<class T>
int Matrix<T>::operator==(const Matrix<T>& ma) const
{
    int comp=0;
    int i=0,j=0;
    if (sz==ma.sz) {
        while ( m[i][j]==ma.m[i][j] && i<sz ) {
            while ( m[i][j]==ma.m[i][j] && j<sz ) {
                j++;
            }
            i++;
        }
        if (i==sz && j==sz) comp=1;
    }
    return comp;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& ma) const
{
    if (sz!=ma.sz) error("Bad Matrix<T> size for addition");
    Matrix<T> msom=ma;
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            msom.m[i][j]=m[i][j]+ma.m[i][j];
    return msom;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& ma) const
{
    if (sz!=ma.sz) error("Bad Matrix<T> size for substraction");
    Matrix<T> msom=ma;
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            msom.m[i][j]=m[i][j]-ma.m[i][j];
    return msom;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& ma) const
{
    if (sz!=ma.sz) error("Bad Matrix<T> size for multiplication\n");
    Matrix<T> mprod(T(0.),sz);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            for (int k=0; k<sz; k++)
                mprod.m[i][j]=mprod.m[i][j]+m[i][k]*ma.m[k][j];
    return mprod;
}

template<class T>
Vector<T> Matrix<T>::operator*(const Vector<T>& va) const
{
    if (sz!=va.GetSize()) error("Bad Vector<T> size for multiplication after Matrix\n");
    if (va.GetGeom()!=1) error("Bad Vector<T> Type for multiplication after Matrix\n");
    Vector<T> mprod(T(0.),sz,1);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            mprod.SetElem(i,mprod.GetElem(i)+m[i][j]*va.GetElem(j));
    return mprod;
}


template<class T>
Matrix<T> Matrix<T>::operator*(const T& a) const
{
    Matrix<T> mprod(T(0.),sz);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            mprod.m[i][j]=a*m[i][j];
    return mprod;
}

template<class T>
Matrix<T> operator*(const T& a,const Matrix<T>& ma)
{
    Matrix<T> mprod(T(0.),ma.GetSize());
    int ssz=ma.GetSize();
    for (int i=0; i<ssz; i++)
        for (int j=0; j<ssz; j++)
            mprod.SetElem(i,j,a*ma.GetElem(i,j));
    return mprod;
}

template<class T>
Matrix<T> operator*(double a,const Matrix<T>& ma)
{
    return ma*Complexe(a);
}

template<class T>
Matrix<T> Matrix<T>::operator/(const Complexe& a) const
{
    if ( a==Complexe(0.) ) error("Division of Matrix<T> by Complexe zero");
    Matrix<T> mprod(T(0.),sz);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            mprod.m[i][j]=m[i][j]/a;
    return mprod;
}


template<class T>
Matrix<T> Matrix<T>::Conj() const
{
    Matrix<T> msom(T(0.),sz);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            msom.m[i][j]=(m[i][j]).Conj();
    return msom;
}

template<class T>
Matrix<T> Matrix<T>::Trans() const
{
    Matrix<T> msom(T(0.),sz);
    for (int i=0; i<sz; i++)
        for (int j=0; j<sz; j++)
            msom.m[i][j]=m[j][i];
    return msom;
}

//    Complexe det();
//    Matrix operator/(const Matrix&);
//    friend Matrix operator/(const Complexe&,const Matrix&);

template<class T>
ostream& operator<<(ostream& f,Matrix<T> mat)
{
    int siz=mat.GetSize();
    for (int i=0; i<siz; i++)
        for (int j=0; j<siz; j++) {
            if (j<siz-1) f<<"\t"<<mat.GetElem(i,j);
            else f<<"\t"<<mat.GetElem(i,j)<<"\n";
        }
    f.clear(ios::goodbit);
    return f;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                     TEMPLATE LORENTZ _ VECTOR
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

template<class T>
Lorentz_Vector<T>::Lorentz_Vector()
{}

template<class T>
Lorentz_Vector<T>::Lorentz_Vector(const Lorentz_Vector<T>& b)
{
    for (int i=0; i<Vector<T>::sz; i++) Vector<T>::v[i]=b.Vector<T>::v[i];
    co=b.co;
    strcpy(ind,b.ind);
}

template<class T>
Lorentz_Vector<T>::Lorentz_Vector(const Index& ty,const T& a,const T& b,const T& c,const T& d):Vector<T>(a,b,c,d),Index(ty)
{}

template<class T>
Lorentz_Vector<T>::Lorentz_Vector(const char *na,int ty,const T& a,const T& b,const T& c,const T& d):Vector<T>(a,b,c,d),Index(na,ty)
{}

template<class T>
Lorentz_Vector<T>::Lorentz_Vector(const char *na,int ty,const Vector<T>& a):Vector<T>(a),Index(na,ty)
{}

template<class T>
Lorentz_Vector<T>::Lorentz_Vector(const Index& ty,const Vector<T>& a):Vector<T>(a),Index(ty)
{}

template<class T>
Lorentz_Vector<T>& Lorentz_Vector<T>::operator=(const Lorentz_Vector<T>& b)
{
    Vector<T>::sz=b.Vector<T>::sz;
    Vector<T>::sh=b.Vector<T>::sh;
    if (Vector<T>::v!=NULL) delete[] Vector<T>::v;
    Vector<T>::v = new T[Vector<T>::sz];
    for (int i=0; i<Vector<T>::sz; i++) Vector<T>::v[i]=b.Vector<T>::v[i];
    co=b.co;
    if (ind!=NULL) delete[] ind;
    ind = new char[10];
    int i=0;
    if ( *(b.ind)=='\0' ) *ind='\0';
    else {
        do {
            ind[i]=b.ind[i];
            i++;
        }
        while ( b.ind[i]!='\0' );
        ind[i]='\0';
    }
    return *this;
}

template<class T>
Lorentz_Vector<T>& Lorentz_Vector<T>::operator=(const Vector<T>& b)
{
    Vector<T>::operator=(b);
    return *this;
}

template<class T>
Lorentz_Vector<T>& Lorentz_Vector<T>::operator=(const Index& b)
{
    Index::operator=(b);
    return *this;
}


template<class T>
int Lorentz_Vector<T>::operator==(const Lorentz_Vector<T>& b) const
{
    int c=( this->Index::operator==(b.GetIndex()) && this->Vector<T>::operator==(b.Get4Vector()) )? 1 : 0 ;
    return c;
}


template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::g(const char* n1,int i1,const char* n2,int i2) const
{
    if (strcmp(n1,GetName())!=0) error("Bad name for contraction with g()");
    if (i1==GetType()) error("Bad index type for contraction with g()");
    Lorentz_Vector<T> c=*this;
    for (int i=1; i<4; i++) c.v[i]=(i1==i2)? Complexe(-1.) * this->v[i] : this->v[i] ;
    c.SetType(i2);
    c.SetName(n2);
    return c;
}

template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::Idxd(const char* n1,int i1) const
{
    if ( strcmp(n1,GetName())==0 && i1==GetType() ) return *this;
    else return this->g(GetName(),!GetType(),n1,i1);
}

template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::Idxd(const Index& idx) const
{
    if ( strcmp(idx.GetName(),GetName())==0 && idx.GetType()==GetType() ) return *this;
    else return this->g(GetName(),!GetType(),idx.GetName(),idx.GetType());
}


template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::operator+(const Lorentz_Vector<T>& b) const
{
    if ( ! this->Index::operator==(b) ) error("Bad Index in Lorentz_Vector addition\n");

    Lorentz_Vector<T> c=b;
    c.Vector<T>::operator=(this->Vector<T>::operator+(b));
    return c;
}

template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::operator-(const Lorentz_Vector<T>& b) const
{
    if ( ! this->Index::operator==(b) ) error("Bad Index in Lorentz_Vector substraction\n");

    Lorentz_Vector<T> c=b;
    c.Vector<T>::operator=(this->Vector<T>::operator-(b));
    return c;
}

template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::operator*(const T& b) const
{
    return Lorentz_Vector<T>(this->GetIndex(),this->Vector<T>::operator*(b));
}

// multiplication Matrixielle non commutative !!
template<class T>
Lorentz_Vector<T> operator*(const T& b,const Lorentz_Vector<T>& a)
{
    return Lorentz_Vector<T>(a.GetIndex(),b*a.Get4Vector());
}


template<class T>
Lorentz_Tensor<T> Lorentz_Vector<T>::operator*(const Lorentz_Vector<T>& b) const
{
    Lorentz_Tensor<T> r;
    r.SetIndex1(GetIndex());
    r.SetIndex2(b.GetIndex());
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            r.SetElem(i,j,GetElem(i)*b.GetElem(j));
    return r;
}

template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::operator/(const Complexe& b) const
{
    if ( b==0.) error("Division of Lorentz_Vector<T> by Complexe zero");
    Vector<T> c=Vector<T>(Vector<T>::v[0],Vector<T>::v[1],Vector<T>::v[2],Vector<T>::v[3]);
    c=c/b;
    return Lorentz_Vector<T>(this->GetIndex(),c);
}

template<class T>
T Lorentz_Vector<T>::C(const Lorentz_Vector<T>& b) const
{
    if ( strcmp(this->GetName(),b.GetName())!=0 ) error("Bad Index name in Lorentz_Vector contraction\n");
    if ( this->GetType()==b.GetType() ) error("Bad Index type in Lorentz_Vector contraction\n");
    T c;
    for (int i=0; i<4; i++) c=c+this->v[i]*b.v[i];
    return c;
}

template<class T>
Lorentz_Vector<T> Lorentz_Vector<T>::C(const Lorentz_Tensor<T>& b) const
{
    if ( strcmp(this->GetName(),b.GetIndex1().GetName())!=0 ) error("Bad Index name in Lorentz_Vector contraction\n");
    if ( this->GetType()==b.GetIndex1().GetType() ) error("Bad Index type in Lorentz_Vector contraction\n");
    Lorentz_Vector<T> c;
    c.SetIndex(b.GetIndex2());
    T sum=T(0.);
    for (int i=0; i<4; i++) {
        sum=T(0.);
        for (int k=0; k<4; k++) sum=sum+this->v[k]*b.GetElem(k,i);
        c.v[i]=sum;
    }
    return c;
}

template<class T>
ostream& operator<<(ostream& f,Lorentz_Vector<T> vec)
{

    f<<"\n\t\tLorentz_vector:\nIndex:\t"<<vec.GetName()<<"\t"<<vec.GetType()<<"\n";
    for (int i=0; i<vec.GetSize(); i++)  f<<"\n"<<vec.GetElem(i);
    f<<"\n\t\tFin Lorentz_Vector.\n";
    f.clear(ios::goodbit);
    return f;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                      TEMPLATE LORENTZ _ TENSOR
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

template<class T>
Lorentz_Tensor<T>::Lorentz_Tensor()
{}

template<class T>
Lorentz_Tensor<T>::Lorentz_Tensor(const Lorentz_Tensor<T>& b)
{

    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Matrix<T>::m[i][j]=b.GetElem(i,j);

    id1=b.id1;
    id2=b.id2;
}

template<class T>
Lorentz_Tensor<T>::Lorentz_Tensor(const Index& i1,const Index& i2,const Matrix<T>& b)
{

    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Matrix<T>::m[i][j]=b.GetElem(i,j);

    id1=i1;
    id2=i2;
}

template<class T>
Lorentz_Tensor<T>::Lorentz_Tensor(const Lorentz_Vector<T>& a,const Lorentz_Vector<T>& b)
{
    id1=a.GetIndex();
    id2=b.GetIndex();
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Matrix<T>::m[i][j]=a.GetElem(i)*b.GetElem(j);
}

template<class T>
Lorentz_Tensor<T>& Lorentz_Tensor<T>::operator=(const Lorentz_Tensor& b)
{
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Matrix<T>::m[i][j]=b.GetElem(i,j);

    id1=b.id1;
    id2=b.id2;
    return *this;
}

template<class T>
Lorentz_Tensor<T>& Lorentz_Tensor<T>::operator=(const Matrix<T>& b)
{
    Matrix<T>::operator=(b);
    return *this;
}

template<class T>
int Lorentz_Tensor<T>::operator==(const Lorentz_Tensor<T>& b) const
{
    int c=( id1.operator==(b.GetIndex1()) && id2.operator==(b.GetIndex2()) && this->Matrix<T>::operator==(b) )? 1 : 0 ;
    return c;
}

template<class T>
Lorentz_Tensor<T> Lorentz_Tensor<T>::operator+(const Lorentz_Tensor<T>& b) const
{
    Lorentz_Tensor<T> som=*this;
    if ( id1==b.GetIndex1() && id2==b.GetIndex2() )
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                som.Matrix<T>::m[i][j]=Matrix<T>::m[i][j]+b.Matrix<T>::m[i][j];
    else  if ( id1==b.GetIndex2() && id2==b.GetIndex1() )
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                som.Matrix<T>::m[i][j]=Matrix<T>::m[i][j]+b.Matrix<T>::m[j][i];
    else error("Bad index in Lorentz_Tensor addition");
    return som;
}

template<class T>
Lorentz_Tensor<T> Lorentz_Tensor<T>::operator-(const Lorentz_Tensor<T>& b) const
{
    Lorentz_Tensor<T> som=*this;
    if ( id1==b.GetIndex1() && id2==b.GetIndex2() )
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                som.Matrix<T>::m[i][j]=Matrix<T>::m[i][j]-b.Matrix<T>::m[i][j];
    else  if ( id1==b.GetIndex2() && id2==b.GetIndex1() )
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                som.Matrix<T>::m[i][j]=Matrix<T>::m[i][j]-b.Matrix<T>::m[j][i];
    else error("Bad index in Lorentz_Tensor substraction");
    return som;
}

template<class T>
Lorentz_Tensor<T> Lorentz_Tensor<T>::operator*(const T& b) const
{
    Lorentz_Tensor<T> som=*this;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            som.Matrix<T>::m[i][j]=Matrix<T>::m[i][j]*b;
    return som;
}

// la multiplication Matrixielle n'est pas commutative!!!!
template<class T>
Lorentz_Tensor<T> operator*(const T& b,const Lorentz_Tensor<T>& a)
{
    Lorentz_Tensor<T> som=a;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            som.SetElem(i,j,b*a.GetElem(i,j));
    return som;

}

template<class T>
Lorentz_Tensor<T> Lorentz_Tensor<T>::operator/(const Complexe& b) const
{
    if (b==Complexe(0.)) error("division by Complexe 0 of LT\n");
    Lorentz_Tensor<T> som=*this;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            som.Matrix<T>::m[i][j]=Matrix<T>::m[i][j]/b;
    return som;
}

template<class T>
Lorentz_Tensor<T> Lorentz_Tensor<T>::C(const Lorentz_Tensor<T>& b) const
{
    if ( strcmp(id2.GetName(),b.id1.GetName())!=0 ) error("Bad Index name in 2 Lorentz_Tensor contraction\n");
    if ( id2.GetType()==b.id1.GetType() ) error("Bad Index type in 2 Lorentz_Tensor contraction\n");

    Lorentz_Tensor<T> c;
    c.id1.SetIndex(id1);
    c.id2.SetIndex(b.id2);

    T sum=T(0.);
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++) {
            sum=T(0.);
            for (int k=0; k<4; k++) sum=sum+this->Matrix<T>::m[i][k]*b.m[k][j];
            c.Matrix<T>::m[i][j]=sum;
        }
    return c;
}

template<class T>
Lorentz_Vector<T> Lorentz_Tensor<T>::C(const Lorentz_Vector<T>& b) const
{
    if ( strcmp(id2.GetName(),b.GetName())!=0 ) error("Bad Index name in <T>LT <T>LV contraction\n");
    if ( id2.GetType()==b.GetType() ) error("Bad Index type in <T>LT <T>LV contraction\n");

    Lorentz_Vector<T> c;
    c.SetIndex(id1);

    T sum=T(0.);
    for (int i=0; i<4; i++) {
        sum=T(0.);
        for (int k=0; k<4; k++) sum=sum+this->Matrix<T>::m[i][k]*b.GetElem(k);
        c.SetElem(i,sum);
    }
    return c;
}


/* template<class T> */
/* Lorentz_Tensor<T> Lorentz_Tensor<T>::C(const Lorentz_Tensor<T>& b) const */
/* { */
/*   if ( strcmp(id2.GetName(),b.GetIndex1().GetName())==0 ) { */
/*    Lorentz_Tensor<T> contr=*this; */
/*    T sum; */
/*    contr.id2=b.GetIndex2(); */
/*    Complexe fact=(id2.GetType()!=b.GetIndex1().GetType() )? 1. :-1. ; */
/*    for (int i=0;i<4;i+)  */
/*      for(int j=0;j<4;j++) {  */
/*       sum=T(0.); */
/*       for (int k=1;k<4;k++) sum=sum+fact*m[i][k]*b.m[k][j]; */
/*       sum=sum+m[i][0]*b.m[0][j]; */
/*       m[i][j]=sum; */
/*      } */
/*  } */
/*  else error("bad index in <T>LT_<T>LT contraction"); */
/*  return contr; */
/* } */

/* template<class T> */
/* Lorentz_Vector<T> Lorentz_Tensor<T>::C(const Lorentz_Vector<T>& b) const */
/* { */

/*  if ( strcmp(id2.GetName(),b.GetName())==0 )  { */
/*    Lorentz_Vector<T> contr=b; */
/*    T sum; */
/*    contr.SetName(GetIndex1().GetName()); */
/*    contr.SetType(GetIndex1().GetType()); */
/*    Complexe fact=(id2.GetType()!=b.GetType() )? 1. : -1.; */
/*    for (int i=0;i<4;i+) { */
/*      for(int j=1;j<4;j++) sum=sum+m[i][j]*b.GetElem(j)*fact; */
/*      sum=sum+m[i][0]*b.GetElem(0) */
/*      contr.SetElem(i,sum); */
/*    } */
/*  } */
/*  else error("bad index in <T>LT_<T>LV contraction"); */
/*  return contr; */
/* } */

template<class T>
ostream& operator<<(ostream& f,const Lorentz_Tensor<T>& ten)
{

    f<<"\n\t\tLorentz_tensor:\nIndex1:\t"<<ten.GetIndex1().GetName()<<"\t"<<ten.GetIndex1().GetType()<<"\n";
    f<<"Index2:\t"<<ten.GetIndex2().GetName()<<"\t"<<ten.GetIndex2().GetType()<<"\n";
    for (int i=0; i<ten.GetSize(); i++) for (int j=0; j<ten.GetSize(); j++)  f<<"\n"<<i<<","<<j<<" "<<ten.GetElem(i,j);
    f<<"\n\t\tFin Lorentz_Tensor.\n";
    f.clear(ios::goodbit);
    return f;
}

#endif
