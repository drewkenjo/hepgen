
#include "myTHEO.hh"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                                COMPLEXE
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

Complexe& Complexe::operator=(const Complexe& c)
{
    re=c.re;
    im=c.im;
    return *this;
}

Complexe operator+(const double& a,const Complexe& b)
{
    return Complexe(a+b.re,b.im);
}

Complexe operator-(const double& a,const Complexe& b)
{
    return Complexe(a-b.re,-b.im);
}

Complexe operator*(const double& a,const Complexe& b)
{
    return Complexe(a*b.re,a*b.im);
}


Complexe Complexe::operator/(Complexe b) const
{
    if ( b==0. ) error("Division of a Complexe by Complexe zero\n") ;

    double rea=(re*b.re+im*b.im)/(b.re*b.re+b.im*b.im);
    double ima=(im*b.re-re*b.im)/(b.re*b.re+b.im*b.im);
    return Complexe(rea,ima);

}

Complexe Complexe::operator/(double b) const
{
    if ( b==0. ) error("Division of a Complexe by double zero\n") ;
    return Complexe(re/b,im/b);
}

Complexe operator/(double a,Complexe b)
{
    if ( b==0. ) error("Division of a double by Complexe zero") ;

    double rea=(a*b.re)/(sq(b.re)+sq(b.im));
    double ima=(-a*b.im)/(sq(b.re)+sq(b.im));
    return Complexe(rea,ima);

}

ostream& operator<<(ostream& f,const Complexe& c)
{
    if (c.imag()==0.) f<<c.real();
    else if (c.real()==0.) f<<c.imag()<<"i";
    else f<<c.real()<<"+"<<c.imag()<<"i";
    return f;
}

Complexe Complexe::sqroot() const
{
    double x=sqrt((re+sqrt(re*re+im*im))/2.);
    double y=(x==0.)? sqrt((-1.)*re) : im/2./x;
    return Complexe(x,y);
}

Complexe Complexe::Conj() const {
    return Complexe(re,-im);
}

//Complexe Conj(const double& a) { return Complexe(a); }


// Complexe cos(const Complexe& a)
//       {return Complexe(cos(a.real()));}
// Complexe sin(const Complexe& a)
//       {return Complexe(sin(a.real()));}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                                INDEX
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

Index::Index(const char* c,int ty)
{
    co=ty;
    ind=new char[10];
    strcpy(ind,c);
}

Index::Index(const Index& b)
{
    co=b.co;
    ind = new char[10];
    strcpy(ind,b.ind);

}

Index::~Index()
{
    if (ind) delete[] ind;
}

Index& Index::operator=(const Index& b)
{
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

int Index::operator==(const Index& b) const
{
    int compar=1;

    int i=0;
    if ( co!=b.co ) compar=0;
    if ( *(b.ind)!='\0' && *ind=='\0' ) compar=0;
    while ( ind[i]!='\0' && b.ind[i]!='\0' && compar==1) {
        if (ind[i]!=b.ind[i]) compar=0;
        i++;
    }
    return compar;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                       C _ LORENTZ _ VECTOR
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

C_Lorentz_Vector::C_Lorentz_Vector():Lorentz_Vector<Complexe>()
{}

C_Lorentz_Vector::C_Lorentz_Vector(const Lorentz_Vector<Complexe>& b)
{
    this->Lorentz_Vector<Complexe>::operator=(b);
}

C_Lorentz_Vector::C_Lorentz_Vector(const char* nm,int typ,const Complexe& a,const Complexe& b,const Complexe& c,const Complexe& d):Lorentz_Vector<Complexe>(nm,typ,a,b,c,d)
{}

C_Lorentz_Vector::C_Lorentz_Vector(const Index& typ,const Complexe& a,const Complexe& b,const Complexe& c,const Complexe& d):Lorentz_Vector<Complexe>(typ,a,b,c,d)
{}

C_Lorentz_Vector::C_Lorentz_Vector(const char* nm,int typ,const Vector<Complexe>& a):Lorentz_Vector<Complexe>(nm,typ,a)
{}

C_Lorentz_Vector::C_Lorentz_Vector(const Index& typ,const Vector<Complexe>& a):Lorentz_Vector<Complexe>(typ,a)
{}

C_Lorentz_Vector::C_Lorentz_Vector(const Complexe& m,const char* nm,int typ,const Complexe& b,const Complexe& c,const Complexe& d):Lorentz_Vector<Complexe>(nm,typ,(m.sq()+b.sq()+c.sq()+d.sq()).sqroot(),b,c,d)
{}

C_Lorentz_Vector::C_Lorentz_Vector(const Complexe& m,const char* nm,int typ,const Vector<Complexe>& v):Lorentz_Vector<Complexe>(nm,typ,(m.sq()+v.vscal(v)).sqroot(),v.GetElem(0),v.GetElem(1),v.GetElem(2))
{}

C_Lorentz_Vector::C_Lorentz_Vector(const Spinor& ub,const CM_Lorentz_Vector& op, const Spinor& u)
{
    if (ub.GetGeom()!=0 || u.GetGeom()!=1 ) error("Bad Geometry of spinor for creation of current");
    Complexe c[4];
    for (int k=0; k<4; k++) c[k]=ub.vscal(op.GetElem(k)*u);
    Lorentz_Vector<Complexe> C(op.GetName(),op.GetType(),c[0],c[1],c[2],c[3]);
    *this=C;
}


Complexe C_Lorentz_Vector::GetTheta() const
{
    if (GetP().real()==0.) return 0.;
    else return acos( ( GetZ()/GetP() ).real() );
}

Complexe C_Lorentz_Vector::GetPhi() const
{
    if (GetP().real()==0. || (GetX().real()==0 && GetY().real()==0) ) return 0.;
    else if (GetX().real() !=0.)
        return (GetX().real()>0)?  atan( (GetY()/GetX()).real() ):
               3.14159265+atan( (GetY()/GetX()).real() );
    else return (GetY().real()>0)? 1.5707963 : -1.5707963 ;
}

C_Lorentz_Vector C_Lorentz_Vector::operator*(double b) const
{
    C_Lorentz_Vector c;
    c.SetIndex(GetIndex());
    for (int k=0; k<4; k++) c.SetElem(k,v[k]*b);
    return c;
}

C_Lorentz_Vector C_Lorentz_Vector::operator*(const Complexe& b) const
{
    C_Lorentz_Vector c;
    c.SetIndex(GetIndex());
    for (int k=0; k<4; k++) c.SetElem(k,b*v[k]);
    return c;
}

CM_Lorentz_Vector C_Lorentz_Vector::operator*(const CMatrix& b) const
{
    CM_Lorentz_Vector c;
    c.SetIndex(GetIndex());
    for (int k=0; k<4; k++) c.SetElem(k,b*v[k]);
    return c;
}

CM_Lorentz_Vector C_Lorentz_Vector::operator+(const CM_Lorentz_Vector& b)
{
    CMatrix unit(Complexe(1.));
    return operator*(unit)+b;
}

CM_Lorentz_Vector C_Lorentz_Vector::operator-(const CM_Lorentz_Vector& b)
{
    CMatrix unit(Complexe(1.));
    return operator*(unit)-b;
}

// Complexe C_Lorentz_Vector::C(const C_Lorentz_Vector& b) const
// {
//  if ( strcmp(this->GetName(),b.GetName())!=0 ) error("Bad Index name in Lorentz_Vector C_C contraction\n");
//  if ( this->GetType()==b.GetType() ) error("Bad Index type in Lorentz_Vector C_C contraction\n");

//  Complexe c;

//  for (int k=0;k<4;k++) c=c + this->GetElem(k) * b.GetElem(k) ;
//  return c;
// }

// CMatrix C_Lorentz_Vector::C(const CM_Lorentz_Vector& b) const
// {
//  if ( strcmp(this->GetName(),b.GetName())!=0 ) error("Bad Index name in Lorentz_Vector C_CM contraction\n");
//  if ( this->GetType()==b.GetType() ) error("Bad Index type in Lorentz_Vector C_CM contraction\n");

//  CMatrix c;

//  for (int i=0;i<4;i++) c=c + this->GetElem(i) * b.GetElem(i) ;
//  return c;
// }


Complexe C_Lorentz_Vector::C(const C_Lorentz_Vector& b) const
{
// En template
    return this->Lorentz_Vector<Complexe>::C(b);
}

CMatrix C_Lorentz_Vector::C(const CM_Lorentz_Vector& b) const
{
    if ( strcmp(this->GetName(),b.GetName())!=0 ) error("Bad Index name in C_LV CM_LV contraction\n");
    if ( this->GetType()==b.GetType() ) error("Bad Index type in C_LV CM_LV contraction\n");

    CMatrix c;

    for (int i=0; i<4; i++) c=c + b.GetElem(i) * (this->GetElem(i)) ;
    return c;
}

C_Lorentz_Vector C_Lorentz_Vector::C(const C_Lorentz_Tensor& b) const
{
// En template
    return C_Lorentz_Vector(this->Lorentz_Vector<Complexe>::C(b));
}

CM_Lorentz_Vector C_Lorentz_Vector::C(const CM_Lorentz_Tensor& b) const
{
    if ( strcmp(this->GetName(),b.GetIndex1().GetName())!=0 ) error("Bad Index name in C_LV CM_LV contraction\n");
    if ( this->GetType()==b.GetIndex1().GetType() ) error("Bad Index type in C_LV CM_LV contraction\n");

    CM_Lorentz_Vector c;
    c.SetIndex( b.GetIndex2() );

    for (int i=0; i<4; i++) {
        CMatrix sum;
        for (int j=0; j<4; j++) sum=sum + b.GetElem(j,i) * GetElem(j) ;
        c.SetElem(i,sum);
    }
    return c;
}

CMatrix C_Lorentz_Vector::Slash() const
{
    // en representation de Dirac:
    CMatrix c;
    c.SetElem(0,0,GetE());
    c.SetElem(1,1,GetE());
    c.SetElem(2,2,(-1.)*GetE());
    c.SetElem(3,3,(-1.)*GetE());

    Complexe f=( GetType()==0 )? 1. : -1. ;
    c.SetElem(0,2,f*GetZ());
    c.SetElem(0,3,f*GetX()-Complexe(0.,1.)*f*GetY());
    c.SetElem(1,2,f*GetX()+Complexe(0.,1.)*f*GetY());
    c.SetElem(1,3,(-1.)*f*GetZ());

    c.SetElem(2,0,(-1.)*f*GetZ());
    c.SetElem(2,1,(-1.)*f*GetX()+Complexe(0.,1.)*f*GetY());
    c.SetElem(3,0,(-1.)*f*GetX()-Complexe(0.,1.)*f*GetY());
    c.SetElem(3,1,f*GetZ());

    return c;
}

CM_Lorentz_Vector C_Lorentz_Vector::Id() const
{
    CMatrix c0(v[0]);
    CMatrix c1(v[1]);
    CMatrix c2(v[2]);
    CMatrix c3(v[3]);

    return CM_Lorentz_Vector(GetName(),GetType(),c0,c1,c2,c3);
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                       CM _ LORENTZ _ VECTOR
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

CM_Lorentz_Vector::CM_Lorentz_Vector():Lorentz_Vector<CMatrix>()
{}

CM_Lorentz_Vector::CM_Lorentz_Vector(const Lorentz_Vector<CMatrix>& b)
{
    this->Lorentz_Vector<CMatrix>::operator=(b);
}

//CM_Lorentz_Vector::CM_Lorentz_Vector(const C_Lorentz_Vector& b)
//{
// this->Index::operator=(b.GetIndex());
// CMatrix unit(Complexe(1.));
// for (int ii=0;ii<4;ii++) this->SetElem(ii,b.GetElem(ii)*unit);
//}

CM_Lorentz_Vector::CM_Lorentz_Vector(const char* nm,int typ,const CMatrix& a,const CMatrix& b,const CMatrix& c,const CMatrix& d):Lorentz_Vector<CMatrix>(nm,typ,a,b,c,d)
{}

CM_Lorentz_Vector::CM_Lorentz_Vector(const Index& typ,const CMatrix& a,const CMatrix& b,const CMatrix& c,const CMatrix& d):Lorentz_Vector<CMatrix>(typ,a,b,c,d)
{}

CM_Lorentz_Vector::CM_Lorentz_Vector(const char* nm,int typ,const Vector<CMatrix>& a):Lorentz_Vector<CMatrix>(nm,typ,a)
{}

CM_Lorentz_Vector::CM_Lorentz_Vector(const Index& typ,const Vector<CMatrix>& a):Lorentz_Vector<CMatrix>(typ,a)
{}

CM_Lorentz_Vector CM_Lorentz_Vector::operator*(double b) const
{
    CM_Lorentz_Vector c;
    c.SetIndex(GetIndex());
    for (int k=0; k<4; k++) c.SetElem(k,v[k]*CMatrix(b));
    return c;
}

CM_Lorentz_Vector CM_Lorentz_Vector::operator*(const Complexe& b) const
{
    CM_Lorentz_Vector c;
    c.SetIndex(GetIndex());
    for (int k=0; k<4; k++) c.SetElem(k,v[k]*CMatrix(b));
    return c;
}

CM_Lorentz_Vector CM_Lorentz_Vector::operator*(const CMatrix& b) const
{
    CM_Lorentz_Vector c;
    c.SetIndex(GetIndex());
    for (int k=0; k<4; k++) c.SetElem(k,v[k]*b);
    return c;
}

CM_Lorentz_Vector CM_Lorentz_Vector::operator+(const C_Lorentz_Vector& b)
{
    CMatrix unit(Complexe(1.));
    return operator+(unit*b);
}

CM_Lorentz_Vector CM_Lorentz_Vector::operator-(const C_Lorentz_Vector& b)
{
    CMatrix unit(Complexe(1.));
    return operator-(unit*b);
}

CMatrix CM_Lorentz_Vector::C(const CM_Lorentz_Vector& b) const
{
    // en template
    return this->Lorentz_Vector<CMatrix>::C(b);
}

CMatrix CM_Lorentz_Vector::C(const C_Lorentz_Vector& b) const
{
    // de meme que C_LV . CM_LV car pas d'indices en sortie
    return b.C(*this);
}

CM_Lorentz_Vector CM_Lorentz_Vector::C(const CM_Lorentz_Tensor& b) const
{
    // en template
    return CM_Lorentz_Vector(this->Lorentz_Vector<CMatrix>::C(b));
}

CM_Lorentz_Vector CM_Lorentz_Vector::C(const C_Lorentz_Tensor& b) const
{
    // de meme que C_LV CM_LT car CM*C=C*CM
    if ( strcmp(this->GetName(),b.GetIndex1().GetName())!=0 ) error("Bad Index name in CM_LV C_LV contraction\n");
    if ( this->GetType()==b.GetIndex1().GetType() ) error("Bad Index type in CM_LV C_LV contraction\n");

    CM_Lorentz_Vector c;
    c.SetIndex( b.GetIndex2() );

    for (int i=0; i<4; i++) {
        CMatrix sum;
        for (int j=0; j<4; j++) sum=sum +GetElem(j) * b.GetElem(j,i) ;
        c.SetElem(i,sum);
    }
    return c;
}


CM_Lorentz_Vector GammaDirac(const char* n,int t)
{
    CMatrix g0,g1,g2,g3;
    Complexe i(0.,1.);

//Matrixes gamma en representation de DIRAC:
    g0.SetElem(0,0,1.);
    g0.SetElem(1,1,1.);
    g0.SetElem(2,2,-1.);
    g0.SetElem(3,3,-1.);
    g1.SetElem(0,3,1.);
    g1.SetElem(1,2,1.);
    g1.SetElem(2,1,-1.);
    g1.SetElem(3,0,-1.);
    g2.SetElem(0,3,-1.*i);
    g2.SetElem(1,2,i);
    g2.SetElem(2,1,i);
    g2.SetElem(3,0,-1.*i);
    g3.SetElem(0,2,1.);
    g3.SetElem(1,3,-1.);
    g3.SetElem(2,0,-1.);
    g3.SetElem(3,1,1.);

    if (t==1) return CM_Lorentz_Vector(n,1,g0,g1,g2,g3);
    else return CM_Lorentz_Vector(n,1,g0,g1,g2,g3).Idxd(n,0);
}

CMatrix Gamma5()
{
    CMatrix g0,g1,g2,g3;
    Complexe i(0.,1.);
//Matrixes gamma en representation de DIRAC:
    g0.SetElem(0,0,1.);
    g0.SetElem(1,1,1.);
    g0.SetElem(2,2,-1.);
    g0.SetElem(3,3,-1.);
    g1.SetElem(0,3,1.);
    g1.SetElem(1,2,1.);
    g1.SetElem(2,1,-1.);
    g1.SetElem(3,0,-1.);
    g2.SetElem(0,3,-1.*i);
    g2.SetElem(1,2,i);
    g2.SetElem(2,1,i);
    g2.SetElem(3,0,-1.*i);
    g3.SetElem(0,2,1.);
    g3.SetElem(1,3,-1.);
    g3.SetElem(2,0,-1.);
    g3.SetElem(3,1,1.);

    return Complexe(0.,1.)*g0*g1*g2*g3;
}

CM_Lorentz_Vector operator*(const CMatrix& a, const CM_Lorentz_Vector& b)
{
    CM_Lorentz_Vector c;
    c.SetIndex(b.GetIndex());
    for (int k=0; k<4; k++) c.SetElem(k,a*b.GetElem(k));
    return c;
}



//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                              C_Lorentz_Tensor
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


C_Lorentz_Tensor::C_Lorentz_Tensor():Lorentz_Tensor<Complexe>()
{}

C_Lorentz_Tensor::C_Lorentz_Tensor(const Lorentz_Tensor<Complexe>& b)
{
    this->Lorentz_Tensor<Complexe>::operator=(b);
}

C_Lorentz_Tensor::C_Lorentz_Tensor(const C_Lorentz_Vector& a,const C_Lorentz_Vector& b)
{
    Lorentz_Tensor<Complexe> t(a,b);
    this->Lorentz_Tensor<Complexe>::operator=(t);
}

C_Lorentz_Tensor::C_Lorentz_Tensor(const char* n1,int t1,const char* n2,int t2)
{
    Lorentz_Tensor<Complexe> t;
    this->Lorentz_Tensor<Complexe>::operator=(t);
    id1=Index(n1,t1);
    id2=Index(n2,t2);
}

C_Lorentz_Tensor::C_Lorentz_Tensor(const char* n1,int t1,const char* n2,int t2,const Complexe& a,const Complexe& b)
{
    Lorentz_Tensor<Complexe> t;
    this->Lorentz_Tensor<Complexe>::operator=(t);
    id1=Index(n1,t1);
    id2=Index(n2,t2);
    Complexe f=( t1==t2 )? -1. : 1;
    m[0][0]=1;
    m[1][1]=f*cos(a)*cos(b);
    m[1][2]=(-1.)*f*cos(a)*sin(b);
    m[1][3]=f*sin(a);
    m[2][1]=f*sin(b);
    m[2][2]=f*cos(b);
    m[2][3]=0.;
    m[3][1]=(-1.)*f*sin(a)*cos(b);
    m[3][2]=f*sin(a)*sin(b);
    m[3][3]=f*cos(a);
}

C_Lorentz_Tensor::C_Lorentz_Tensor(const Spinor& ub,const CM_Lorentz_Tensor& op, const Spinor& u)
{
    if (ub.GetGeom()!=0 || u.GetGeom()!=1 ) error("Bad Geometry of spinor for creation of Tensor current");
    Complexe c;
    Lorentz_Tensor<Complexe> t;
    this->Lorentz_Tensor<Complexe>::operator=(t);
    SetIndex1( op.GetIndex1() );
    SetIndex2( op.GetIndex2() );
    for (int k1=0; k1<4; k1++)
        for (int k2=0; k2<4; k2++) {
            c=ub.vscal(op.GetElem(k1,k2)*u);
            SetElem(k1,k2,c);
        }
}


CM_Lorentz_Tensor C_Lorentz_Tensor::operator*(const CMatrix& b) const
{
    CM_Lorentz_Tensor c;
    c.SetIndex1( GetIndex1() );
    c.SetIndex2( GetIndex2() );
    for (int k1=0; k1<4; k1++)
        for (int k2=0; k2<4; k2++) c.SetElem(k1,k2,GetElem(k1,k2)*b);
    return c;
}

C_Lorentz_Tensor C_Lorentz_Tensor::operator*(const Complexe& b) const
{
    return C_Lorentz_Tensor(Lorentz_Tensor<Complexe>::operator*(b)) ;
}

C_Lorentz_Tensor C_Lorentz_Tensor::operator*(double b) const
{
    return C_Lorentz_Tensor(Lorentz_Tensor<Complexe>::operator*(Complexe(b))) ;
}

CM_Lorentz_Tensor C_Lorentz_Tensor::operator+(const CM_Lorentz_Tensor& b)
{
    CMatrix unit(Complexe(1.));
    return operator*(unit)+b;
}

CM_Lorentz_Tensor C_Lorentz_Tensor::operator-(const CM_Lorentz_Tensor& b)
{
    CMatrix unit(Complexe(1.));
    return operator*(unit)-b;
}

C_Lorentz_Tensor C_Lorentz_Tensor::C(const C_Lorentz_Tensor& b) const
{
    // en template
    return C_Lorentz_Tensor(this->Lorentz_Tensor<Complexe>::C(b));
}

CM_Lorentz_Tensor C_Lorentz_Tensor::C(const CM_Lorentz_Tensor& b) const
{
    if ( strcmp(id2.GetName(),b.GetIndex1().GetName())!=0 ) error("Bad index name in C_LT CM_LT contraction/n");
    if ( id2.GetType()==b.GetIndex1().GetType() ) error("Bad index name in C_LT CM_LT contraction/n");

    CM_Lorentz_Tensor contr;
    contr.SetIndex1( GetIndex1() );
    contr.SetIndex2( b.GetIndex2() );
    for (int i=0; i<4; i++)
        for(int j=0; j<4; j++) {
            CMatrix sum(0.);
            for (int k=0; k<4; k++) sum=sum+b.GetElem(k,j)*m[i][k];
            contr.SetElem(i,j,sum);
        }
    return contr;
}

C_Lorentz_Vector C_Lorentz_Tensor::C(const C_Lorentz_Vector& b) const
{
    // en template
    return C_Lorentz_Vector(this->Lorentz_Tensor<Complexe>::C(b));
}

CM_Lorentz_Vector C_Lorentz_Tensor::C(const CM_Lorentz_Vector& b) const
{
    if ( strcmp(id2.GetName(),b.GetName())!=0 ) error("Bad index name in C_LT CM_LT contraction/n");
    if ( id2.GetType()==b.GetType() ) error("Bad index name in C_LT CM_LT contraction/n");

    CM_Lorentz_Vector contr;
    contr.SetIndex( GetIndex1() );
    for (int i=0; i<4; i++) {
        CMatrix sum(0.);
        for (int k=0; k<4; k++) sum=sum+b.GetElem(k)*m[i][k];
        contr.SetElem(i,sum);
    }
    return contr;
}

C_Lorentz_Tensor gT(const char* n1,int t1,const char* n2,int t2)
{
    C_Lorentz_Tensor gg;
    gg.id1=Index(n1,t1);
    gg.id2=Index(n2,t2);
    Complexe f=(t1==t2)? -1 : 1. ;
    gg.m[0][0]=1.;
    gg.m[1][1]=f;
    gg.m[2][2]=f;
    gg.m[3][3]=f;
    return gg;
}

C_Lorentz_Tensor TLorentz(const C_Lorentz_Vector& b,const char* n1,int t1,const char* n2,int t2 ,const char* dir)
{
    Complexe fact=(strcmp(dir,"inverse")==0)? -1. : 1.;
    Complexe mb=b.Mass2().sqroot();

    C_Lorentz_Vector br=(b.GetIndex().GetType()==0)? b.Idxd("ukwn",1) : b;

    Complexe E=br.GetE();
    Complexe g=E/mb;

    Complexe bv[3]= {fact*br.GetX()/E,fact*br.GetY()/E,fact*br.GetZ()/E};
    Complexe bmod=br.GetP()/E;

    C_Lorentz_Tensor res(n1,1,n2,0);

    res.SetElem(0,0,g);
    res.SetElem(0,1,-1.*g*bv[0]);
    res.SetElem(0,2,-1.*g*bv[1]);
    res.SetElem(0,3,-1.*g*bv[2]);
    res.SetElem(1,0,-1.*g*bv[0]);
    res.SetElem(2,0,-1.*g*bv[1]);
    res.SetElem(3,0,-1.*g*bv[2]);
    for (int k=0; k<3; k++)
        for (int l=0; l<3; l++)
            if (k!=l) res.SetElem(k+1,l+1,(g-1.)*bv[k]*bv[l]/bmod.sq());
            else  res.SetElem(k+1,l+1,(g-1.)*bv[k]*bv[l]/bmod.sq()+1.);

    if (t1==0) res=gT(n1,0,n1,0).C(res);
    if (t2==1) res=res.C(gT(n2,1,n2,1));

    return res;
}

C_Lorentz_Tensor LeviCita(const char* n1,int t1,const char* n2,int t2,const C_Lorentz_Tensor& b)
{
    //convention indices en bas 0123=-1
    Complexe coef=-1.;

    C_Lorentz_Tensor c,bcor=b;
    c.SetIndex1( Index(n1,0) );
    c.SetIndex2( Index(n2,0));

    if (b.GetIndex1().GetType()==0) bcor=gT("ukwn",1,b.GetIndex1().GetName(),1).C(b);
    if (b.GetIndex2().GetType()==0) bcor=bcor.C(gT(b.GetIndex2().GetName(),1,"ukwn",1));

    static int Levi_deja_cree=0;
    static double levi[4][4][4][4];
    Complexe sgn;

    if (Levi_deja_cree==0) {
        epsil(levi);
        Levi_deja_cree=1;
    }

    for (int mu=0; mu<4; mu++)
        for (int nu=0; nu<4; nu++) {
            Complexe cadd;
            for (int k=0; k<4; k++)
                for (int l=0; l<4; l++) {
                    sgn=levi[mu][nu][k][l];
                    cadd=cadd+sgn*bcor.GetElem(k,l);
                    c.SetElem(mu,nu,cadd);
                }
        }
    if (t1==1) c=gT(n1,1,n1,1).C(c);
    if (t2==1) c=c.C(gT(n2,1,n2,1));
    return -1.*c;
}

void epsil(double e[4][4][4][4])
{
    for (int i1=0; i1<4; i1++)
        for (int i2=0; i2<4; i2++)
            for (int i3=0; i3<4; i3++)
                for (int i4=0; i4<4; i4++)
                    if ( i1==i2 || i1==i3 || i1==i4 || i2==i3 || i2==i4 || i3==i4 )
                        e[i1][i2][i3][i4]=0.;
                    else {
                        int a1=i1;
                        int a2=i2;
                        int a3=i3;
                        int a4=i4;
                        int t;
                        double sgn=1.;
                        while ( !( a1<a2 && a2<a3 && a3<a4 ) ) {
                            if (a3>a4) {
                                t=a4;
                                a4=a3;
                                a3=t;
                                sgn=-1.*sgn;
                            }
                            if (a2>a3) {
                                t=a3;
                                a3=a2;
                                a2=t;
                                sgn=-1.*sgn;
                            }
                            if (a1>a2) {
                                t=a2;
                                a2=a1;
                                a1=t;
                                sgn=-1.*sgn;
                            }
                        }
                        e[i1][i2][i3][i4]=sgn;
                    }
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                              CM_Lorentz_Tensor
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


CM_Lorentz_Tensor::CM_Lorentz_Tensor():Lorentz_Tensor<CMatrix>()
{}

CM_Lorentz_Tensor::CM_Lorentz_Tensor(const Lorentz_Tensor<CMatrix>& b)
{
    this->Lorentz_Tensor<CMatrix>::operator=(b);
}

//CM_Lorentz_Tensor::CM_Lorentz_Tensor(const C_Lorentz_Tensor& b)
//{
// id1=b.GetIndex1();
// id2=b.GetIndex2();
// CMatrix unit(Complexe(1.));
// for (int ii=0;ii<4;ii++)
//   for (int jj=0;jj<4;jj++) this->SetElem(ii,jj,b.GetElem(ii,jj)*unit);
//}

CM_Lorentz_Tensor::CM_Lorentz_Tensor(const CM_Lorentz_Vector& a,const CM_Lorentz_Vector& b)
{
    Lorentz_Tensor<CMatrix> t(a,b);
    this->Lorentz_Tensor<CMatrix>::operator=(t);
}

CM_Lorentz_Tensor::CM_Lorentz_Tensor(const char* n1,int t1,const char* n2,int t2)
{
    Lorentz_Tensor<CMatrix> t;
    this->Lorentz_Tensor<CMatrix>::operator=(t);
    id1=Index(n1,t1);
    id2=Index(n2,t2);
}

CM_Lorentz_Tensor::CM_Lorentz_Tensor(const char* n1,int t1,const char* n2,int t2,const char* n)
{
    CM_Lorentz_Vector gam=GammaDirac(n1,t1);
    this->Lorentz_Tensor<CMatrix>::operator=( Complexe(0.,0.5)*( CM_Lorentz_Tensor(gam,gam.Idxd(n2,t2))-
            CM_Lorentz_Tensor(gam.Idxd(n2,t2),gam) ));
}

CM_Lorentz_Tensor CM_Lorentz_Tensor::operator*(double b) const
{
    CM_Lorentz_Tensor c;
    c.SetIndex1(id1);
    c.SetIndex2(id2);

    for (int k=0; k<4; k++)
        for (int l=0; l<4; l++) c.SetElem(k,l,m[k][l]*CMatrix(b));
    return c;
}

CM_Lorentz_Tensor CM_Lorentz_Tensor::operator*(const Complexe& b) const
{
    CM_Lorentz_Tensor c;
    c.SetIndex1(id1);
    c.SetIndex2(id2);

    for (int k=0; k<4; k++)
        for (int l=0; l<4; l++) c.SetElem(k,l,m[k][l]*CMatrix(b));
    return c;
}

CM_Lorentz_Tensor CM_Lorentz_Tensor::operator*(const CMatrix& b) const
{
    CM_Lorentz_Tensor c;
    c.SetIndex1(id1);
    c.SetIndex2(id2);

    for (int k=0; k<4; k++)
        for (int l=0; l<4; l++) c.SetElem(k,l,m[k][l]*b);
    return c;
}

CM_Lorentz_Tensor CM_Lorentz_Tensor::operator+(const C_Lorentz_Tensor& b)
{
    CMatrix unit(Complexe(1.));
    return operator+(unit*b);
}

CM_Lorentz_Tensor CM_Lorentz_Tensor::operator-(const C_Lorentz_Tensor& b)
{
    CMatrix unit(Complexe(1.));
    return operator-(unit*b);
}


CM_Lorentz_Tensor CM_Lorentz_Tensor::C(const CM_Lorentz_Tensor& b) const
{
    // en template
    return CM_Lorentz_Tensor(this->Lorentz_Tensor<CMatrix>::C(b));
}

CM_Lorentz_Tensor CM_Lorentz_Tensor::C(const C_Lorentz_Tensor& b) const
{
    if ( strcmp(id2.GetName(),b.GetIndex1().GetName())!=0 ) error("Bad index name in CM_LT C_LT contraction/n");
    if ( id2.GetType()==b.GetIndex1().GetType() ) error("Bad index name in CM_LT C_LT contraction/n");

    CM_Lorentz_Tensor contr=*this;
    contr.SetIndex2( b.GetIndex2() );
    for (int i=0; i<4; i++)
        for(int j=0; j<4; j++) {
            CMatrix sum(0.);
            for (int k=0; k<4; k++) sum=sum+m[i][k]*b.GetElem(k,j);
            contr.m[i][j]=sum;
        }
    return contr;
}

CM_Lorentz_Vector CM_Lorentz_Tensor::C(const CM_Lorentz_Vector& b) const
{
    // en template
    return CM_Lorentz_Vector(this->Lorentz_Tensor<CMatrix>::C(b));
}


CM_Lorentz_Vector CM_Lorentz_Tensor::C(const C_Lorentz_Vector& b) const
{
    if ( strcmp(id2.GetName(),b.GetName())!=0 ) error("Bad index name in CM_LT C_LV contraction/n");
    if ( id2.GetType()==b.GetType() ) error("Bad index name in CM_LT C_LV contraction/n");

    CM_Lorentz_Vector contr;
    contr.SetIndex( GetIndex1() );
    for (int i=0; i<4; i++) {
        CMatrix sum(0.);
        for (int k=0; k<4; k++) sum=sum+m[i][k]*b.GetElem(k);
        contr.SetElem(i,sum);
    }
    return contr;
}



C_Lorentz_Tensor operator*(double a, const C_Lorentz_Tensor& CLT)
// il faut que ce soit explicite pour que ca prevalle a double*Matrix<Complexe>
{
    return CLT.Lorentz_Tensor<Complexe>::operator*(Complexe(a));
}

C_Lorentz_Tensor operator*(const Complexe& a, const C_Lorentz_Tensor& CLT)
// il faut que ce soit explicite pour que ca prevalle a double*Matrix<Complexe>
{
    return CLT.Lorentz_Tensor<Complexe>::operator*(a);
}

CM_Lorentz_Tensor operator*(const CMatrix& a, const C_Lorentz_Tensor& b)
{
    return CM_Lorentz_Tensor(b*a);
}


CM_Lorentz_Tensor operator*(double a, const CM_Lorentz_Tensor& b)
{
    return CM_Lorentz_Tensor(b*CMatrix(a));
}

CM_Lorentz_Tensor operator*(const Complexe& a, const CM_Lorentz_Tensor& b)
{
    return CM_Lorentz_Tensor(b*CMatrix(a));
}

CM_Lorentz_Tensor operator*(const CMatrix& a, const CM_Lorentz_Tensor& b)
{
    CM_Lorentz_Tensor c;
    c.SetIndex1(b.GetIndex1());
    c.SetIndex2(b.GetIndex2());

    for (int k=0; k<4; k++)
        for (int l=0; l<4; l++) c.SetElem(k,l,a*b.GetElem(k,l));
    return c;
}

CM_Lorentz_Tensor TSigmaDirac(const char* n1,int t1,const char* n2,int t2)
{
    CM_Lorentz_Vector gam=GammaDirac(n1,t1);
    return Complexe(0.,0.5)*( CM_Lorentz_Tensor(gam,gam.Idxd(n2,t2))-
                              CM_Lorentz_Tensor(gam.Idxd(n2,t2),gam) );
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                              BISPINOR
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


Spinor::Spinor():Vector<Complexe>()
{
    h=pos;
    pa=part;
// v[0]=0.;v[1]=0.;v[2]=0.;v[3]=0.;
}

Spinor::Spinor(double):Vector<Complexe>()
{
    h=pos;
    pa=part;
//v[0]=0.;v[1]=0.;v[2]=0.;v[3]=0.;
}

Spinor::Spinor(helicity he,particule par,invar in,const Vector<Complexe>& vec):Vector<Complexe>(vec)
{
    h=he;
    pa=par;
    inv=in;
}



Spinor::Spinor(particule par,const C_Lorentz_Vector& Cv,helicity hel,
               Complexe ma,invar inv):Vector<Complexe>()
// Construction d'un spineur de Dirac pour une particule par
//          de masse ma
//          de vecteur de Lorentz Cv
//          d'helicite hel
// spineur normalise a 2m

{
    Complexe i(0.,1.);
    if (inv==spin) {

        Vector<Complexe> X1(0.,2,1);
//  if (!( ma.real()==Cv.Mass().real() )) printf("****Warning: particule not on mass shell for spinor creation\n");

        if ( (hel==pos && par==part) || (hel==neg && par==anti) ) {
            X1.SetElem(0,1.);
            X1.SetElem(1,0.);
        }
        else {
            X1.SetElem(0,0.);
            X1.SetElem(1,1.);
        }

        Matrix<Complexe> SP(0.,2);
        SP.SetElem(0,0,Cv.GetElem(3));
        SP.SetElem(0,1,Cv.GetElem(1)+(-1.)*i*Cv.GetElem(2));
        SP.SetElem(1,0,Cv.GetElem(1)+i*Cv.GetElem(2));
        SP.SetElem(1,1,(-1.)*Cv.GetElem(3));
        SP=SP/(Cv.GetElem(0)+ma);

        Vector<Complexe> X2=SP*X1;

        Complexe N=(Cv.GetElem(0)+ma).sqroot();

        h=hel;
        pa=par;
        if (par==part) {
            v[0]=N*X1.GetElem(0);
            v[1]=N*X1.GetElem(1);
            v[2]=N*X2.GetElem(0);
            v[3]=N*X2.GetElem(1);
        }
        else {
            v[0]=N*X2.GetElem(0);
            v[1]=N*X2.GetElem(1);
            v[2]=N*X1.GetElem(0);
            v[3]=N*X1.GetElem(1);
        }
    }
    else {
        Vector<Complexe> X1(0.,2,1);

        Complexe coef=1.;
        if (hel==pos) {
            X1.SetElem(0,cos(Cv.GetTheta()/2.));
            X1.SetElem(1,exp(i*Cv.GetPhi())*sin(Cv.GetTheta()/2.));
        }
        else {
            coef=-1.;
            X1.SetElem(0,-1.*exp(-1.*i*Cv.GetPhi())*sin(Cv.GetTheta()/2.));
            X1.SetElem(1,cos(Cv.GetTheta()/2.));
        }
        h=hel;
        pa=par;
        v[0]=(Cv.GetE()+ma).sqroot()*X1.GetElem(0);
        v[1]=(Cv.GetE()+ma).sqroot()*X1.GetElem(1);
        v[2]=coef*(Cv.GetE()-ma).sqroot()*X1.GetElem(0);
        v[3]=coef*(Cv.GetE()-ma).sqroot()*X1.GetElem(1);
    }
}

Spinor::Spinor(particule par,const C_Lorentz_Vector& Cv,double hel,
               Complexe ma,invar inv):Vector<Complexe>()
// Construction d'un spineur de Dirac pour une particule par
//          de masse ma
//          de vecteur de Lorentz Cv
//          d'helicite hel
// spineur normalise a 2m

{
    Complexe i(0.,1.);
    if (hel!=-0.5 && hel!=0.5 ) {
        v[0]=0.;
        v[1]=0.;
        v[2]=0.;
        v[3]=0.;
    }
    else {

        if (inv==spin) {

            Vector<Complexe> X1(0.,2,1);
//  if (!( ma.real()==Cv.Mass().real() )) printf("****Warning: particule not on mass shell for spinor creation\n");

            if ( (hel==.5 && par==part) || (hel==-.5 && par==anti) ) {
                X1.SetElem(0,1.);
                X1.SetElem(1,0.);
            }
            else {
                X1.SetElem(0,0.);
                X1.SetElem(1,1.);
            }

            Matrix<Complexe> SP(0.,2);
            SP.SetElem(0,0,Cv.GetElem(3));
            SP.SetElem(0,1,Cv.GetElem(1)+(-1.)*i*Cv.GetElem(2));
            SP.SetElem(1,0,Cv.GetElem(1)+i*Cv.GetElem(2));
            SP.SetElem(1,1,(-1.)*Cv.GetElem(3));
            SP=SP/(Cv.GetElem(0)+ma);

            Vector<Complexe> X2=SP*X1;

            Complexe N=(Cv.GetElem(0)+ma).sqroot();

            h=(hel==.5)? pos : neg ;
            pa=par;
            if (par==part) {
                v[0]=N*X1.GetElem(0);
                v[1]=N*X1.GetElem(1);
                v[2]=N*X2.GetElem(0);
                v[3]=N*X2.GetElem(1);
            }
            else {
                v[0]=N*X2.GetElem(0);
                v[1]=N*X2.GetElem(1);
                v[2]=N*X1.GetElem(0);
                v[3]=N*X1.GetElem(1);
            }
        }
        else {
            Vector<Complexe> X1(0.,2,1);

            Complexe coef=1.;
            if (hel==.5) {
                X1.SetElem(0,cos(Cv.GetTheta()/2.));
                X1.SetElem(1,exp(i*Cv.GetPhi())*sin(Cv.GetTheta()/2.));
                h=pos;
            }
            else {
                coef=-1.;
                X1.SetElem(0,-1.*exp(-1.*i*Cv.GetPhi())*sin(Cv.GetTheta()/2.));
                X1.SetElem(1,cos(Cv.GetTheta()/2.));
                h=neg;
            }
            pa=par;
            v[0]=(Cv.GetE()+ma).sqroot()*X1.GetElem(0);
            v[1]=(Cv.GetE()+ma).sqroot()*X1.GetElem(1);
            v[2]=coef*(Cv.GetE()-ma).sqroot()*X1.GetElem(0);
            v[3]=coef*(Cv.GetE()-ma).sqroot()*X1.GetElem(1);
        }

    } // hel not defined
}

Spinor Spinor::bar() const
{
    Spinor bs;
    Vector<Complexe> bsh;

    bsh=Trans().Conj();
    bs.v[0]=bsh.GetElem(0);
    bs.v[1]=bsh.GetElem(1);
    bs.v[2]=(-1.)*bsh.GetElem(2);
    bs.v[3]=(-1.)*bsh.GetElem(3);
    bs.sh=bsh.GetGeom();
    bs.h=h;
    bs.pa=pa;
    return bs;
}


C_Lorentz_Vector Spinor::mult(const CM_Lorentz_Vector& CMLV,const Spinor& s)
{
    if (GetGeom()==1 || s.GetGeom()==0) error("pb in ubar O u: geometry of spinor");
    C_Lorentz_Vector c;
    c.SetIndex(CMLV.GetIndex());
    for (int i=0; i<4; i++)
        c.SetElem(i,vscal(CMLV.GetElem(i)*s));
    return c;
}

C_Lorentz_Tensor Spinor::mult(const CM_Lorentz_Tensor& CMLT,const Spinor& s)
{
    if (GetGeom()==1 || s.GetGeom()==0) error("pb in ubar O u: geometry of spinor");
    C_Lorentz_Tensor c;
    c.SetIndex1(CMLT.GetIndex1());
    c.SetIndex2(CMLT.GetIndex2());
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            c.SetElem(i,j,vscal(CMLT.GetElem(i,j)*s));
    return c;
}

ostream& operator<<(ostream& f,const Spinor& bs)
{
    if (bs.GetGeom()==0) f<<"line";
    else f<<"column"<<"\n";
    for (int i=0; i<bs.GetSize(); i++)  f<<"\n"<<bs.GetElem(i);
    f<<"\n";
    f.clear(ios::goodbit);
    return f;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//
//                              SPINOR FOR DELTA
//
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

Dspinor::Dspinor(const char* n,int ty,const C_Lorentz_Vector& p,double Sz,Complexe mD)
{
    SetIndex(Index(n,ty));
// covariant or contravariant:
    double coef=(ty==1)? 1. : -1. ;

// coef de Clebsh-Gordon < 1/2 sz 1 la | 3/2 Sz >
    double CSz=( Sz==1.5 || Sz==-1.5 )? 6. : 2. ;

    Vector<Complexe> El[3];
    El[0]=Vector<Complexe>( Complexe(1./sqrt(2.)),-1.*Complexe(0.,1./sqrt(2.)),Complexe(0.) );
    El[1]=Vector<Complexe>( Complexe(0.),Complexe(0.),Complexe(1.) );
    El[2]=Vector<Complexe>( Complexe(-1./sqrt(2.)),-1.*Complexe(0.,1./sqrt(2.)),Complexe(0.) );


    Complexe pla[3];
    for (int ix=1; ix<=3; ix++)
        pla[ix]=p.Get3Vector().vscal(El[ix]);

    Spinor sum[4];

    for (int isz=0; isz<2; isz++) {
        helicity sz=(isz==0)? neg : pos ;
        for (int la=-1; la<=1; la++) {
            double Cl=( la== 0 )? 1. : 2. ;

            sum[0]=sum[0]+Spinor(part,p,sz,mD,spin)*Complexe(sqrt(CSz)/sqrt(3.*Cl)*pla[la+1]/mD);

            for (int ix=1; ix<=3; ix++)
                sum[ix]=sum[ix]+Spinor(part,p,sz,mD,spin)*Complexe(sqrt(CSz)/sqrt(3.*Cl)*(
                            El[la+1].GetElem(ix-1)+pla[la+1]/mD/(p.GetE()+mD)*p.GetElem(ix) )*coef);
        }
    }

    for (int ix=0; ix<=3; ix++)
        SetElem(ix,sum[ix]);
}

Dspinor::Dspinor(const Lorentz_Vector<Spinor>& Ds)
{
    SetIndex(Ds.GetIndex());
    for (int ix=0; ix<=3; ix++)
        SetElem(ix,Ds.GetElem(ix));
}


Dspinor Dspinor::bar()
{
    return Lorentz_Vector<Spinor>(GetIndex(), GetElem(0).bar(), GetElem(1).bar(), GetElem(1).bar(), GetElem(3).bar() );
}

Spinor Dspinor::C(const C_Lorentz_Vector& n)
{
    if ( strcmp(GetName(),n.GetName())!=0 ) error("Bad Index name in Dspinor Lorentz_Vector contraction\n");
    if ( GetType()==n.GetType() ) error("Bad Index type in Dspinor Lorentz_Vector contraction\n");
    return GetElem(0)*n.GetElem(0) + GetElem(1)*n.GetElem(1) +  GetElem(2)*n.GetElem(2) + GetElem(3)*n.GetElem(3);
}

