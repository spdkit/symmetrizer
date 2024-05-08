#include "matrix3x3.h"
#include <iostream>

using namespace std;


void matrix3x3::SetupReflectiontMatrix(vector3  vec0)
{
    vec0.normalize();

    double  x = vec0[0],
            y = vec0[1],
            z = vec0[2];

    ele[0][0] = 1.0-2.0*x*x ;
    ele[0][1] = -2.0*x*y;
    ele[0][2] = -2.0*x*z;

    ele[1][0] = -2.0*x*y;
    ele[1][1] = 1.0-2.0*y*y;
    ele[1][2] = -2.0*y*z;

    ele[2][0] = -2.0*x*z;
    ele[2][1] = -2.0*y*z;
    ele[2][2] = 1.0-2.0*z*z;

    /*
      reflection = np.array([
          [1.0-2.0*x*x,    -2.0*x*y,      -2.0*x*z    ],
          [-2.0*x*y,      1.0-2.0*y*y,    -2.0*y*z    ],
          [-2.0*x*z,      -2.0*y*z,       1.0-2.0*z*z ]
      ])
      */
}


void matrix3x3::SetupRotMatrix(vector3  vec0, unsigned int  n /*order of axis*/)
{
    double angle= 360.0/double(n);
    SetupRotMatrix(vec0,angle);
}

void matrix3x3::SetupRotMatrix(vector3  vec0, double angle)
{
    double ang = angle*DEG_TO_RAD;

    double c = cos(ang), t=1.0-c,
            s = sin(ang);

    vector3 vec = vec0.normalize();

    //cout << angle << endl;
    //cout << vec <<endl;


    double  x = vec[0],
            y = vec[1],
            z = vec[2];

    ele[0][0] = x*x*t+c ;
    ele[0][1] = x*y*t-z*s;
    ele[0][2] = x*z*t+y*s;

    ele[1][0] = x*y*t+z*s;
    ele[1][1] = y*y*t+c;
    ele[1][2] = y*z*t-x*s;

    ele[2][0] = x*z*t-y*s;
    ele[2][1] = y*z*t+x*s;
    ele[2][2] = z*z*t+c;

    /*
      r = np.array([
          [x*x*(1.0-c)+c  , x*y*(1.0-c)-z*s, x*z*(1.0-c)+y*s],
          [x*y*(1.0-c)+z*s, y*y*(1.0-c)+c  , y*z*(1.0-c)-x*s],
          [x*z*(1.0-c)-y*s, y*z*(1.0-c)+x*s, z*z*(1.0-c)+c  ]
      ]);
  */

}


void matrix3x3::SetupRotMatrix(double phi,double theta,double psi)
{
    double p  = phi * DEG_TO_RAD;
    double h  = theta * DEG_TO_RAD;
    double b  = psi * DEG_TO_RAD;

    double cx = cos(p);
    double sx = sin(p);
    double cy = cos(h);
    double sy = sin(h);
    double cz = cos(b);
    double sz = sin(b);

    ele[0][0] = cy*cz;
    ele[0][1] = cy*sz;
    ele[0][2] = -sy;

    ele[1][0] = sx*sy*cz-cx*sz;
    ele[1][1] = sx*sy*sz+cx*cz;
    ele[1][2] = sx*cy;

    ele[2][0] = cx*sy*cz+sx*sz;
    ele[2][1] = cx*sy*sz-sx*cz;
    ele[2][2] = cx*cy;
}




void matrix3x3::SetColumn(int col, const vector3 &v)
{
    ele[0][col] = v.x();
    ele[1][col] = v.y();
    ele[2][col] = v.z();
}

void matrix3x3::SetRow(int row, const vector3 &v)
{
    ele[row][0] = v.x();
    ele[row][1] = v.y();
    ele[row][2] = v.z();
}

vector3 matrix3x3::GetColumn(unsigned int col) const
{
    return vector3(ele[0][col], ele[1][col], ele[2][col]);
}



vector3 matrix3x3::GetRow(unsigned int row) const
{
    return vector3(ele[row][0], ele[row][1], ele[row][2]);
}



//! Calculates the product m*v of the matrix m and the column vector represented by v
vector3 operator *(const matrix3x3 &m,const vector3 &v)
{
    vector3 vv;

    vv.x() = v.x()*m.ele[0][0] + v.y()*m.ele[0][1] + v.z()*m.ele[0][2];
    vv.y() = v.x()*m.ele[1][0] + v.y()*m.ele[1][1] + v.z()*m.ele[1][2];
    vv.z() = v.x()*m.ele[2][0] + v.y()*m.ele[2][1] + v.z()*m.ele[2][2];

    return(vv);
}



vector3 operator *(const vector3 &v, const matrix3x3 &m)
{
    vector3 vv;

    vv._vx = v._vx*m.ele[0][0] + v._vy*m.ele[1][0] + v._vz*m.ele[2][0];
    vv._vy = v._vx*m.ele[0][1] + v._vy*m.ele[1][1] + v._vz*m.ele[2][1];
    vv._vz = v._vx*m.ele[0][2] + v._vy*m.ele[1][2] + v._vz*m.ele[2][2];

    return(vv);
}



matrix3x3 operator *(const matrix3x3 &A,const matrix3x3 &B)
{
    matrix3x3 result;

    result.ele[0][0] = A.ele[0][0]*B.ele[0][0] + A.ele[0][1]*B.ele[1][0] + A.ele[0][2]*B.ele[2][0];
    result.ele[0][1] = A.ele[0][0]*B.ele[0][1] + A.ele[0][1]*B.ele[1][1] + A.ele[0][2]*B.ele[2][1];
    result.ele[0][2] = A.ele[0][0]*B.ele[0][2] + A.ele[0][1]*B.ele[1][2] + A.ele[0][2]*B.ele[2][2];

    result.ele[1][0] = A.ele[1][0]*B.ele[0][0] + A.ele[1][1]*B.ele[1][0] + A.ele[1][2]*B.ele[2][0];
    result.ele[1][1] = A.ele[1][0]*B.ele[0][1] + A.ele[1][1]*B.ele[1][1] + A.ele[1][2]*B.ele[2][1];
    result.ele[1][2] = A.ele[1][0]*B.ele[0][2] + A.ele[1][1]*B.ele[1][2] + A.ele[1][2]*B.ele[2][2];

    result.ele[2][0] = A.ele[2][0]*B.ele[0][0] + A.ele[2][1]*B.ele[1][0] + A.ele[2][2]*B.ele[2][0];
    result.ele[2][1] = A.ele[2][0]*B.ele[0][1] + A.ele[2][1]*B.ele[1][1] + A.ele[2][2]*B.ele[2][1];
    result.ele[2][2] = A.ele[2][0]*B.ele[0][2] + A.ele[2][1]*B.ele[1][2] + A.ele[2][2]*B.ele[2][2];

    return(result);
}



//! calculates the product m*(*this) of the matrix m and the column vector represented by *this
vector3 &vector3::operator *= (const matrix3x3 &m)
{
    vector3 vv;

    vv.SetX(x()*m.Get(0,0) + y()*m.Get(0,1) + z()*m.Get(0,2));
    vv.SetY(x()*m.Get(1,0) + y()*m.Get(1,1) + z()*m.Get(1,2));
    vv.SetZ(x()*m.Get(2,0) + y()*m.Get(2,1) + z()*m.Get(2,2));
    x() = vv.x();
    y() = vv.y();
    z() = vv.z();

    return(*this);
}


matrix3x3 matrix3x3::inverse(void) const
{
    double det = determinant();
    matrix3x3 inverse;
    if (fabs(det) <= 1e-6) return inverse;

    inverse.ele[0][0] = ele[1][1]*ele[2][2] - ele[1][2]*ele[2][1];
    inverse.ele[1][0] = ele[1][2]*ele[2][0] - ele[1][0]*ele[2][2];
    inverse.ele[2][0] = ele[1][0]*ele[2][1] - ele[1][1]*ele[2][0];
    inverse.ele[0][1] = ele[2][1]*ele[0][2] - ele[2][2]*ele[0][1];
    inverse.ele[1][1] = ele[2][2]*ele[0][0] - ele[2][0]*ele[0][2];
    inverse.ele[2][1] = ele[2][0]*ele[0][1] - ele[2][1]*ele[0][0];
    inverse.ele[0][2] = ele[0][1]*ele[1][2] - ele[0][2]*ele[1][1];
    inverse.ele[1][2] = ele[0][2]*ele[1][0] - ele[0][0]*ele[1][2];
    inverse.ele[2][2] = ele[0][0]*ele[1][1] - ele[0][1]*ele[1][0];

    inverse /= det;

    return(inverse);
}



matrix3x3 matrix3x3::transpose(void) const
{
    matrix3x3 transpose;

    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            transpose.ele[i][j] = ele[j][i];

    return(transpose);
}


double matrix3x3::determinant(void) const
{
    return( ele[0][0] * (ele[1][1] * ele[2][2] - ele[1][2] * ele[2][1])
            + ele[0][1] * (ele[1][2] * ele[2][0] - ele[1][0] * ele[2][2])
            + ele[0][2] * (ele[1][0] * ele[2][1] - ele[1][1] * ele[2][0]) );
}


//!    return False if there are indices i,j such that
//!    fabs(*this[i][j]-*this[j][i]) > 1e-6. Otherwise, it returns true
bool matrix3x3::isSymmetric(void) const
{
    return( IsApprox( ele[0][1], ele[1][0], 1e-6 )
            && IsApprox( ele[0][2], ele[2][0], 1e-6 )
            && IsApprox( ele[1][2], ele[2][1], 1e-6 ) );
}



bool matrix3x3::isDiagonal(void) const
{
    return( IsNegligible( ele[1][0], ele[0][0], 1e-6 )
            && IsNegligible( ele[2][0], ele[0][0], 1e-6 )
            && IsNegligible( ele[0][1], ele[1][1], 1e-6 )
            && IsNegligible( ele[2][1], ele[1][1], 1e-6 )
            && IsNegligible( ele[0][2], ele[2][2], 1e-6 )
            && IsNegligible( ele[1][2], ele[2][2], 1e-6 ) );
}


bool matrix3x3::isUnitMatrix(void) const
{
    return ( isDiagonal()
             && IsApprox( ele[0][0], 1.0, 1e-6 )
            && IsApprox( ele[1][1], 1.0, 1e-6 )
            && IsApprox( ele[2][2], 1.0, 1e-6 ) );
}





static inline double SQUARE( double x ) { return x*x; }

void matrix3x3::FillOrth(double Alpha,double Beta, double Gamma,
                         double A, double B, double C)
{
    double V;

    Alpha *= DEG_TO_RAD;
    Beta  *= DEG_TO_RAD;
    Gamma *= DEG_TO_RAD;

    // from the PDB specification:
    //  http://www.rcsb.org/pdb/docs/format/pdbguide2.2/part_75.html


    // since we'll ultimately divide by (a * b), we've factored those out here
    V = C * sqrt(1 - SQUARE(cos(Alpha)) - SQUARE(cos(Beta)) - SQUARE(cos(Gamma))
                 + 2 * cos(Alpha) * cos(Beta) * cos(Gamma));

    ele[0][0] = A;
    ele[0][1] = B*cos(Gamma);
    ele[0][2] = C*cos(Beta);

    ele[1][0] = 0.0;
    ele[1][1] = B*sin(Gamma);
    ele[1][2] = C*(cos(Alpha)-cos(Beta)*cos(Gamma))/sin(Gamma);

    ele[2][0] = 0.0;
    ele[2][1] = 0.0;
    ele[2][2] = V / (sin(Gamma)); // again, we factored out A * B when defining V
}




void matrix3x3::rotateXAngle(const double angle)
{
    double Angle=angle*DEG_TO_RAD;
    double cosAngle=cos(Angle),sinAngle=sin(Angle);

    ele[0][0] = 1.0;      ele[0][1] = 0.0;       	   ele[0][2] = 0.0;
    ele[1][0] = 0.0;      ele[1][1] = cosAngle;       ele[1][2] = -sinAngle;
    ele[2][0] = 0.0;      ele[2][1] = sinAngle;       ele[2][2] = cosAngle;
}

void matrix3x3::rotateYAngle(const double angle)
{
    double Angle=angle*DEG_TO_RAD;
    double cosAngle=cos(Angle),sinAngle=sin(Angle);

    ele[0][0] = cosAngle;     	 ele[0][1] = 0.0;       ele[0][2] = sinAngle;
    ele[1][0] = 0.0;      	     ele[1][1] = 1.0;       ele[1][2] = 0.0;
    ele[2][0] = -sinAngle;     	 ele[2][1] = 0.0;       ele[2][2] = cosAngle;

}

void matrix3x3::rotateZAngle(const double angle)
{
    double Angle=angle*DEG_TO_RAD;
    double cosAngle=cos(Angle),sinAngle=sin(Angle);

    ele[0][0] = cosAngle;     ele[0][1] = -sinAngle;      ele[0][2] = 0.0;
    ele[1][0] = sinAngle;     ele[1][1] = cosAngle;       ele[1][2] = 0.0;
    ele[2][0] = 0.0;    	  ele[2][1] = 0.0;            ele[2][2] = 1.0;
}


//----------------------------------------------------------
vector3 rotateVector(vector3 from, double angle, string Axis)
{
    vector3 ResultVector;
    ResultVector.Set(0.0f,0.0f,0.0f);

    matrix3x3 OperateMatrix;

    if(Axis.find("-"))   angle=-angle;

    transform(Axis.begin(),Axis.end(),Axis.begin(),::tolower);

    if (Axis.find("x"))
        OperateMatrix.rotateXAngle(angle);

    if (Axis.find("y"))
        OperateMatrix.rotateYAngle(angle);

    if (Axis.find("z"))
        OperateMatrix.rotateZAngle(angle);

    ResultVector=OperateMatrix*from;

    return ResultVector;
}


vector3 rotateVector(vector3 from, double angle, vector3 axis)
{
    vector3 resultVector;resultVector.Set(0.0f,0.0f,0.0f);
    matrix3x3 OperateMatrix;

    OperateMatrix.SetupRotMatrix(axis,angle);
    resultVector=OperateMatrix*from;
    return resultVector;
}

//get rotation from vec0->vec1
void matrix3x3::SetupRotMatrix(vector3 vec0, vector3 vec1)
{
    vec0.normalize();
    vec1.normalize();

    if(isEqual(vec0,vec1)) {
        this->UnitMatrix();
        return;
    }

    if(isNegativeEqual(vec0,vec1)) {
        this->UnitMatrix();
        Set(0,0,-1.0);
        Set(1,1,-1.0);
        Set(2,2,-1.0);
        return;
    }

    double angle=vectorAngle(vec0,vec1);

    if(angle<2.0) {
        this->UnitMatrix();
        return;
    }

    vector3 axis=cross(vec0,vec1);
    this->SetupRotMatrix(axis,angle);
}






void matrix3x3::reflectionXOY()
{
    ele[0][0] = 1.0; ele[0][1] = 0.0; ele[0][2] = 0.0;
    ele[1][0] = 0.0; ele[1][1] = 1.0; ele[1][2] = 0.0;
    ele[2][0] = 0.0; ele[2][1] = 0.0; ele[2][2] = -1.0;
}

void matrix3x3::reflectionYOZ()
{
    ele[0][0] = -1.0;  ele[0][1] = 0.0;  ele[0][2] = 0.0;
    ele[1][0] =  0.0;  ele[1][1] = 1.0;  ele[1][2] = 0.0;
    ele[2][0] =  0.0;  ele[2][1] = 0.0;  ele[2][2] = 1.0;
}

void matrix3x3::reflectionXOZ()
{
    ele[0][0] =  1.0;  ele[0][1] =  0.0;  ele[0][2] = 0.0;
    ele[1][0] =  0.0;  ele[1][1] = -1.0;  ele[1][2] = 0.0;
    ele[2][0] =  0.0;  ele[2][1] =  0.0;  ele[2][2] = 1.0;
}



// rotate around Z axis and reflection by XOY
void matrix3x3::rotateZReflection(double & angle)
{
    double Angle=angle*DEG_TO_RAD;
    double cosAngle=cos(Angle),sinAngle=sin(Angle);

    ele[0][0] = cosAngle; ele[0][1] = -sinAngle; ele[0][2] = 0.0;
    ele[1][0] = sinAngle; ele[1][1] =  cosAngle; ele[1][2] = 0.0;
    ele[2][0] = 0.0;      ele[2][1] =  0.0;      ele[2][2] =-1.0;
}


void matrix3x3::Cartesian2Fractional(double a0,double b0, double c0,
                                     double alpha, double beta, double gamma)
{
    double pi=3.141592654;
    double cosA=cos(alpha*DEG_TO_RAD);
    double cosB=cos(beta*DEG_TO_RAD);
    double cosG=cos(gamma*DEG_TO_RAD);

    double sinA=sqrt(1.0-cosA*cosA);
    double sinB=sqrt(1.0-cosB*cosB);
    double sinG=sqrt(1.0-cosG*cosG);


    double V=a0*b0*c0*sqrt(1.0 - cosA*cosA - cosB*cosB - cosG*cosG + 2.0*cosA*cosB*cosG);

    ele[0][0] = 1.0/a0;    ele[0][1] = -cosG/(a0*sinG);     ele[0][2] = b0*c0*(cosA*cosG-cosB)/(V*sinG);
    ele[1][0] = 0.0;       ele[1][1] = 1.0/(b0*sinG);       ele[1][2] = a0*c0*(cosB *cosG-cosA)/(V*sinG);
    ele[2][0] = 0.0;       ele[2][1] = 0.0;                 ele[2][2] = a0*b0* sinG/V;
}


void matrix3x3::Fractional2Cartesian(double a0,double b0, double c0,
                                     double alpha, double beta, double gamma)
{
    double pi=3.141592654;
    double cosA=cos(alpha*pi/180.0);
    double cosB=cos(beta*pi/180.0);
    double cosG=cos(gamma*pi/180.0);

    //double sinA=sqrt(1.0-cosA*cosA);
    //double sinB=sqrt(1.0-cosB*cosB);
    double sinG=sqrt(1.0-cosG*cosG);


    double V=a0*b0*c0*sqrt(1.0 -cosA*cosA -cosB*cosB -cosG*cosG + 2.0*cosA*cosB*cosG);

    ele[0][0] = a0;    ele[0][1] = b0*cosG;     ele[0][2] = c0*cosB;
    ele[1][0] = 0.0;   ele[1][1] = b0*sinG;     ele[1][2] = c0*(cosA- cosB*cosG)/sinG;
    ele[2][0] = 0.0;   ele[2][1] = 0.0;         ele[2][2] = V/(a0*b0*sinG);
}


ostream& operator<< ( ostream& co, const matrix3x3& m )
{
    co << "[ "
       << m.ele[0][0]
            << ", "
            << m.ele[0][1]
            << ", "
            << m.ele[0][2]
            << " ]" << endl;

    co << "[ "
       << m.ele[1][0]
            << ", "
            << m.ele[1][1]
            << ", "
            << m.ele[1][2]
            << " ]" << endl;

    co << "[ "
       << m.ele[2][0]
            << ", "
            << m.ele[2][1]
            << ", "
            << m.ele[2][2]
            << " ]" << endl;

    return co ;
}



matrix3x3 matrix3x3::findEigenvectorsIfSymmetric(vector3 &eigenvals) const
{
    matrix3x3 result;

    double d[3];
    matrix3x3 copyOfThis = *this;

    jacobi(3, copyOfThis.ele[0], d, result.ele[0]);
    eigenvals.Set(d);

    return result;
}



void matrix3x3::jacobi(unsigned int n, double *a, double *d, double *v)
{
    double onorm, dnorm;
    double b, dma, q, t, c, s;
    double  atemp, vtemp, dtemp;
    int i, j, k, l;
    int nrot;

    int MAX_SWEEPS=50;

    // Set v to the identity matrix, set the vector d to contain the
    // diagonal elements of the matrix a
    for (j = 0; j < static_cast<int>(n); j++)
    {
        for (i = 0; i < static_cast<int>(n); i++)
            v[n*i+j] = 0.0;
        v[n*j+j] = 1.0;
        d[j] = a[n*j+j];
    }

    nrot = MAX_SWEEPS;
    for (l = 1; l <= nrot; l++)
    {
        // Set dnorm to be the maximum norm of the diagonal elements, set
        // onorm to the maximum norm of the off-diagonal elements
        dnorm = 0.0;
        onorm = 0.0;
        for (j = 0; j < static_cast<int>(n); j++)
        {
            dnorm += (double)fabs(d[j]);
            for (i = 0; i < j; i++)
                onorm += (double)fabs(a[n*i+j]);
        }
        // Normal end point of this algorithm.
        if((onorm/dnorm) <= 1.0e-12)
            goto Exit_now;

        for (j = 1; j < static_cast<int>(n); j++)
        {
            for (i = 0; i <= j - 1; i++)
            {

                b = a[n*i+j];
                if(fabs(b) > 0.0)
                {
                    dma = d[j] - d[i];
                    if((fabs(dma) + fabs(b)) <=  fabs(dma))
                        t = b / dma;
                    else
                    {
                        q = 0.5 * dma / b;
                        t = 1.0/((double)fabs(q) + (double)sqrt(1.0+q*q));
                        if (q < 0.0)
                            t = -t;
                    }

                    c = 1.0/(double)sqrt(t*t + 1.0);
                    s = t * c;
                    a[n*i+j] = 0.0;

                    for (k = 0; k <= i-1; k++)
                    {
                        atemp = c * a[n*k+i] - s * a[n*k+j];
                        a[n*k+j] = s * a[n*k+i] + c * a[n*k+j];
                        a[n*k+i] = atemp;
                    }

                    for (k = i+1; k <= j-1; k++)
                    {
                        atemp = c * a[n*i+k] - s * a[n*k+j];
                        a[n*k+j] = s * a[n*i+k] + c * a[n*k+j];
                        a[n*i+k] = atemp;
                    }

                    for (k = j+1; k < static_cast<int>(n); k++)
                    {
                        atemp = c * a[n*i+k] - s * a[n*j+k];
                        a[n*j+k] = s * a[n*i+k] + c * a[n*j+k];
                        a[n*i+k] = atemp;
                    }

                    for (k = 0; k < static_cast<int>(n); k++)
                    {
                        vtemp = c * v[n*k+i] - s * v[n*k+j];
                        v[n*k+j] = s * v[n*k+i] + c * v[n*k+j];
                        v[n*k+i] = vtemp;
                    }

                    dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
                    d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
                    d[i] = dtemp;
                } /* end if */
            } /* end for i */
        } /* end for j */
    } /* end for l */

Exit_now:

    // Now sort the eigenvalues (and the eigenvectors) so that the
    // smallest eigenvalues come first.
    nrot = l;

    for (j = 0; j < static_cast<int>(n)-1; j++)
    {
        k = j;
        dtemp = d[k];
        for (i = j+1; i < static_cast<int>(n); i++)
            if(d[i] < dtemp)
            {
                k = i;
                dtemp = d[k];
            }

        if(k > j)
        {
            d[k] = d[j];
            d[j] = dtemp;
            for (i = 0; i < static_cast<int>(n); i++)
            {
                dtemp = v[n*i+k];
                v[n*i+k] = v[n*i+j];
                v[n*i+j] = dtemp;
            }
        }
    }
}




bool IsApprox(const double & a, const double & b,  const double precision)
{
    return( fabs(a - b) <= precision * std::min( fabs(a), fabs(b) ) );
}
bool IsNegligible(const double & a, const double & b,  const double precision)
{
    return( fabs(a) <= precision * fabs(b) );
}

bool CanBeSquared(const double &a)
{
    if( a == 0 ) return true;
    const double max_squarable_double = 1e150;
    const double min_squarable_double = 1e-150;
    double abs_a = fabs(a);
    return(abs_a < max_squarable_double && abs_a > min_squarable_double);
}

matrix3x3 SetupRotationReflectiontMatrix(vector3 vec0, unsigned int n)
{
    //be careful here!!!
    double angle = 360.0/double(n);
    return SetupRotationReflectiontMatrix (vec0,angle);
}

matrix3x3 SetupRotationReflectiontMatrix(vector3 vec0, double angle)
{
    matrix3x3 reflection;
    reflection.SetupReflectiontMatrix(vec0);

    matrix3x3 rotation;
    rotation.SetupRotMatrix(vec0,angle);

    return  reflection*rotation;
}

// reflection by Mirror which parallel to XOZ,
// with angle from XOZ (anti_clock)
void matrix3x3::reflectionXOZwithAngle(double angle)
{
    double ang = angle*DEG_TO_RAD;
    vector3 axis,temp;
    axis.Set(cos(ang),sin(ang),0.0);
    temp=cross(axis,VZ);
    matrix3x3 result;
    result.SetupReflectiontMatrix(temp);
    *this =  result;
}




