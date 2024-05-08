/**********************************************************************
vector3.h - Handle 3D coordinates.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2006 by Benoit Jacob

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef VECTOR_H
#define VECTOR_H
#include <ostream>
#include <math.h>


#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180.0/M_PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI/180.0)
#endif



class matrix3x3;

// class introduction in vector3.cpp
class	 vector3
{
private :
    double _vx, _vy, _vz ;

public :
    vector3 (const double inX=0.0, const double inY=0.0, const double inZ=0.0):
        _vx(inX), _vy(inY), _vz(inZ)
    {}

    //! Copy Constructor
    vector3 (const vector3& v): _vx(v._vx), _vy(v._vy), _vz(v._vz)
    {}

    vector3 (float c[3]):  _vx(c[0]), _vy(c[1]), _vz(c[2])
    {}



    //! Set x,y and z-component of a vector
    void Set(const double x, const double y, const double z)
    {
        _vx = x;
        _vy = y;
        _vz = z;
    }

    void Set(unsigned int idx, double v)
    {
        if(idx==0) _vx =v;
        if(idx==1) _vy =v;
        if(idx==2) _vz =v;
    }

    //! Set x,y and z-component of a vector from c[0]..c[2]
    void Set(const double *c)
    {
        _vx = c[0];
        _vy = c[1];
        _vz = c[2];
    }

    //! Access function to get the x-coordinate of the vector
    void SetX(const double x)
    {
        _vx = x;
    }

    //! Access function to get the y-coordinate of the vector
    void SetY(const double y)
    {
        _vy = y;
    }

    //! Access function to get the z-coordinate of the vector
    void SetZ(const double z)
    {
        _vz = z;
    }


    //! \brief Set c[0]..c[2] to the components of the vector
    //! \warning No error checking is performed
    void Get(float *c)
    {
        c[0]=_vx;
        c[1]=_vy;
        c[2]=_vz;
    }



    //! Assignment
    vector3& operator= ( const vector3& v)
    {
        _vx = v._vx;
        _vy = v._vy;
        _vz = v._vz;
        return *this;
    }

    //! \return the vector as a const double *
    const double *AsArray() const
    {
        return &_vx;
    }


    //! \return *this + v
    vector3& operator+= ( const vector3& v)
    {
        _vx += v._vx;
        _vy += v._vy;
        _vz += v._vz;
        return *this;
    }


    //! \return *this - v
    vector3& operator-= ( const vector3& v)
    {
        _vx -= v._vx;
        _vy -= v._vy;
        _vz -= v._vz;
        return *this;
    };

    //! \return *this + f
    vector3& operator+= ( const double* f)
    {
        _vx += f[0];
        _vy += f[1];
        _vz += f[2];
        return *this;
    };

    //! \return *this - f
    vector3& operator-= ( const double* f)
    {
        _vx -= f[0];
        _vy -= f[1];
        _vz -= f[2];
        return *this;
    };


    //! \return *this * c
    vector3& operator*= ( const double& c)
    {
        _vx *= c;
        _vy *= c;
        _vz *= c;
        return *this;
    };


    //! \return *this divided by c
    vector3& operator/= ( const double& c)
    {
        double inv = 1.0 / c;
        return( (*this) *= inv );
    };

    //! Multiplication of matrix and vector
    //! \return the result (i.e., the updated vector)
    vector3& operator*= ( const matrix3x3 &);

    //! Create a random unit vector
    //   void randomUnitVector(OBRandom *oeRand= NULL);


    //! Scales a vector to give it length one.
    //! \return the result (i.e., the normalized vector)
    vector3 & normalize () ;
    double norm () ;


    bool CanBeNormalized () const;


    double length_2 () const
    {
        return _vx*_vx + _vy*_vy + _vz*_vz;
    };


    double length () const
    {
        return sqrt( length_2() );
    };


    const double & x () const
    {
        return _vx ;
    } ;


    const double & y () const
    {
        return _vy ;
    } ;


    const double & z () const
    {
        return _vz ;
    } ;


    double & x ()
    {
        return _vx ;
    }


    double & y ()
    {
        return _vy ;
    }


    double & z ()
    {
        return _vz ;
    }



    double operator[] ( unsigned int i) const
    {
        if (i > 2)
            return 0.0;


        if (i == 0)    return _vx;
        if (i == 1)    return _vy;
        else return _vz;
    }




    int operator== ( const vector3& ) const;


    //!    Use vector3::IsApprox()
    int operator!= ( const vector3& other ) const
    {
        return ! ( (*this) == other );
    }



    friend bool isEqual (const vector3& v1, const vector3& v2) {
        if ( (v1-v2).length() < 0.1 ) return true;
        return false;
    }

    friend bool isNegativeEqual (const vector3& v1, const vector3& v2) {
        if ( (v1+v2).length() < 0.1) return true;
        return false;
    }




    friend  vector3 operator+ ( const vector3& v1, const vector3& v2)
    {
        return vector3(v1.x()+v2.x(), v1.y()+v2.y(), v1.z()+v2.z());
    };




    friend  vector3 operator- ( const vector3& v1, const vector3& v2)
    {
        return vector3(v1.x()-v2.x(), v1.y()-v2.y(), v1.z()-v2.z());
    };


    friend  vector3 operator- ( const vector3& v)
    {
        return vector3(-v.x(), -v.y(), -v.z());
    };


    friend  vector3 operator * ( const double& c, const vector3& v)
    {
        return vector3( c*v.x(), c*v.y(), c*v.z());
    };


    friend  vector3 operator* ( const vector3& v, const double& c)
    {
        return vector3( c*v.x(), c*v.y(), c*v.z());
    };


    friend  vector3 operator/ ( const vector3& v, const double& c)
    {
        return vector3( v.x()/c, v.y()/c, v.z()/c);
    };


    friend  double dot ( const vector3& v1, const vector3& v2 )
    {
        return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z() ;
    }

    //! A random access iterator over x, y, z
    typedef double* iterator;

    //! A random access iterator over const x, y, z
    typedef const double* const_iterator;

    //! A signed integral type for differences between two iterators
    typedef std::ptrdiff_t difference_type;
    //! \return iterator to beginning
    iterator begin() { return &_vx; }

    //! \return iterator to end
    iterator end() { return &_vx + 3; }

    //! /return const_iterator to beginning
    const_iterator begin() const { return &_vx; }

    //! /return const_iterator to end
    const_iterator end() const { return &_vx + 3; }


    //!vector and matrix operation
    friend vector3 operator *(const vector3 &v,const matrix3x3 &m);
    friend vector3 operator *(const matrix3x3 &m,const vector3 &v);

    //! Cross product of two vectors
    friend vector3 cross ( const vector3&, const vector3& );


    friend double vectorDistance ( const vector3& v1, const vector3& v2 );

    //! Calculate the angle between vectors (in degrees)
    friend double vectorAngle ( const vector3& v1, const vector3& v2 );
    friend double vectorAngle ( const vector3& v1, const vector3& v2, const vector3& v3 );

    friend double vectorAngleRadian ( const vector3& v1, const vector3& v2 );



    //! Calculate the torsion angle between vectors (in degrees)
    friend double TorsionAngle(const vector3 &a, const vector3 &b,const vector3 &c, const vector3 &d);

    //! Calculate the distance of point a to the plane determined by b,c,d
    friend double Point2Plane(vector3 a, vector3 b, vector3 c, vector3 d);

    //! Calculate the angle between point a and the plane determined by b,c,d
    friend double Point2PlaneAngle(const vector3 a, const vector3 b, const vector3 c, const vector3 d);


    //! Safe comparison for doubleing-point vector3
    bool IsApprox( const vector3 & other, const double & precision ) const;


    //! \return square of the distance between *this and vv
    //! equivalent to length_2(*this-vv)

    double distSq(const vector3 &vv) const
    {
        double dx = x() - vv.x();
        double dy = y() - vv.y();
        double dz = z() - vv.z();
        return( dx*dx + dy*dy + dz*dz );
    }

    double distance(const vector3 &vv) const
    {
        return sqrt(distSq(vv));
    }

    //! Creates a vector of length one, orthogonal to *this.
    //! \return Whether the method was successful
    bool createOrthoVector(vector3 &v) const;



    //a, b, c of cell from scalar parameters
    void abcFromCellParameters(vector3 & a, vector3 & b, vector3 & c,
                               double a0,float b0, double c0,
                               double alpha, double beta, double gamma)
    {
        double pi=3.141592654;
        double cosA=cos(alpha*pi/180.0);
        double cosB=cos(beta*pi/180.0);
        double cosG=cos(gamma*pi/180.0);

        double sinA=sqrt(1.0-cosA*cosA);
        double sinB=sqrt(1.0-cosB*cosB);
        double sinG=sqrt(1.0-cosG*cosG);


        double V=a0*b0*c0*sqrt(1.0 -cosA*cosA -cosB*cosB -cosG*cosG +2.0*cosA*cosB*cosG);

        a.Set(a0, 0.0, 0.0);
        b.Set(b0*cosG, b0*sinG, 0.0);
        c.Set(c0*cosB, c0*(cosA-cosB*cosG)/sinG,  V/(a0*b0*sinG));

    }

    /*
    double areaVector3(vector3  a, vector3  b, vector3  c )
    {
        double l1 = (b-a).length();

        double l2 = (c-b).length();

        double l3 = (a-c).length();

        double p = (l1 + l2 + l3) / 2.0;

        return sqrt(p * (p - l1) * (p - l2) * (p - l3));
    }

*/


};





//! Prints a representation of the vector as a row vector of the form "<0.1,1,2>"
std::ostream& operator<< ( std::ostream&, const vector3& );

//! The zero vector: <0.0, 0.0, 0.0>
extern  const vector3 VZero;

//! The x unit vector: <1.0, 0.0, 0.0>
extern  const vector3 VX;

//! The y unit vector: <0.0, 1.0, 0.0>
extern  const vector3 VY;

//! The z unit vector: <0.0, 0.0, 1.0>
extern  const vector3 VZ;

#endif

















