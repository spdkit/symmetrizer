/**********************************************************************
vector3.cpp - Handle 3D coordinates.
 
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

#include <iostream>

#include "vector3.h"
#include "matrix3x3.h"



using namespace std;






int vector3::operator== ( const vector3& other ) const
{
    return ( ( x() == other.x() ) &&
             ( y() == other.y() ) &&
             ( z() == other.z() ) );
}


bool vector3::IsApprox(const vector3 & other, const double & precision) const
{
    return( distSq( other )
            <= precision * precision
            * std::min( length_2(), other.length_2() ) );
}



bool vector3::CanBeNormalized () const
{
    if( _vx == 0.0 && _vy == 0.0 && _vz == 0.0 ) return false;
    return( CanBeSquared(_vx)
            && CanBeSquared(_vy)
            && CanBeSquared(_vz) );
}




vector3 & vector3 :: normalize ()
{
    if( CanBeNormalized() )
        (*this) /= length();
    else
        (*this) /= length();

    return(*this);
}

double vector3 ::norm()
{
 return length();
}


vector3 cross ( const vector3& v1, const vector3& v2 )
{
    vector3 vv ;

    vv.x() =   v1.y()*v2.z() - v1.z()*v2.y() ;
    vv.y() = - v1.x()*v2.z() + v1.z()*v2.x() ;
    vv.z() =   v1.x()*v2.y() - v1.y()*v2.x() ;

    return ( vv ) ;
}



double vectorDistance ( const vector3& v1, const vector3& v2 )
{
    vector3 vv = v2-v1;
    return vv.length();
}

double vectorAngleRadian ( const vector3& v1, const vector3& v2 )
{
    double dp;

    dp = dot(v1,v2)/ ( v1.length() * v2.length() );


    if (dp < -0.999999)
      dp = -0.9999999;

    if (dp > 0.9999999)
      dp = 0.9999999;

     return acos(dp);
}

double vectorAngle ( const vector3& v1, const vector3& v2 )
{
    double dp;

    dp = dot(v1,v2)/ ( v1.length() * v2.length() );


    if (dp < -0.999999)
      dp = -0.9999999;

    if (dp > 0.9999999)
      dp = 0.9999999;

    return((RAD_TO_DEG * acos(dp)));
}

double vectorAngle ( const vector3& v1, const vector3& v2, const vector3& v3 )
{
	vector3 vv1 = v2-v1,vv2=v3-v2;
	return vectorAngle(vv1,vv2);
}


double CalcTorsionAngle(const vector3 &a, const vector3 &b,
                       const vector3 &c, const vector3 &d)
{
    double torsion;
    vector3 b1,b2,b3,c1,c2,c3;

    b1 = a - b;
    b2 = b - c;
    b3 = c - d;

    c1 = cross(b1,b2);
    c2 = cross(b2,b3);
    c3 = cross(c1,c2);


    if (c1.length() * c2.length() < 0.001)
    {
      torsion = 0.0;
      return torsion;
    }


    torsion = vectorAngle(c1,c2);
    if (dot(b2,c3) > 0.0)
      torsion = -torsion;

    return(torsion);
  }




bool vector3::createOrthoVector(vector3 &res) const
{

    // sanity check
    if( ! CanBeNormalized() )  {
        res=VZero;
        return false;
    }

    if( ! IsNegligible( _vx, _vz ) || ! IsNegligible( _vy, _vz ) )
    {
      double norm = sqrt( _vx*_vx + _vy*_vy );
      res._vx = -_vy / norm;
      res._vy = _vx / norm;
      res._vz = 0.0;
    }

    else
    {
      double norm = sqrt( _vy*_vy + _vz*_vz );
      res._vx = 0.0;
      res._vy = -_vz / norm;
      res._vz = _vy / norm;
    }

    return true;
  }






//! Calculate the distance of point a to the plane determined by b,c,d */
double Point2Plane(vector3 a, vector3 b, vector3 c, vector3 d)
{
    vector3 v_ba = a-b;
    vector3 v_normal = cross(c-b, d-b);
    return fabs( dot( v_normal, v_ba ) / v_normal.length() );
}

//! Calculate the angle between point a and the plane determined by b,c,d */
double Point2PlaneAngle(const vector3 a, const vector3 b, const vector3 c, const vector3 d)
{
    vector3 ac, bc, cd, normal;
    double angle;

    ac = a - c;
    bc = b - c;
    cd = c - d;
 
    normal = cross(bc, cd);
    angle = 90.0 - vectorAngle(normal, ac);

    return angle;
}


ostream& operator<< ( ostream& co, const vector3& v )
{
  co << "  " << v.x() << "  " << v.y() << "  " << v.z() << " " ;
  return co ;
}

const vector3 VZero ( 0.0, 0.0, 0.0 ) ;
const vector3 VX    ( 1.0, 0.0, 0.0 ) ;
const vector3 VY    ( 0.0, 1.0, 0.0 ) ;
const vector3 VZ    ( 0.0, 0.0, 1.0 ) ;








