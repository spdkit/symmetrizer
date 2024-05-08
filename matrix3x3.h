#pragma once

#include <ostream>
#include "vector3.h"

#include <algorithm>
#include <vector>


#ifndef M_PI
#define M_PI 3.141592654
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180.0/M_PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI/180.0)
#endif

#include <string>

using namespace std;
class matrix3x3
{
    double ele[3][3];

public:
      matrix3x3(void)
      {
          ele[0][0] = 1.0;      ele[0][1] = 0.0;       ele[0][2] = 0.0;
          ele[1][0] = 0.0;      ele[1][1] = 1.0;       ele[1][2] = 0.0;
          ele[2][0] = 0.0;      ele[2][1] = 0.0;       ele[2][2] = 1.0;
      }

      //! Constructs s times the unit matrix
      matrix3x3(double s)
      {
          ele[0][0] = s;            ele[0][1] = 0.0;        ele[0][2] = 0.0;
          ele[1][0] = 0.0;          ele[1][1] = s;          ele[1][2] = 0.0;
          ele[2][0] = 0.0;          ele[2][1] = 0.0;        ele[2][2] = s;
      }

      //! Constructs a matrix from row/column vectors, and type = 0 for column and 1 for row
      matrix3x3(vector3 v1,vector3 v2,vector3 v3, int type)
      {
          //from row
          if(type == 1) {
              ele[0][0] = v1.x();          ele[0][1] = v1.y();          ele[0][2] = v1.z();
              ele[1][0] = v2.x();          ele[1][1] = v2.y();          ele[1][2] = v2.z();
              ele[2][0] = v3.x();          ele[2][1] = v3.y();          ele[2][2] = v3.z();
          }

          // from column
          else {
              ele[0][0] = v1.x();          ele[0][1] = v1.x();          ele[0][2] = v1.x();
              ele[1][0] = v1.y();          ele[1][1] = v2.y();          ele[1][2] = v2.y();
              ele[2][0] = v1.z();          ele[2][1] = v3.z();          ele[2][2] = v3.z();
          }
      }

      matrix3x3(double d[3][3])
      {
          ele[0][0] = d[0][0];          ele[0][1] = d[0][1];          ele[0][2] = d[0][2];
          ele[1][0] = d[1][0];          ele[1][1] = d[1][1];          ele[1][2] = d[1][2];
          ele[2][0] = d[2][0];          ele[2][1] = d[2][1];          ele[2][2] = d[2][2];
      }

      matrix3x3(double d[9])
      {
        for (unsigned int i = 0; i < 3; ++ i)
            for (unsigned int j = 0; j < 3; ++ j)
                    ele[i][j]=d[i*3 + j];
      }


      matrix3x3(double m11,double m12,double m13,
                double m21,double m22,double m23,
                double m31,double m32,double m33)
      {
          ele[0][0] = m11;          ele[0][1] = m12;          ele[0][2] = m13;
          ele[1][0] = m21;          ele[1][1] = m22;          ele[1][2] = m23;
          ele[2][0] = m31;          ele[2][1] = m32;          ele[2][2] = m33;
      }


      matrix3x3(vector3 row1,vector3 row2,vector3 row3)
      {
          ele[0][0] = row1.x();          ele[0][1] = row1.y();          ele[0][2] = row1.z();

          ele[1][0] = row2.x();          ele[1][1] = row2.y();          ele[1][2] = row2.z();

          ele[2][0] = row3.x();          ele[2][1] = row3.y();          ele[2][2] = row3.z();
      }


      void GetArray(double *m)
      {
          m[0] = ele[0][0];          m[1] = ele[0][1];          m[2] = ele[0][2];
          m[3] = ele[1][0];          m[4] = ele[1][1];          m[5] = ele[1][2];
          m[6] = ele[2][0];          m[7] = ele[2][1];          m[8] = ele[2][2];
      }

     void UnitMatrix()
     {
          ele[0][0] = 1.0;          ele[0][1] = 0.0;        ele[0][2] = 0.0;
          ele[1][0] = 0.0;          ele[1][1] = 1.0;        ele[1][2] = 0.0;
          ele[2][0] = 0.0;          ele[2][1] = 0.0;        ele[2][2] = 1.0;
     }


      //!   return a constant reference to an element of the matrix.
      //!   row and column must be between 0 and 2. No check is done.
      const double & operator() (int row, int column ) const
      {
        return ele[row][column];
      }

      //!    return a non-constant reference to an element of the matrix.
      //!    row and column must be between 0 and 2. No check is done.
      double & operator() (int row, int column )
      {
          return ele[row][column];
      }

      //! Calculates the inverse of a matrix.
      matrix3x3 inverse(void) const;


      //! Calculates the transpose of a matrix.
      matrix3x3 transpose(void) const;

      //! Generates a matrix for a random rotation
//      void randomRotation(OBRandom &rnd);

      //! \return The determinant of the matrix
      double determinant() const;


      //! Checks if a matrix is orthogonal
      //! This method checks if a matrix is orthogonal, i.e.
      //!  if all column vectors are normalized and are mutually orthogonal. A matrix is orthogonal if,
      //!  and only if the transformation it describes is orthonormal.
      //!  An orthonormal transformation is a transformation that preserves length and angle.

      bool isOrthogonal(void) const
      {
          return (*this * transpose()).isUnitMatrix();
      }

      bool isSymmetric(void) const;
      bool isDiagonal(void) const;
      bool isUnitMatrix(void) const;


      double Get(int row,int column) const
      {
          return(ele[row][column]);
      }
      void Set(int row,int column, double v)
      {
         ele[row][column]= v;
      }


      void SetColumn(int column, const vector3 &v) ;
      void SetRow(int row, const vector3 &v);

      vector3 GetColumn(unsigned int col) const  ;
      vector3 GetRow(unsigned int row) const;

      //! Multiplies all entries of the matrix by a scalar c
      matrix3x3 &operator*=(const double &c)
      {
          for( int i = 0; i < 3; i++ )
              for( int j = 0; j < 3; j++ )
                  ele[i][j] *= c;
          return *this;
      }

      //! Divides all entries of the matrix by a scalar c
      matrix3x3 &operator/=(const double &c)
      {
         return( (*this) *= ( 1.0 / c ) );
      }



      //! Calculate an orthogonalisation matrix for a unit cell
      //! specified by the parameters alpha, beta, gamma, a, b, c
      //! where alpha, beta, and gamma are the cell angles (in degrees)
      //! and a, b, and c are the cell vector lengths
      //! Used by OBUnitCell
      void FillOrth(double alpha, double beta, double gamma,
                    double a, double b, double c);

      //! Find the eigenvalues and -vectors of a symmetric matrix
      matrix3x3 findEigenvectorsIfSymmetric(vector3 &eigenvals) const ;

      //! Matrix-vector multiplication
      friend  vector3 operator *(const matrix3x3 &,const vector3 &);
      friend  vector3 operator *(const vector3 &, const matrix3x3 &);

      //! Matrix-matrix multiplication
      friend  matrix3x3 operator *(const matrix3x3 &,const matrix3x3 &);


      //! Calculates a rotation matrix, rotating around the specified axis by the specified angle (in degrees)
      void rotateXAngle(const double );
      void rotateYAngle(const double );
      void rotateZAngle(const double );


      //! \brief Calculate a rotation matrix for rotation about the x, y, and z
      //! axes by the angles specified (in degrees)
      void SetupRotMatrix(double x, double y, double z);
      void SetupRotMatrix(vector3  vec, double  angle);
      void SetupRotMatrix(vector3  vec, unsigned int  n);
      void SetupRotMatrix(vector3  vec0, vector3 vec1);

      void SetupRotateMatrix(double x, double y, double z) {SetupRotMatrix(x,y,z);};
      void SetupRotateMatrix(vector3  vec, double  angle) {SetupRotMatrix(vec,angle);};
      void SetupRotateMatrix(vector3  vec, unsigned int  n) {SetupRotMatrix(vec,n);};

      void SetupReflectiontMatrix(vector3  );


      //! Calculates a matrix that represents reflection on a plane
      void reflectionXOY();
      void reflectionYOZ();
      void reflectionXOZ();
      void rotateZReflection(double & );
      void reflectionXOZwithAngle(double);



      //transforma matrix  from cartesian  to fractional
      void Cartesian2Fractional(double a0,double b0, double c0,
                                double alpha, double beta, double gamma);

      //transforma matrix  from fractional to cartesian
      void Fractional2Cartesian(double a0,double b0, double c0,
                                double alpha, double beta, double gamma);



      //! Output a text representation of a matrix
      friend  std::ostream& operator<< ( std::ostream&, const matrix3x3 & ) ;

      //! Eigenvalue calculation
      static void jacobi(unsigned int n, double *a, double *d, double *v);
};


class  matrix4
{
    double array1[16];

public:
      matrix4()
      {
         for (int i=0;i<16;i++) array1[i] =0.0;
      }
     ~matrix4(){};


    double & operator () ( int row, int col )
    {
        return array1[row + 4* col ]; //opengl right multiply
    }

    double * array()
    {
      return array1;
    }

};


bool IsApprox(const double & a, const double & b,
                             const double precision = 1e-11);
bool IsNegligible(const double & a, const double & b,
                                 const double precision = 1e-11);
bool CanBeSquared(const double &a);

vector3 rotateVector(vector3, double, string);
vector3 rotateVector(vector3, double, vector3);

matrix3x3 SetupRotationReflectiontMatrix(vector3 , unsigned int  );
matrix3x3 SetupRotationReflectiontMatrix(vector3 , double);













