#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <cstddef>
#include <math.h>

typedef unsigned int uint;
typedef unsigned long int ulong;

#define UNUSED(x) (void)(x)

//single XYZ for use as a temporary and return type
struct XYZ
{
   double x, y, z;

   XYZ() : x(0.0), y(0.0), z(0.0) {}
   XYZ(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}

   XYZ& operator=(XYZ const& rhs) 
   { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }
   XYZ& operator+=(XYZ const& rhs)
   { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
   XYZ& operator-=(XYZ const& rhs)
   { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
   XYZ& operator*=(XYZ const& rhs)
   { x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this; }
   XYZ& operator/=(XYZ const& rhs)
   { x /= rhs.x; y /= rhs.y; z /= rhs.z; return *this; }

   XYZ& operator*=(const double a)
   { x *= a; y *= a; z *= a; return *this; }

   XYZ& operator-=(const double a)
   { x -= a; y -= a; z -= a; return *this; }

  XYZ& operator+=(const double a)
   { x += a; y += a; z += a; return *this; }

   XYZ operator+(XYZ const& rhs) const
   { return XYZ(*this) += rhs; }
   XYZ operator-(XYZ const& rhs) const
   { return XYZ(*this) -= rhs; }
   XYZ operator*(XYZ const& rhs) const
   { return XYZ(*this) *= rhs; }
  XYZ operator/(XYZ const& rhs) const
   { return XYZ(*this) /= rhs; }


   XYZ operator*(const double a) const
   { return XYZ(*this) *= a; }

  XYZ operator-(const double a) const
   { return XYZ(*this) -= a; }

  XYZ operator+(const double a) const
   { return XYZ(*this) += a; }

   XYZ operator-() const { return XYZ(*this) * -1.0; }

   void Inverse()
  { 
     x = 1.0 / x; 
     y = 1.0 / y; 
     z = 1.0 / z; 
  }

   double Length() const { return sqrt(LengthSq()); }
   double LengthSq() const { return x * x + y * y + z * z; }
   XYZ& Normalize()
   {
      *this *= (1 / Length());
      return *this;
   }

  double DotProduct(XYZ const& a)
  {
    double sum = 0.0;
    sum += x * a.x;
    sum += y * a.y;
    sum += z * a.z;
    return sum;
  }

  //Calc AxB product 
  XYZ CrossProduct(const XYZ &A, const XYZ &B) const 
  {
    XYZ temp;
    temp.x = A.y * B.z - A.z * B.y;
    temp.y = A.z * B.x - A.x * B.z;
    temp.z = A.x * B.y - A.y * B.x;
  
return temp;
  }
};

#endif /*BASIC_TYPES_H*/
