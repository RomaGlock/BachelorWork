//
#include "math_our_cu.h"
namespace math_our_cu{
  __device__ Point_cu::Point_cu():x(0),y(0),z(0){}
  __device__ Point_cu::Point_cu(double a,double b,double c):x(a),y(b),z(c){}
  __device__ void Point_cu::set(double a,double b,double c){
    this->x=a;
    this->y=b;
    this->z=c;
  }
  
  __device__ double Point_cu::module() const{
    return norm3d(x,y,z);
  }
  __device__ Point_cu Point_cu::operator+(const Point_cu& p) const{
    return Point_cu(x+p.x, y+p.y, z+p.z);
  }
  __device__ Point_cu Point_cu::operator^(const Point_cu& p) const{
    return Point_cu(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);
  }
  __device__ Point_cu Point_cu::operator*(double p) const{
    return Point_cu(x*p, y*p, z*p);
  }
  __device__ Point_cu Point_cu::operator/(double p) const{
    return Point_cu(x/p, y/p, z/p);
  }
  __device__ const Point_cu & Point_cu::operator+=(const Point_cu& p){
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
  }
  __device__ const Point_cu & Point_cu::operator-=(const Point_cu& p){
    x -= p.x;
    y -= p.y;
    z -= p.z;
    return *this;
  }
  __device__ const Point_cu & Point_cu::operator/=(double p){
    x /= p;
    y /= p;
    z /= p;
    return *this;
  }
  __device__ const Point_cu & Point_cu::operator*=(double p){
    x *= p;
    y *= p;
    z *= p;
    return *this;
  }
  __device__ const Point_cu& Point_cu::normalize()
  {
    (*this) /= module();
    return *this;
  }
  __device__ Point_cu Point_cu::operator-(const Point_cu& p) const{
    return Point_cu(x-p.x, y-p.y, z-p.z);
  }
  __device__ double Point_cu::operator*(const Point_cu& p) const{
    return x*p.x+y*p.y+z*p.z;
  }
}
