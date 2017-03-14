#include <math.h>
#include <vector>
#include<stdio.h>
#include <stdlib.h>
#include "math_our.h"


namespace math_our{
  double TOL = 1e-7;
  Point::Point(double x,double y,double z) : x(x), y(y), z(z) {
  }
  
  void Point::set(double a,double b,double c){
    this->x=a;
    this->y=b;
    this->z=c;
  }
  
  double Point::module() const{
    return sqrt(x*x+y*y+z*z);
  }
  double Point::value() const{
    return x*x+y*y+z*z+x+y+z+2;
  }
  bool Point::operator<(const Point& p) const{
    if((*this)==p) return false;    
    else{
      if(value()<p.value()) return true;
      else return false;
    }
  }
  bool Point::operator<(const double p) const{
    if(value()<p) return true;
    else return false;
  }
  bool cmpDouble(double a,double b){
    long c;
    c=0;
    if(fabs(a)<TOL && fabs(b)<TOL) c=1;
    if(fabs(1-a/b)<TOL && fabs(1-b/a)<TOL) c=1;
    if(c){
      return true;
    }else{
      return false;
    }
  }
  bool Point::operator==(const Point& p) const{
    long c1,c2,c3;
    c1=0;
    c2=0;
    c3=0;
    if(fabs(p.x)<TOL && fabs(x)<TOL) c1=1;
    if(fabs(p.y)<TOL && fabs(y)<TOL) c2=1;
    if(fabs(p.z)<TOL && fabs(z)<TOL) c3=1;
    if(fabs(1-p.x/x)<TOL) c1=1;
    if(fabs(1-p.y/y)<TOL) c2=1;
    if(fabs(1-p.z/z)<TOL) c3=1;
    if(c1*c2*c3>0)
      return true;
    else 
      return false;
  }
  Point Point::operator+(const Point& p) const{
    return Point(x+p.x, y+p.y, z+p.z);
  }
  Point Point::operator^(const Point& p) const{
    return Point(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);
  }
  Point Point::operator*(double p) const{
    return Point(x*p, y*p, z*p);
  }
  Point Point::operator/(double p) const{
    return Point(x/p, y/p, z/p);
  }
  const Point & Point::operator+=(const Point& p){
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
  }
  const Point & Point::operator-=(const Point& p){
    x -= p.x;
    y -= p.y;
    z -= p.z;
    return *this;
  }
  const Point & Point::operator/=(double p){
    x /= p;
    y /= p;
    z /= p;
    return *this;
  }
  const Point & Point::operator*=(double p){
    x *= p;
    y *= p;
    z *= p;
    return *this;
  }
  const Point& Point::normalize()
  {
    (*this) /= module();
    return *this;
  }
  Point Point::operator-(const Point& p) const{
    return Point(x-p.x, y-p.y, z-p.z);
  }
  double Point::operator*(const Point& p) const{
    return x*p.x+y*p.y+z*p.z;
  }
  double mix_mul(const Point& a,const Point& b,const Point& c){
    return a.x*(b.y*c.z-c.y*b.z)-a.y*(b.x*c.z-c.x*b.z)+a.z*(b.x*c.y-c.x*b.y);
  }
  std::ostream& operator << (std::ostream& os, const Point& p)
  {
    return os  << p.x << " " << p.y << " " << p.z;
  }  
  struct PointTreeCompare{
    bool operator()(const math_our::Point& p1, const math_our::Point& p2){
      if(p1==p2){
	return false;
      }
      bool result = true;
      if(math_our::cmpDouble(p1.x,p2.x)){
	if(math_our::cmpDouble(p1.y,p2.y)){
	  if(math_our::cmpDouble(p1.z,p2.z)){
	    result=false;
	  }else{
	    if(p1.z > p2.z){
	      result = false;
	    }
	  }
	}else{
	  if(p1.y > p2.y){
	    result = false;
	  }
	}
      }else{
	if(p1.x > p2.x){
	  result = false;
	}
      }
      return result;
    }
  };
  Map::~Map(){
    if(!root) return;
    std::vector<element*> stack;
    stack.reserve(size_);
    stack.push_back(root);
    ulong i=0;
    while(i<stack.size()){
      if(stack[i]->left) stack.push_back(stack[i]->left);
      if(stack[i]->right) stack.push_back(stack[i]->right);
      free(stack[i]);
      ++i;
    }
  }
  std::pair<Map::iterator,bool>Map::insert(const std::pair<Point,long>& p){
    if(!root){
      root=add_element(p);
      return std::pair<Map::iterator,bool>(&(root->d),true);
    }
    element* c;
    c=root;
    while(c){
      if(PointTreeCompare()(p.first,c->d.first)){
	if(c->left){
	  c=c->left;
	}else{
	  c->left=add_element(p);
	  return std::pair<Map::iterator,bool>(&(c->left->d),true);
	}
      }else if(PointTreeCompare()(c->d.first,p.first)){
	if(c->right){
	  c=c->right;
	}else{
	  c->right=add_element(p);
	  return std::pair<Map::iterator,bool>(&(c->right->d),true);
	}
      }else return std::pair<Map::iterator,bool>(&(c->d),false);
    }
    return std::pair<Map::iterator,bool>(NULL,false);
  }
  Map::element *Map::add_element(const std::pair<Point,long>& p){
    element* c;
    c=(element*)malloc(sizeof(element));
    c->d.first=p.first;
    c->d.second=p.second;
    c->left=NULL;
    c->right=NULL;
    ++size_;
    return c;
  }
}