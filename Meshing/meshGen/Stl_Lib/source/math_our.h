#pragma once

#include <iostream>

#define ulong unsigned long

namespace math_our{
  extern double TOL;
  bool cmpDouble(double a,double b);
  class Point{
  public:
    Point():x(0),y(0),z(0){}
    Point(double,double,double);
    double x;
    double y;
    double z;
    bool operator<(const Point& p) const;
    bool operator<(const double p) const;
    bool operator==(const Point& p) const;
    Point operator+(const Point& p) const;
    Point operator^(const Point& p) const;
    const Point& operator+=(const Point& p);
    const Point& operator-=(const Point& p);
    const Point& operator/=(double p);
    const Point& operator*=(double p);
    Point operator-(const Point& p) const;
    double operator*(const Point& p) const;
    Point operator*(double p) const;
    Point operator/(double p) const;
    double module() const;
    double value() const;
    const Point& normalize();
    void set(double,double,double);
  };
  class Map{
  private:
    class element{
    public:
      struct data{
	Point first;
	long second;
      };
      data d;
      element* left;
      element* right;
    };
  public:
    typedef element::data *iterator;
    std::pair<iterator,bool>insert(const std::pair<Point,long>&);
  private:
      ulong size_;
      element* root; 
      element *add_element(const std::pair<Point,long>&);
  public:
    ulong size()const;
    Map():size_((ulong)0),root(NULL){};
    ~Map();
  };
  double mix_mul(const Point&,const Point&,const Point&);
  std::ostream& operator << (std::ostream& os, const Point& p);
}