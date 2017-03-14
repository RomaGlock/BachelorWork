#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>


namespace stl
{
  using namespace std;
  bool Triangle::operator<(const Triangle&  p)const{
    long a1=a;
    long a2=p.a;
    long b1=b;
    long b2=p.b;
    long c1=c;
    long c2=p.c;
    long buff;
    if(a==b || a==c || b==c)return false;
    //---------------this
    if(a1>b1){
      buff=a1;
      a1=b1;
      b1=buff;
    }
    if(b1>c1){
      buff=b1;
      b1=c1;
      c1=buff;
    }
    if(a1>b1){
      buff=a1;
      a1=b1;
      b1=buff;
    }
    if(b1>c1){
      buff=b1;
      b1=c1;
      c1=buff;
    }
    //---------------that
    if(a2>b2){
      buff=a2;
      a2=b2;
      b2=buff;
    }
    if(b2>c2){
      buff=b2;
      b2=c2;
      c2=buff;
    }
    if(a2>b2){
      buff=a2;
      a2=b2;
      b2=buff;
    }
    if(b2>c2){
      buff=b2;
      b2=c2;
      c2=buff;
    }
    //-----------------------
    if(a1<a2)return true;
    else if(a1>a2) return false;
    else if(a1==a2 && b1<b2) return true;
    else if(b1==b2 && c1<c2) return true;
    else return false;
  }
  Triangle::Triangle(long a,long b,long c, const math_our::Point& n,long gr):a(a),b(b),c(c),n(n),gr(gr){  
    
  }
  
} 
