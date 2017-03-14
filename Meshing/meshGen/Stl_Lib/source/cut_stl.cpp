#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>


namespace stl
{
  using namespace std;
  
  
  math_our::Point projection(const math_our::Point& n,const math_our::Point& n1,
                            double D,const math_our::Point& r0);
  math_our::Point findNormal(long a,long b,Triangle* tri,unsigned long size);
  Stl_io Stl_io::Cut_stl(const char* fileName){
    CUT_arg cut_arg(fileName);
    math_our::Point n;
    double D;
    n=cut_arg.n/cut_arg.n.module();
    D=cut_arg.D;
    Group::TriangleArray::iterator jd;
    //--------------------------normal list-----------------------------------
    math_our::Point *nList;
    nList=(math_our::Point*)malloc(sizeof(math_our::Point)*nPoints);
//     long nNList[nPoints];
    for(long i=0;i<nPoints;i++) {
      nList[i]=math_our::Point(0.,0.,0.);
//       nNList[i]=0;
    }
    vector<Triangle> onEdgeTri;
    for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::iterator j = i->tri.begin(); j != i->tri.end(); j=jd) {
        jd=j;
        ++jd;
        if(pointArray[j->a]*n+D>0 || pointArray[j->b]*n+D>0 || pointArray[j->c]*n+D>0){
          if(pointArray[j->a]*n+D>0 && pointArray[j->b]*n+D>0 && pointArray[j->c]*n+D>0){
          }else{
            onEdgeTri.push_back(*j);
          }
          i->tri.erase(j);
        }
      }
    }
    set<edge> edgeTree;
    long k=0;
    long buff;
    long nNeighbours=0;
    long nOwners=0;
    for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        addTriInEdgeTree(k,*j,edgeTree);
        k++;
      }
    }
    
    //------------------------------make owner<neighbour--------------------------
    for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
      if(i->neig!=-1){
        if(i->neig<i->owner){
          buff=i->neig;
          i->neig=i->owner;
          i->owner=buff;
          buff=i->neig3rd;
          i->neig3rd=i->owner3rd;
          i->owner3rd=buff;	  
        }
        nNeighbours++;
      }
      nOwners++;
    }    
    vector<long> trans;
    trans.resize(pointArray.size());
    long buffer;
    for(unsigned long i=0;i<pointArray.size();i++)trans[i]=-1;
    pointArray.resize((nOwners-nNeighbours)*2+pointArray.size());    
    double diag1,diag2;
    math_our::Point newPoint;
    Triangle newTriangle1;
    Triangle newTriangle2;
    math_our::Point p1,p2,p3;
    D-=cut_arg.cs;
    for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
      if(i->neig==-1){
        newPoint=findNormal(i->a,i->b,&(onEdgeTri[0]),onEdgeTri.size());
        if(newPoint.module()>0.5){
          nList[i->a]=newPoint;
          nList[i->b]=newPoint;
        }
      }
    }
    for(long i=0;i<nPoints;i++) {
      nList[i].normalize();
    }
    for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
      if(i->neig==-1){
        if(trans[i->a]==-1){
          newPoint=projection(n,nList[i->a],D,pointArray[i->a]);
          p1=find_image(newPoint,n,D);
          trans[i->a]=pointArray.size();
          newPoint=(newPoint+p1)/2.0;
          pointArray.push_back(newPoint);
        }
        if(trans[i->b]==-1){
          newPoint=projection(n,nList[i->b],D,pointArray[i->b]);
          p1=find_image(newPoint,n,D);
          trans[i->b]=pointArray.size();
          newPoint=(newPoint+p1)/2.0;
          pointArray.push_back(newPoint);
        }
        diag1=(pointArray[i->a]-pointArray[trans[i->b]]).module();
        diag2=(pointArray[i->b]-pointArray[trans[i->a]]).module();
        if (diag1<diag2){
          newTriangle1.a=i->a;
          newTriangle1.b=trans[i->b];
          newTriangle1.c=trans[i->a];
          newTriangle1.n=i->n;
          newTriangle1.gr=i->group;
          newTriangle2.a=i->a;
          newTriangle2.b=trans[i->b];
          newTriangle2.c=i->b;
          newTriangle2.n=i->n;
          newTriangle2.gr=i->group;
        }else{
          newTriangle1.a=i->a;
          newTriangle1.b=i->b;
          newTriangle1.c=trans[i->a];
          newTriangle1.n=i->n;
          newTriangle1.gr=i->group;
          newTriangle2.a=i->b;
          newTriangle2.b=trans[i->a];
          newTriangle2.c=trans[i->b];
          newTriangle2.n=i->n;
          newTriangle2.gr=i->group;
        }
        p1=(pointArray[newTriangle1.b]-pointArray[newTriangle1.a]);
        p2=(pointArray[newTriangle1.c]-pointArray[newTriangle1.a]);
        p3=(p1^p2).normalize();
        if(p3*newTriangle1.n<0){
          newTriangle1.n=p3*(-1.);
          buffer=newTriangle1.a;
          newTriangle1.a=newTriangle1.b;
          newTriangle1.b=buffer; 
        }else{
          newTriangle1.n=p3*(1.);	  
        }
        p1=(pointArray[newTriangle2.b]-pointArray[newTriangle2.a]);
        p2=(pointArray[newTriangle2.c]-pointArray[newTriangle2.a]);
        p3=(p1^p2).normalize();
        if(p3*newTriangle2.n<0){
          newTriangle2.n=p3*(-1.);
          buffer=newTriangle2.a;
          newTriangle2.a=newTriangle2.b;
          newTriangle2.b=buffer;
        }else{
          newTriangle2.n=p3*(1.);	  
        }
        groupArray[i->group].tri.push_back(newTriangle1);
        groupArray[i->group].tri.push_back(newTriangle2);
      }
    }
    //------------deleting unused points--------------------
    
    std::vector<long>newName;
    newName.resize(pointArray.size());
    for(unsigned long i=0;i<pointArray.size();i++) newName[i]=-1;
    vector<math_our::Point> unique;
    unique.reserve(pointArray.size());
    for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        if(newName[j->a]==-1){
          newName[j->a]=unique.size();
          unique.push_back(pointArray[j->a]);
        }
        if(newName[j->b]==-1){
          newName[j->b]=unique.size();
          unique.push_back(pointArray[j->b]);
        }
        if(newName[j->c]==-1){
          newName[j->c]=unique.size();
          unique.push_back(pointArray[j->c]);
        }
        j->a=newName[j->a];
        j->b=newName[j->b];
        j->c=newName[j->c];
      }
    }
    for(vector<long>::iterator i=trans.begin();i!=trans.end();++i) 
      if((*i)!=-1) (*i)=newName[*i];
    vector<long> redPoint;
    redPoint.resize(unique.size());
    for(unsigned long i=0;i<unique.size();++i) {
      redPoint[i]=-2;
    }
    long nRedPoints=0;
    for(vector<long>::iterator i=trans.begin();i!=trans.end();++i) {
      if((*i)!=-1) {
        redPoint[*i]=-1;
        nRedPoints++;
      }
    }
    //redPoint: -1 = point on cut edge;
    //		-2 = common point;
    
    swap(unique,pointArray);
    //------------phantom points and elements--------------------
    vector<math_our::Point> phantomPoints;
    typedef std::list<Triangle> TriangleArray;
    TriangleArray phantomElements;
    phantomPoints.resize(2*pointArray.size()-nRedPoints);
    long pNewPoint=pointArray.size();
    for(unsigned long i=0;i<pointArray.size();i++){
      phantomPoints[i]=pointArray[i];
    }
    for(unsigned long i=0;i<pointArray.size();i++){
      if(redPoint[i]==-2){
        newPoint=find_image(pointArray[i],n,D);
        phantomPoints[pNewPoint]=newPoint;
        redPoint[i]=pNewPoint;
        pNewPoint++;
      }
      if(redPoint[i]==-1) redPoint[i]=i;
    }
    math_our::Point phantomZero=find_image(math_our::Point(0,0,0),n,D);
    for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        newTriangle1.a=j->a;
        newTriangle1.b=j->b;
        newTriangle1.c=j->c;
        newTriangle1.n=j->n;
        newTriangle1.gr=0;
        phantomElements.push_back(newTriangle1);
      }
    }
    for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        newTriangle1.a=redPoint[j->a];
        newTriangle1.b=redPoint[j->c];
        newTriangle1.c=redPoint[j->b];
        newTriangle1.n=find_image(j->n,n,D)-phantomZero;
        newTriangle1.gr=0;
        phantomElements.push_back(newTriangle1);
      }
    }
    is_Cutted=true;
    typedef std::vector<Group> GroupArray;
    GroupArray phantomGroupArray;
    phantomGroupArray.resize(1);
    phantomGroupArray[0].tri.assign(phantomElements.begin(),phantomElements.end());
    phantomGroupArray[0].name="phantom";
    Stl_io phantomStl(phantomPoints,phantomGroupArray);
    nPoints=(long)pointArray.size();
    return phantomStl;
  }
  math_our::Point findNormal(long a,long b,Triangle* tri,unsigned long size){
    math_our::Point normal;
    for(unsigned long i=0;i<size ;++i){
      if( (a==tri[i].a && b==tri[i].b) ||
          (a==tri[i].b && b==tri[i].a) ||
          (a==tri[i].b && b==tri[i].c) ||
          (a==tri[i].c && b==tri[i].b) ||
          (a==tri[i].a && b==tri[i].c) ||
          (a==tri[i].c && b==tri[i].a)
      ){
        normal=tri[i].n;
        normal.normalize();
        return normal;
      }
    }
    printf("%s(%d):edge not found in cut triangles\n",__FILE__,__LINE__);
    return math_our::Point(0.,0.,0.);
  }
  math_our::Point projection(const math_our::Point& n,const math_our::Point& n1,
                            double D,const math_our::Point& r0){
    math_our::Point newPoint,a,b;
    double tau,alpha,beta,D1;
    D1=-(r0*n1);
    a=n1^n;
    beta=(D-D1*(n1*n))/(pow(n1*n,2.)-1);
    alpha=-beta*(n1*n)-D1;
    b=n1*alpha+n*beta;
    
    tau=-D-(r0*n);
    newPoint=r0+n*tau;
    
    tau=((a*newPoint)-(a*b))/(a*a);
    newPoint=a*tau+b;
    return newPoint;
  }
  math_our::Point find_image(math_our::Point p,math_our::Point n,double D){
    double tau;
    math_our::Point newPoint;
    tau=-D-(p*n);
    newPoint=p+n*(2.*tau);
    return newPoint;
  }
}
