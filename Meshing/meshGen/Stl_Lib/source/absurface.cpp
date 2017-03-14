#include "stl_io.h"
#include <stdio.h>
namespace stl
{
  using namespace std;
  pair<long,double> solve(double*,ulong);
  Stl_io Stl_io::absurface(const char* fileName)
  {
    vector<double>phi;
    double tmpDouble;
    math_our::Point tmpPoint;
    phi.reserve(pointArray.size());
    FILE *file;    
    file=fopen(fileName,"rb");
    for(unsigned long i=0;file && i<pointArray.size();i++){
      if(!fscanf(file,"%20lf",&tmpDouble)){
	printf("Aborted!\nnumber of lines is not equal to namber of points\n");
	fclose(file);
	return *this;
      }
      phi.push_back(tmpDouble);
    }
    fclose(file);
    typedef std::vector<Group> GroupArray;
    typedef std::vector<Point> PointArray;
    GroupArray ablationGroupArray;
    PointArray ablationPointArray;
    
    ablationGroupArray=groupArray;
    ablationPointArray.reserve(nPoints);
    vector<double>line;
    line.resize(nLayers+1);
    pair<long,double>solution;
    double currl;
    file=fopen("meshGen/labla","wb");
    long nwrongs=0;
    for(long  i=0;i<nPoints;i++){
      currl=0;
      for(long j=0;j<=nLayers;j++){
	line[j]=phi[j*nPoints+i];
      }
      if(line[0]<0 || line[nLayers]>0){
        nwrongs++;
      }
      solution=solve(&(line[0]),line.size());
      tmpPoint=pointArray[solution.first*nPoints+i];
      if(solution.first){
        for(long j=1;j<=solution.first;j++){
          currl+=(pointArray[j*nPoints+i]-pointArray[(j-1)*nPoints+i]).module();
        }
      }
      if(solution.first!=nLayers){
	tmpPoint+=(pointArray[(solution.first+1)*nPoints+i]-tmpPoint)*solution.second;
        currl+=(tmpPoint-pointArray[solution.first*nPoints+i]).module();
      }
      ablationPointArray.push_back(tmpPoint);
      fprintf(file,"%14lf\n",currl);
    }
    fclose(file);
    printf("nWrong %ld\n",nwrongs);
    Stl_io ablationStl(ablationPointArray,ablationGroupArray);
    return ablationStl;
  }
  pair<long,double> solve(double*a,ulong n){
    if(a[0]<0)return pair<long,double>((long)0,(double)0.);
    else if(a[n-1]>0) return pair<long,double>((long)(n-1),(double)0.);
    ulong i;
    for(i=1;a[i]>0;++i){
    }
    return pair<long,double>((long)(i-1),(double)(0.-a[i-1])/(a[i]-a[i-1]));
  }
} 
