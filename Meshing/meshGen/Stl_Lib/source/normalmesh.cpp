#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>
namespace stl
{
  using namespace std;
  void Stl_io::nmesh(const char* fileName){
    MMLE_arg arg(fileName);
    long FirststLayerOri;
    FirststLayerOri=(arg.nCudaIter>0)?1:-1;
    arg.nCudaIter=abs(arg.nCudaIter);
    //--------------------------internal stl information-------------------------
    nPoints=pointArray.size()/(nLayers+1);
    //--------------------------Tmp Variables-------------------------------------
    double tmpDouble;
    math_our::Point tmpPoint;
    long offset;
    //--------------------------normal list-----------------------------------
    math_our::Point *nList;
    long *nNList;    
    //--------------------------length list-----------------------------------
    double *lList;
    lList=get_lengths(fileName);
    //--------------------------profile-----------------------------------------
    vector<double> profile;
    get_profile(profile,fileName);
    //--------------------------meshing------------------------------------------
    arg.nCudaIter=profile.size();
    printf("%ld layers will be created\n",arg.nCudaIter);
    for(long ii=0;ii<arg.nCudaIter;ii++){
      get_normals_av(&nList,&nNList);
    //--------------------------internal stl information-------------------------
      nPoints=pointArray.size()/(nLayers+1);
      pointArray.resize((nLayers+2)*nPoints);
      offset=nLayers*nPoints;
    //--------------------------first layer--------------------------------------
      for(long i=0;i<nPoints;i++){
        tmpDouble=nList[i].module();
        tmpPoint=nList[i]*lList[i]*profile[ii]/tmpDouble*((double)FirststLayerOri);
        pointArray[offset+nPoints+i]=pointArray[offset+i]+tmpPoint;
      }
      nLayers++;
   //---------------------------relaxation--------------------------------------
      llrelax_cu_fast(arg.nStep,arg.elas,arg.elas2);
    //--------------------------after meshing actions----------------------------
      is_3Dmesh=(true);
      free(nList);
      free(nNList);
    }
    free(lList);
  }
}
