#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>
namespace stl
{
  using namespace std;

  void CudaSetDevice(int dev)
  {
    cudaSetDevice(dev);
  }

  __global__ void makeAllNewThreads_kernel(long start,double pow_E,long nPoints, double *dev_pointArray, double *dev_newPoint, unsigned long phantom_pointArray_size, double *dev_charges);

  __global__ void makeAllNewThreads_kernel(long start,double pow_E,long nPoints, double *dev_pointArray, double *dev_newPoint, unsigned long phantom_pointArray_size, double *dev_charges)
  {
    long i = start+threadIdx.x + blockIdx.x*blockDim.x;
    if (i < nPoints)
    {
      //local var
      double x=dev_pointArray[SSD_1*i+0];
      double y=dev_pointArray[SSD_1*i+1];
      double z=dev_pointArray[SSD_1*i+2];
      double l=dev_pointArray[SSD_1*i+3];
      double nlT0;
      double nlT1;
      double nlT2;
      double Ex,Ey,Ez;
      double E2;
      double X,Y,Z,Q;
      double Rx,Ry,Rz;
      double r2,r;
      unsigned long k;
      Ex=0;
      Ey=0;
      Ez=0;
      for(k=0;k<phantom_pointArray_size;k++){
        X=dev_charges[k*SSD_3+0];
        Y=dev_charges[k*SSD_3+1];
        Z=dev_charges[k*SSD_3+2];
        Q=dev_charges[k*SSD_3+3];
        Rx=x-X;
        Ry=y-Y;
        Rz=z-Z;
        r2=Rx*Rx+Ry*Ry+Rz*Rz;
        r=pow(r2,pow_E);
        Ex+=(Q*Rx/r);
        Ey+=(Q*Ry/r);
        Ez+=(Q*Rz/r);	      
      }
      E2=Ex*Ex+Ey*Ey+Ez*Ez;
      E2=sqrt(E2);
      Ex=Ex/E2*l;
      Ey=Ey/E2*l;
      Ez=Ez/E2*l;
      nlT0=x+Ex;
      nlT1=y+Ey;
      nlT2=z+Ez;
      dev_newPoint[SSD_1*i+0] = nlT0;
      dev_newPoint[SSD_1*i+1] = nlT1;
      dev_newPoint[SSD_1*i+2] = nlT2;
    }
  }  
  void Stl_io::MMPoi(const char* fileName,Stl_io& phantom){
    
    MMLE_arg arg(fileName);
    long FirststLayerOri;
    FirststLayerOri=(arg.nCudaIter>0)?1:-1;
    arg.nCudaIter=abs(arg.nCudaIter);
    //--------------------------internal stl information-------------------------
    nPoints=pointArray.size()/(nLayers+1);
    //nLayers
    //--------------------------Tmp Variables-------------------------------------
    math_our::Point e1,e2,e3;
    double tmpDouble;
    math_our::Point tmpPoint;
    long offset=nLayers*nPoints;
    //----------------------------------------------------------------------------
    typedef std::vector<Group> GroupArray;
    typedef std::list<Triangle> TriangleArray;
    //--------------------------normal list-----------------------------------
    math_our::Point *nList;
    long *nNList;
    get_normals_av(&nList,&nNList);
    //--------------------------length list-----------------------------------
    double *lList;
    lList=get_lengths(fileName);
//     return;
    //--------------------------profile-----------------------------------------
    vector<double> profile;
    get_profile(profile,fileName);
    //-------------------------charges----------------------------------------
    double *charges;
    ulong nCharges;
    phantom_charge charge_pair;
    charge_pair=phantom.get_charges(arg.pow_q,FirststLayerOri,fileName);
    charges=charge_pair.charge;
    nCharges=charge_pair.nCharge;
    {
      long negChC=0;
      for(unsigned long i=0;i<nCharges;i++){
        if(charges[SSD_3*i+3]>0){
        }else{
          negChC++;
        }
      }
      if(negChC)printf("negative charge  %ld |%lu %lu\n",negChC,nCharges,(ulong)(phantom.pointArray.size()));
    }
    //--------------------------statistic----------------------------------------
    unsigned int start_time =  clock();
    //--------------------------meshing------------------------------------------
    /* unique data:
    * nList[] 	- array of normal
    * lList[] 	- array of length
    * arg_n 	- number of point in each tread
    * arg.m	- mesh maker parameter
    * phantom.pointArray.size() number of cHarged points
    * pointArray[] - array of coordinates
    * shared data
    * charges[] array of x,y,z,q
    */
    vector<math_our::Point> layerBuffer;
    double *newPoint1D;
    newPoint1D=(double*)malloc(sizeof(double)*SSD_1*nPoints);
    double *pointsCuda;
    pointsCuda=(double*)malloc(sizeof(double)*SSD_2*nPoints);
    arg.nCudaIter=profile.size();
    for(long ii=0;ii<arg.nCudaIter;ii++){
      start_time =  clock();
      long arg_n=1;
//     //--------------------------internal stl information-------------------------
      nPoints=pointArray.size()/(nLayers+1);
      pointArray.resize((nLayers+arg_n+1)*nPoints);
      offset=nLayers*nPoints;
    //--------------------------first layer--------------------------------------
      bool flag=false;
      if(!is_3Dmesh){
        if(arg.firstLayer){
          for(long i=0;i<nPoints;i++) {
            tmpDouble=nList[i].module();
            tmpPoint=nList[i]*lList[i]*profile[ii]/tmpDouble*((double)FirststLayerOri);
            pointArray[offset+nPoints+i]=pointArray[offset+i]+tmpPoint;
          }
          nLayers++;
          offset+=nPoints;
        }else{
          flag=true;
          layerBuffer.reserve(nPoints);
          for(long i=0;i<nPoints;i++) {
            tmpDouble=nList[i].module();
            tmpPoint=nList[i]*lList[i]*profile[ii]/tmpDouble*0.05*((double)FirststLayerOri);
            layerBuffer.push_back(pointArray[offset+i]);
//             pointArray[offset+i]+=tmpPoint;
          }
          arg_n++;
        }
      }
      else{
        arg_n++;
      }
      //-------------------------other layers--------------------------------------
      
      if(arg_n>1){
        for(long i=0;i<nPoints;i++){
          pointsCuda[SSD_2*i+0]=pointArray[offset+i].x;
          pointsCuda[SSD_2*i+1]=pointArray[offset+i].y;
          pointsCuda[SSD_2*i+2]=pointArray[offset+i].z;
          pointsCuda[SSD_2*i+3]=lList[i]*profile[ii];
        }
      //--------------------------cuda place---------------------------------------
        
        double *dev_charges, *dev_newPoint;
        double *dev_pointArray;

        cudaMalloc( (void**)&dev_charges, SSD_3*nCharges*sizeof(double) );
        cudaMalloc( (void**)&dev_pointArray, nPoints*SSD_2*sizeof(double) );
        cudaMalloc( (void**)&dev_newPoint, SSD_1*nPoints*sizeof(double) );

        cudaGetErrorString (cudaMemcpy(dev_charges, charges, SSD_3*nCharges*sizeof(double), cudaMemcpyHostToDevice));
        cudaGetErrorString (cudaMemcpy(dev_pointArray, pointsCuda, nPoints*SSD_2*sizeof(double), cudaMemcpyHostToDevice));

        for(long internal_i=0;internal_i<nPoints/(long)2048+(long)1;++internal_i){
          long threads = 256;
          long blocks = 8;//(nPoints)/threads+1;
          long start = internal_i*2048;
          makeAllNewThreads_kernel<<<blocks, threads>>>(
            start,
            arg.pow_E,nPoints,dev_pointArray,
            dev_newPoint,
            nCharges,dev_charges);
        }
        cudaGetErrorString (cudaMemcpy(newPoint1D, dev_newPoint, SSD_1*nPoints*sizeof(double), cudaMemcpyDeviceToHost));
        cudaFree(dev_charges);
        cudaFree(dev_pointArray);
        cudaFree(dev_newPoint);
        for(long i=0;i<nPoints;i++){
          for(long j=1;j<arg_n;j++){
            pointArray[offset+j*nPoints+i].x=newPoint1D[SSD_1*i+0];
            pointArray[offset+j*nPoints+i].y=newPoint1D[SSD_1*i+1];
            pointArray[offset+j*nPoints+i].z=newPoint1D[SSD_1*i+2];
          }
        }
      }     
      nLayers+=arg_n;
      nLayers--;
      if(flag){
        for(long i=0;i<nPoints;i++) {
          pointArray[offset+i]=layerBuffer[i];
        }
        layerBuffer.clear();
      }
      //---------------------------relaxation--------------------------------------
      llrelax_cu_fast(arg.nStep,arg.elas,arg.elas2);
      //--------------------------after meshing actions----------------------------
      unsigned int end_time = clock();
      unsigned int search_time = end_time - start_time;
      fprintf(stdout,"Time (s) %ld,%06ld, Time per point %ld.%06ld (s). %4ld Layers\n",
              search_time/1000000,search_time%1000000,search_time/nPoints/1000000,search_time/nPoints%1000000,nLayers);
      fflush(stdout);
      is_3Dmesh=(true);
    }
    free(nList);
    free(nNList);
    free(lList);
    free(charges);
    free(newPoint1D);
    free(pointsCuda);
  }
}
