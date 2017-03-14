#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>
#include "math_our_cu.h"
#define DXRELAX 0.001
namespace stl
{
  using namespace std;
  __global__ void makelrelax_fast(
          long start,
          long nPoints,
          long* dev_nnei,
          math_our_cu::Point_cu* dev_neigth_points,
          math_our_cu::Point_cu* dev_next_layer,
          double alpha,
          long mem_struct_size,
          long* dev_exept
        ){
    long i = start+threadIdx.x + blockIdx.x*blockDim.x;
    if (i >= nPoints) return;
    long n=dev_nnei[i];
    long index;
    math_our_cu::Point_cu e,downP,bufP,x0,e0;
    x0=dev_neigth_points[mem_struct_size*i+2*n];
    downP=dev_neigth_points[mem_struct_size*i+2*n+1];
    bufP=dev_neigth_points[mem_struct_size*i+2*n+2];
    e=bufP-downP;
    e0=e;
    double tmpDouble;
    double N=double(n);
    double cs=e.module();
    e0/=cs;
    math_our_cu::Point_cu f(0.,0.,0.);
    double f_,f_elas,dfdt,dfdt_elas,f_elas_;
    math_our_cu::Point_cu tmpPoint;
    math_our_cu::Point_cu tmpPoint_;
    double tmpPointModule,tmpPoint_e;
    double tmpPointModule_,tmpPoint_e_;
    for(long j=0;j<n;++j){
      index=mem_struct_size*i+j;
      f+=dev_neigth_points[index];
    }
    for(long iter=0;iter<10;++iter){
      e=(x0-downP).normalize();
//       tmpPoint=bufP-x0;
      tmpPoint=downP+e*cs-x0;
      tmpPointModule=tmpPoint.module();
      tmpPoint_e=tmpPoint*e;
//       tmpPoint_=bufP-x0-e*DXRELAX*cs;      //tmpPoint_=bufP-x0-e*DXRELAX;
      tmpPoint_=downP+e*cs-x0-e*DXRELAX*cs;
      tmpPointModule_=tmpPoint_.module();
      tmpPoint_e_=tmpPoint_*e;
      f_=(f-x0*N)*e;
      f_elas=tmpPoint_e*tmpPointModule;
      f_elas_=tmpPoint_e_*tmpPointModule_;
      dfdt=-N;
      dfdt_elas=(f_elas_-f_elas)/DXRELAX/cs; //dfdt_elas=(f_elas_-f_elas)/DXRELAX
      f_elas/=cs;
      dfdt_elas/=cs;
      f_elas*=alpha;
      dfdt_elas*=alpha;
      tmpDouble=(f_+f_elas);
      tmpDouble/=(dfdt+dfdt_elas);
      if((((x0-e*tmpDouble)-downP)*e)>0.5*cs){
        x0-=e*tmpDouble;
      }else{
        x0=downP+e*cs*0.5;
      }
    }    
    dev_next_layer[i]=x0;
  }
  __global__ void refreshnei(
          long start,
          long nPoints,
          long* dev_nei,long* dev_nnei,long* dev_offset_nei,
          math_our_cu::Point_cu* dev_neigth_points,
          math_our_cu::Point_cu* dev_next_layer,
          long mem_struct_size
        ){
    long i = start+threadIdx.x + blockIdx.x*blockDim.x;
    if (i >= nPoints) return;
    long n=dev_nnei[i];
    long offset=dev_offset_nei[i];
    for(ulong j=0;j<n;j++){
      dev_neigth_points[mem_struct_size*i+j]=dev_next_layer[dev_nei[offset+j]];
    }
    dev_neigth_points[mem_struct_size*i+2*n]=dev_next_layer[i];
  }
  __device__ math_our_cu::Point_cu find_force1(math_our_cu::Point_cu* nei,
                                              long n,
                                              math_our_cu::Point_cu x0,
                                              math_our_cu::Point_cu downP,
                                              math_our_cu::Point_cu bufP,
                                              math_our_cu::Point_cu normal0,
                                              double cs,
                                              double alpha,long i
                                             ){
    math_our_cu::Point_cu normal,currEdge,f(0.,0.,0.);
    double tmpDouble;
    normal=(x0-downP).normalize();
    for(ulong j=0;j<n;j++){
      currEdge=x0-nei[j];
      tmpDouble=currEdge.x*currEdge.x+currEdge.y*currEdge.y+currEdge.z*currEdge.z;
      f+=currEdge/tmpDouble*cs;
    }
    tmpDouble=(bufP-x0).module();
    if(tmpDouble>cs*0.01){
      f+=(bufP-x0)/tmpDouble*sqrt(fabs(1.-normal0*normal))*alpha;
      f+=(bufP-x0)/cs*alpha;
    }
    return f;
  }
  __device__ math_our_cu::Point_cu find_force(math_our_cu::Point_cu* nei,
                                              long nn,
                                              math_our_cu::Point_cu x0,
                                              math_our_cu::Point_cu downP,
                                              math_our_cu::Point_cu bufP,
                                              math_our_cu::Point_cu normal0,
                                              double cs,
                                              double alpha,long i, long out
                                             ){
    math_our_cu::Point_cu normal,n,currEdge,f(0.,0.,0.),sf(0.,0.,0.),uf(0.,0.,0.),a,b;
    double tmpDouble,D,l;
    normal=x0-downP;
    l=normal.module();
    normal/=l;
    for(ulong j=0;j<nn;j++){
      currEdge=nei[j]-x0;
//       tmpDouble=(nei[j-n]-downP).module();
      f+=currEdge;// /tmpDouble;
    }
    for(ulong j=0;j<(nn-1);j++){
      a=(nei[j+nn]-nei[j]);
      b=(nei[j+1+nn]-nei[j]);
      n=a^b;      
      tmpDouble=n.module();
      n/=tmpDouble;
      D=-(nei[j+nn]*n);
      if((downP-nei[j+nn])*n<0){
        n*=(-1.);
        D=-D;
      }
      currEdge=nei[j]-nei[j+1];
      sf+=n*(currEdge.module()*0.86-((x0*n)+D));
      //force side//
      a=(nei[j]-nei[j+1+nn]);
      b=(nei[j+1]-nei[j+1+nn]);
      n=a^b;
      tmpDouble=n.module();
      n/=tmpDouble;
      D=-(nei[j+1+nn]*n);
      if((downP-nei[j+1+nn])*n<0){
        n*=(-1.);
        D=-D;
      }
      sf+=n*(currEdge.module()*0.86-((x0*n)+D));
      //force side//
      a=(nei[j+nn]-downP);
      b=(nei[j+1+nn]-downP);
      n=a^b;
      tmpDouble=n.module();
      n/=tmpDouble;
      D=-(downP*n);
      if((bufP-downP)*n<0){
        n*=(-1.);
        D=-D;
      }
      uf+=n*(l-((x0*n)+D));
      //force up//
    }
    if(i){
      a=(nei[nn-1+nn]-nei[nn-1]);
      b=(nei[0+nn]-nei[nn-1]);
      n=a^b;
      tmpDouble=n.module();
      n/=tmpDouble;
      D=-(nei[nn-1+nn]*n);
      if((downP-nei[nn-1+nn])*n<0){
        n*=(-1.);
        D=-D;
      }
      currEdge=nei[nn-1]-nei[0];
      sf+=n*(currEdge.module()*0.86-((x0*n)+D));
      //force side//   
      a=(nei[nn-1]-nei[0+nn]);
      b=(nei[0]-nei[0+nn]);
      n=a^b;
      tmpDouble=n.module();
      n/=tmpDouble;
      D=-(nei[0+nn]*n);
      if((downP-nei[0+nn])*n<0){
        n*=(-1.);
        D=-D;
      }
      sf+=n*(currEdge.module()*0.86-((x0*n)+D));
      //force side//
      a=(nei[nn-1+nn]-downP);
      b=(nei[0+nn]-downP);
      n=a^b;
      tmpDouble=n.module();
      n/=tmpDouble;
      D=-(downP*n);
      if((bufP-downP)*n<0){
        n*=(-1.);
        D=-D;
      }
      uf+=n*(l-((x0*n)+D));      
      //force up//
    }
//     tmpDouble=(bufP-x0).module();
//     if(tmpDouble>cs*0.01){
//       f+=(bufP-x0)/tmpDouble*sqrt(fabs(1.-normal0*normal))*alpha*cs;
//       f+=(bufP-x0)*alpha;
//     }
//     return f+sf*alpha+uf*alpha;
    return f+uf*alpha+sf*alpha;
  }
  __device__ math_our_cu::Point_cu spring_relax_solve(math_our_cu::Point_cu x0,
                                                      math_our_cu::Point_cu dx,
                                                      math_our_cu::Point_cu downP,
                                                      double D,
                                                      math_our_cu::Point_cu n,
                                                      double min_cos,long i){
    double f1,f2,fmid;
    math_our_cu::Point_cu x1,x2,xmid;
    x1=x0;
    x2=x1+dx;
    f2=(x2-downP).normalize()*n-min_cos;
    if(f2>0)return dx;
    f1=(x1-downP).normalize()*n-min_cos;
    if(f1<0)return math_our_cu::Point_cu(0.,0.,0.); //no idea how to solve
    for(long iter=0;iter<10;++iter){
      xmid=(x1+x2)*0.5;
      fmid=(xmid-downP).normalize()*n-min_cos;
      if(fmid>0){
        x1=xmid;
      }else{
        x2=xmid;
      }
    }
    return (x1+x2)*0.5-x0;
  }
  __device__ math_our_cu::Point_cu spring_relax_check(math_our_cu::Point_cu x0,
                                                math_our_cu::Point_cu dx,
                                                math_our_cu::Point_cu bufP,
                                                math_our_cu::Point_cu downP,
                                                math_our_cu::Point_cu* nei,
                                                long nnei,
                                                double min_cos,long i
                                             ){
    double D;
    math_our_cu::Point_cu n;
    for(ulong j=0;j<(nnei-1);j++){
      n=(nei[j]-downP)^(nei[j+1]-downP);
      n.normalize();
      D=-(downP*n);
      if((bufP-downP)*n<0){
        n*=(-1.);
        D=-D;
      }
      
      dx=spring_relax_solve(x0,dx,downP,D,n,min_cos,i);
    }
    if(i){
      n=(nei[nnei-1]-downP)^(nei[0]-downP);
      n.normalize();
      D=-(downP*n);
      if((bufP-downP)*n<0){
        n*=(-1.);
        D=-D;
      }
      dx=spring_relax_solve(x0,dx,downP,D,n,min_cos,i);
    }
    return dx;
  }
  __global__ void spring_relax(
          long start,
          long nPoints,
          long* dev_nnei,
          math_our_cu::Point_cu* dev_neigth_points,
          math_our_cu::Point_cu* dev_next_layer,
          double alpha,
          long mem_struct_size,
          long* dev_exept
        ){
    long i = start+threadIdx.x + blockIdx.x*blockDim.x;
    if (i >= nPoints) return;
    long out=0;
    long n=dev_nnei[i];
    long offset_curr=mem_struct_size*i;
    long offset_prev=offset_curr+n;
    math_our_cu::Point_cu x0=dev_neigth_points[offset_prev+n];
    if(!dev_exept[i]){
//       dev_next_layer[i]=x0;
      return;
    }
    math_our_cu::Point_cu downP=dev_neigth_points[offset_prev+n+1];
    math_our_cu::Point_cu bufP=dev_neigth_points[offset_prev+n+2];
    math_our_cu::Point_cu f,f_dx,f_dy;
    math_our_cu::Point_cu force_dir_x,force_dir_y;
    math_our_cu::Point_cu x_dx,x_dy;
    double cs,f_n,tmpDouble;
    double J11,J12,J21,J22;
    double J_11,J_21;
    double dx,dy;
    math_our_cu::Point_cu normal0=bufP-downP;
    cs=normal0.module();
    normal0/=cs;
    math_our_cu::Point_cu normal;
    for(long iter=0;iter<1;++iter){
      normal=(x0-downP).normalize();
      f=find_force(dev_neigth_points+offset_curr,n,x0,downP,bufP,normal0,cs,alpha,1,out);
      if( (isnan(f.x)) && (isnan(f.y)) && (isnan(f.z))){
        printf("f   %ld\n",i);
      }
      f=f-(normal*(f*normal));
      f_n=f.module();     //force module
      if(f_n<1e-7){
        dev_next_layer[i]=x0; 
        return;
      }
      force_dir_x=f/f_n;  //force direction
      force_dir_y=normal^force_dir_x;
      x_dx=x0+force_dir_x*DXRELAX*cs;
      x_dy=x0+force_dir_y*DXRELAX*cs;
      f_dx=find_force(dev_neigth_points+offset_curr,n,x_dx,downP,bufP,normal0,cs,alpha,1,out);
      if( (isnan(f_dx.x)) && (isnan(f_dx.y)) && (isnan(f_dx.z))){
        printf("fdx %ld\n",i);
      }
      f_dy=find_force(dev_neigth_points+offset_curr,n,x_dy,downP,bufP,normal0,cs,alpha,1,out);
      if( (isnan(f_dy.x)) && (isnan(f_dy.y)) && (isnan(f_dy.z))){
        printf("fdy %ld\n",i);
      }
      J11=(f_dx*force_dir_x-f_n)/(DXRELAX*cs);
      J12=(f_dy*force_dir_x-f_n)/(DXRELAX*cs);
      J21=(f_dx*force_dir_y -0.)/(DXRELAX*cs);
      J22=(f_dy*force_dir_y -0.)/(DXRELAX*cs);
      tmpDouble=J11*J22-J12*J21;
      J_11=J22/tmpDouble;
      J_21=-J21/tmpDouble;
      dx=J_11*f_n;
      dy=J_21*f_n;
      x_dx=force_dir_x*dx*(-1)-force_dir_y*dy;
      x_dx=spring_relax_check(x0,x_dx,bufP,downP,dev_neigth_points+offset_prev,n,0.3,1);
      if( (!isnan(x_dx.x)) && (!isnan(x_dx.y)) && (!isnan(x_dx.z))){
        x0+=x_dx;
      }
    }
    dev_next_layer[i]=x0;    
  }
  __global__ void spring_relax_edge(
          long start,
          long nPoints,
          long* dev_nnei,
          math_our_cu::Point_cu* dev_neigth_points,
          math_our_cu::Point_cu* dev_next_layer,
          double alpha,
          long mem_struct_size,
          long* dev_exept
        ){
    long i = start+threadIdx.x + blockIdx.x*blockDim.x;
    if (i >= nPoints) return;
    long out=0;
    if(i==13737){
      out=1;
    }
    long n=dev_nnei[i];
    long offset_curr=mem_struct_size*i;
    long offset_prev=offset_curr+n;
    math_our_cu::Point_cu x0=dev_neigth_points[offset_prev+n];
    if(dev_exept[i]){
      return;
    }
    math_our_cu::Point_cu downP=dev_neigth_points[offset_prev+n+1];
    math_our_cu::Point_cu bufP=dev_neigth_points[offset_prev+n+2];
    math_our_cu::Point_cu f,f_dx;
    math_our_cu::Point_cu force_dir_x;
    math_our_cu::Point_cu x_dx;
    math_our_cu::Point_cu edge_n;
    double cs,f_n;
    double J11;
    double dx;
    math_our_cu::Point_cu normal0=bufP-downP;
    cs=normal0.module();
    normal0/=cs;
    edge_n=(dev_neigth_points[offset_prev]-bufP)^(dev_neigth_points[offset_prev+n-1]-bufP);
    edge_n.normalize();
//     force_dir_x=(edge_n^normal0).normalize();
    math_our_cu::Point_cu normal;
    for(long iter=0;iter<2;++iter){
      normal=(x0-downP).normalize();
      force_dir_x=(edge_n^normal).normalize();
      f=find_force(dev_neigth_points+offset_curr,n,x0,downP,bufP,normal0,cs,alpha,0,out);
      if( (isnan(f.x)) && (isnan(f.y)) && (isnan(f.z))){
        printf("f   %ld\n",i);
      }
      f_n=f*force_dir_x;
      if(fabs(f_n)<1e-10){
        dev_next_layer[i]=x0; 
        return;
      }
      x_dx=x0+force_dir_x*DXRELAX*cs;
      f_dx=find_force(dev_neigth_points+offset_curr,n,x_dx,downP,bufP,normal0,cs,alpha,0,out);
      if( (isnan(f_dx.x)) && (isnan(f_dx.y)) && (isnan(f_dx.z))){
        printf("fdx   %ld\n",i);
      }
      J11=(f_dx*force_dir_x-f_n)/(DXRELAX*cs);
      dx=-f_n/J11;
      x_dx=force_dir_x*dx;
      x_dx=spring_relax_check(x0,x_dx,bufP,downP,dev_neigth_points+offset_prev,n,0.3,0);
      if(out){
        printf("f    %lf %lf %lf\n",f.x,f.y,f.z);
        printf("f_dx %lf %lf %lf\n",f_dx.x,f_dx.y,f_dx.z);
        printf("n %lf %lf %lf\n",normal.x,normal.y,normal.z);
        printf("f_n %lf J11 %lf dx %lf x+dx %lf\n",f_n,J11,dx,(x0+x_dx-downP).module());
        printf("xdx   %lf %lf %lf\n",x_dx.x,x_dx.y,x_dx.z);
      }
      if( (!isnan(x_dx.x)) && (!isnan(x_dx.y)) && (!isnan(x_dx.z))){
        x0+=x_dx;
      }
    }
    dev_next_layer[i]=x0;
  }
  void Stl_io::llrelax_cu_fast(long nStep,double elas,double elas2){
    vector<math_our::Point> layerBuffer;
    math_our::Point tmpPoint;
//     double tmpDouble;
    long tmpLong;
    //------------------------exeption list--------------------------------------------
    ulong k=0;
    set<edge> edgeTree;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        addTriInEdgeTree(k,*j,edgeTree);
        k++;
      }
    }
    //------------------------------make owner<neighbour--------------------------
    long* exept;
    exept=(long*)malloc(nPoints*sizeof(long));
    for(long i=0;i<nPoints;i++){
      exept[i]=1;
    }
    for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
      if(i->neig==-1){
        exept[i->a]=0;
        exept[i->b]=0;
      }
    }  
    //--------------------------neighbours------------------------------------
    vector<vector<long> > neiList;
    vector<long> pointNeiPoint;
    vector<long> nei_offset;
    vector<long> pointNNeiPoint;
    pointNNeiPoint.reserve(nPoints);
    neiList.resize(nPoints);
    nei_offset.reserve(nPoints);
    
    findPnP(neiList);
    long pNP=0;
    nei_offset.push_back(0);
    ulong max_nnei=0;
    for(long i=0;i<nPoints;++i){
      pointNNeiPoint.push_back(neiList[i].size());
      pNP+=neiList[i].size();
      max_nnei=max_nnei>neiList[i].size()?max_nnei:neiList[i].size();
      if(i>0){
        nei_offset.push_back(nei_offset[i-1]+pointNNeiPoint[i-1]);
      }
    }
    pointNeiPoint.reserve(pNP);
    for(long i=0;i<nPoints;++i){
      pointNeiPoint.insert(pointNeiPoint.end(),neiList[i].begin(),neiList[i].end());
    }
    layerBuffer.reserve(nPoints);
    layerBuffer.insert<math_our::Point*>(layerBuffer.end(),&(pointArray[nLayers*nPoints]),&(pointArray[nPoints*(nLayers+1)]));
    math_our::Point* layer_result_buffer;
    math_our::Point *neigth_points;
    layer_result_buffer=(math_our::Point*)malloc(nPoints*sizeof(math_our::Point));
    long mem_struct_size=(long)(2*max_nnei+4);
    neigth_points=(math_our::Point*)malloc(mem_struct_size*nPoints*sizeof(math_our::Point));
    ulong tmpUlong;
    for(long i=0;i<nPoints;i++){
      tmpUlong=neiList[i].size();
      for(ulong j=0;j<tmpUlong;j++){
        tmpLong=pointNeiPoint[nei_offset[i]+j];
        tmpPoint=layerBuffer[tmpLong];
        neigth_points[mem_struct_size*i+j]=tmpPoint;
        tmpPoint=pointArray[(nLayers-1)*nPoints+tmpLong];
        neigth_points[mem_struct_size*i+j+tmpUlong]=tmpPoint;
      }
      tmpPoint=layerBuffer[i];
      neigth_points[mem_struct_size*i+2*tmpUlong]=tmpPoint;
      neigth_points[mem_struct_size*i+2*tmpUlong+2]=tmpPoint;
      tmpPoint=pointArray[(nLayers-1)*nPoints+i];
      neigth_points[mem_struct_size*i+2*tmpUlong+1]=tmpPoint;
      
    }      
    //---------------------------------------------------------------------------------
    long* dev_nei;
    long* dev_nnei;
    long* dev_offset_nei;
    long* dev_exept;
    math_our_cu::Point_cu* dev_neigth_points;
    math_our_cu::Point_cu* dev_next_layer;
    cudaMalloc( (void**)&dev_nei, pNP*sizeof(long) );
    cudaMalloc( (void**)&dev_nnei, nPoints*sizeof(long) );
    cudaMalloc( (void**)&dev_exept, nPoints*sizeof(long) );
    cudaMalloc( (void**)&dev_offset_nei, nPoints*sizeof(long) );
    cudaMalloc( (void**)&dev_next_layer, nPoints*sizeof(math_our::Point) );
    cudaMalloc( (void**)&dev_neigth_points,mem_struct_size*nPoints*sizeof(math_our_cu::Point_cu) );
    cudaMemcpy((void*)dev_nei, (void*)&(pointNeiPoint[0]), pNP*sizeof(long), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_offset_nei, (void*)&(nei_offset[0]), nPoints*sizeof(long), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_nnei, (void*)&(pointNNeiPoint[0]), nPoints*sizeof(long), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_exept, (void*)&(exept[0]), nPoints*sizeof(long), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_next_layer, (void*)&(layerBuffer[0]), nPoints*sizeof(math_our::Point), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_neigth_points, (void*)&(neigth_points[0]), mem_struct_size*nPoints*sizeof(math_our_cu::Point_cu), cudaMemcpyHostToDevice);
    ulong mem_used=mem_struct_size*nPoints*sizeof(math_our_cu::Point_cu)+nPoints*sizeof(math_our::Point)+nPoints*sizeof(long);    
    printf("mem used %10lu KB\n",mem_used/1024);
    ulong start_time =  clock();
    for(long n_elas=0;n_elas<nStep;++n_elas){
      for(long internal_i=0;internal_i<nPoints/(long)2048+(long)1;++internal_i){
        long threads = 256;
        long blocks = 8;//(nPoints)/threads+1;
        long start = internal_i*2048;
        makelrelax_fast<<<blocks, threads>>>(
          start,
          nPoints,
          dev_nnei,
          dev_neigth_points,
          dev_next_layer,
          elas,
          mem_struct_size,
          dev_exept
        );
      }      
      for(long internal_i=0;internal_i<nPoints/(long)2048+(long)1;++internal_i){
        long threads = 256;
        long blocks = 8;//(nPoints)/threads+1;
        long start = internal_i*2048;
        refreshnei<<<blocks, threads>>>(
          start,
          nPoints,
          dev_nei,dev_nnei,dev_offset_nei,
          dev_neigth_points,
          dev_next_layer,
          mem_struct_size
        );
      }
      for(long internal_i=0;internal_i<nPoints/(long)2048+(long)1;++internal_i){
        long threads = 256;
        long blocks = 8;//(nPoints)/threads+1;
        long start = internal_i*2048;
        spring_relax<<<blocks, threads>>>(
          start,
          nPoints,
          dev_nnei,
          dev_neigth_points,
          dev_next_layer,
          elas2,
          mem_struct_size,
          dev_exept
        );
      }
      for(long internal_i=0;internal_i<nPoints/(long)2048+(long)1;++internal_i){
        long threads = 256;
        long blocks = 8;//(nPoints)/threads+1;
        long start = internal_i*2048;
        spring_relax_edge<<<blocks, threads>>>(
          start,
          nPoints,
          dev_nnei,
          dev_neigth_points,
          dev_next_layer,
          elas2,
          mem_struct_size,
          dev_exept
        );
      }
      for(long internal_i=0;internal_i<nPoints/(long)2048+(long)1;++internal_i){
        long threads = 256;
        long blocks = 8;//(nPoints)/threads+1;
        long start = internal_i*2048;
        refreshnei<<<blocks, threads>>>(
          start,
          nPoints,
          dev_nei,dev_nnei,dev_offset_nei,
          dev_neigth_points,
          dev_next_layer,
          mem_struct_size
        );
      }
      if(nStep>40){
        long xy=nStep/40;
        if(!(n_elas%xy)){
          printf("#");
          fflush(stdout);
        }
      }
    }
    ulong end_time = clock();
    ulong search_time = end_time - start_time;
    fprintf(stdout,"\nRelax time (s) %ld,%06ld\n",
            search_time/1000000,search_time%1000000);
    fflush(stdout);
    cudaMemcpy((void*)&(pointArray[nLayers*nPoints]), (void*)dev_next_layer, nPoints*sizeof(math_our::Point), cudaMemcpyDeviceToHost);
//     tmpLong=0;
//     for(long i=0;i<nPoints;i++){
//       math_our::Point p1,p2;
//       p1=pointArray[nLayers*nPoints+i]-pointArray[(nLayers-1)*nPoints+i];
//       p2=layerBuffer[i]-pointArray[(nLayers-1)*nPoints+i];
//       p1.normalize();
//       p2.normalize();
//       if(p1*p2<0.999){
//         tmpLong++;
//       }
//     }
//     printf("%ld\n",tmpLong);
    cudaFree(dev_nnei);
    cudaFree(dev_next_layer);
    cudaFree(dev_neigth_points);
    cudaFree(dev_offset_nei);
    cudaFree(dev_nei);
    cudaFree(dev_exept);
    free(layer_result_buffer);
    free(neigth_points);
    free(exept);
  }
}

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