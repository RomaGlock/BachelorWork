#include "stl_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
namespace stl
{
  using namespace std;
  void Stl_io::merge_out_in(const Stl_io& in,const char* fileName){
    if(nPoints!=in.nPoints){
      printf("!error nPoints does not match!\n");
      return;
    }
    adjust_numbers(in,fileName);
    double tmpDouble;
    vector<double> lList,lList2;
    lList.resize(nPoints);
    lList2.resize(nPoints);
    for(long i=0;i<nPoints;i++){
      lList[i]=0;
      for(long j=(nLayers-1);j>=0;j--){
        tmpDouble=(pointArray[(j+1)*nPoints+i]-pointArray[j*nPoints+i]).module();
        lList[i]+=tmpDouble;
      }
    }
    reverse_layers();
    pointArray.resize((nLayers+in.nLayers+1)*nPoints);
    for(long i=1;i<=in.nLayers;i++){
      pointArray.insert<vector<math_our::Point>::const_iterator>(
        pointArray.begin()+nPoints*(nLayers+i),
        in.pointArray.begin()+in.nPoints*i,
        in.pointArray.begin()+in.nPoints*(i+1)
      );
    }
    PChange_arg arg(fileName);
    FILE* file;
    file=fopen(arg.lengthFile,"rb");
    if(file){      
      read_length_file(&(lList2[0]),file);
      fclose(file);
    }else{
      printf("%s: Info: no lengthFile; current Length is used\n",__FILE__);
    }
    vector<double> profile;
    file=fopen(arg.profileFile,"rb");
    if(file){      
      profile_read(profile,file);
      fclose(file);
    }else{
      printf("%s: Aborted!: no profileFile;\n",__FILE__);
      return;
    }
    file=fopen(arg.profileFile,"wb");
    for(long i=(profile.size()-1);i>=0;--i){
      fprintf(file,"%14lf\n",profile[i]);
    }
    fclose(file);
//     void profile_read(std::vector<double>& profile,FILE* file){
    for(long i=0;i<nPoints;i++){
      lList[i]+=lList2[i];
    }
    file=fopen(arg.lengthFile,"wb");
    for(long i=0;i<nPoints;i++){
     fprintf(file,"%14lf\n",lList[i]);
    }
    fclose(file);
    changeProfile(fileName);
    reverse_layers();
    file=fopen(arg.lengthFile,"wb");
    for(long i=0;i<nPoints;i++){
     fprintf(file,"%14lf\n",lList2[i]);
    }
    fclose(file);
    file=fopen(arg.profileFile,"wb");
    for(ulong i=0;i<profile.size();++i){
      fprintf(file,"%14lf\n",profile[i]);
    }
    fclose(file);
  }
}
