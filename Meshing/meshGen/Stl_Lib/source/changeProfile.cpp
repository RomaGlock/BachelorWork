#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
namespace stl
{
  using namespace std; 
  void Stl_io::changeProfile(const char* fileName){
    PChange_arg arg(fileName);
    vector<double> profile;
    long tmpLong;
    double tmpDouble;
    double *lList;
    lList=(double*)malloc(sizeof(double)*pointArray.size());
    FILE* file;
    file=fopen(arg.profileFile,"rb");
    if(!file){
      printf("%s does not exist\n",arg.profileFile);
      free(lList);      
      return;
    }else{
      profile_read(profile,file);
      fclose(file);
    }
    file=fopen(arg.lengthFile,"rb");
    for(long i=0;i<nPoints;i++){
      lList[i]=0;
      for(long j=(nLayers-1);j>=0;j--){
	tmpDouble=(pointArray[(j+1)*nPoints+i]-pointArray[j*nPoints+i]).module();
	lList[i]+=tmpDouble;
      }
    }
    if(file){      
      read_length_file(lList,file);
      fclose(file);
    }else{
      printf("%s: Info: no lengthFile; current Length is used\n",__FILE__);
    }
    if(arg.createProfile){
      file=fopen(arg.profileFile,"wb");
      if(file){
	for(unsigned long i=0;i<profile.size();i++){
	  fprintf(file,"%lf\n",profile[i]);
	}
	fclose(file);
      }
    }
    vector<math_our::Point> newPointArray;
    newPointArray.resize(nPoints*(profile.size()+1));
    for(unsigned long j=1;j<profile.size();j++){
      profile[j]+=profile[j-1];
    }
    double credit;
    double pos1,pos2;
    long next,prev;
    for(long i=0;i<nPoints;i++){
      newPointArray[i]=pointArray[i];
    }
    for(long i=0;i<nPoints;i++){
      tmpLong=0;
      pos1=0;
      pos2=0;
      for(unsigned long j=0;j<profile.size();j++){
	credit=profile[j]*lList[i];
	while(pos2<credit && tmpLong<nLayers){
	  pos1=pos2;
	  pos2+=(pointArray[i+nPoints*(tmpLong+1)]-pointArray[i+nPoints*tmpLong]).module();
	  ++tmpLong;
	}
	tmpDouble=(credit-pos1)/(pos2-pos1);
	prev=i+nPoints*(tmpLong-1);
	next=i+nPoints*tmpLong;
	newPointArray[i+(j+1)*nPoints]=pointArray[prev]+(pointArray[next]-pointArray[prev])*tmpDouble;
      }     
    }
    pointArray.clear();
    pointArray.resize(nPoints*(profile.size()+1));
    swap(newPointArray,pointArray);
    nLayers=profile.size();
    free(lList);
  }
}