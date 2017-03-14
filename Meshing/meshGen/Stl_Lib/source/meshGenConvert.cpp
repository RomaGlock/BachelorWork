#include "stl_io.h" 
#include "string.h"
using namespace stl;
using namespace std;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  if(argc>2){
    Stl_io stl_;
    char* arg1;
    char* arg2;
    arg1=stl::get_suff(argv[1]);
    arg2=stl::get_suff(argv[2]);
    if(!strcmp(arg1,"stl")){
      if(argc>3){
	double TOL;
	sscanf(argv[3],"%20lf",&TOL);
	stl_=Stl_io(argv[1],TOL,1);
      }else{
	stl_=Stl_io(argv[1],1e-7,1);
      }
    }else if(!strcmp(arg1,"off")){
      stl_.ReadOFF(argv[1]);
    }else if(!strcmp(arg1,"bin")){
      stl_.ReadBin(argv[1]);
    }else{
      printf("Aborted:\nstl,off,bin input formats are supported\n");
      return 0;
    }
    printf("input mesh has %ld points\n",stl_.nPoints);
    if(!strcmp(arg2,"stl")){
      stl_.Write(argv[2]);
    }else if(!strcmp(arg2,"off")){
      stl_.WriteOFF(argv[2]);
    }else if(!strcmp(arg2,"bin")){
      stl_.WriteBin(argv[2]);
    }else if(!strcmp(arg2,"lilu")){
      stl_.WriteLiLu(argv[2]);
    }else{
      printf("Aborted:\nstl,off,bin,lilu output formats are supported\n");
      return 0;
    }
  }else{
    printf("2 args are expected\n");
  }
  return 0;
}