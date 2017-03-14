#include <stdio.h>
#include "stl_io.h" 
using namespace stl;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io stl_;
  FILE *file;
  stl_.ReadOFF(argv[1]);
  file=fopen(argv[2],"wb");
  fprintf(file,"//\n");
  for(unsigned long i=0;i<stl_.pointArray.size();i++){
    fprintf(file,"Point(%lu)={ %14lf , %14lf , %14lf };\n",
	    i+1,stl_.pointArray[i].x,stl_.pointArray[i].y,stl_.pointArray[i].z);
  }
  fclose(file);
  return 0;
}