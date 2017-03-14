#include <stdio.h>
#include "stl_io.h"  
using namespace stl;
using namespace std;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io stl_;
  Stl_io stl_1;
  stl_.ReadBin(argv[1]);
  stl_1=stl_.absurface(argv[2]);
  stl_1.WriteOFF(argv[3]);
  stl_1.WriteBin("2d.bin");
  return 0;
} 
