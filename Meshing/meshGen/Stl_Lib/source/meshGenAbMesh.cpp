#include <stdio.h>
#include "stl_io.h"  
using namespace stl;
using namespace std;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  if(argc == 5) {
    CudaSetDevice(1);
  }
  Stl_io stl_;
  stl_.ReadOFF(argv[1]);
  stl_.abmesh_cu(argv[2]);
//   stl_.reverse_layers();
  stl_.WriteBin(argv[3]);
  stl_.MakeOFMesh(argv[2]);
  stl_.ab_field(argv[2]);
  return 0;
}