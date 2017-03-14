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
  Stl_io phantom_stl;
  stl_.ReadBin(argv[1]);  
  phantom_stl.ReadBin("meshGen/phantom.bin");
  stl_.MMPoi(argv[2],phantom_stl);
  stl_.WriteBin(argv[3]);
  stl_.MakeOFMesh(argv[2]);
  return 0;
} 
