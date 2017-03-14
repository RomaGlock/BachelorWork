#include <stdio.h>
#include "stl_io.h"  
using namespace stl;
using namespace std;
//3d settings
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io stl_;
  stl_.ReadBin(argv[1]);
  stl_.changeProfile(argv[2]);
  stl_.WriteBin(argv[1]);
  stl_.MakeOFMesh(argv[2]);
  return 0;
}  
