#include <stdio.h>
#include "stl_io.h"  
using namespace stl;
using namespace std;
//2d 3d settings field
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io phantom_stl;
  phantom_stl.ReadBin(argv[1]);
  phantom_stl.findWave(argv[2]);
  return 0;
} 