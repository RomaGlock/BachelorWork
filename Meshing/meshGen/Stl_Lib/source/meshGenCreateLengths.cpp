#include <stdio.h>
#include "stl_io.h"  
using namespace stl;
using namespace std;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io stl_;
  stl_.ReadBin(argv[1]);
  stl_.createLengthFile(argv[2]);
  return 0;
} 