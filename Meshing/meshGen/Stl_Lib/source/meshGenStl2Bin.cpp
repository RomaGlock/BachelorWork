#include <stdio.h>
#include "stl_io.h"  
using namespace stl;
using namespace std;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io stl_;
  Stl_io phantom_stl;
  stl_=Stl_io(argv[1],1e-7);
  phantom_stl=stl_.Cut_stl(argv[2]);
  stl_.WriteBin(argv[3]);
  stl_.Write("meshGen/current.stl");
  phantom_stl.WriteBin("meshGen/phantom.bin");
  phantom_stl.Write("meshGen/phantom.stl");
  return 0;
}