#include "stl_io.h" 
using namespace stl;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io stl_;
  stl_=Stl_io(argv[1],1e-7,1);
  stl_.WriteOFF(argv[2]);
  return 0;
}