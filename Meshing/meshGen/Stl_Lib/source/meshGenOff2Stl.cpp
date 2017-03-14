#include "stl_io.h" 
using namespace stl;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io stl_;
  stl_.ReadOFF(argv[1]);
  stl_.Write(argv[2]);
  return 0;
}