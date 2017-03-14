#include "stl_io.h"  
using namespace stl;
using namespace std;
int main(int argc, char** argv){
  printf("%s\n",argv[0]);
  printf("compilation from %s ",__DATE__);
  printf("%s\n",__TIME__);
  Stl_io out,in;
  out.ReadBin(argv[1]);
  in.ReadBin(argv[2]);
  out.merge_out_in(in,argv[3]);
  out.WriteBin(argv[4]);
  out.MakeOFMesh(argv[3]);
  return 0;
}  
 
