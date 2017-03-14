#include "stl_io.h"
using namespace stl;
using namespace std;
int main(int argc, char** argv){
  Stl_io stl_;
  stl_.ReadBin(argv[1]);
  stl_.MakeOFMesh(argv[2]);
  stl_.WriteOFF("meshGen/of.off");
  return 0;
}