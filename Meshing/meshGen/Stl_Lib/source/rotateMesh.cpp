#include "stl_io.h"
#include <fstream>
#include <sstream>
#include <math.h>
#define M_PI       3.14159265358979323846

using namespace stl;
using namespace std;
using namespace math_our;

void rotateMesh(vector<math_our::Point> &points, const Point& a, const Point& v, double o)
{
	o* = M_PI/180.;
	const double xx = cos(o)+(1-cos(o))*v.x*v.x;
	const double xy = (1-cos(o))*v.x*v.y-sin(o)*v.z;
	const double xz = (1-cos(o))*v.x*v.z+(sin(o))*v.y;
	const double yx = (1-cos(o))*v.y*v.x+sin(o)*v.z;
	const double yy = cos(o)+(1-cos(o))*v.y*v.y;
	const double yz = (1-cos(o))*v.y*v.z-sin(o)*v.x;
	const double zx = (1-cos(o))*v.z*v.x-sin(o)*v.y;
	const double zy = (1-cos(o))*v.z*v.y+sin(o)*v.x;
	const double zz = cos(o)+(1-cos(o))*v.z*v.z;

	for(vector<Point>::iterator it=points.begin(); it!=points.end(); ++it)
	{
		*it -= a;
		Point buf1 = *it;
		Point buf2;
		buf2.x = xx*buf1.x + xy*buf1.y + xz*buf1.z;
		buf2.y = yx*buf1.x + yy*buf1.y + yz*buf1.z;
		buf2.z = zx*buf1.x + zy*buf1.y + zz*buf1.z;
		*it = buf2;
		*it += a;
	}
}

inline double readDouble(ifstream &os)
{
	double a;
	char str[256];
	os>>str;
	os>>str;
	istringstream s(str);
	s>>a;
	return a;
}

inline Point readPoint(ifstream &os)
{
	Point a;
	a.x = readDouble(os);
	a.y = readDouble(os);
	a.z = readDouble(os);
	return a;
}

bool readFile(char *nameFile, Point &a, Point &v, double &angle)
{
	bool flag = false;
	ifstream file(nameFile, ios::in);
	if(file.is_open())
	{
		a = readPoint(file);
		v = readPoint(file);
		angle = readDouble(file);
		flag = true;
		file.close();
	}
	return flag;
}

int main(int argc, char** argv)
{
  Stl_io stl_;
  Point a,v;
  double angle;
  if(readFile(argv[2], a, v, angle))
  {
	  stl_.ReadBin(argv[1]);
	  v.normalize();
	  rotateMesh(stl_.pointArray, a, v, angle);
	  if(argc > 3)
	  {
	      stl_.WriteBin(argv[3]);
	      cout<<"mesh rotate on "<<angle<<" and write to "<<argv[3]<<endl;
	  }
	  else
	  {
	      stl_.WriteBin(argv[1]);
	      cout<<"mesh rotate on "<<angle<<" and write to "<<argv[1]<<endl;
	  }
	  return 0;
  }
  cout<<"not read file: "<<argv[2]<<endl;
  return 0;
}