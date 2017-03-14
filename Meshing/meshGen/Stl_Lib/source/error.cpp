#include "stl_io.h"
#include<stdio.h>
namespace stl
{
  using namespace std;
  long error(long i,const char *s,int d){
    switch (i){
      case 1:
	fprintf(stderr,"Can not open file %s. (%d)\n",s,d);
	return -1;
      case 2:
	fprintf(stderr,"Wrong stl: %s. (%d)\n",s,d);
	return -1;
      case 3:
	fprintf(stderr,"Internal error: %s. (%d)\n",s,d);
	return -1;
    }
    return -1;
  } 
}