#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <string.h>

#define POINTS "constant/polyMesh/points"

namespace stl
{
  using namespace std;
  
  struct data{
    math_our::Point p;
    mutable double v[30];
    bool operator<(const data& a) const{
      if((*this)==a) return false;    
      else{
	if(this->p<a.p) return true;
	else return false;
      }
    }
    bool operator==(const data& a) const{
      if(this->p==a.p) 
	return true;    
      else
	return false;
    }
    bool operator<(const double& a) const{      
      if(this->p<a) return true;
      else return false;
    }
  };
  typedef std::multiset<data> DataTree;
  double* find(data* t,const data& p,unsigned long nPoints){
    long min;
    long max;
    long n,m;
    long  i;
    long middle;
    data *point_array;
    point_array=t;
    double a1=0.99999*p.p.value();
    double a2=1.00001*p.p.value();
    n=0;
    m=nPoints;
    
    while((m-n)>1){
      middle=n+(m-n)/2;
      if(point_array[middle]<a1)
	n=middle;
      else 
	m=middle;  
    }
    min=n;
    
    n=0;
    m=nPoints;
    while((m-n)>1){
      middle=n+(m-n)/2;
      if(point_array[middle]<a2)
	n=middle;
      else 
	m=middle;  
    }
    max=m;
    if(min>max){
      long buff;
      buff=max;
      max=min;
      min=buff;
    }
    for(i=min;i<=max;i++){
      if((p.p-point_array[i].p).module()<math_our::TOL) {
	return &(point_array[i].v[0]);
      }
    }
    return NULL;
  }
  
  double solve_r_sch(math_our::Point x0,math_our::Point e,double t,math_our::Point a);
  double FINTRPQ(
    double T0,
    double T1,
    double T2,
    double T3,
    double F0,
    double F1,
    double F2,
    double F3,
    double T
  );
  long Stl_io::findWave(const char* fileName)const{
    printf("info : adaptation | version from 24 nov 2016\n");
    fflush(stdout);
    findWave_arg arg(fileName);
    double eps=fabs(arg.eps);    
    double tmpDouble,limitL;
    double *tmpDoubleP;
    double *field,*currentLength,*wavePos;
    long *startLayer;    
    bool global_flag=false;
    limitL=arg.minL;
    FILE *file;    
    file=fopen(arg.fieldFile,"rb");
    field=(double*)malloc(pointArray.size()*sizeof(double));
    currentLength=(double*)malloc(nPoints*sizeof(double));
    wavePos=(double*)malloc(nPoints*sizeof(double));
    startLayer=(long*)malloc(nPoints*sizeof(long));
    long tmpLong;
    math_our::Point tmpPoint;
    double credit;
    double pos1,pos2;
    long next,prev;
    if(!file){
      global_flag=true;
    }else{
      if(arg.made==0){
	for(unsigned long i=0;i<pointArray.size();i++){
	  fscanf(file,"%20lf",field+i);
	}
	fclose(file);
      }else if(arg.made==1){
	long nVertex=0;
	char string[300];
	data tmpData;
	FILE *pfile;
	
	pfile=fopen(POINTS,"rb");
	while(nVertex==0 && fgets(string,300,pfile)){
	  sscanf(string,"%20ld",&nVertex);
	}
	fgets(string,300,pfile);
	if((unsigned long)nVertex!=pointArray.size()){
	  printf("OFmesh nPoints error\n");		
	  free(field);
	  free(currentLength);
	  free(wavePos);
	  free(startLayer);
	  fclose(pfile);
	  return -1;
	}
	DataTree dataTree;
	math_our::TOL=arg.TOL;
	for(long i=0;i<nVertex;i++){
	  fscanf(pfile,"(%20lf%20lf%20lf)\n",&(tmpData.p.x),&(tmpData.p.y),&(tmpData.p.z));
	  fscanf(file,"%20lf",&(tmpData.v[3]));
	  dataTree.insert(tmpData);
	}
	fclose(file);
	fclose(pfile);
	vector<data> dataVector;
	dataVector.reserve(dataTree.size());
	dataVector.assign(dataTree.begin(),dataTree.end());
	for(long i=0;i<(nLayers+1)*nPoints;i++){
	  tmpData.p=pointArray[i];
	  tmpDoubleP=find(&(dataVector[0]),tmpData,dataTree.size());
	  if(!tmpDoubleP){
	    printf("Error: points do not match\n");
	    free(field);
	    free(currentLength);
	    free(wavePos);
	    free(startLayer);
	    return 0;
	  }
	  field[i]=tmpDoubleP[3];
	}
      }else if(arg.made==2){
	long nVertex=pointArray.size();
	data tmpData;
	DataTree dataTree;
	math_our::TOL=arg.TOL;
	for(long i=0;i<nVertex;i++){
	  fscanf(file,"%20lf%20lf%20lf%20lf",&(tmpData.p.x),&(tmpData.p.y),&(tmpData.p.z),&(tmpData.v[3]));
	  dataTree.insert(tmpData);
	}
	fclose(file);
	vector<data> dataVector;
	dataVector.reserve(dataTree.size());
	dataVector.assign(dataTree.begin(),dataTree.end());
	for(long i=0;i<nVertex;i++){
	  tmpData.p=pointArray[i];
	  tmpDoubleP=find(&(dataVector[0]),tmpData,dataTree.size());
	  if(!tmpDoubleP){
	    printf("Error: points do not match\n# %ld :%.13le %.13le %.13le\n",i,tmpData.p.x,tmpData.p.y,tmpData.p.z);	    
	    free(field);
	    free(currentLength);
	    free(wavePos);
	    free(startLayer);
	    return 0;
	  }
	  field[i]=tmpDoubleP[3];
	}
      }else if(arg.made==3){
	fclose(file);
      }
      bool flag;
      double *curAbs;
      curAbs=(double*)malloc(sizeof(double)*(nLayers+1));
      for(long i=0;i<nPoints;i++){
	currentLength[i]=0;
	for(long j=(nLayers-1);j>=0;j--){
	  tmpDouble=(pointArray[(j+1)*nPoints+i]-pointArray[j*nPoints+i]).module();
	  currentLength[i]+=tmpDouble;
	}
      }
      long info_1=0;
      long info_2=0;
      long info_3=0;
      double maxlvl=0;
      for(long i=0;i<nPoints;i++){
	curAbs[0]=0.;
	for(long j=0;j<nLayers;++j){
	  tmpDouble=(pointArray[(j+1)*nPoints+i]-pointArray[j*nPoints+i]).module();
	  curAbs[j+1]=curAbs[j]+tmpDouble/currentLength[i];
	}
	wavePos[i]=0;
	startLayer[i]=-1;
	flag=false;
	for(long j=(nLayers-1);j>=0 && !flag;j--){
          if(j){
            tmpDouble=field[j*nPoints+i]/field[nLayers*nPoints+i];
            if(tmpDouble>1.+eps /*|| tmpDouble<1.-eps*/){
              flag=true;
              startLayer[i]=j;
              double T0,T1,T2,T3,F0,F1,F2,F3,T;
              T0=field[(j-1)*nPoints+i];
              T1=field[j*nPoints+i];
              T2=field[(j+1)*nPoints+i];
              T3=((j+2)>nLayers)?field[nLayers*nPoints+i]*0.999:field[(j+2)*nPoints+i];
              F0=curAbs[j-1];
              F1=curAbs[j];
              F2=curAbs[j+1];
              if(maxlvl<curAbs[j])maxlvl=curAbs[j];
              F3=((j+2)>nLayers)?(2.-curAbs[j]):curAbs[j+2];
              T=field[nLayers*nPoints+i]*(1.+arg.eps);
              tmpDouble=FINTRPQ(T0,T1,T2,T3,F0,F1,F2,F3,T);
              if(tmpDouble<F1 || tmpDouble>F2 ){
                tmpDouble=F1+(F2-F1)*(T-T1)/(T2-T1);
                wavePos[i]=currentLength[i]*tmpDouble;
                info_1++;
              }else{ 
                if(tmpDouble>0){
                  wavePos[i]=currentLength[i]*tmpDouble;
                }else{	      
                  tmpDouble=F1+(F2-F1)*(T-T1)/(T2-T1);
                  wavePos[i]=currentLength[i]*tmpDouble;
                  info_2++;
                }
              }
            }
          }else{
            wavePos[i]=currentLength[i];
          }
        }
      }
      printf("nPoints      :%8ld\n",nPoints);
      printf("ENO -> linear:%8ld\n",info_1);
      printf("ENO       out:%8ld\n",info_2);
      for(long i=0;i<nPoints;i++){
	if(wavePos[i]>0 && wavePos[i]<currentLength[i]){
	}else{
	  wavePos[i]=currentLength[i]*maxlvl;
	  info_3++;
	}
      }
      printf("ENO     error:%8ld\n",info_3);
      //-------------------------------x----------
      std::vector<math_our::Point> currentX;
      currentX.resize(nPoints);
      {
        for(long i=0;i<nPoints;i++){
          tmpLong=0;
          pos1=0;
          pos2=0;
          for(unsigned long j=0;j<1;j++){
            credit=wavePos[i];
            while(pos2<credit && tmpLong<nLayers){
              pos1=pos2;
              pos2+=(pointArray[i+nPoints*(tmpLong+1)]-pointArray[i+nPoints*tmpLong]).module();
              ++tmpLong;
            }
            tmpDouble=(credit-pos1)/(pos2-pos1);
            prev=i+nPoints*(tmpLong-1);
            next=i+nPoints*tmpLong;
            tmpPoint=(pointArray[next]-pointArray[prev]);
            currentX[i]=pointArray[prev]+tmpPoint*tmpDouble;
          }
        }
      }
      //--------------------------------------------------
      credit=1e100;
      for(long i=0;i<nPoints;i++){
        if(wavePos[i]<credit){
          credit=wavePos[i];
        }
      }
      CUT_arg cut_arg(fileName);
      double minX=1e100;
      double maxX=-1e100;
      for(long i=0;i<nPoints;i++) {
        if(currentX[i]*cut_arg.n<minX)minX=currentX[i]*cut_arg.n;
        if(currentX[i]*cut_arg.n>maxX)maxX=currentX[i]*cut_arg.n;
      }
      for(long i=0;i<nPoints;i++){
        tmpDouble=(currentX[i]*cut_arg.n-minX)/(maxX-minX);
        wavePos[i]+=credit*(arg.mul-1.)+tmpDouble*credit*(arg.minL-arg.mul+1.);
// 	if(((arg.mul-1.)*wavePos[i])<(tmpDouble*arg.minL)){
// 	  wavePos[i]*=arg.mul;
// 	}else{
// 	  wavePos[i]+=arg.minL*tmpDouble;
// 	}
      }
      arg.minL=credit;
      free(curAbs);
    }
    if(global_flag){
      printf("field file does not exist\n");
      tmpDouble=1e100;      
      for(long i=0;i<nPoints;i++){
	currentLength[i]=0;
	for(long j=(nLayers-1);j>=0;j--){
	  tmpDouble=(pointArray[(j+1)*nPoints+i]-pointArray[j*nPoints+i]).module();
	  currentLength[i]+=tmpDouble;
	}
	if(currentLength[i]<tmpDouble){
	  tmpDouble=currentLength[i];
	}
      }
      for(long i=0;i<nPoints;i++){
	wavePos[i]=currentLength[i]*((arg.mul+1.)/2.);
      }
//       if(tmpDouble<arg.minL)
      arg.minL=tmpDouble;
    }
    //--------------------------------------neighbour list---------------------------------------
    std::vector<std::list<long> > neiList;
    neiList.resize(nPoints);
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	neiList[j->a].push_back(j->b);
	neiList[j->a].push_back(j->c);
	neiList[j->b].push_back(j->c);
	neiList[j->b].push_back(j->a);
	neiList[j->c].push_back(j->a);
	neiList[j->c].push_back(j->b);
      }
    }
    for(long i=0;i<nPoints;++i){
      neiList[i].sort();
      neiList[i].unique();
    }
    std::vector<std::vector<long> >neiVec;
    neiVec.resize(nPoints);
    for(long i=0;i<nPoints;++i){
      neiVec[i].reserve(neiList[i].size());
      neiVec[i].assign(neiList[i].begin(),neiList[i].end());
    }
    //--------------------------------------relaxation-------------------------------------------
    std::vector<math_our::Point> currentPos;
    std::vector<math_our::Point> startPos;
    std::vector<math_our::Point> currentDir;
    std::vector<math_our::Point> currentNeig;
    currentPos.reserve(nPoints);
    startPos.reserve(nPoints);
    currentDir.reserve(nPoints);
    long info_1=0;
    for(long i_rel=0;i_rel<arg.nRel;++i_rel){
      for(long i=0;i<nPoints;i++){
	tmpLong=0;
	pos1=0;
	pos2=0;
	for(unsigned long j=0;j<1;j++){
	  credit=wavePos[i];
	  while(pos2<credit && tmpLong<nLayers){
	    pos1=pos2;
	    pos2+=(pointArray[i+nPoints*(tmpLong+1)]-pointArray[i+nPoints*tmpLong]).module();
	    ++tmpLong;
	  }
	  tmpDouble=(credit-pos1)/(pos2-pos1);
	  prev=i+nPoints*(tmpLong-1);
	  next=i+nPoints*tmpLong;
	  tmpPoint=(pointArray[next]-pointArray[prev]);
	  currentPos[i]=pointArray[prev]+tmpPoint*tmpDouble;
	  currentDir[i]=tmpPoint;
	}
      }
      if(!i_rel){
	for(long i=0;i<nPoints;i++){
	  startPos.push_back(currentPos[i]);
	}
      }
      for(long i=0;i<nPoints;i++){
	currentNeig.reserve(neiList[i].size()+1);
	for(std::list<long>::iterator j=neiList[i].begin();j!=neiList[i].end();++j){
	  currentNeig.push_back(currentPos[*j]);
	}
	currentNeig.push_back(startPos[i]);
	tmpDouble=solve(currentPos[i],currentDir[i],&(currentNeig[0]),currentNeig.size(),arg.elas);
	currentNeig.clear();
	tmpDouble*=arg.iterLimit;
	if((wavePos[i]+tmpDouble)>arg.minL 
	  && (wavePos[i]+tmpDouble)<(currentLength[i]*arg.mul) && 
	  (wavePos[i]+tmpDouble)<(currentLength[i]+arg.minL*limitL)){
	  wavePos[i]+=tmpDouble;	  
	}else{
	  info_1++;
	}
      }
      if(arg.nRel>40){
	long xy=arg.nRel/40;
	if(!(i_rel%xy)){
	  printf("#");
	  fflush(stdout);
	}
      }
    }
    printf("\n");
    printf("nRelax out   :%8ld\n",info_1);
    //-------------------------------------------------------------------------------------------
   
    
    file=fopen(arg.outputFile,"wb");
    for(long i=0;i<nPoints;i++){
     fprintf(file,"%14lf\n",wavePos[i]);
    }
    fclose(file);
    free(field);
    free(currentLength);
    free(wavePos);
    free(startLayer);
    return 1;
  }
  double solve(
    const math_our::Point& x0,math_our::Point& e,
	       const math_our::Point* a,unsigned long n,
	       double alpha
	      )
  {
    double t=0;
    double cs=e.module();
    e.normalize();
    double tmpDouble;
    double N=double(n-1);
    double f_e=0.;
    math_our::Point tmpPoint;
    math_our::Point tau;
    math_our::Point f;
    f=math_our::Point(0.,0.,0.);
    for(unsigned long i=0;i<n-1;++i){
	f+=a[i];  
    }
    f-=x0*N;
    f_e=f*e;
    double f_;
    double f_elas=0.;
    double dfdt=0.;
    double dfdt_elas=0.;
    for(long iter=0;iter<10;++iter){
      tmpPoint=a[n-1]-x0-e*t;
      f_=f_e-N*t;
      f_elas=tmpPoint*e*(tmpPoint).module();
      dfdt=-N;
      dfdt_elas=solve_r_sch(x0,e,t,a[n-1])*(tmpPoint*e)-(tmpPoint).module();
      f_elas/=cs;
      dfdt_elas/=cs;
      f_elas*=alpha;
      dfdt_elas*=alpha;
      tmpDouble=(f_+f_elas);
      tmpDouble/=(dfdt+dfdt_elas);
      t-=tmpDouble;
    }
    return t;
  }
  double solve_r_sch(math_our::Point x0,math_our::Point e,double t,math_our::Point a){
    return 0.5*pow((a-x0)*(a-x0)+e*(a-x0)*(-2.)*t+t*t,0.5)*(e*(a-x0)*(-2.)+2.*t);
  }
  double FINTRPQ(double T0,double T1,double T2,double T3,double F0,double F1,double F2,double F3,double T){
    double kt,T10,T21,T32,T20,T31,A1,A2;
    T10 = T1 - T0;
    T21 = T2 - T1;
    T32 = T3 - T2;
    T20 = T2 - T0;
    T31 = T3 - T1;
    A1 = F0/(T10*T20) - F1/(T10*T21) + F2/(T21*T20);
    A2 = F1/(T21*T31) - F2/(T21*T32) + F3/(T32*T31);

    if(A1*A2<=0.){
	  kt = (T-T1)/(T2-T1);
	  return (1.-kt)*F1+kt*F2;
    }else if(fabs(A1)<=fabs(A2)){
	  return F0*(T-T1)*(T-T2)/(T10*T20) - F1*(T-T0)*(T-T2)/(T10*T21) + F2*(T-T0)*(T-T1)/(T21*T20);
    }else{
	  return F1*(T-T2)*(T-T3)/(T21*T31) - F2*(T-T1)*(T-T3)/(T21*T32) + F3*(T-T1)*(T-T2)/(T32*T31);
    }
  }
}

