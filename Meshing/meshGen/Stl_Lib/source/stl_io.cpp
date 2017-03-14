#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>


namespace stl
{
  using namespace std;
  
  void Stl_io::ReadStl(char* fileName)
  {
    FILE* file;
    long i;
    long i2=0;
    math_our::Point point;
    char string_[300];
    char string__buf[300];
    long inGroupFlag=0;
    GroupTree::iterator currentGroup;
//     PointTree pointTree;
    math_our::Map pointTree;
    GroupTree groupTree;
    
    
    long v[3];
    math_our::Point nn;
    math_our::Point vv[3];
    
    if(!(file=fopen(fileName,"rb"))){
      error(1,fileName,__LINE__);
    }
    i=0;
    while(fgets(string_,300,file)){
      i++;
    }
    rewind(file);
    pointArray.reserve(i);
    while(fgets(string_,300,file)){
      sscanf(string_," %40s",string__buf);    
      if(!strcmp(string__buf,"solid")){
	i=6;
	while(string_[i]!='\n') {
	  string__buf[i-6]=string_[i];
	  i++;
	}
	string__buf[i-6]='\0';
	if(inGroupFlag) error(2,"Group is not closed",__LINE__);
	currentGroup = groupTree.find(string__buf);
	if(currentGroup == groupTree.end()) {
	  Group g(string__buf);
	  pair<GroupTree::iterator, bool> groupInserted = groupTree.insert(pair<string,Group>(string__buf,g));
	  currentGroup = groupInserted.first;
	}
	if(currentGroup != groupTree.end())inGroupFlag=1;
      }else{
	if(!strcmp(string__buf,"endsolid")){
	  if(!inGroupFlag) error(2,"Group is not opened",__LINE__);
	  inGroupFlag=0;
	}else{
	  if(!strcmp(string__buf,"facet")){
	    if(sscanf(string_,"%*s normal %20lf %20lf %20lf",&(nn.x),&(nn.y),&(nn.z))==4)error(2,"strange facet line",__LINE__);
	    i2=0;
	    while(i2!=3){
	      fgets(string_,300,file);
	      if(sscanf(string_," vertex %20lf %20lf %20lf",&(vv[i2].x),&(vv[i2].y),&(vv[i2].z))==3){
// 		pair<PointTree::iterator, bool> pt = pointTree.insert(pair<Point, long>(vv[i2], pointArray.size()));
		pair<math_our::Map::iterator, bool> pt = pointTree.insert(pair<Point, long>(vv[i2], pointArray.size()));
		if(pt.second) {
		  v[i2]=pointArray.size();
		  pointArray.push_back(vv[i2]);
		}else{
		  v[i2]=pt.first->second;
		}
		i2++;
	      }else{
		i2=0;
	      }
	    }
	    currentGroup->second.tri.push_back(Triangle(v[0],v[1],v[2],nn,0));
	  }
	}
      }
    }
    fclose(file);
    groupArray.reserve(groupTree.size());
    for(GroupTree::iterator i = groupTree.begin(); i != groupTree.end(); ++i) {
      groupArray.push_back(i->second);
    }
    long groupArraySize=groupArray.size();
    for(i=0;i<groupArraySize;i++){
      for(Group::TriangleArray::iterator j = groupArray[i].tri.begin(); j != groupArray[i].tri.end(); ++j) {
	j->gr=i;
      }
    }
    
  }
  Stl_io::Stl_io():is_3Dmesh(false),is_Cutted(false),nLayers(0),nPoints(0){}
  Stl_io::Stl_io(char* fileName):is_3Dmesh(false),is_Cutted(false),nLayers(0),nPoints(0){
    ReadStl(fileName);
    nPoints=pointArray.size();
    Delete_duplicate();
  }
  Stl_io::Stl_io(char* fileName, double TOL):is_3Dmesh(false),is_Cutted(false),nLayers(0),nPoints(0){
    math_our::TOL=TOL;
    ReadStl(fileName);
    nPoints=pointArray.size();
    Delete_duplicate();
  }
  Stl_io::Stl_io(char* fileName, double TOL,int aaa):is_3Dmesh(false),is_Cutted(false),nLayers(0),nPoints(0){
    math_our::TOL=TOL;
    ReadStl(fileName);
    nPoints=pointArray.size();
  }
  Stl_io::Stl_io(const Stl_io::PointArray& a,const Stl_io::GroupArray& b):
		pointArray(a),groupArray(b),is_3Dmesh(false),is_Cutted(false),nLayers(0),nPoints(pointArray.size())
  {
//     Delete_duplicate();
  }
  void Stl_io::Delete_duplicate(){
    std::set<Triangle> triTree;
    std::pair<std::set<Triangle>::iterator,bool> p;
    bool flag=false;
    Group::TriangleArray::iterator k;
    if(groupArray.begin()!=groupArray.end()){
      k=groupArray.begin()->tri.begin();
    }
    for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	if(flag){
          i->tri.erase(k);
	}
	p=triTree.insert((*j));
	if(!(p.second)){
	  k=j;
	  flag=true;
	}else{
	  flag=false;
	}
      }
    }
//     for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
//       i->tri.clear();
//     }
//     for(std::set<Triangle>::iterator i=triTree.begin();i!=triTree.end();++i){
//       groupArray[i->gr].tri.push_back(*i);
//     }
  }
  
  void Stl_io::Write(const char* fileName)const{
    FILE *file;
    file=fopen(fileName,"wb");
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      fprintf(file,"solid %s\n",i->name.c_str());
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	fprintf(file,"  facet normal %12le %15le %15le\n",j->n.x,j->n.y,j->n.z);
	fprintf(file,"    outer loop\n");
	fprintf(file,"       vertex %.12le %.12le %.12le\n",
		pointArray[j->a].x,pointArray[j->a].y,pointArray[j->a].z);
	fprintf(file,"       vertex %.12le %.12le %.12le\n",
		pointArray[j->b].x,pointArray[j->b].y,pointArray[j->b].z);
	fprintf(file,"       vertex %.12le %.12le %.12le\n",
		pointArray[j->c].x,pointArray[j->c].y,pointArray[j->c].z);
	fprintf(file,"    endloop\n");
	fprintf(file,"  endfacet\n");
      }
      fprintf(file,"endsolid %s\n",i->name.c_str());
    }
    fclose(file);
  }
  void Stl_io::WriteBin(const char* fileName)const{
    long pas;
    long gas;
    FILE *file;
    file=fopen(fileName,"wb");
    //pointArray
    pas=pointArray.size();
    math_our::Point* pointArrayD;
    fwrite((void*)&pas,sizeof(long),(size_t)1,file);
    pointArrayD=(math_our::Point*)malloc(sizeof(math_our::Point)*pas);
    for(long i=0;i<pas;i++)pointArrayD[i]=pointArray[i];
    fwrite((void*)pointArrayD,sizeof(math_our::Point),(size_t)pas,file);
    free(pointArrayD);
    //groupArray
    long ss;
    long ts;
    char name[300];
    gas=groupArray.size();
    fwrite((void*)&gas,sizeof(long),(size_t)1,file);
    Triangle* tri;
    for(long i=0;i<gas;i++){
      ss=groupArray[i].name.size();
      for(long j=0;j<ss;j++){
	name[j]=groupArray[i].name[j];
      }
      fwrite((void*)&ss,sizeof(long),(size_t)1,file);
      fwrite((void*)&name,sizeof(char),(size_t)ss,file);
      ts=groupArray[i].tri.size();
      tri=static_cast<Triangle*>(malloc(sizeof(Triangle)*ts));
      long k=0;
      for(list<Triangle>::const_iterator j=groupArray[i].tri.begin();j!=groupArray[i].tri.end();++j){
	tri[k]=*j;
	k++;
      }
      fwrite((void*)&ts,sizeof(long),(size_t)1,file);
      fwrite((void*)tri,sizeof(Triangle),(size_t)ts,file);
      free(tri);
    }
    fwrite((void*)&is_3Dmesh,sizeof(bool),(size_t)1,file);
    fwrite((void*)&is_Cutted,sizeof(bool),(size_t)1,file);
    fwrite((void*)&nLayers,sizeof(long),(size_t)1,file);
    fwrite((void*)&nPoints,sizeof(long),(size_t)1,file);
    fclose(file);
  }
  long Stl_io::ReadBin(const char* fileName){
    long pas;
    long gas;
    
    FILE *file;
    file=fopen(fileName,"rb");
    if(!file)return 0;
    if(pointArray.size())pointArray.clear();
    if(groupArray.size())groupArray.clear();
    //pointArray
    fread((void*)&pas,sizeof(long),(size_t)1,file);
    pointArray.resize(pas);
    math_our::Point *pointArrayD;
    pointArrayD=(math_our::Point*)malloc(sizeof(math_our::Point)*pas);
    fread((void*)pointArrayD,sizeof(math_our::Point),(size_t)pas,file);
    for(long i=0;i<pas;i++)pointArray[i]=pointArrayD[i];
    free(pointArrayD);
    //groupArray
    long ss;
    long ts;
    char name[300];
    fread((void*)&gas,sizeof(long),(size_t)1,file);
    groupArray.resize(gas);
    Triangle* tri;
    for(long i=0;i<gas;i++){
      fread((void*)&ss,sizeof(long),(size_t)1,file);
      fread((void*)&name,sizeof(char),(size_t)ss,file);
      fread((void*)&ts,sizeof(long),(size_t)1,file);     
      tri=static_cast<Triangle*>(malloc(sizeof(Triangle)*(ts)));
      fread((void*)(tri/*+hmr*/),sizeof(Triangle),(size_t)ts,file);
      groupArray[i].name.resize(ss);
      for(long j=0;j<ss;j++){
	groupArray[i].name[j]=name[j];
      }
      for(long j=0;j<ts;j++){
	groupArray[i].tri.push_back(tri[j]);
      }
      free(tri);
    }
    fread((void*)&is_3Dmesh,sizeof(bool),(size_t)1,file);
    fread((void*)&is_Cutted,sizeof(bool),(size_t)1,file);
    fread((void*)&nLayers,sizeof(long),(size_t)1,file);
    fread((void*)&nPoints,sizeof(long),(size_t)1,file);   
    fclose(file);
    return 1;
  }
  void Stl_io::WriteOFF(const char* fileName)const{
    FILE *file;
    file=fopen(fileName,"wb");
    long k=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	k++;
      }
    }
    fprintf(file,"OFF\n");
    fprintf(file,"%ld %ld %ld\n",nPoints,k,(long)0);
    for(long i=0;i<nPoints;i++){
      fprintf(file,"%14le %14le %14le\n",pointArray[i].x,pointArray[i].y,pointArray[i].z);
    }
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	fprintf(file,"%d %ld %ld %ld\n",3,j->a,j->b,j->c);
      }
    }
    fclose(file);
  }
  long Stl_io::ReadOFF(const char* fileName){
    long pas;
    long ts;
    long e_n;
    long m;
    long a,b,c;
    double x,y,z;
    FILE *fin;
    char line[256];
    Triangle tri;
    if(pointArray.size())pointArray.clear();
    if(groupArray.size())groupArray.clear();
    //pointArray
     if ((fin = fopen(fileName, "r"))==NULL){
	printf("Read error. File \'%s\' could not be read.\n", fileName);
	return 0;
    };
  
    while (fgets(line, 256, fin) != NULL) {
      if (line[0]=='O' && line[1]=='F' && line[2]=='F')
	break;
    }    
    if (fscanf(fin,"%20ld %20ld %20ld\n", &pas, &ts, &e_n) != 3){
      printf("Read error. Expected 3 integer number for the number of vertices "
	    "and simplices.\n");
      fclose(fin);
      return 0;
    }
    
    printf("   vertices: %ld --- simplices: %ld \n", pas, ts);
  
    // Assign memory
    pointArray.resize(pas);
    groupArray.resize(1);
    groupArray[0].name="surface";
    // Read vertex coordinates
    for (long i = 0; i < pas; i++) {
      if (fscanf(fin, "%20lf %20lf %20lf\n", &x, &y, &z) != 3){
	printf("Read error. Expected 3 floats for the coordinates of vertex: %ld.\n", i);
    fclose(fin);
	return 0;
      }
      pointArray[i].x = x;
      pointArray[i].y = y;
      pointArray[i].z = z;
    }
    
    // Read the number of vertices per simplex from the first line of simplices
    fscanf(fin, "%20ld", &m);
    
  // Check format of the read simplices
    if (m != 3 && m != 4) {
      printf("Read error. Expected a 3 or 4 for the first value in the first simplex line.\n");
      fclose(fin);
      return 0;
    }
  
  // Input is surface mesh
    else if (m == 3) {
      printf("   Input is surface mesh.\n");    
      // Read the rest of the first line
      if (fscanf(fin, "%20ld %20ld %20ld", &a, &b, &c) != 3){
	printf("Read error. Expected 3 integers for the first simplex.\n");
    fclose(fin);
	return 0;
      }
      tri.a=a;
      tri.b=b;
      tri.c=c;
      tri.n=(pointArray[b]-pointArray[a])^(pointArray[c]-pointArray[a]).normalize();
      tri.gr=0;
      groupArray[0].tri.push_back(tri);
      for (long i = 1; i < ts; i++) {
	if (fscanf(fin, "%20ld %20ld %20ld %20ld", &m, &a, &b, &c) != 4){
	  printf("Read error. Expected 4 integers for simplex %ld.\n", i);
      fclose(fin);
	  return 0;
	}
	tri.a=a;
	tri.b=b;
	tri.c=c;
	tri.n=(pointArray[b]-pointArray[a])^(pointArray[c]-pointArray[a]).normalize();
	tri.gr=0;
	groupArray[0].tri.push_back(tri);
	// Skip any additional character on this line
	while (fgetc(fin) != '\n'){}
      }
    }
    nPoints=pointArray.size();
    // Close file
    fclose(fin);
    return 1;
  }
  void Stl_io::WriteLiLu(const char* fileName)const{
    printf("sizeof(int)     %5d\nsizeof(double)  %5d\n",(int)sizeof(int),(int)sizeof(double));
    FILE *file;
    file=fopen(fileName,"wb");
    math_our::Point tmpPoint;
    int pas=(int)nPoints;
    int tas=0;
    vector<double> faceCenter; 
    vector<double> faceNormals;
    vector<int> faceNNeiFaces;
    vector<int> pointNNeiFaces;
    vector<int> pointNNeiPoint;
    vector<int> faceNeiFaces;
    vector<int> pointNeiFaces;
    vector<int> pointNeiPoint;
    vector<int> facesArray;
    pointNNeiFaces.resize(nPoints);
    pointNNeiPoint.reserve(nPoints);
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	++tas;
      }
    }
    faceCenter.resize(tas*3);
    faceNormals.resize(tas*3);
    facesArray.reserve(tas*3);
    tas=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	tmpPoint=((pointArray[j->a]+pointArray[j->b]+pointArray[j->c])/3.);
	faceCenter[3*tas+0]=tmpPoint.x;
	faceCenter[3*tas+1]=tmpPoint.y;
	faceCenter[3*tas+2]=tmpPoint.z;
	tmpPoint=j->n/j->n.module();
	faceNormals[3*tas+0]=tmpPoint.x;
	faceNormals[3*tas+1]=tmpPoint.y;
	faceNormals[3*tas+2]=tmpPoint.z;
	facesArray.push_back((int)(j->a));
	facesArray.push_back((int)(j->b));
	facesArray.push_back((int)(j->c));
	++tas;
      }
    }
    printf("nPoints           %10d\n",pas);
    printf("nElements         %10d\n\n",tas);
    fflush(stdout);
    //-------------------neighbours--------------------------
    std::vector <std::list<int> >neiList;
    neiList.resize(nPoints);
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	neiList[j->a].push_back((int)(j->b));
	neiList[j->a].push_back((int)(j->c));
	neiList[j->b].push_back((int)(j->c));
	neiList[j->b].push_back((int)(j->a));
	neiList[j->c].push_back((int)(j->a));
	neiList[j->c].push_back((int)(j->b));
      }
    }
    for(long i=0;i<nPoints;++i){
      neiList[i].sort();
      neiList[i].unique();
    }
    std::vector< std::vector<int> >neiVec;
    neiVec.resize(nPoints);
    for(long i=0;i<nPoints;++i){
      neiVec[i].reserve(neiList[i].size());
      neiVec[i].assign(neiList[i].begin(),neiList[i].end());
    }
    long pNP=0;
    for(long i=0;i<nPoints;++i){
      pointNNeiPoint.push_back(neiList[i].size());
      pNP+=neiList[i].size();
    }
    pointNeiPoint.reserve(pNP);
    for(long i=0;i<nPoints;++i){
      pointNeiPoint.insert(pointNeiPoint.end(),neiList[i].begin(),neiList[i].end());
    }
    printf("pnp/p             %10lf\n",(double)pointNeiPoint.size()/(double)pas);
    fflush(stdout);
    set<edge> edgeTree;
    long k=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	addTriInEdgeTree(k,*j,edgeTree);
	k++;
      }
    }
    faceNNeiFaces.resize(k);
    long fNF=0;
    vector<long>offset_fnf;
    for(set<edge>::iterator face=edgeTree.begin();face!=edgeTree.end();++face){
      if(face->owner!=-1 && face->neig!=-1){
	faceNNeiFaces[face->owner]++;
	faceNNeiFaces[face->neig]++;
	fNF+=2;
      }
    }
    faceNeiFaces.resize(fNF,-1);
    offset_fnf.resize(k);
    for(long i=1;i<k;++i){
      offset_fnf[i]=offset_fnf[i-1]+faceNNeiFaces[i-1];
    }
    fflush(stdout);
    for(set<edge>::iterator face=edgeTree.begin();face!=edgeTree.end();++face){
      if(face->owner!=-1 && face->neig!=-1){
	faceNeiFaces[offset_fnf[face->owner]]=face->neig;
	faceNeiFaces[offset_fnf[face->neig]]=face->owner;
	offset_fnf[face->owner]++;
	offset_fnf[face->neig]++;
      }
    }
    printf("fnf/f             %10lf\n",(double)faceNeiFaces.size()/(double)tas);
    fflush(stdout);
    offset_fnf.clear();
    offset_fnf.resize(nPoints);
    long pnf=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	pointNNeiFaces[j->a]++;
	pointNNeiFaces[j->b]++;
	pointNNeiFaces[j->c]++;
	pnf+=3;
      }
    }
    for(long i=1;i<nPoints;++i){
      offset_fnf[i]=offset_fnf[i-1]+pointNNeiFaces[i-1];
    }
    k=0;
    pointNeiFaces.resize(pnf,-1);
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	pointNeiFaces[offset_fnf[j->a]]=k;
	pointNeiFaces[offset_fnf[j->b]]=k;
	pointNeiFaces[offset_fnf[j->c]]=k;
	offset_fnf[j->a]++;
	offset_fnf[j->b]++;
	offset_fnf[j->c]++;
	k++;
      }
    }
    printf("pnf/p             %10lf\n",(double)pointNeiFaces.size()/(double)pas);
    fflush(stdout);
    long info_err=0;
    for(long i=0;i<fNF;i++){
      if(faceNeiFaces[i]==-1)info_err++;
    }
    if(info_err)fprintf(stderr,"face-face connection error  %6ld\n",info_err);
    info_err=0;
    for(long i=0;i<pnf;i++){
      if(pointNeiFaces[i]==-1)info_err++;
    }
    if(info_err)fprintf(stderr,"point-face connection error %6ld\n",info_err);
    
    offset_fnf.clear();
    offset_fnf.resize(nPoints);
    int *one;
    int *three;
    one=(int*)malloc(sizeof(int)*tas);
    three=(int*)malloc(sizeof(int)*tas);
    for(int i=0;i<tas;i++){
      one[i]=1;
      three[i]=3;
    }
    fwrite((void*)&pas,sizeof(int),(size_t)1,file);
    fwrite((void*)&tas,sizeof(int),(size_t)1,file);
    fwrite((void*)&(pointArray[0]),sizeof(math_our::Point),(size_t)nPoints,file);
    fwrite((void*)three,sizeof(int),(size_t)(tas),file);
    fwrite((void*)&(facesArray[0]),sizeof(int),(size_t)(tas*3),file);
    fwrite((void*)&(faceNormals[0]),sizeof(double),(size_t)(tas*3),file);
    fwrite((void*)&(faceCenter[0]),sizeof(double),(size_t)(tas*3),file);
    fwrite((void*)one,sizeof(int),(size_t)(tas),file);
    fwrite((void*)&(pointNNeiFaces[0]),sizeof(int),(size_t)(pas),file);
    fwrite((void*)&(pointNNeiPoint[0]),sizeof(int),(size_t)(pas),file);
    fwrite((void*)&(faceNNeiFaces[0]),sizeof(int),(size_t)(tas),file);
    fwrite((void*)&(pointNeiFaces[0]),sizeof(int),(size_t)(pnf),file);
    fwrite((void*)&(pointNeiPoint[0]),sizeof(int),(size_t)(pNP),file);
    fwrite((void*)&(faceNeiFaces[0]),sizeof(int),(size_t)(fNF),file);
    free(one);
    free(three);
    fclose(file);
  }
  long Stl_io::find_point(const math_our::Point& p)const{
    long min;
    long max;
    long n,m;
    long  i;
    long middle;
    const math_our::Point *point_array;
    double a1=0.99*p.value();
    double a2=1.01*p.value();
    long nPoints;
    point_array=&pointArray[0];
    nPoints=pointArray.size();
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
    for(i=min;i<=max;i++)if((p-point_array[i]).module()<math_our::TOL) return i;
    return -1;
  }
//--------------------------------------------------------------------------------------------------------------------------
  void STLtoMESH(char *filename)
  {
	  char str[255];
	  sprintf(str, "gmsh %s.stl -2 -o %s.mesh", filename, filename);
	  system(str);
  }
  
}