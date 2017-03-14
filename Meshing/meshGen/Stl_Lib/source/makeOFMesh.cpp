#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>

#define POINTFILE "constant/polyMesh/points"
#define OWNERFILE "constant/polyMesh/owner"
#define NEIGHFILE "constant/polyMesh/neighbour"
#define FACESFILE "constant/polyMesh/faces"
#define BOUNDFILE "constant/polyMesh/boundary" 
namespace stl
{
  using namespace std;
  long Stl_io::addTriInEdgeTree(long nTri,const Triangle& p,std::set<edge>& Tree)const{
    long a1,b1,c1,buff;
    a1=p.a;
    b1=p.b;
    c1=p.c;
    if(a1>b1){
      buff=a1;
      a1=b1;
      b1=buff;
    }
    if(b1>c1){
      buff=b1;
      b1=c1;
      c1=buff;
    }
    if(a1>b1){
      buff=a1;
      a1=b1;
      b1=buff;
    }
    if(b1>c1){
      buff=b1;
      b1=c1;
      c1=buff;
    }
    math_our::Point zero(0,0,0);
    std::set<edge>::iterator founded_edge;
    founded_edge=Tree.find(edge(a1,b1,-1,-1,-1,-1,zero,0));
    if(founded_edge==Tree.end()){
      Tree.insert(edge(a1,b1,nTri,-1,c1,-1,p.n,p.gr));
    }else if(founded_edge->neig!=-1) {
      error(3,"Error 3rd neighbour",__LINE__);
    }else {
      founded_edge->neig=nTri;
      founded_edge->neig3rd=c1;
      founded_edge->n=p.n;
      founded_edge->group=p.gr;
    }
    founded_edge=Tree.find(edge(a1,c1,-1,-1,-1,-1,zero,0));
    if(founded_edge==Tree.end()){
      Tree.insert(edge(a1,c1,nTri,-1,b1,-1,p.n,p.gr));
    }else if(founded_edge->neig!=-1) {
      error(3,"Error 3rd neighbour",__LINE__);
    }else {
      founded_edge->neig=nTri;
      founded_edge->neig3rd=b1;
      founded_edge->n=p.n;
      founded_edge->group=p.gr;
    }
    founded_edge=Tree.find(edge(b1,c1,-1,-1,-1,-1,zero,0));
    if(founded_edge==Tree.end()){
      Tree.insert(edge(b1,c1,nTri,-1,a1,-1,p.n,p.gr));
    }else if(founded_edge->neig!=-1) {
      error(3,"Error 3rd neighbour",__LINE__);      
    }else {
      founded_edge->neig=nTri;
      founded_edge->neig3rd=a1;
      founded_edge->n=p.n;
      founded_edge->group=p.gr;
    }
    return 0;
  }
  
  long Stl_io::MakeOFMesh(math_our::Point Mach,const char* fileName)const{
    if(!is_3Dmesh){
      fprintf(stderr,"Can not create OFMesh: stl is not 3d mesh\n");
      return -1;
    }
    OFMesh_arg arg(fileName);
//-------------------Make own/neig array
    set<edge> edgeTree;
    long k=0;
    long nSurfaceElements;
    long buff;
    long nNeighbours=0;
    long nOwners=0;
    long nCells=0;
    long nFaces=0;
    long nInternalFaces=0;
    FILE* pfile;
    FILE* ofile;
    FILE* nfile;
    FILE* ffile;
    FILE* bfile;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	addTriInEdgeTree(k,*j,edgeTree);
	k++;
      }
    }
    nSurfaceElements=k;
    //------------------------------make owner<neighbour--------------------------
    for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
      if(i->neig!=-1){
	if(i->neig<i->owner){
	  buff=i->neig;
	  i->neig=i->owner;
	  i->owner=buff;
	  buff=i->neig3rd;
	  i->neig3rd=i->owner3rd;
	  i->owner3rd=buff;	  
	}
	nNeighbours++;
      }
      nOwners++;
    }
    if(nOwners!=nNeighbours){
      pfile=fopen(POINTFILE,"wb");
      ofile=fopen(OWNERFILE,"wb");
      nfile=fopen(NEIGHFILE,"wb");
      ffile=fopen(FACESFILE,"wb");
      bfile=fopen(BOUNDFILE,"wb");
      //-------------------------------statistic------------------------------------------------------------------------------
      //nPoints:14  nCells:24  nFaces:60  nInternalFaces:36
      nCells=nSurfaceElements*nLayers;
      nInternalFaces=nLayers*nNeighbours+nSurfaceElements*(nLayers-1);
      nFaces=nInternalFaces+2*nSurfaceElements+(nOwners-nNeighbours)*nLayers;
      printf("nPoints:%lu  nCells:%ld  nFaces:%ld  nInternalFaces:%ld\n",(ulong)(pointArray.size()),nCells,nFaces,nInternalFaces);
      //-------------------------------points---------------------------------------------------------------------------------
      fprintf(pfile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n");
      fprintf(pfile,"    class       vectorField;\n    location    \"constant/polyMesh\";\n    object      points;\n}\n//\n");
      fprintf(pfile,"%lu\n(\n",(ulong)(pointArray.size()));
      for(unsigned long i=0;i<pointArray.size();i++)fprintf(pfile,"(%.12le %.12le %.12le)\n",pointArray[i].x,pointArray[i].y,pointArray[i].z);
      fprintf(pfile,")\n\n\n//\n");
      fclose(pfile);
      //-------------------------------owners----------------------------------------------------------------------------------
      fprintf(ofile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       labelList;\n");
      fprintf(ofile,"    note        \"nPoints:%lu  nCells:%ld  nFaces:%ld  nInternalFaces:%ld\";\n",(ulong)(pointArray.size()),nCells,nFaces,nInternalFaces);
      fprintf(ofile,"    location    \"constant/polyMesh\";\n    object      owner;\n}\n//\n\n");
      fprintf(ofile,"%ld\n(\n",nFaces);
      
      //-------------------------------neighbour-------------------------------------------------------------------------------
      fprintf(nfile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       labelList;\n");
      fprintf(nfile,"    note        \"nPoints:%lu  nCells:%ld  nFaces:%ld  nInternalFaces:%ld\";\n",(ulong)(pointArray.size()),nCells,nFaces,nInternalFaces);
      fprintf(nfile,"    location    \"constant/polyMesh\";\n    object      neighbour;\n}\n//\n\n");
      fprintf(nfile,"%ld\n(\n",nInternalFaces);
      //-------------------------------faces-----------------------------------------------------------------------------------
      fprintf(ffile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n");
      fprintf(ffile,"    class       faceList;\n    location    \"constant/polyMesh\";\n    object      faces;\n}\n\n//\n\n");
      fprintf(ffile,"%ld\n(\n",nFaces);
      //-----------------------------------------------------------------------------------------------------------------------
      bool orientation;
      math_our::Point e1,e2,e3,eb;
      for(set<edge>::iterator face=edgeTree.begin();face!=edgeTree.end();++face){
	if(face->neig!=-1){
	  e1=pointArray[face->b]-pointArray[face->a];
	  e2=pointArray[face->a+nPoints]-pointArray[face->a];
	  e3=pointArray[face->neig3rd]-pointArray[face->a];
	  orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	  if(orientation){
	    for(long i=0;i<nLayers;i++){
	      fprintf(ofile,"%ld\n",(face->owner)+i*nSurfaceElements);
	      fprintf(nfile,"%ld\n",(face->neig)+i*nSurfaceElements);
	      fprintf(ffile,"4(%ld %ld %ld %ld)\n",
		      (face->a)+(nPoints*i),face->b+(nPoints*i),face->b+nPoints+(nPoints*i),face->a+nPoints+(nPoints*i));
	    }
	  }else{
	    for(long i=0;i<nLayers;i++){
	      fprintf(ofile,"%ld\n",(face->owner)+i*nSurfaceElements);
	      fprintf(nfile,"%ld\n",(face->neig)+i*nSurfaceElements);
	      fprintf(ffile,"4(%ld %ld %ld %ld)\n",
		      face->b+(nPoints*i),face->a+(nPoints*i),face->a+nPoints+(nPoints*i),face->b+nPoints+(nPoints*i));
	    }
	  }
	}
      }
      //-------------Internal---------------------------------------------------------------
      k=0;
      for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
	for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	  e1=pointArray[j->b]-pointArray[j->a];
	  e2=pointArray[j->c]-pointArray[j->a];
	  e3=pointArray[j->a+nPoints]-pointArray[j->a];
	  orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	  if(orientation){
	    for(long i2=0;i2<(nLayers-1);i2++){
	      fprintf(ofile,"%ld\n",k+i2*nSurfaceElements);
	      fprintf(nfile,"%ld\n",k+(i2+1)*nSurfaceElements);
	      fprintf(ffile,"3(%ld %ld %ld)\n",
		      j->a+(nPoints*(i2+1)),j->b+(nPoints*(i2+1)),j->c+(nPoints*(i2+1)));
	    }
	  }else{
	    for(long i2=0;i2<(nLayers-1);i2++){
	      fprintf(ofile,"%ld\n",k+i2*nSurfaceElements);
	      fprintf(nfile,"%ld\n",k+(i2+1)*nSurfaceElements);
	      fprintf(ffile,"3(%ld %ld %ld)\n",
		      j->a+(nPoints*(i2+1)),j->c+(nPoints*(i2+1)),j->b+(nPoints*(i2+1)));
	    }
	  }
	  k++;
	}
      }
      //---------------Wall------------------------------------------------------
      k=0;
      for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
	for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	  e1=pointArray[j->b]-pointArray[j->a];
	  e2=pointArray[j->c]-pointArray[j->a];
	  e3=pointArray[j->a+nPoints]-pointArray[j->a];
	  orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	  if(!orientation){
	    fprintf(ofile,"%ld\n",k);
	    fprintf(ffile,"3(%ld %ld %ld)\n",
		      j->a,j->b,j->c);
	  }else{
	    fprintf(ofile,"%ld\n",k);
	    fprintf(ffile,"3(%ld %ld %ld)\n",
		      j->a,j->c,j->b);
	  }
	  k++;
	}
      }
      //---------------Outer-------------------------------------------------
      k=0;
      long inlet_counter=0;
      long *inlet_info;
      inlet_info=(long*)malloc(sizeof(long)*4*nSurfaceElements);
      long owner;
      long p_1,p_2,p_3;
      for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
	for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	  e1=pointArray[j->b]-pointArray[j->a];
	  e2=pointArray[j->c]-pointArray[j->a];
	  e3=pointArray[j->a+nPoints]-pointArray[j->a];
	  orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	  if(orientation){
	    //
	    owner=k+(nLayers-1)*nSurfaceElements;
	    p_1=j->a+nLayers*nPoints;
	    p_2=j->b+nLayers*nPoints;
	    p_3=j->c+nLayers*nPoints;
	  }else{
	    //
	    owner=k+(nLayers-1)*nSurfaceElements;
	    p_1=j->a+nLayers*nPoints;
	    p_2=j->c+nLayers*nPoints;
	    p_3=j->b+nLayers*nPoints;
	  }
	  e1=pointArray[p_2]-pointArray[p_1];
	  e2=pointArray[p_3]-pointArray[p_1];
	  e3=(e1^e2);
	  e3=e3.normalize();
	  inlet_info[inlet_counter*4  ]=owner;
	  inlet_info[inlet_counter*4+1]=p_1;
	  inlet_info[inlet_counter*4+2]=p_2;
	  inlet_info[inlet_counter*4+3]=p_3;
	  inlet_counter++;	  
	  k++;
	}
      }
      for(long i=0;i<inlet_counter;i++){
	fprintf(ofile,"%ld\n",inlet_info[i*4]);
	fprintf(ffile,"3(%ld %ld %ld)\n",
		      inlet_info[i*4+1],inlet_info[i*4+2],inlet_info[i*4+3]);
      }
      //-------------------------------Outlet---------------------------------------------------------------------------------
      for(set<edge>::iterator face=edgeTree.begin();face!=edgeTree.end();++face){
	if(face->neig==-1){
	  e1=pointArray[face->b]-pointArray[face->a];
	  e2=pointArray[face->a+nPoints]-pointArray[face->a];
	  e3=pointArray[face->owner3rd]-pointArray[face->a];
	  orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	  if(!orientation){
	    for(long i=0;i<nLayers;i++){
	      fprintf(ofile,"%ld\n",(face->owner)+i*nSurfaceElements);
// 	      fprintf(nfile,"%ld\n",(face->neig)+i*nSurfaceElements);
	      fprintf(ffile,"4(%ld %ld %ld %ld)\n",
		      (face->a)+(nPoints*i),face->b+(nPoints*i),face->b+nPoints+(nPoints*i),face->a+nPoints+(nPoints*i));
	    }
	  }else{
	    for(long i=0;i<nLayers;i++){
	      fprintf(ofile,"%ld\n",(face->owner)+i*nSurfaceElements);
// 	      fprintf(nfile,"%ld\n",(face->neig)+i*nSurfaceElements);
	      fprintf(ffile,"4(%ld %ld %ld %ld)\n",
		      face->b+(nPoints*i),face->a+(nPoints*i),face->a+nPoints+(nPoints*i),face->b+nPoints+(nPoints*i));
	    }
	  }
	}
      }
      //-------------------------------boundary-------------------------------------------------------------------------------
      fprintf(bfile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n");
      fprintf(bfile,"    class       polyBoundaryMesh;\n    location    \"constant/polyMesh\";\n    object      boundary;\n}\n//\n");
      fprintf(bfile,"%ld\n(\n",arg.def?(groupArray.size()+2):3);
      long boundary_base=nInternalFaces;
      char str[300];
      unsigned long j;
      if(arg.def){
	for(vector<Group>::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
	  for(j=0;j<i->name.size() && i->name[j]!=' ';j++)str[j]=i->name[j];
	  str[j]='\0';
	    fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",str);
	    fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
		    i->tri.size(),boundary_base);
	    boundary_base+=i->tri.size();
	}
      }else{
	fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",arg.wall);
	fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
		nSurfaceElements,boundary_base);
	  boundary_base+=nSurfaceElements;
      }
      fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",arg.inlet);
      fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
	      nSurfaceElements,boundary_base);
      boundary_base+=nSurfaceElements;
      fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",arg.outlet);
      fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
	      (nOwners-nNeighbours)*nLayers,boundary_base);
      //-------------------------------close files----------------------------------------------------------------------------
      fprintf(ofile,")\n\n//\n");
      fclose(ofile);
      fprintf(nfile,")\n\n//\n");
      fclose(nfile);
      fprintf(ffile,")\n\n//\n");
      fclose(ffile);
      fprintf(bfile,")\n\n//\n");
      fclose(bfile);
      free(inlet_info);
      return -1;
    }
    pfile=fopen(POINTFILE,"wb");
    ofile=fopen(OWNERFILE,"wb");
    nfile=fopen(NEIGHFILE,"wb");
    ffile=fopen(FACESFILE,"wb");
    bfile=fopen(BOUNDFILE,"wb");
    //-------------------------------statistic------------------------------------------------------------------------------
    //nPoints:14  nCells:24  nFaces:60  nInternalFaces:36
    nCells=nSurfaceElements*nLayers;
    nInternalFaces=nLayers*nOwners+nSurfaceElements*(nLayers-1);
    nFaces=nInternalFaces+2*nSurfaceElements;
    printf("nPoints:%lu  nCells:%ld  nFaces:%ld  nInternalFaces:%ld\n",(ulong)(pointArray.size()),nCells,nFaces,nInternalFaces);
    //-------------------------------points---------------------------------------------------------------------------------
    fprintf(pfile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n");
    fprintf(pfile,"    class       vectorField;\n    location    \"constant/polyMesh\";\n    object      points;\n}\n//\n");
    fprintf(pfile,"%lu\n(\n",(ulong)(pointArray.size()));
    for(unsigned long i=0;i<pointArray.size();i++)fprintf(pfile,"(%.12le %.12le %.12le)\n",pointArray[i].x,pointArray[i].y,pointArray[i].z);
    fprintf(pfile,")\n\n\n//\n");
    fclose(pfile);
    //-------------------------------owners----------------------------------------------------------------------------------
    fprintf(ofile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       labelList;\n");
    fprintf(ofile,"    note        \"nPoints:%lu  nCells:%ld  nFaces:%ld  nInternalFaces:%ld\";\n",(ulong)(pointArray.size()),nCells,nFaces,nInternalFaces);
    fprintf(ofile,"    location    \"constant/polyMesh\";\n    object      owner;\n}\n//\n\n");
    fprintf(ofile,"%ld\n(\n",nFaces);
    
    //-------------------------------neighbour-------------------------------------------------------------------------------
    fprintf(nfile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       labelList;\n");
    fprintf(nfile,"    note        \"nPoints:%lu  nCells:%ld  nFaces:%ld  nInternalFaces:%ld\";\n",(ulong)(pointArray.size()),nCells,nFaces,nInternalFaces);
    fprintf(nfile,"    location    \"constant/polyMesh\";\n    object      neighbour;\n}\n//\n\n");
    fprintf(nfile,"%ld\n(\n",nInternalFaces);
    //-------------------------------faces-----------------------------------------------------------------------------------
    fprintf(ffile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n");
    fprintf(ffile,"    class       faceList;\n    location    \"constant/polyMesh\";\n    object      faces;\n}\n\n//\n\n");
    fprintf(ffile,"%ld\n(\n",nFaces);
    //-----------------------------------------------------------------------------------------------------------------------
    bool orientation;
    math_our::Point e1,e2,e3,eb;
    for(set<edge>::iterator face=edgeTree.begin();face!=edgeTree.end();++face){
      e1=pointArray[face->b]-pointArray[face->a];
      e2=pointArray[face->a+nPoints]-pointArray[face->a];
      e3=pointArray[face->neig3rd]-pointArray[face->a];
      orientation=(mix_mul(e1,e2,e3)>0)?true:false;
      if(orientation){
	for(long i=0;i<nLayers;i++){
	  fprintf(ofile,"%ld\n",(face->owner)+i*nSurfaceElements);
	  fprintf(nfile,"%ld\n",(face->neig)+i*nSurfaceElements);
	  fprintf(ffile,"4(%ld %ld %ld %ld)\n",
		  (face->a)+(nPoints*i),face->b+(nPoints*i),face->b+nPoints+(nPoints*i),face->a+nPoints+(nPoints*i));
	}
      }else{
	for(long i=0;i<nLayers;i++){
	  fprintf(ofile,"%ld\n",(face->owner)+i*nSurfaceElements);
	  fprintf(nfile,"%ld\n",(face->neig)+i*nSurfaceElements);
	  fprintf(ffile,"4(%ld %ld %ld %ld)\n",
		  face->b+(nPoints*i),face->a+(nPoints*i),face->a+nPoints+(nPoints*i),face->b+nPoints+(nPoints*i));
	}
      }
    }
    //-------------Internal---------------------------------------------------------------
    k=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	e1=pointArray[j->b]-pointArray[j->a];
	e2=pointArray[j->c]-pointArray[j->a];
	e3=pointArray[j->a+nPoints]-pointArray[j->a];
	orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	if(orientation){
	  for(long i2=0;i2<(nLayers-1);i2++){
	    fprintf(ofile,"%ld\n",k+i2*nSurfaceElements);
	    fprintf(nfile,"%ld\n",k+(i2+1)*nSurfaceElements);
	    fprintf(ffile,"3(%ld %ld %ld)\n",
		    j->a+(nPoints*(i2+1)),j->b+(nPoints*(i2+1)),j->c+(nPoints*(i2+1)));
	  }
	}else{
	  for(long i2=0;i2<(nLayers-1);i2++){
	    fprintf(ofile,"%ld\n",k+i2*nSurfaceElements);
	    fprintf(nfile,"%ld\n",k+(i2+1)*nSurfaceElements);
	    fprintf(ffile,"3(%ld %ld %ld)\n",
		    j->a+(nPoints*(i2+1)),j->c+(nPoints*(i2+1)),j->b+(nPoints*(i2+1)));
	  }
	}
	k++;
      }
    }
    //---------------Wall------------------------------------------------------
    k=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	e1=pointArray[j->b]-pointArray[j->a];
	e2=pointArray[j->c]-pointArray[j->a];
	e3=pointArray[j->a+nPoints]-pointArray[j->a];
	orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	if(!orientation){
	  fprintf(ofile,"%ld\n",k);
	  fprintf(ffile,"3(%ld %ld %ld)\n",
		    j->a,j->b,j->c);
	}else{
	  fprintf(ofile,"%ld\n",k);
	  fprintf(ffile,"3(%ld %ld %ld)\n",
		    j->a,j->c,j->b);
	}
	k++;
      }
    }
    //---------------Outer-------------------------------------------------
    k=0;
    long inlet_counter=0;
    long outlet_counter=0;
    long *inlet_info;
    inlet_info=(long*)malloc(sizeof(long)*4*nSurfaceElements);
    long *outlet_info;
    outlet_info=(long*)malloc(sizeof(long)*4*nSurfaceElements);
    long owner;
    long p_1,p_2,p_3;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	e1=pointArray[j->b]-pointArray[j->a];
	e2=pointArray[j->c]-pointArray[j->a];
	e3=pointArray[j->a+nPoints]-pointArray[j->a];
	orientation=(mix_mul(e1,e2,e3)>0)?true:false;
	if(orientation){
	  //
	  owner=k+(nLayers-1)*nSurfaceElements;
	  p_1=j->a+nLayers*nPoints;
	  p_2=j->b+nLayers*nPoints;
	  p_3=j->c+nLayers*nPoints;
	  //
// 	  fprintf(ofile,"%ld\n",k+(nLayers-1)*nSurfaceElements);
// 	  fprintf(ffile,"3(%ld %ld %ld)\n",
// 		    j->a+nLayers*nPoints,j->b+nLayers*nPoints,j->c+nLayers*nPoints);
	}else{
	  //
	  owner=k+(nLayers-1)*nSurfaceElements;
	  p_1=j->a+nLayers*nPoints;
	  p_2=j->c+nLayers*nPoints;
	  p_3=j->b+nLayers*nPoints;
	  //
// 	  fprintf(ofile,"%ld\n",k+(nLayers-1)*nSurfaceElements);
// 	  fprintf(ffile,"3(%ld %ld %ld)\n",
// 		    j->a+nLayers*nPoints,j->c+nLayers*nPoints,j->b+nLayers*nPoints);
	}
	e1=pointArray[p_2]-pointArray[p_1];
	e2=pointArray[p_3]-pointArray[p_1];
	e3=(e1^e2);
	e3=e3.normalize();
	if ((Mach*e3)>1.1){
	  outlet_info[outlet_counter*4  ]=owner;
	  outlet_info[outlet_counter*4+1]=p_1;
	  outlet_info[outlet_counter*4+2]=p_2;
	  outlet_info[outlet_counter*4+3]=p_3;
	  outlet_counter++;
	}else{
	  inlet_info[inlet_counter*4  ]=owner;
	  inlet_info[inlet_counter*4+1]=p_1;
	  inlet_info[inlet_counter*4+2]=p_2;
	  inlet_info[inlet_counter*4+3]=p_3;
	  inlet_counter++;
	}
	k++;
      }
    }
    for(long i=0;i<inlet_counter;i++){
      fprintf(ofile,"%ld\n",inlet_info[i*4]);
      fprintf(ffile,"3(%ld %ld %ld)\n",
		    inlet_info[i*4+1],inlet_info[i*4+2],inlet_info[i*4+3]);
    }
    for(long i=0;i<outlet_counter;i++){
      fprintf(ofile,"%ld\n",outlet_info[i*4]);
      fprintf(ffile,"3(%ld %ld %ld)\n",
		    outlet_info[i*4+1],outlet_info[i*4+2],outlet_info[i*4+3]);
    }
    //-------------------------------boundary-------------------------------------------------------------------------------
    fprintf(bfile,"//\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n");
    fprintf(bfile,"    class       polyBoundaryMesh;\n    location    \"constant/polyMesh\";\n    object      boundary;\n}\n//\n");
    long nbc=2;
    if(!strcmp(arg.inlet,arg.outlet)){
      nbc=1;
    }
    fprintf(bfile,"%ld\n(\n",arg.def?(groupArray.size()+nbc):(nbc+1));
    long boundary_base=nInternalFaces;
    char str[300];
    unsigned long j;
    if(arg.def){
      for(vector<Group>::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
	for(j=0;j<i->name.size();j++)str[j]=i->name[j];
	str[j]='\0';
	  fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",str);
	  fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
		  i->tri.size(),boundary_base);
	  boundary_base+=i->tri.size();
      }
    }else{
      fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",arg.wall);
      fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
	      nSurfaceElements,boundary_base);
	boundary_base+=nSurfaceElements;
    }
    if(strcmp(arg.inlet,arg.outlet)){
      fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",arg.inlet);
      fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
	      inlet_counter,boundary_base);
      boundary_base+=inlet_counter;
      fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",arg.outlet);
      fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
	      outlet_counter,boundary_base);
    }else{
      fprintf(bfile,"    %s\n    {\n        type            patch;\n        physicalType    patch;\n",arg.inlet);
      fprintf(bfile,"        nFaces          %ld;\n        startFace       %ld;\n    }\n",
	      inlet_counter+outlet_counter,boundary_base);
      //boundary_base+=inlet_counter;
      //boundary_base+=outlet_counter;
    }
    //-------------------------------close files----------------------------------------------------------------------------
    fprintf(ofile,")\n\n//\n");
    fclose(ofile);
    fprintf(nfile,")\n\n//\n");
    fclose(nfile);
    fprintf(ffile,")\n\n//\n");
    fclose(ffile);
    fprintf(bfile,")\n\n//\n");
    fclose(bfile);
    free(inlet_info);
    free(outlet_info);
    //----------------------------------------------------------------------------------------------------------------------
//-------------------
    return 0;
  }
}
