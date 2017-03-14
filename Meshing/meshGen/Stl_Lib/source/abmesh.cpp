#include "stl_io.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include <math.h>
// PointArray >= 4
#define SSD_1 4
// NewPointArray >= 4
#define SSD_2 4
// ChargeArray >= 4
#define SSD_3 4

#define MAX_N 100
namespace stl
{
  using namespace std;
  void Stl_io::abmesh_cu(const char* fileName){
    Stl_io phantom;
    FILE* file;
    if(is_3Dmesh){
      printf("aborted: it is 3D mesh\n");
      return;
    }
    const char lfile[]="meshGen/length_abl";
    const char pfile[]="meshGen/profile_abl";
    const char sfile[]="meshGen/settings_abl";
    
    abmesh_arg arg_abl(fileName);
    MMLE_arg arg_mmle(fileName);
    CUT_arg arg_cut(fileName);
    
    sprintf(arg_mmle.profileFile,"%s",pfile);
    sprintf(arg_mmle.lengthFile,"%s",lfile);
    
    file=fopen(sfile,"wb");
    if(!file){
      error(1,sfile,__LINE__);
      return;
    }
    arg_mmle.addInFile(file);
    arg_cut.addInFile(file);
    fflush(file);
    fclose(file);
      
    if(arg_abl.n_u){
      file=fopen(pfile,"wb");
      if(!file){
	error(1,pfile,__LINE__);
	return;
      }
      fprintf(file,"1 %d\n",arg_abl.n_u);
      fflush(file);
      fclose(file);
      file=fopen(lfile,"wb");
      if(!file){
	error(1,lfile,__LINE__);
	return;
      }
      fprintf(file,"%lf %ld\n",arg_abl.l_u,nPoints);
      fflush(file);
      fclose(file);
      
      phantom=extract_layer(0);
      MMPoi(sfile,phantom);
      reverse_layers();
    }
    for(int i=0;i<arg_abl.n_d;++i){      
      file=fopen(pfile,"wb");
      if(!file){
	error(1,pfile,__LINE__);
	return;
      }
      fprintf(file,"1\n");
      fflush(file);
      fclose(file);
      file=fopen(lfile,"wb");
      if(!file){
	error(1,lfile,__LINE__);
	return;
      }
      fprintf(file,"%lf %ld\n",-arg_abl.l_d/(double)arg_abl.n_d,nPoints);
      fflush(file);
      fclose(file);
      phantom=extract_layer(nLayers);
      MMPoi(sfile,phantom);
    }
  }
  void Stl_io::ab_field(const char* fileName)const{
    abmesh_arg arg_abl(fileName);
    FILE *file,*tmpFile;
    file=fopen(arg_abl.inFile,"rb");
    if(!file){
      printf("error open %s\n",arg_abl.inFile);
      return;
    }
    long nElements=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	nElements++;
      }
    }
    vector<pair<char*,double*> >fields;
    char string[300];
    char* tmpString;
    fgets(string,300,file);
    tmpFile=fopen("tmpFile0112358","wb");
    fprintf(tmpFile,"%s",string);
    fflush(tmpFile);
    fclose(tmpFile);
    tmpFile=fopen("tmpFile0112358","rb");
    long tmpLong=0;
    double *tmpP;
    do{
      tmpString=(char*)malloc(sizeof(char)*30);
      tmpLong=fscanf(tmpFile," %30s",tmpString);
      if(tmpLong==1){
	tmpP=(double*)malloc(sizeof(double)*nElements);
	fields.push_back(pair<char*,double*>(tmpString,tmpP));
      }
    }while(tmpLong==1);
    fclose(tmpFile);
    system("rm tmpFile0112358");
    tmpLong=0;
    for(long i=0;i<nElements;i++){
      for(ulong j=0;j<fields.size();j++){
	tmpLong=fscanf(file," %20lf",&(fields[j].second[i]));
	if(tmpLong!=1){
	  printf("error fields file has %ld of %ld strings\n",i,nElements);
	  fclose(file);
	  return;
	}
      }
    }
    fclose(file);
    system("mkdir -p 0");
    for(ulong j=0;j<fields.size();j++){
//      if(!strcmp(fields[j].first,"SU")){
      tmpString[0]='0';
      tmpString[1]='/';
      tmpString[2]='\0';
      strcat(tmpString,fields[j].first);
      file=fopen(tmpString,"wb");
      if(!file){
	printf("error open %s\n",fields[j].first);
	return;
      }
      fprintf(file,"//\nFoamFile \n\
{\n\
    version     2.0;\n\
    format      ascii;\n\
    class       volScalarField;\n\
    location    \"0\";\n\
    object      %s;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\
\n\
dimensions      [0 0 0 0 0 0 0];\n\
\n\
internalField   uniform 0; \n\
boundaryField \n\
{ \n\
    border\n\
    {\n\
    type          fixedValue;\n\
    value         nonuniform List<scalar>\n\
%ld\n\
(\n",fields[j].first,nElements);
	for(long i=0;i<nElements;i++){
	  fprintf(file,"%lf\n",fields[j].second[i]);
	}
    fprintf(file,");\n\
    }\n\
    in\n\
    {\n\
        type            zeroGradient;\n\
    }\n\
    outlet\n\
    {\n\
        type            zeroGradient;\n\
    }\n\
}\n");
      fclose(file);
//      }
    }
  }
}