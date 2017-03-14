#include "stl_io.h"
#include<string.h>
namespace stl
{
  using namespace std;
  void Settings::Write(const char* fileName){
    FILE *file;
    file=fopen(fileName,"wb");
    if(!file){
      printf("could not write settings to \"%s\"\n",fileName);
      return;
    }
    pch.addInFile(file);
    mml.addInFile(file);
    cut.addInFile(file);
    ofm.addInFile(file);
    ada.addInFile(file);
    abm.addInFile(file);
    fclose(file);
  }
  Settings::Settings(const char* fileName){
    pch=PChange_arg(fileName);
    mml=MMLE_arg(fileName);
    cut=CUT_arg(fileName);
    ofm=OFMesh_arg(fileName);
    ada=findWave_arg(fileName);
    abm=abmesh_arg(fileName);
  }
  void Settings::modi_par(const char* parameter){
    char sec[30];
    char par[30];
    char val[30];
    char type[30];
    if(sscanf(parameter,"%30s %30s %30s %30s",sec,par,val,type)!=4){
      printf("bad parameter");
    }
    if(!strcmp(sec,"pch")){
      if(!strcmp(par,"len")){
        if(val[0]=='~'){
          pch.lengthFile[0]='\0';
        }else{
          strcpy(pch.lengthFile,val);
        }
      }else if(!strcmp(par,"pro")){
        if(val[0]=='~'){
          pch.profileFile[0]='\0';
        }else{
          strcpy(pch.profileFile,val);
        }
      }else if(!strcmp(par,"cre")){
        
      }
    }else if(!strcmp(sec,"mml")){
    }else if(!strcmp(sec,"cut")){
    }else if(!strcmp(sec,"ofm")){
    }else if(!strcmp(sec,"ada")){
    }else if(!strcmp(sec,"abm")){
    }
  }
  MMLE_arg::MMLE_arg():	l(1.),l1(0.),l2(0.),la(0.),
			lBL(0.1),nBL(5),mBL(1.2),
			lSW(0.1),nSW(5),mSW(1.2),
			pow_q(1.),pow_E(1.),
			nStep(0),elas(0.01),elas2(2),nCudaIter(20),
			createProfile(0),firstLayer(0){
    lengthFile=(char*)malloc(300*sizeof(char));
    profileFile=(char*)malloc(300*sizeof(char)); 
    profileFile[0]='a';
    lengthFile[0]='a';
    lengthFile[1]='\0';
    profileFile[1]='\0';
  }
  MMLE_arg::MMLE_arg(const MMLE_arg& a){
    l=a.l;
    l1=a.l1;
    l2=a.l2;
    la=a.la;
    lBL=a.lBL;
    nBL=a.nBL;
    mBL=a.mBL;
    lSW=a.lSW;
    nSW=a.nSW;
    mSW=a.mSW;
    pow_q=a.pow_q;
    pow_E=a.pow_E;
    nStep=a.nStep;
    elas=a.elas;
    elas2=a.elas2;
    nCudaIter=a.nCudaIter;
    createProfile=a.createProfile;
    firstLayer=a.firstLayer;
    lengthFile=(char*)malloc(300*sizeof(char));
    profileFile=(char*)malloc(300*sizeof(char));
    for(long i=0;i<300;i++){
      profileFile[i]=a.profileFile[i];
      lengthFile[i]=a.lengthFile[i];
    }
  }
  MMLE_arg& MMLE_arg::operator=(const MMLE_arg& a){
    l=a.l;
    l1=a.l1;
    l2=a.l2;
    la=a.la;
    lBL=a.lBL;
    nBL=a.nBL;
    mBL=a.mBL;
    lSW=a.lSW;
    nSW=a.nSW;
    mSW=a.mSW;
    pow_q=a.pow_q;
    pow_E=a.pow_E;
    nStep=a.nStep;
    elas=a.elas;
    elas2=a.elas2;
    nCudaIter=a.nCudaIter;
    createProfile=a.createProfile;
    firstLayer=a.firstLayer;
    for(long i=0;i<300;i++){
      profileFile[i]=a.profileFile[i];
      lengthFile[i]=a.lengthFile[i];
    }
    return *this;
  }
  CUT_arg::CUT_arg():n(math_our::Point(1,0,0)),D(-100),cs(1){
  }
  MMLE_arg::MMLE_arg(const char* fileName){
    FILE* file;
    char str[300];
    char str_a[300];
    char* q=(char*)1;
    file=fopen(fileName,"rb");
    strcpy(str_a,"");
    lengthFile=(char*)malloc(300*sizeof(char));
    profileFile=(char*)malloc(300*sizeof(char));
    profileFile[0]='\0';
    lengthFile[0]='\0';
    do{
      q=fgets(str,300,file);
      sscanf(str,"%200s",str_a);
    }while(strcmp(str_a,"MMPE_options") && q);
    if(q){
      fgets(str,300,file);
      sscanf(str," %*s %20lf %20lf %20lf %20lf",&l,&l1,&l2,&la);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&lBL);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&nBL);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&mBL);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&lSW);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&nSW);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&mSW);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&pow_q);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&pow_E);
      fgets(str,300,file);
      sscanf(str," %*s %20ld %20lf %20lf",&nStep,&elas,&elas2);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&nCudaIter);
      fgets(str,300,file);
      sscanf(str," %*s %300s",lengthFile);
      fgets(str,300,file);
      sscanf(str," %*s %300s",profileFile);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&createProfile);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&firstLayer);
    }
    fclose(file);
  }
  void MMLE_arg::addInFile(FILE* file){
    fprintf(file,"MMPE_options\n");
    fprintf(file,"            l %20lf %20lf %20lf %20lf\n",l,l1,l2,la);
    fprintf(file,"          lBL %20lf\n",lBL);
    fprintf(file,"          nBL %20ld\n",nBL);
    fprintf(file,"          mBL %20lf\n",mBL);
    fprintf(file,"          lSW %20lf\n",lSW);
    fprintf(file,"          nSW %20ld\n",nSW);
    fprintf(file,"          mSW %20lf\n",mSW);
    fprintf(file,"        pow_q %20lf\n",pow_q);
    fprintf(file,"        pow_E %20lf\n",pow_E);
    fprintf(file,"        nStep %20ld %20lf %20lf\n",nStep,elas,elas2);
    fprintf(file,"      nLayers %20ld\n",nCudaIter);
    fprintf(file,"   lengthFile %20s%s\n","",lengthFile);
    fprintf(file,"  profileFile %20s%s\n","",profileFile);
    fprintf(file,"createProfile %20ld\n",createProfile);
    fprintf(file,"   firstLayer %20ld\n",firstLayer);
  }
  MMLE_arg::~MMLE_arg(){
    free(lengthFile);
    free(profileFile);
  }
  void CUT_arg::addInFile(FILE* file){
    math_our::Point a=n/n.module();
    a=a*(-D);
    fprintf(file,"cut_options\n");
    fprintf(file,"       normal %20lf %20lf %20lf\n",n.x,n.y,n.z);
    fprintf(file," pnt_on_plane %20lf %20lf %20lf\n",a.x,a.y,a.z);
    fprintf(file,"           cs %20lf\n",cs);
  }
  abmesh_arg::abmesh_arg(const char* fileName){
    FILE* file;
    char str[300];
    char str_a[300];
    char* q=(char*)1;
    inFile=(char*)malloc(300*sizeof(char));
    inFile[0]='\0';
    file=fopen(fileName,"rb");
    strcpy(str_a,"");
    do{
      q=fgets(str,300,file);
      sscanf(str,"%200s",str_a);
    }while(strcmp(str_a,"abmesh_options") && q);
    if(q){
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&l_u);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&l_d);
      fgets(str,300,file);
      sscanf(str," %*s %20d",&n_u);
      fgets(str,300,file);
      sscanf(str," %*s %20d",&n_d);
      fgets(str,300,file);
      sscanf(str," %*s %200s",inFile);
    }
    fclose(file);
  }

  void abmesh_arg::addInFile(FILE* file)
  {
      fprintf(file,"abmesh_options\n");
      fprintf(file," l_u %20lf\n",l_u);
      fprintf(file," l_d %20lf\n",l_d);
      fprintf(file," n_u %20d\n",n_u);
      fprintf(file," n_d %20d\n",n_d);
      fprintf(file,"file %20s%s\n","",inFile);
  }

  abmesh_arg::abmesh_arg():l_u(0.),l_d(0.),n_u(0),n_d(0){
    inFile=(char*)malloc(300*sizeof(char));
    inFile[0]='\0';
  }
  abmesh_arg::~abmesh_arg(){
    free(inFile);
  }
  abmesh_arg::abmesh_arg(const abmesh_arg& a){
    l_u=a.l_u;
    l_d=a.l_d;
    n_u=a.n_u;
    n_d=a.n_d;
    inFile=(char*)malloc(300*sizeof(char));
    for(long i=0;i<300;i++){
      inFile[i]=a.inFile[i];
    }
  }
  abmesh_arg& abmesh_arg::operator=(const abmesh_arg& a){
    l_u=a.l_u;
    l_d=a.l_d;
    n_u=a.n_u;
    n_d=a.n_d;
    for(long i=0;i<300;i++){
      inFile[i]=a.inFile[i];
    }
    return *this;
  }
  CUT_arg::CUT_arg(const char* fileName){
    FILE* file;
    char str[300];
    char str_a[300];
    char* q=(char*) 1;
    math_our::Point r0,n_normal;
    file=fopen(fileName,"rb");
    strcpy(str_a,"");
    do{
      q=fgets(str,300,file);
      sscanf(str,"%40s",str_a);
    }while(strcmp(str_a,"cut_options") && q);
    if(q){
      fgets(str,300,file);
      sscanf(str," %*s %20lf %20lf %20lf",&n.x,&n.y,&n.z);
      fgets(str,300,file);
      sscanf(str," %*s %20lf %20lf %20lf",&r0.x,&r0.y,&r0.z);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&cs);
      double moduleN;
      moduleN=n.module();
      n_normal.x=n.x/moduleN;
      n_normal.y=n.y/moduleN;
      n_normal.z=n.z/moduleN;
      D=-(n_normal*r0);
    }
    fclose(file);
  }
  findWave_arg::findWave_arg(const char* fileName){
    FILE* file;
    char str[300];
    char str_a[300];
    char* q=(char*) 1;
    file=fopen(fileName,"rb");
    strcpy(str_a,"");
    fieldFile=(char*)malloc(300*sizeof(char));
    outputFile=(char*)malloc(300*sizeof(char));
    do{
      q=fgets(str,300,file);
      sscanf(str,"%40s",str_a);
    }while(strcmp(str_a,"findWave_options") && q);
    if(q){
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&mul);
      fgets(str,300,file);
      sscanf(str," %*s %20s",fieldFile);
      fgets(str,300,file);
      sscanf(str," %*s %20s",outputFile);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&eps);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&made);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&TOL);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&nRel);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&elas);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&minL);
      fgets(str,300,file);
      sscanf(str," %*s %20lf",&iterLimit);
    }
    fclose(file);
  }
  void findWave_arg::addInFile(FILE* file)
  {
      fprintf(file,"findWave_options\n");
      fprintf(file,"  extraCellsMul %20lf\n",mul);
      fprintf(file,"  fieldFileName %20s%s\n","",fieldFile);
      fprintf(file," outputFileName %20s%s\n","",outputFile);
      fprintf(file,"eps_field/inlet %20lf\n",eps);
      fprintf(file,"           made %20ld\n",made);
      fprintf(file,"            TOL %20lf\n",TOL);
      fprintf(file,"           nRel %20ld\n",nRel);
      fprintf(file,"           elas %20lf\n",elas);
      fprintf(file,"      maxLength %20lf\n",minL);
      fprintf(file,"      iterLimit %20lf\n",iterLimit);
  }
  findWave_arg::~findWave_arg(){
    free(fieldFile);
    free(outputFile);
  }
  findWave_arg::findWave_arg():mul(0.),eps(0.),made(0),
                                TOL(0.),nRel(0),elas(0.),
                                minL(0.),iterLimit(0.){
    fieldFile=(char*)malloc(300*sizeof(char));
    outputFile=(char*)malloc(300*sizeof(char));
    fieldFile[0]='\0';
    outputFile[0]='\0';
  }
  findWave_arg::findWave_arg(const findWave_arg& a){
    mul=a.mul;
    eps=a.eps;
    made=a.made;
    TOL=a.TOL;
    nRel=a.nRel;
    elas=a.elas;
    minL=a.minL;
    iterLimit=a.iterLimit;
    fieldFile=(char*)malloc(300*sizeof(char));
    outputFile=(char*)malloc(300*sizeof(char));
    for(long i=0;i<300;i++){
      fieldFile[i]=a.fieldFile[i];
      outputFile[i]=a.outputFile[i];
    }
  }
  findWave_arg& findWave_arg::operator=(const findWave_arg& a){
    mul=a.mul;
    eps=a.eps;
    made=a.made;
    TOL=a.TOL;
    nRel=a.nRel;
    elas=a.elas;
    minL=a.minL;
    iterLimit=a.iterLimit;
    for(long i=0;i<300;i++){
      fieldFile[i]=a.fieldFile[i];
      outputFile[i]=a.outputFile[i];
    }
    return *this;
  }
  OFMesh_arg::OFMesh_arg(const char* fileName){
    FILE* file;
    char str[300];
    char str_a[300];
    char* q=(char*)1;
    file=fopen(fileName,"rb");
    strcpy(str_a,"");
    wall=(char*)malloc(300*sizeof(char));
    inlet=(char*)malloc(300*sizeof(char));
    outlet=(char*)malloc(300*sizeof(char));
    do{
      q=fgets(str,300,file);
      sscanf(str,"%200s",str_a);
    }while(strcmp(str_a,"OFMesh_options") && q);
    if(q){
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&def);
      fgets(str,300,file);
      sscanf(str," %*s %300s",wall);
      fgets(str,300,file);
      sscanf(str," %*s %300s",inlet);
      fgets(str,300,file);
      sscanf(str," %*s %300s",outlet);
    }
    fclose(file);
  }
  void OFMesh_arg::addInFile(FILE* file)
  {
      fprintf(file,"OFMesh_options\n");
      fprintf(file,"default %20ld\n",def);
      fprintf(file,"   wall %20s%s\n","",wall);
      fprintf(file,"  inlet %20s%s\n","",inlet);
      fprintf(file," outlet %20s%s\n","",outlet);
  }
  OFMesh_arg::~OFMesh_arg(){
    free(wall);
    free(inlet);
    free(outlet);
  }
  OFMesh_arg::OFMesh_arg():def(0){
    wall=(char*)malloc(300*sizeof(char));
    inlet=(char*)malloc(300*sizeof(char));
    outlet=(char*)malloc(300*sizeof(char));
  }
  OFMesh_arg::OFMesh_arg(const OFMesh_arg& a){
    def=a.def;
    wall=(char*)malloc(300*sizeof(char));
    inlet=(char*)malloc(300*sizeof(char));
    outlet=(char*)malloc(300*sizeof(char));
    for(long i=0;i<300;i++){
      wall[i]=a.wall[i];
      inlet[i]=a.inlet[i];
      outlet[i]=a.outlet[i];
    }
  }
  OFMesh_arg& OFMesh_arg::operator=(const OFMesh_arg& a){
    def=a.def;
    for(long i=0;i<300;i++){
      wall[i]=a.wall[i];
      inlet[i]=a.inlet[i];
      outlet[i]=a.outlet[i];
    }
    return *this;
  }
  PChange_arg::PChange_arg(const char* fileName){
    FILE* file;
    char str[300];
    char str_a[300];
    char* q=(char*)1;
    file=fopen(fileName,"rb");
    strcpy(str_a,"");
    lengthFile=(char*)malloc(300*sizeof(char));
    profileFile=(char*)malloc(300*sizeof(char));
    do{
      q=fgets(str,300,file);
      sscanf(str,"%40s",str_a);
    }while(strcmp(str_a,"changeProfile_options") && q);
    if(q){
      fgets(str,300,file);
      sscanf(str," %*s %300s",lengthFile);
      fgets(str,300,file);
      sscanf(str," %*s %20ld",&nStep);
      fgets(str,300,file);
      sscanf(str," %*s %300s",profileFile);
      fgets(str,300,file);
      sscanf(str," %*s %10ld",&createProfile);
    }
    fclose(file);
  }
  void PChange_arg::addInFile(FILE* file)
  {
      fprintf(file,"changeProfile_options\n");
      fprintf(file,"   lengthFile %20s%s\n","",lengthFile);
      fprintf(file,"   nLengthRel %20ld\n",nStep);
      fprintf(file,"  profileFile %20s%s\n","",profileFile);
      fprintf(file,"createProfile %20ld\n",createProfile);
  }
  PChange_arg::~PChange_arg(){
    free(lengthFile);
    free(profileFile);    
  }
  PChange_arg::PChange_arg():nStep(0),createProfile(0){
    lengthFile=(char*)malloc(300*sizeof(char));
    profileFile=(char*)malloc(300*sizeof(char));
  }
  PChange_arg::PChange_arg(const PChange_arg& a){
    nStep=a.nStep;
    createProfile=a.createProfile;
    lengthFile=(char*)malloc(300*sizeof(char));
    profileFile=(char*)malloc(300*sizeof(char));
    for(long i=0;i<300;i++){
      profileFile[i]=a.profileFile[i];
      lengthFile[i]=a.lengthFile[i];
    }
  }
  PChange_arg& PChange_arg::operator=(const PChange_arg& a){
    nStep=a.nStep;
    createProfile=a.createProfile;
    for(long i=0;i<300;i++){
      profileFile[i]=a.profileFile[i];
      lengthFile[i]=a.lengthFile[i];
    }
    return *this;
  }
}