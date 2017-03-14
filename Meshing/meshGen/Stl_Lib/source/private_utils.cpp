#include "stl_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
namespace stl
{
  using namespace std;
  void Stl_io::adjust_numbers(const Stl_io& stl2,const char* fileName){
    if(nPoints!=stl2.nPoints){
      printf("Aborted! nPoints does not match\n");
      return;
    }
    findWave_arg arg(fileName);
    math_our::TOL=arg.TOL;
    printf("TOL=%le\nrenumbering...\n",arg.TOL);
    typedef map<math_our::Point,long> mapPoint;
    mapPoint pmap;
    mapPoint::iterator pmapiter;
    for(long i=0;i<stl2.nPoints;i++){
      pmap.insert(pair<math_our::Point,long>(stl2.pointArray[i],i));
    }
    vector<long> newname;
    newname.resize(nPoints);
    for(long i=0;i<stl2.nPoints;i++){
      newname[i]=-1;
    }
    for(long i=0;i<stl2.nPoints;i++){
      pmapiter=pmap.find(pointArray[i]);
      if(pmapiter!=pmap.end()){
        if(newname[i]==-1){
          newname[i]=pmapiter->second;
        }else{
          printf("Aborted! Doublicated point\n");
          return;
        }
      }else{
        printf("Aborted! Point # %ld of %ld not found\n",i,nPoints);
        return;
      }
    }
    for(GroupArray::iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        j->a=newname[j->a];
        j->b=newname[j->b];
        j->c=newname[j->c];
      }
    }
    vector<math_our::Point> layerBuffer;
    layerBuffer.resize(pointArray.size());
    for(long j=0;j<=nLayers;j++){
      for(long i=0;i<nPoints;i++){
        layerBuffer[nPoints*j+newname[i]]=pointArray[nPoints*j+i];
      } 
    }
    swap(pointArray,layerBuffer);
    layerBuffer.clear();
    printf("success\n");
  }
  void Stl_io::reverse_layers(){
    vector<math_our::Point> layerBuffer;
    layerBuffer.reserve(pointArray.size());
    for(long i=nLayers;i>=0;--i){
//       layerBuffer.insert<math_our::Point*>(layerBuffer.end(),&(pointArray[nPoints*i]),&(pointArray[nPoints*(i+1)]));      
      layerBuffer.insert<vector<math_our::Point>::iterator>(layerBuffer.end(),pointArray.begin()+nPoints*i,pointArray.begin()+nPoints*(i+1));      
    }
    swap(pointArray,layerBuffer);
    layerBuffer.clear();
  }
  Stl_io Stl_io::extract_layer(long n)const{
    typedef vector<Group> GroupArray;
    GroupArray layerGroup;
    vector<math_our::Point> layerPoint;
    layerPoint.reserve(nPoints);
    layerGroup=groupArray;
//     layerPoint.insert<const math_our::Point*>(layerPoint.end(),&(pointArray[nPoints*n]),&(pointArray[nPoints*(n+1)]));
    layerPoint.insert<vector<math_our::Point>::const_iterator>(layerPoint.end(),pointArray.begin()+nPoints*n,pointArray.begin()+nPoints*(n+1)); 
    return Stl_io(layerPoint,layerGroup);
  }
  void Stl_io::findPnP(vector<vector<long> > & neiList)const{
    long tmpLong;
    //------------------------exeption list--------------------------------------------
    ulong k=0;
    set<edge> edgeTree;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        addTriInEdgeTree(k,*j,edgeTree);
        k++;
      }
    }
    //------------------------------make owner<neighbour--------------------------
    long* exept;
    exept=(long*)malloc(nPoints*sizeof(long));
    for(long i=0;i<nPoints;i++){
      exept[i]=1;
    }
    for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
      if(i->neig==-1){
        exept[i->a]=0;
        exept[i->b]=0;
      }
    }  
    typedef vector<vector<pair<long,long> > > NeiFaces;
    NeiFaces neiFaces;
    neiFaces.resize(nPoints);
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
        neiFaces[j->a].push_back(pair<long,long>(j->b,j->c));
        neiFaces[j->b].push_back(pair<long,long>(j->c,j->a));
        neiFaces[j->c].push_back(pair<long,long>(j->a,j->b));
      }
    }
    {
      bool flag;
      pair<long,long> buff;
      list<long> pts;
      list<long>::iterator jj,jjj;
      for(long i=0;i<nPoints;++i){
        if(exept[i]){
          for(ulong j=0;j<(neiFaces[i].size()-1);++j){
            flag=true;
            for(k=j+1;k<neiFaces[i].size() && flag;++k){
              if(neiFaces[i][k].first==neiFaces[i][j].second){
                flag=false;
                buff=neiFaces[i][k];
              }else if(neiFaces[i][k].second==neiFaces[i][j].second){
                flag=false;
                buff.first=neiFaces[i][k].second;
                buff.second=neiFaces[i][k].first;
              }
            }
            if(flag){
              printf("error connections %d\n",__LINE__);
            }else{
              neiFaces[i][k-1]=neiFaces[i][j+1];
              neiFaces[i][j+1]=buff;
            }
          }
          for(ulong j=0;j<(neiFaces[i].size()-1);++j){
            if(neiFaces[i][j].second!=neiFaces[i][j+1].first){
              printf("error connections %d\n",__LINE__);
            }
          }
          if(neiFaces[i][neiFaces[i].size()-1].second!=neiFaces[i][0].first){
            printf("error connections %d\n",__LINE__);
          }
        }else{
          pts.clear();
          for(ulong j=0;j<neiFaces[i].size();++j){
            pts.push_back(neiFaces[i][j].first);
            pts.push_back(neiFaces[i][j].second);
          }
          pts.sort();
          tmpLong=-1;
          jj=pts.begin();
          ++jj;
          jjj=jj;
          ++jjj;
          for(list<long>::iterator j=pts.begin();tmpLong ==-1 && j!=pts.end() && jj!=pts.end();++j){
            jj=j;
            ++jj;
            jjj=jj;
            ++jjj;
            if(j==pts.begin()){
              if(*jj!=*j){
                tmpLong=*j;
              }
            }else if(jjj==pts.end()){
              if(*jj!=*j){
                tmpLong=*jj;
              }
            }else if(*jj!=*j && *jj!=*jjj){
              tmpLong=*jj;
            }
          }
          if(tmpLong==-1){
            printf("error connections %d\n",__LINE__);
          }
          flag=true;
          for(k=0;flag && k<neiFaces[i].size();++k){
            if(neiFaces[i][k].first==tmpLong){
              flag=false;
              buff=neiFaces[i][k];
            }else if(neiFaces[i][k].second==tmpLong){
              flag=false;
              buff.first=neiFaces[i][k].second;
              buff.second=neiFaces[i][k].first;
            }
          }
          if(flag){
            printf("error connections %d\n",__LINE__);
          }else{
            neiFaces[i][k-1]=neiFaces[i][0];
            neiFaces[i][0]=buff;
          }
          for(ulong j=0;j<(neiFaces[i].size()-1);++j){
            flag=true;
            for(k=j+1;k<neiFaces[i].size() && flag;++k){
              if(neiFaces[i][k].first==neiFaces[i][j].second){
                flag=false;
                buff=neiFaces[i][k];
              }else if(neiFaces[i][k].second==neiFaces[i][j].second){
                flag=false;
                buff.first=neiFaces[i][k].second;
                buff.second=neiFaces[i][k].first;
              }
            }
            if(flag){
              printf("error connections %d\n",__LINE__);
            }else{
              neiFaces[i][k-1]=neiFaces[i][j+1];
              neiFaces[i][j+1]=buff;
            }
          }
          for(ulong j=0;j<(neiFaces[i].size()-1);++j){
            if(neiFaces[i][j].second!=neiFaces[i][j+1].first){
              printf("error connections %d\n",__LINE__);
            }
          }
          neiList[i].push_back(neiFaces[i][0].first);          
        }
        for(ulong j=0;j<neiFaces[i].size();++j){
          neiList[i].push_back(neiFaces[i][j].second);
        }
      }
    }
    free(exept);
  }
  void Stl_io::writeMSHfield(const char* fileName,double* field)const{
    typedef std::list<Triangle> TriangleArray;
    typedef std::vector<Group> GroupArray;
    long k=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(TriangleArray::const_iterator j=i->tri.begin();j!=i->tri.end();++j){
	k++;
      }
    }
    FILE *out;
    int N_dif = (int)nPoints;
    int N_Triangle = (int)k;
    out = fopen (fileName, "wb");
    fprintf (out, "$MeshFormat\n");
    fprintf (out, "2.2 0 8\n");          //0 - ASCII, 1 - binary
    fprintf (out, "$EndMeshFormat\n");

//#############################   запись кол-во и все точки различные
    fprintf (out, "$Nodes\n");
    fprintf (out, "%d\n", N_dif);
    for (int i=0; i<N_dif; i++){
	    fprintf (out, "%d  %le  %le  %le\n",
		     i+1, pointArray[i].x, pointArray[i].y, pointArray[i].z);
    }
    fprintf (out, "$EndNodes\n");
	
//#############################   запись кол-ва и все частицы
    fprintf (out, "$Elements\n");
    fprintf (out, "%d\n", N_Triangle);
    k=1;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(TriangleArray::const_iterator j=i->tri.begin();j!=i->tri.end();++j){
	  fprintf (out, "%ld  %ld  2  6  100  ", k, (long)2);
	  fprintf (out, "%ld  %ld  %ld  ", (j->a)+(long)1,(j->b)+(long)1,(j->c)+(long)1);
	  fprintf (out, "\n");
	  k++;
      }
    }
    fprintf (out, "$EndElements\n");

//#############################   запись данных и количества
    fprintf (out, "$ElementData\n");

    fprintf(out,"%i\n",1);//one string tag:                                                  //??? ??? ????? ??!!
    fprintf(out,"A scalar view\n");//the name of the view ("A scalar view")
    fprintf(out,"%i\n",1);//one real tag
    fprintf(out,"%lf\n",0.0);//the time value 
    fprintf(out,"%i\n",3);//three integer tags
    fprintf(out,"%i\n",0);//the time step (0; time steps always start at 0)
    fprintf(out,"%i\n",1);//1-component (scalar) field
    fprintf (out, "%d\n", N_Triangle);
    double* efield;
    efield=(double*)malloc(sizeof(double)*N_Triangle);
    math_our::Point e1,e2,e3;
    k=0;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(TriangleArray::const_iterator j=i->tri.begin();j!=i->tri.end();++j){
	e1=pointArray[j->b]-pointArray[j->a];
	e2=pointArray[j->c]-pointArray[j->a];
	e3=e1^e2;
	efield[k]=(field[j->a]+field[j->b]+field[j->c])/fabs(e3.module());
	k++;
      }
    }
    for (int i=0; i<N_Triangle; i++)
	    fprintf (out, "%d  %le\n", i+1, efield[i]);
    fprintf (out, "$EndElementData\n");
    fclose(out);
    free(efield);
  }
  char* get_suff(char* a){
    long i=0;
    char* r=a;
    while(i<300 && a[i]!='\0'){
      if(a[i]=='.'){
	r=&(a[i+1]);
      }
      ++i;
    }
    return r;
  }
  bool Stl_io::get_normals_av(math_our::Point** np,long ** cp)const{
    typedef std::vector<Group> GroupArray;
    typedef std::list<Triangle> TriangleArray;
    math_our::Point *nList;
    math_our::Point zeroPoint(0,0,0);
    math_our::Point tmpPoint;
    nList=(math_our::Point*)malloc(sizeof(math_our::Point)*nPoints);
    long *nNList;
    nNList=(long*)malloc(sizeof(long)*nPoints);
    *np=nList;
    *cp=nNList;
    for(long i=0;i<nPoints;i++) {
      nList[i]=zeroPoint;
      nNList[i]=0;
    }
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(TriangleArray::const_iterator j=i->tri.begin();j!=i->tri.end();++j){
	tmpPoint=((pointArray[nLayers*nPoints+j->b]-pointArray[nLayers*nPoints+j->a])^(pointArray[nLayers*nPoints+j->c]-pointArray[nLayers*nPoints+j->a])).normalize();
// 	tmpPoint=j->n;	
// 	tmpPoint.normalize();
	nList[j->a]+=tmpPoint;
	nList[j->b]+=tmpPoint;
	nList[j->c]+=tmpPoint;
	nNList[j->a]++;
	nNList[j->b]++;
	nNList[j->c]++;
      }
    }
    for(long i=0;i<nPoints;i++) {
      if(nNList[i]==0 || nList[i]==zeroPoint){
	return false;
      }else{
	nList[i]/=nNList[i];
      }
    }
    return true;
  }
  void Stl_io::read_length_file(double* lList,FILE* file)const{
    long tmpLong;
    double tmpDouble;
    long n;
    char string[300];
    char* a=string;
    long i=0;
    long lust=-1;
    while(a && i<nPoints){
      a=fgets(string,300,file);
      if(a){
	tmpLong=sscanf(string,"%20lf %20ld",&tmpDouble,&n);
	if(tmpLong==1){
	  
	  lList[i]=tmpDouble;
	  lust=i;
	  ++i;
	}else if(tmpLong==2){
	  for(long j=0;j<n;j++){
	    if((i+j)<nPoints){
	      lList[i+j]=tmpDouble;
	      lust=i+j;
	    }
	  }
	  i+=n;
	}
      }
    }
    ++lust;
    if(lust!=nPoints){
      printf("length File is incorrect\
      \n %ld of %ld lines only\
      \n. Results is unpredictable!\n",lust,nPoints);
    }
  }
  void Stl_io::createLengthFile(const char* fileName)const{
    MMLE_arg arg(fileName);
    FILE* file;    
    char str[300];    
    sprintf(str,"rm -f %s",arg.lengthFile);    
    system(str);    
    double* lengths=get_lengths(fileName);    
    file=fopen(arg.lengthFile,"wb");
    if(!file){
      printf("createLengthFile error: cannot open %s\n",arg.lengthFile);
    }
    for(long i=0;i<nPoints;i++){
      fprintf(file,"%15lf\n",lengths[i]);
    }
    free(lengths);
    fclose(file);
  }
  double* Stl_io::get_lengths(const char* fileName)const{
    double tmpDouble;
    FILE *file;
    MMLE_arg arg(fileName);
    CUT_arg cut_arg(fileName);
    double minX=1e100;
    double maxX=-1e100;
    double* lList;
    double f_;    
    lList=(double*)malloc(sizeof(double)*nPoints);    
    for(long i=0;i<nPoints;i++) {
      if(pointArray[i]*cut_arg.n<minX)minX=pointArray[i]*cut_arg.n;
      if(pointArray[i]*cut_arg.n>maxX)maxX=pointArray[i]*cut_arg.n;
    }    
    for(long i=0;i<nPoints;i++) {
      tmpDouble=(pointArray[i]*cut_arg.n-minX)/(maxX-minX);
      f_=pow((pow(arg.l2,1./arg.la)-pow(arg.l1,1./arg.la))*tmpDouble+pow(arg.l1,1./arg.la),arg.la);
      lList[i]=f_*(maxX-minX)/8.;
    }    
    file=fopen(arg.lengthFile,"rb");
    if(file){
      read_length_file(lList,file);
      fclose(file);
    }else{
      printf("%s: information: no lengthFile \"%s\"; approxLength is used\n",__FILE__,arg.lengthFile);
      for(long ii=0;ii<nPoints;ii++) {
	lList[ii]*=arg.l;
      }
    }
    return lList;
  }
  phantom_charge Stl_io::get_charges(double pow_q,long FirststLayerOri,const char* fileName)const{
    typedef std::vector<Group> GroupArray;
    typedef std::list<Triangle> TriangleArray;
    double *charges;
    double tmpDouble;
    math_our::Point e1,e2,e3,tmpPoint;
    charges=(double*)malloc(sizeof(double)*SSD_3*(pointArray.size()*2)); //x,y,z,q
    math_our::Point *phnList;
    long *phnNList;
    get_normals_av(&phnList,&phnNList);
    CUT_arg cut_arg(fileName);
    double minX=1e100;
    double maxX=-1e100;
    for(long i=0;i<nPoints;i++) {
      if(pointArray[i]*cut_arg.n<minX)minX=pointArray[i]*cut_arg.n;
      if(pointArray[i]*cut_arg.n>maxX)maxX=pointArray[i]*cut_arg.n;
    }
    for(unsigned long i=0;i<pointArray.size();i++){
      tmpDouble=phnList[i].module();
      tmpPoint=phnList[i]*(maxX-minX)/300./phnNList[i]*1./tmpDouble*(-(double)FirststLayerOri);
      charges[SSD_3*i+0]=pointArray[i].x+tmpPoint.x;
      charges[SSD_3*i+1]=pointArray[i].y+tmpPoint.y;
      charges[SSD_3*i+2]=pointArray[i].z+tmpPoint.z;
      charges[SSD_3*i+3]=0.;      
    }
    free(phnNList);
    free(phnList);
    double minS=1e100;
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(TriangleArray::const_iterator j=i->tri.begin();j!=i->tri.end();++j){
	e1=pointArray[j->b]-pointArray[j->a];
	e2=pointArray[j->c]-pointArray[j->a];
	e3=e1^e2;
	minS=(e3.module()<minS)?e3.module():minS;
      }
    }
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(TriangleArray::const_iterator j=i->tri.begin();j!=i->tri.end();++j){
	e1=pointArray[j->b]-pointArray[j->a];
	e2=pointArray[j->c]-pointArray[j->a];
	e3=e1^e2;
	tmpDouble=pow(e3.module()/minS,pow_q);
	charges[SSD_3*j->a+3]+=tmpDouble;
	charges[SSD_3*j->b+3]+=tmpDouble;
	charges[SSD_3*j->c+3]+=tmpDouble;
	tmpDouble=pow(e3.module(),0.5);
      }
    }
    //----------------------mix charge-----------------------------------------
    {
      double* tmpCharges;
      tmpCharges=(double*)malloc(sizeof(double)*nPoints);
      for(long i=0;i<nPoints;i++){
	tmpCharges[i]=charges[SSD_3*i+3];
      }
      mix(tmpCharges,0);
      for(long i=0;i<nPoints;i++){
	charges[SSD_3*i+3]=tmpCharges[i];
      }
      free(tmpCharges);
    }
    //--------------------------big experimental charge------------------------
    long k=0;
    long nNeighbours=0;
    long nOwners=0;
    set<edge> edgeTree;
    vector<long> exept;
    long nExept=0;
    exept.reserve(nPoints);
    for(long i=0;i<nPoints;i++){
      exept.push_back(1);
    }
    for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
      for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
	addTriInEdgeTree(k,*j,edgeTree);
	k++;
      }
    }
    //------------------------------make owner<neighbour--------------------------
    for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
      if(i->neig!=-1){
	if(i->neig<i->owner){
	  long buff;
	  buff=i->neig;
	  i->neig=i->owner;
	  i->owner=buff;
	  buff=i->neig3rd;
	  i->neig3rd=i->owner3rd;
	  i->owner3rd=buff;	  
	}
	nNeighbours++;
      }else if(i->neig==-1){
        exept[i->a]=0;
        exept[i->b]=0;
      }
      nOwners++;
    }
    for(long i=0;i<nPoints;i++){
      if(!(exept[i])){
        ++nExept;
      }
    }
    if(nNeighbours==nOwners){
      phantom_charge ret(charges,pointArray.size(),math_our::Point(0.,0.,0.),0.);
      return ret;
    }else{
      math_our::Point tmpPoint;
      math_our::Point mass_centre(0.,0.,0.);
      vector<math_our::Point> circle;
      circle.reserve(edgeTree.size());
      k=0;
      for(set<edge>::iterator i=edgeTree.begin();i!=edgeTree.end();++i){
	if(i->neig==-1){
	  tmpPoint=(pointArray[i->a]+pointArray[i->b])/2.;
	  circle.push_back(tmpPoint);
	  mass_centre+=tmpPoint;
	  k++;
	}
      }
      mass_centre/=circle.size();
      tmpDouble=0.;
      long farest=-1;
      for(ulong i=0;i<circle.size();++i){
	if((circle[i]-mass_centre).module()>tmpDouble){
	  tmpDouble=(circle[i]-mass_centre).module();
	  farest=i;
	}
      }
      if(farest<0){
	printf("%s(%d) phantom big charge error \n",__FILE__,__LINE__);
	fflush(stdout);
      }
      tmpDouble=0.;
      long tmpLong=farest;
      for(ulong i=0;i<circle.size();++i){
	double mul;
	mul=((circle[farest]-mass_centre)^(circle[i]-mass_centre)).module();
	if(mul>tmpDouble){
	  tmpDouble=mul;
	  tmpLong=(long)i;
	}
      }
      if(tmpLong==farest){
	printf("%s(%d) phantom big charge error \n",__FILE__,__LINE__);
	fflush(stdout);
      }
      math_our::Point cp[3];
      cp[0]=mass_centre;
      cp[1]=circle[tmpLong];
      cp[2]=circle[farest];
      math_our::Point charge_centre(0.,0.,0.);
      double big_charge=0.;
      for(unsigned long i=0;i<pointArray.size();i++){
	tmpPoint.x=charges[SSD_3*i+0];
	tmpPoint.y=charges[SSD_3*i+1];
	tmpPoint.z=charges[SSD_3*i+2];
	tmpDouble =charges[SSD_3*i+3];
	charge_centre+=((tmpPoint)*tmpDouble);
	big_charge+=tmpDouble;
      }
      charge_centre/=big_charge;
      math_our::Point n=((cp[1]-cp[0])^(cp[2]-cp[0])).normalize();
      double D=-(n*cp[0]);
      k=0;
      for(long i=0;i<nPoints;i++){
        if(exept[i]){
          tmpPoint.x=charges[SSD_3*i+0];
          tmpPoint.y=charges[SSD_3*i+1];
          tmpPoint.z=charges[SSD_3*i+2];
          tmpPoint=find_image(tmpPoint,n,D);
          charges[SSD_3*(nPoints+k)+0]=tmpPoint.x;
          charges[SSD_3*(nPoints+k)+1]=tmpPoint.y;
          charges[SSD_3*(nPoints+k)+2]=tmpPoint.z;
          charges[SSD_3*(nPoints+k)+3]=charges[SSD_3*(i)+3];
          ++k;
        }else{
          tmpPoint.x=charges[SSD_3*i+0];
          tmpPoint.y=charges[SSD_3*i+1];
          tmpPoint.z=charges[SSD_3*i+2];
          e1=find_image(tmpPoint,n,D);
          tmpPoint=(tmpPoint+e1)*0.5;
          charges[SSD_3*i+0]=tmpPoint.x;
          charges[SSD_3*i+1]=tmpPoint.y;
          charges[SSD_3*i+2]=tmpPoint.z;
        }
      }
      phantom_charge ret(charges,nPoints*2-nExept,n,D);
      return ret;
    }
    //-----------------------------------------------------------------------------
  }
  void profile_read(std::vector<double>& profile,FILE* file){
    double tmpDouble;
    long tmpLong;
    long n;
    char string[300];
    char* a=string;
    while(a){
      a=fgets(string,300,file);
      if(a){
	tmpLong=sscanf(string,"%20lf %20ld",&tmpDouble,&n);
	if(tmpLong==1){	  
	  profile.push_back(tmpDouble);
	}else if(tmpLong==2){
	  for(long j=0;j<n;j++){
	      profile.push_back(tmpDouble);
	  }
	}
      }
    }
    printf("Info: %lu layers will be created\n",(ulong)(profile.size()));
    tmpDouble=0;
    for(ulong i=0;i<profile.size();i++){
      tmpDouble+=profile[i];
    }
    for(ulong i=0;i<profile.size();i++){
      profile[i]/=tmpDouble;
    }
  }
  void Stl_io::get_profile(std::vector<double> & profile,const char* fileName)const{
    //S(n)=b*(1-q^n) / (1-q)
    MMLE_arg arg(fileName);
    FILE *file;
    file=fopen(arg.profileFile,"rb");
    if(file){
      profile_read(profile,file);
      fclose(file);
    }else{
      printf("%s does not exist\n",arg.profileFile);
      arg.nCudaIter=abs(arg.nCudaIter);
      profile.resize(arg.nCudaIter);
      if((arg.nBL+arg.nSW)<=arg.nCudaIter){
	double Sn=arg.lBL;
	double q=arg.mBL;
	double n=arg.nBL;
	double b1=(1-q)/(1-pow(q,n))*Sn;
	for(long i=0;i<n;i++){
	  profile[i]=b1*pow(q,i);
	}
	long nIso=arg.nCudaIter-arg.nSW-arg.nBL;
	double lIso=(1-arg.lSW-arg.lBL)/nIso;
	for(long i=0;i<nIso;i++){
	  profile[i+arg.nBL]=lIso;
	}
	Sn=arg.lSW;
	q=arg.mSW;
	n=arg.nSW;
	b1=(1-q)/(1-pow(q,n))*Sn;
	for(long i=0;i<n;i++){
	  profile[arg.nCudaIter-i-1]=b1*pow(q,i);
	}
      }else{
	printf("Warning! nBL+nSW>nLayers : \n\
	thickness of layers will be constant\n");
	for(long i=0;i<arg.nCudaIter;i++){
	  profile[i]=1./(double)arg.nCudaIter;
	}
      }
    }
    if(arg.createProfile){
      file=fopen(arg.profileFile,"wb");
      printf("creating profile file: %lu layers\n",(unsigned long)profile.size());
      if(file){
	for(ulong i=0;i<profile.size();i++){
	  fprintf(file,"%lf\n",profile[i]);
	}
	fclose(file);
      }
    }
  }
}
