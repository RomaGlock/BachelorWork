#pragma once


#include <set>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "math_our.h"
#include "args.h"
#include "structs.h"
// PointArray >= 4
#define SSD_1 4
// NewPointArray >= 4
#define SSD_2 4
// ChargeArray >= 4
#define SSD_3 4



namespace stl
{
    void CudaSetDevice(int dev);

    double solve(
        const math_our::Point& x0,math_our::Point& e,
        const math_our::Point* a,unsigned long n,
        double alpha
        );
    long error(long,const char*,int);
    void profile_read(std::vector<double>& profile,FILE* file);
    char* get_suff(char*);
    math_our::Point find_image(math_our::Point p,math_our::Point n,double D);
    long addTriInEdgeTree(long nTri,const Triangle& p,std::set<edge>& Tree);
    long error(long i,const char *s,int d);
    class Stl_io
    {
        typedef math_our::Point Point;
        typedef std::vector<Point> PointArray;
        typedef std::map<std::string,Group> GroupTree;
        typedef std::vector<Group> GroupArray;
    public:    
        Stl_io();
        explicit Stl_io(char* fileName);
        Stl_io(char* fileName,double TOL);
        Stl_io(char* fileName,double TOL,int);
        Stl_io(const PointArray&,const GroupArray&);
        PointArray pointArray;
        GroupArray groupArray;
        bool is_3Dmesh;
        bool is_Cutted;//not used
        long nLayers;
        long nPoints;
        void Write(const char* fileName)const;
        void MMLE(MMLE_arg& p);
        void MMPoi(const char* fileName,Stl_io& phantom);
        void nmesh(const char* fileName);
        Stl_io extract_layer(long n)const;
        void abmesh_cu(const char* fileName);
        void ab_field(const char* fileName)const;
        Stl_io absurface(const char* fileName);
        void reverse_layers();
        long MakeOFMesh(math_our::Point Mach,const char* fileName)const;
        long MakeOFMesh(const char* fileName)const{
            return this->MakeOFMesh(math_our::Point (100,0,0),fileName);
        }
        void WriteBin(const char* fileName)const;
        void WriteOFF(const char* fileName)const;
        void WriteLiLu(const char* fileName)const;
        long ReadBin(const char* fileName);
        long ReadOFF(const char* fileName);
        Stl_io Cut_stl(const char* fileName);
        long findWave(const char*)const;
        template<class T> 
        void mix(T* a)const;
        template<class T> 
        void mix(T* a,long n)const;
        void changeProfile(const char* fileName);
        bool get_normals_av(math_our::Point **,long **)const;
        void findPnP(std::vector<std::vector<long> > & neiList)const;
        phantom_charge get_charges(double,long,const char*)const;
        void createLengthFile(const char* fileName)const;
        void writeMSHfield(const char* fileName,double* field)const;
        void llrelax(long nStep,double elas);
        void llrelax_cu(long nStep,double elas);
        void llrelax_cu_fast(long nStep,double elas,double elas2);
        void merge_out_in(const Stl_io& in,const char* fileName);
        void adjust_numbers(const Stl_io&,const char* fileName);
        //     Stl_io(const Stl_io&);
        //     const Stl_io& operator= (const Stl_io&);
    private:
        long addTriInEdgeTree(long nTri,const Triangle& p,std::set<edge>& Tree)const;
        void Delete_duplicate();
        void ReadStl(char* fileName);
        long find_point(const Point& p)const;
        void get_profile(std::vector<double>& profile,const char* fileName)const;
        double* get_lengths(const char* fileName)const;
        void read_length_file(double* l,FILE* file)const;
    };
    template<class T>
    void Stl_io::mix(T* a,long n)const{
        for(long i=0;i<abs(n);i++){
            mix<T>(a);
        }
    }
    template<class T>
    void Stl_io::mix(T* a)const{
        long nSurfaceElements=0;
        std::vector<long> nOwnerTriangle;
        std::vector<T> triValue;
        std::vector<T> pointValue;
        for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
            for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
                nSurfaceElements++;
            }
        }
        triValue.resize(nSurfaceElements);
        pointValue.resize(nPoints);
        nOwnerTriangle.resize(nPoints);
        long k=0;
        for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
            for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
                triValue[k]=(a[j->a]+a[j->b]+a[j->c])/3.;
                k++;
            }
        }
        k=0;
        for(GroupArray::const_iterator i=groupArray.begin();i!=groupArray.end();++i){
            for(Group::TriangleArray::const_iterator j = i->tri.begin(); j != i->tri.end(); ++j) {
                pointValue[j->a]+=triValue[k];
                pointValue[j->b]+=triValue[k];
                pointValue[j->c]+=triValue[k];
                nOwnerTriangle[j->a]++;
                nOwnerTriangle[j->b]++;
                nOwnerTriangle[j->c]++;
                k++;
            }
        }
        for(long i=0;i<nPoints;i++){
            pointValue[i]/=(double)nOwnerTriangle[i];
            a[i]=pointValue[i];
        }
    }

    void STLtoMESH(char *filename); //save stl file as mesh(yams) using gmsh
}