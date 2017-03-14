namespace stl
{
  struct abmesh_arg{
        double l_u;
        double l_d;
        int n_u;
        int n_d;
        char* inFile;
        void addInFile(FILE* file);
        ~abmesh_arg();
        abmesh_arg();
        abmesh_arg(const abmesh_arg&);
        abmesh_arg& operator=(const abmesh_arg& );
        explicit abmesh_arg(const char* fileName);
    };
    struct MMLE_arg{
        double l;
        double l1,l2,la;
        double lBL;
        long nBL;
        double mBL;
        double lSW;
        long nSW;
        double mSW;
        double pow_q;
        double pow_E;
        long nStep;
        double elas;
        double elas2;
        long nCudaIter;
        char* lengthFile;
        char* profileFile;
        long createProfile;
        long firstLayer;
        void addInFile(FILE* file);
        MMLE_arg();
        ~MMLE_arg();
        MMLE_arg(const MMLE_arg&);
        MMLE_arg& operator=(const MMLE_arg&);
        explicit MMLE_arg(const char* fileName);
    };
    struct CUT_arg{
        math_our::Point n;
        double D;
        double cs; //minimal size
        void addInFile(FILE* file);
        CUT_arg();
        explicit CUT_arg(const char* fileName);
    };
    struct findWave_arg{
        double mul;
        char* fieldFile;
        char* outputFile;
        double eps;
        long made;//0 - OF, 1 - OF renumbered, 2 - Slava bin, 3 - Slava txt
        double TOL;
        long nRel;
        double elas;
        double minL;
        double iterLimit;
        void addInFile(FILE* file);
        ~findWave_arg();
        findWave_arg(const findWave_arg&);
        findWave_arg();
        findWave_arg& operator=(const findWave_arg&);
        explicit findWave_arg(const char* fileName);
    };
    struct OFMesh_arg{
        long def;
        char* wall;
        char* inlet;
        char* outlet;
        void addInFile(FILE* file);
        ~OFMesh_arg();
        OFMesh_arg();
        OFMesh_arg(const OFMesh_arg&);
        OFMesh_arg& operator=(const OFMesh_arg&);
        explicit OFMesh_arg(const char* fileName);
    };
    struct PChange_arg{
        char* lengthFile;
        long nStep;
        char* profileFile;
        long createProfile;
        void addInFile(FILE* file);
        ~PChange_arg();
        PChange_arg();
        PChange_arg(const PChange_arg&);
        PChange_arg& operator=(const PChange_arg&);
        explicit PChange_arg(const char* fileName);
    };
    struct Settings{
      PChange_arg pch;
      MMLE_arg mml;
      CUT_arg cut;
      OFMesh_arg ofm;
      findWave_arg ada;
      abmesh_arg abm;
      explicit Settings(const char* fileName);
      void Write(const char* fileName);
      void modi_par(const char*);
    };
}
