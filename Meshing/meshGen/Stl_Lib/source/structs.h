namespace stl
{
  struct Triangle{
        long a;
        long b;
        long c;
        math_our::Point n;
        long gr;
        Triangle(long,long,long,const math_our::Point&,long);
        Triangle():a(-1),b(-1),c(-1),n(math_our::Point(0,0,0)),gr(-1){};
        bool operator<(const Triangle& )const;
    };    
    struct Group{
        typedef std::list<Triangle> TriangleArray;
        std::string name;
        TriangleArray tri;
        explicit Group(const std::string& name):name(name){}
        explicit Group(const std::pair<std::string, Group>& p);
        const Group& operator= (const std::pair<std::string, Group>& p);
        Group(){};
    };
    struct edge{
        long a;
        long b;
        mutable long owner;
        mutable long neig;
        mutable long owner3rd;
        mutable long neig3rd;
        mutable math_our::Point n;
        mutable long group;
        bool operator<(const edge& p)const{
            if(a<p.a)return true;
            else if(a==p.a && b<p.b) return true;
            else return false;
        }    
        edge():a(-1),b(-1),owner(-1),neig(-1),owner3rd(-1),neig3rd(-1),n(math_our::Point(0,0,0)),group(-1){};
        edge(long a,long b,long owner,long neig,long owner3rd,long neig3rd,math_our::Point n,long group):
            a(a),b(b),owner(owner),neig(neig),owner3rd(owner3rd),neig3rd(neig3rd),n(n),group(group){};
        long addTriInEdgeTree(long nTri,const Triangle&,std::set<edge>&);
    };
    struct phantom_charge{
        double *charge;
        ulong nCharge;
        math_our::Point n;
        double D;
        phantom_charge():
            charge(NULL),nCharge(0),n(math_our::Point(0.,0.,0.)),D(0){
        }
        phantom_charge(double *a,ulong b,math_our::Point c,double d):
            charge(a),nCharge(b),n(c),D(d){
        }
    };
}
