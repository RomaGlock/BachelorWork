namespace math_our_cu{
  class Point_cu{
  public:
    __device__ Point_cu();
    __device__ Point_cu(double,double,double);    
    double x;
    double y;
    double z;
    __device__ Point_cu operator+(const Point_cu& p) const;
    __device__ Point_cu operator^(const Point_cu& p) const;
    __device__ const Point_cu& operator+=(const Point_cu& p);
    __device__ const Point_cu& operator-=(const Point_cu& p);
    __device__ const Point_cu& operator/=(double p);
    __device__ const Point_cu& operator*=(double p);
    __device__ Point_cu operator-(const Point_cu& p) const;
    __device__ double operator*(const Point_cu& p) const;
    __device__ Point_cu operator*(double p) const;
    __device__ Point_cu operator/(double p) const;
    __device__ double module() const;
    __device__ const Point_cu& normalize();
    __device__ void set(double,double,double);
  }; 
}