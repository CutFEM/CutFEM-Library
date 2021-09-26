#ifndef _DA_HPP_
#define _DA_HPP_


#include <stdlib.h>
#include <cmath>
#include <iostream>


template <class R,int N=1> 
struct Diff 
{
  R val; 
  R d[N]; 
  Diff() : val(0)
  {
    for(int i=0; i<N;++i)
      d[i]=0;
  }
  Diff(const R &a) : val(a) 
  {
    for(int i=0; i<N;++i)
      d[i]=0;
  }
  // Diff(const R &a,const R &da) : val(a)
  // {   
  //   for(int i=0; i<N;++i)
  //     d[i]=0;
  //   d[0]=da; 
  // }
  Diff(const Diff &a,const int k) : val(a.val) 
  { 
    for(int i=0; i<N;++i)
      d[i]=(i==k);
  }
  Diff(const double &a,const int k) : val(a)
  {   
    for(int i=0; i<N;++i)
      d[i]=(i==k);
  }
  Diff& operator  + (){return *this;};
  Diff& operator  = (double);
  Diff& operator -= (double);
  Diff& operator -= (const Diff&);
  Diff& operator += (double) ;
  Diff& operator += (const Diff&) ;
  Diff& operator *= (double);
  Diff& operator *= (const Diff&);
  Diff& operator /= (double) ;
  Diff& operator /= (const Diff&) ;
  
};

template <class R,int N>  
std::ostream& operator<<(std::ostream& f, const Diff<R,N>& a)
{
  f << "val = "  << a.val << std::endl; 
  f << "[ ";
  for(int i=0;i<N;++i)
    f << a.d[i] << "  ";
  f << "] " << std::endl;
  return f;}


template <class R,int N>  
Diff<R,N> operator + (double x, const Diff<R,N>& y)
{   Diff<R,N> r(y); r.val += x;  return r;}

template <class R,int N>  
Diff<R,N> operator + (const Diff<R,N>& y, double x)
{   Diff<R,N> r(y);  r.val += x;  return r;}

template <class R,int N>  
Diff<R,N> operator-(const Diff<R,N>  & x,const Diff<R,N>  & y)
{
  Diff<R,N> r(x.val-y.val);
  for(int i=0; i<N;++i)	     
    r.d[i] = x.d[i]-y.d[i]; 
  return r; 
}
template <class R,int N>  
Diff<R,N> operator - (double x, const Diff<R,N>& y)
{  Diff<R,N> r;
  for(int i=0; i<N;++i) r.d[i]  = (-1)* y.d[i];
  r.val = x - y.val;   return r;}

template <class R,int N>  
Diff<R,N> operator - (const Diff<R,N>& x, double y)
{  Diff<R,N>  r(x);  r.val -= y;  return r;        }



template <class R,int N>  
Diff<R,N> operator+(const Diff<R,N>  & x,const Diff<R,N>  & y)
{
  Diff<R,N> r(x.val+y.val);
  for(int i=0; i<N;++i)	     
    r.d[i] = x.d[i]+y.d[i]; 
  return r; 
}



template <class R,int N>
Diff<R,N> operator*(const Diff<R,N>  & x,const Diff<R,N>  & y)
{
  Diff<R,N> r(x.val*y.val);
  for(int i=0; i<N;++i)	     
    r.d[i] = x.d[i]*y.val + x.val*y.d[i]; 
  return r; 
}

template <class R,int N>
Diff<R,N> operator*(const Diff<R,N>  & x,const double &y)
{
  Diff<R,N> r(x.val*y);
  for(int i=0; i<N;++i)	     
    r.d[i] = x.d[i]*y ; 
  return r; 
}

template <class R,int N>
Diff<R,N> operator*(const double & y, const Diff<R,N>  & x)
{
  Diff<R,N> r(x.val*y);
  for(int i=0; i<N;++i)	     
    r.d[i] = y*x.d[i] ; 
  return r; 
}

template <class R,int N>
Diff<R,N> sqrt(const Diff<R,N>  & x)
{
  Diff<R,N> r(sqrt(x.val));
  for(int i=0; i<N;++i)	     
    r.d[i] = x.d[i]*0.5/r.val; 
  return r; 
}

template <class R,int N>
Diff<R,N> operator/(const Diff<R,N> & x, const Diff<R,N>  & y)
{
  Diff<R,N> r(x.val/y.val);
  for(int i=0; i<N;++i)	     
    r.d[i] = (x.d[i]*y.val-x.val*y.d[i])/(y.val*y.val); 
  return r; 
}

template <class R,int N>
Diff<R,N> operator/(const Diff<R,N> & x, const double  & y)
{
  Diff<R,N> r(x.val/y);
  for(int i=0; i<N;++i)	     
    r.d[i] = x.d[i] / y; 
  return r; 
}

template <class R,int N>
Diff<R,N> operator/(const double & x, const Diff<R,N>  & y)
{
  Diff<R,N> r(x/y.val);
  for(int i=0; i<N;++i)	     
    r.d[i] = (-x*y.d[i])/(y.val*y.val); 
  return r; 
}


template <class R,int N>
Diff<R,N>& Diff<R,N>::operator = (double y)
{  val =y;
  for(int i=0; i<N;++i) d[i] = 0;
  return *this;
}



template <class R,int N>
Diff<R,N>& Diff<R,N>::operator += (double y)
{  val += y; return *this; }

template <class R,int N>
Diff<R,N>& Diff<R,N>::operator += (const Diff<R,N>& y)
{   val += y.val;
  for(int i=0; i<N;++i) d[i] += y.d[i];
  return *this; }

template <class R,int N>
Diff<R,N>& Diff<R,N>::operator -= (double y)
{  val -= y; return *this; }

template <class R,int N>
Diff<R,N>& Diff<R,N>::operator -= (const Diff<R,N>& y)
{   val -= y.val;
  for(int i=0; i<N;++i) d[i] -= y.d[i];
  return *this; }


template <class R,int N>
Diff<R,N>& Diff<R,N>::operator *= (double y)
{  val*=y;
  for(int i=0; i<N;++i) d[i]*=y;
  return *this;
}
template <class R,int N>
Diff<R,N>& Diff<R,N>::operator *= (const Diff<R,N>& y)
{ return *this = *this * y;}


template <class R,int N>
Diff<R,N>& Diff<R,N>::operator /= (const Diff<R,N>& y)
{  return *this = *this / y;}

template <class R,int N>
Diff<R,N>& Diff<R,N>::operator /= (double y)
{ const double inv = 1.0 / y;
    val *= inv;
    for(int i=0; i<N;++i) d[i] *= inv;
    return *this;
}

template <class R,int N>
Diff<R,N> pow (const Diff<R,N>& x,const int y)
  {
    Diff<R,N> r(1);    
    if(y>=0) for(int i=0;i<y;i++) r*=x;
    else for(int i=0; i<-y;i++) r/=x;
    return r;
  }


template<int N>
R laplacian(Diff< Diff<double,N>,N> & x) {
  R l = 0;
  for(int i=0; i<N;++i) l += x.d[i].d[i];

  return l;
}



template <class R,int N>
Diff<R,N> sin (const Diff<R,N>& x)
{   Diff<R,N> r; r.val=sin(x.val);
  for(int i=0; i<N;++i) r.d[i]=x.d[i]*cos(x.val);
  return r;
}

template <class R,int N>
Diff<R,N> cos (const Diff<R,N>& x)
{   Diff<R,N> r; r.val=cos(x.val);
  for(int i=0; i<N;++i) r.d[i]=(-1)*x.d[i]*sin(x.val);
  return r;
}


template <class R,int N>
Diff<R,N> atan (const Diff<R,N>& x)
{   Diff<R,N> r; r.val=atan(x.val);
  for(int i=0; i<N;++i) r.d[i] = 1 / (1 + x.val*x.val) * x.d[i];
  return r;
}


#endif
