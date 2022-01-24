#ifndef HASH_TABLE_HPP
#define HASH_TABLE_HPP

#include "../util/assertion.hpp"

template<typename T,int N>
struct SortArray {

};

template<typename T>
struct SortArray<T,1> {
  T v[1];
  SortArray(T *a)
  {
    v[0]=a[0];
  }
  SortArray(const T& a0)
  {
    v[0]=a0;
  }
  SortArray(){}
  bool operator == (const SortArray<T,1> & t)  const
  {  return v[0] == t.v[0]  ;}
  bool operator<(const SortArray<T,1> & t)  const
  {  return v[0] < t.v[0]  ;}
  size_t hash() const {return (size_t) v[0];}

};


template<typename T>
struct SortArray<T,2> {
  //  using std::std::swap;
  T v[2];
  SortArray(T *a)
  {
    v[0]=a[0];
    v[1]=a[1];
    if(v[0]>v[1]) std::swap(v[0],v[1]);
  }
  SortArray(const T& a0,const T &a1)
  {
    v[0]=a0;
    v[1]=a1;
    if(v[0]>v[1]) std::swap(v[0],v[1]);
  }
  SortArray(){}
  bool operator == (const SortArray<T,2> & t)  const
  {  return v[0] == t.v[0] && v[1] == t.v[1] ;}
  bool operator<(const SortArray<T,2> & t)  const
    {  return v[0] != t.v[0] ? v[0] < t.v[0] : v[1] < t.v[1] ;}
  T operator[](const int i)  const {  assert(i < 2);return v[i];}
  size_t hash() const {return (size_t) v[0];}
};


template<typename T>
struct SortArray<T,3> {
  T v[3];
  SortArray(T *a)
  {
    v[0]=a[0];
    v[1]=a[1];
    v[2]=a[2];
    if(v[0]>v[1]) std::swap(v[0],v[1]);
    if(v[1]>v[2]) {
      std::swap(v[1],v[2]);
      if(v[0]>v[1]) std::swap(v[0],v[1]);
    ASSERTION(v[0] <= v[1] && v[1] <= v[2] );
    }
  }

  SortArray(const T& a0,const T &a1, const T &a2)
  {
    v[0]=a0;
    v[1]=a1;
    v[2]=a2;
    if(v[0]>v[1]) std::swap(v[0],v[1]);
    if(v[1]>v[2]) {
      std::swap(v[1],v[2]);
      if(v[0]>v[1]) std::swap(v[0],v[1]);
      ASSERTION(v[0] <= v[1] && v[1] <= v[2] );
    }

  }


  SortArray(){}
  bool operator == (const SortArray<T,3> & t)  const
  {  return v[0] == t.v[0] && v[1] == t.v[1]  && v[2] == t.v[2] ;}

  bool operator<(const SortArray<T,3> & t)  const
    {  return v[0] != t.v[0] ? v[0] < t.v[0] :
           ( v[1] != t.v[1] ? v[1] < t.v[1] :  v[2] < t.v[2] );}
  T operator[](const int i)  const {  assert(i < 3);return v[i];}
  size_t hash() const {return (size_t) v[0];}
};



template<typename T>
struct SortArray<T,4> {
  T v[4];
  SortArray(T *a)
  {
    v[0]=a[0];
    v[1]=a[1];
    v[2]=a[2];
    v[3]=a[3];
    if(v[0]>v[1]) std::swap(v[0],v[1]);
    if(v[1]>v[2]) {
      std::swap(v[1],v[2]);
      if(v[0]>v[1]) std::swap(v[0],v[1]);
    }
    if(v[2]>v[3]) {
      std::swap(v[2],v[3]);
      if(v[1]>v[2]) {std::swap(v[1],v[2]);
	if(v[0] > v[1]) std::swap(v[0],v[1]);
      }
      ASSERTION(v[0] <= v[1] && v[1] <= v[2] && v[2] <= v[3] );
    }


  }

  SortArray(const T& a0,const T &a1, const T &a2, const T &a3)
  {
    v[0]=a0;
    v[1]=a1;
    v[2]=a2;
    v[3]=a3;
    if(v[0]>v[1]) std::swap(v[0],v[1]);
    if(v[1]>v[2]) {
      std::swap(v[1],v[2]);
      if(v[0]>v[1]) std::swap(v[0],v[1]);
    }
    if(v[2]>v[3]) {
      std::swap(v[2],v[3]);
      if(v[1]>v[2]) {std::swap(v[1],v[2]);
	if(v[0] > v[1]) std::swap(v[0],v[1]);
      }
      ASSERTION(v[0] <= v[1] && v[1] <= v[2] && v[2] <= v[3] );
    }

  }


  SortArray(){}
  bool operator == (const SortArray<T,4> & t)  const
  {  return v[0] == t.v[0] && v[1] == t.v[1]  && v[2] == t.v[2] && v[3] == t.v[3];}

  bool operator<(const SortArray<T,3> & t)  const
    {  return v[0] != t.v[0] ? v[0] < t.v[0] :
           ( v[1] != t.v[1] ? v[1] < t.v[1] :
	     ((v[2] != t.v[2]) ? v[2] < t.v[2] : v[3] < t.v[3]) );}
  T operator[](const int i)  const {  assert(i < 4);return v[i];}
  size_t hash() const {return (size_t) v[0];}
};


template<typename T,int N>
std::ostream & operator<<(std::ostream & f,const SortArray<T,N> & item)
{
    for (int i=0;i<N;++i) f << " " << item.v[i];
    return f;
}



template<class K,class V>
class HashTable {
public:
  struct   nKV { size_t next; K k; V v;
    nKV(){} };
  typedef nKV *iterator;
  size_t n,nx,nk,ncol,nfind;
  size_t * head;
  nKV * t;
  static const  size_t endhash= (size_t) -1;

  HashTable(size_t nnx,size_t nnk)
    :    n(0),nx(nnx),nk(nnk),ncol(0),nfind(0),
	 head(new size_t[nk]),t(new nKV[nx])
  {  reset();}

  void reset()
  {
    n=0;
    ncol=0;
    for (size_t j=0;j<nk;++j)
      head[j]=endhash;
  }

  nKV *  find(const K & key)
  {
    nfind++;
    for (size_t k=head[key.hash() %nk];k!=endhash;k=t[k].next)
      {
	++ncol;
	if(key == t[k].k) return t+k;
      }
    return 0;
  }
    // add FH  21 avril 2009
  size_t  operator()(nKV * p) { return p ? p-t : n;}

  iterator end(){ return t+n;}
  iterator begin(){ return t;}

  nKV *add(const K & key,const V & v)
  {
    size_t k =key.hash()%nk;
    assert(n<nx);
    t[n].v = v;
    t[n].k = key;
    t[n].next=head[k];
    head[k]=n;
    return t+ n++;
  }

  V & operator[](const K & key)
  {
    nKV *p = find(key);
    if(p) return p->v;
    else return t[add(key,V())].v;
  }
  ~HashTable()
  {
    // if(nfind && verbosity>4)
    //   cout << "    ~HashTable:   Cas moyen : " << (double) ncol/ nfind << endl;
    delete [] head;
    delete [] t;
  }

  //  pas de copie ....
private:
  HashTable(const HashTable&);
  void operator=(const HashTable&);
};


#endif
