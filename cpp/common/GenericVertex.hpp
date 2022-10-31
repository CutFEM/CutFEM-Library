#ifndef _GENERIC_VERTEX_HPP
#define _GENERIC_VERTEX_HPP

#include "Label.hpp"


/*
 *    Class of the generic vertex
 *
 */
template<typename Rn>
class GenericVertex : public Rn,public Label
{
  template<typename T,typename B,typename V>
  friend class GenericMesh;                      

  friend inline std::ostream& operator <<(std::ostream& f, const GenericVertex & v )
  { f << (const Rn &) v << ' ' << (const Label &) v   ; return f; }
  friend inline std::istream& operator >> (std::istream& f,  GenericVertex & v )
        { f >> (Rn &) v >> (Label &) v ; return f; }
    
//   Rn *normal; // pointeur sur la normal exterieur pour filtre des points de departs
   
public:
  typedef Rn Rd;
  static const int d=Rd::d;  
  GenericVertex() : Rd(),Label() {};//,normal(0) {};
  GenericVertex(const Rd & P,int r=0): Rd(P),Label(r) {};//,normal(0){}
    
//   void SetNormal(Rd *&n,const Rd & N)
//   { if (normal) { 
//       Rd NN=*normal+N; 
//       *normal= NN/NN.norme(); }
//     else *(normal=n++)=N;

//     std::cout << *normal << std::endl;
//     getchar();
//   }
    
//   Rd Ne() const {return normal ? *normal: Rd();}
//   bool ninside(const Rd & P) const
//   {
//     return normal? (Rd(*this,P),*normal)<=0: true;
//   }

    
private: // pas de copie pour ne pas prendre l'adresse
  GenericVertex(const GenericVertex &);
  void operator=(const GenericVertex &); 
  
};



#endif
