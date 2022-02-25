#ifndef FRACTURE_HPP_
#define FRACTURE_HPP_

#include "GenericInterface.hpp"

// A partition of a cut element should be like a little mesh
// array points
// array elements etc
class Element2 {
  typedef typename Mesh2::Rd Rd;
  typedef SortArray<int, 3> TriaIdx;


  public :
  static const int nv = 3;
  Rd  vertices[nv];
  int lab;

  Element2() {}
  Element2(const std::vector<Rd>& v0, TriaIdx& iv,int r=0) {
    for(int i=0;i<nv;++i) vertices[i]=v0[iv[i]];
    lab=r;
  }
  const R2& operator[](int i) const {return vertices[i];}
  // R2& operator[](int i) {return vertices[i];}

  void print() const  {
    for(int i=0;i<nv;++i) {
      std::cout << vertices[i] << std::endl;
    }
  }
};

class Local_Partition{

  typedef Element2 Element;
  typedef typename Mesh2::Rd Rd;
  typedef SortArray<int, 3> TriaIdx;

public:

private:
  std::vector<Rd>       vertices_;
  std::vector<TriaIdx>  elements_idx;

public:


  Local_Partition() {}
  void reset() {
    vertices_.resize(0);
    elements_idx.resize(0);
  }
  void add_node(const Rd P) { vertices_.push_back(P);}
  void add_triangle(const TriaIdx& iv) {
    elements_idx.push_back(iv);
  }
  void erase_triangle(int k) {
    elements_idx.erase(elements_idx.begin()+k);
  }
  int nb_element() const {return elements_idx.size();}
  int nb_node() const {return vertices_.size();}

  Element get_element(int k) { return Element(vertices_,elements_idx[k]);}


};

class Cut_Base {

  typedef SortArray<int, 3> TriaIdx;

  public :
  int label;

  Cut_Base(int lab) : label(lab) {}
private:
  std::map<int, int> map_element_to_face;
  virtual double f(const R2 A) const = 0;

public:
  bool is_negative(const R2 P) const {
    return (f(P) <=0);
  }
  bool is_positive(const R2 P) const {
    return (f(P) > 0);
  }
  bool cut_edge(const R2 A, const R2 B) const {
    return (f(A)*f(B)< 0);
  }
  R2 get_intersection_node(const R2 A, const R2 B) const {
    double t = -f(A)/(f(B)-f(A));
    return (1-t) * A + t * B;
  }
  double eval(const R2 A) const {return f(A);}

  void add_element_and_face(int k, int iface) { map_element_to_face[k] = iface;}
  bool is_cut_element(int k) const {
    return (map_element_to_face.find(k) != map_element_to_face.end());
  }

  void cut_partition(Local_Partition& local_partition, vector<TriaIdx>& new_element_idx, std::list<int>& erased_element)const {
    new_element_idx.resize(0);
    erased_element.resize(0);
    byte ls[3];
    for(int k=0; k<local_partition.nb_element();++k){

      // BUILD ELEMENT
      const Element2 K = local_partition.get_element(k);
      for(int i=0;i<3;++i) ls[i] = util::sign(eval(K[i]));

      //COMPUTE THE PARTITION
      const RefPartition2& patch(RefPartition2::instance(ls));

      // THE FRACTURE i CUT THIS LOCAL ELEMENT k
      if(patch.is_cut()) {
        // add k to the list to element to remove
        erased_element.push_back(k);

        // LOOP OVER ELEMENTS IN PATCH
        for(auto it = patch.element_begin(); it != patch.element_end(); ++it) {
          // create the Nodes
          std::cout << " index node to create " << std::endl;
          int idx_begin = local_partition.nb_node();
          TriaIdx idx_array(idx_begin, idx_begin+1, idx_begin+2);
          for(int i=0; i<3;++i) {
            Uint idx = (*it)[i];
            std::cout << idx << std::endl;

            if(idx < 3) { local_partition.add_node(K[idx]);}
            else{
              int i0 = Mesh2::Element::nvedge[idx - 3][0];
              int i1 = Mesh2::Element::nvedge[idx - 3][1];
              local_partition.add_node(get_intersection_node(K[i0], K[i1]));
            }
          }
          // ADD THE INDICES
          new_element_idx.push_back(idx_array);
        }

        std::cout << " local element " << k << " is cut" << std::endl;
      }

      else {
        std::cout << " local element " << k << " is not cut" << std::endl;


      }

    }


  }



};

class Cut_LevelSet : public Cut_Base {
  double (*fun_levelSet)(const R2 P);

public:
  Cut_LevelSet(double (*ff)(const R2 P), int lab = 0) : Cut_Base(lab), fun_levelSet(ff) {}

private:
  double f(const R2 A) const {
    return fun_levelSet(A);
  }

};

class Cut_Line : public Cut_Base {

  // a side is the part that satisfies ax+by+c < 0
  double a, b, c;
public:
  Cut_Line(double aa, double bb, double cc, int lab = 0) : Cut_Base(lab), a(aa), b(bb), c(cc) {}

private:
  double f(const R2 A) const {
    return (a*A.x+b*A.y+c);
  }

};

class Cut_Parametric_Line : public Cut_Base {
  Interval t;
  Droite x;
  Droite y;

public:
  Cut_Parametric_Line(Interval tt, Droite xx, Droite yy, int lab = 0) : Cut_Base(lab), t(tt), x(xx), y(yy){}

};





class Fracture {

  typedef FaceInterface<2> Face;
  typedef SortArray<int, 3> TriaIdx;
  typedef typename Mesh2::Element Element;
  typedef typename Mesh2::RefPatch RefPatch;
  typedef typename Mesh2::Rd Rd;

  std::vector<const Cut_Base*> cut_lines_;

  const Mesh2& Th;
  std::vector<Face>        faces_;
  std::vector<Rd>       vertices_;
  std::vector<Rd> outward_normal_;


public:

  Fracture(const Mesh2 & MM) : Th(MM){
    faces_.resize( 0);
    vertices_.resize(0);
    outward_normal_.resize(0);
  }
  void add(Cut_Base& L1) {
    cut_lines_.push_back(&L1);
    // int nline = cut_lines_.size();
    build(L1);
  }


  Rd operator()(const int k, const int i) const {return vertices_[faces_[k][i]];}
  const Rd& operator()(const int i) const {return vertices_[CheckV(i)];}
  const Rd& get_node(const int i) const {return vertices_[CheckV(i)];}
  const Face& operator[](const int k) const {return faces_[CheckT(k)];}
  const Face& get_face(const int k)    const {return faces_[CheckT(k)];}
  int CheckV(int i) const { ASSERTION(i>=0 && i < vertices_.size()); return i;}
  int CheckT(int i) const { ASSERTION(i>=0 && i < faces_.size()); return i;}

  Uint nb_element  () const { return faces_.size(); }
  Rd normal(const int k) const { return outward_normal_[k];}

  bool is_cut_element(const int k) const {
    for(int i=0; i<cut_lines_.size(); ++i) {
      if(cut_lines_[i]->is_cut_element(k)) return true;
    }
    return false;
  }

  void build_local_partition(const int k, Local_Partition& local_partition) const {


    // GET THE ELEMENT TO BE CUT
    const Element& K(Th[k]);
    const int nv = K.nv;
    int iv[nv] = {0,1,2};

    // INITIALIZE THE LOCAL_PARTITION WITH K
    local_partition.reset();
    for(int i=0; i<3; ++i) {
      local_partition.add_node(K[i]);
    }
    local_partition.add_triangle(TriaIdx(iv));

    if(!is_cut_element(k)) return;

    // START CUTTING PROCEDURE
    typedef SortArray<int, nv> TriaIdx;
    std::list<int> erased_element;
    vector<TriaIdx> new_element_idx;

    std::cout << " GLOBAL ELEMENT \t" << k << std::endl;
    for(int i=0; i<cut_lines_.size(); ++i) {

      // FRACTURE i DOES NOT CUT THIS ELEMENT => NEXT
      if(!cut_lines_[i]->is_cut_element(k)) continue;
      std::cout << " is cut by fracture \t" << i << std::endl;
      // WE PERFORM THE CUT
      cut_lines_[i]->cut_partition(local_partition, new_element_idx, erased_element);

      // SORT THE LIST FROM HIGHER TO LOWER INDEX
      // TO NOT DESTROYED INDICES WHEN ERASING
      erased_element.sort(std::greater<int>());
      std::cout << " nb K to erase " << erased_element.size() << std::endl;
      // ERASING THE CUT ELEMENTS
      for(auto it = erased_element.begin(); it != erased_element.end(); ++it){
        std::cout << " erased element " << *it << std::endl;
        local_partition.erase_triangle(*it);
      }

      // CREATING THE NEW ELEMENTS
      for(auto it = new_element_idx.begin(); it != new_element_idx.end(); ++it){
        local_partition.add_triangle(*it);
      }
    }
  }



private:
  void build(Cut_Base& L1) {

    const Uint nv = Element::nv;
    double loc_ls[nv];
    byte   loc_ls_sign[nv];

    for (int k=0; k<Th.nbElmts(); k++) {
      const Element& K(Th[k]);

      // EVALUATE LOCAL VALUES AND SIGNS
      // std::vector<RemumberVertexPairT> zero_vertex_uses; // used to renumber the zero_vertexes
      for (Uint i= 0; i < nv; ++i) {
        loc_ls     [i] = L1.eval(K[i]);
        loc_ls_sign[i] = util::sign(loc_ls[i]);
      }

      // BUILD THE CUT PATCH FOR INTERFACE
      const RefPatch& cut =  RefPatch::instance( loc_ls_sign);


      for (typename RefPatch::const_face_iterator it= cut.face_begin(), end= cut.face_end();
      it != end; ++it) {
        // face_of_element_[k] = element_of_face_.size();
        // faces_.push_back( make_face(*it, K, loc_ls, zero_vertex_uses, label));

        L1.add_element_and_face(k, faces_.size());
        faces_.push_back( make_face(*it, K, loc_ls, L1.label));
        // element_of_face_.push_back(k);
        outward_normal_.push_back(make_normal(K, loc_ls));
      }

    }

  }
  Face make_face (const typename RefPatch::FaceIdx& ref_tri, const Element& K, const double lset[Element::nv], int label) {

    const Uint nv = Element::nv;
    const Uint nve = Element::nva;

    Uint loc_vert_num;
    Uint triIdx[nve];

    // int idxK = (*backMesh)(K);

    for (Uint j= 0; j < nve; ++j) {
      loc_vert_num= ref_tri[j];
      if (loc_vert_num < nv) {                            // zero vertex

        // const Uint idx = (*backMesh)(K[loc_vert_num]);
        // Rd Q = (*backMesh)(idx);
        // vertices_.push_back(Q);
        // triIdx[j] = vertices_.size() - 1;
        assert(0);
      }
      else { // genuine edge vertex

        const Ubyte i0 = Element::nvedge[loc_vert_num - nv][0],
        	i1 = Element::nvedge[loc_vert_num - nv][1];
        const double t = lset[i0]/(lset[i0] - lset[i1]);
        Rd Q = (1.0 - t) * ((Rd) K[i0]) + t * ((Rd) K[i1]); // linear interpolation
        vertices_.push_back(Q);
        triIdx[j] = vertices_.size() - 1;

        // edge_of_node_.push_back(loc_vert_num - K.nv);
      }
    }
    return Face(triIdx, label);
  }
  Rd make_normal (const Element& K, const double lset[Element::nv]) {

    Rd grad[Element::nv];
    K.Gradlambda(grad);
    Rd normal_ls;
    for(int i=0; i<Element::nv;++i) {

      normal_ls += grad[i]*lset[i];
    }
    normal_ls /= normal_ls.norm();
    return normal_ls;
  }

};





































#endif
