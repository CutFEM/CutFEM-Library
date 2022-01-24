#ifndef PARAVIEW_HPP
#define PARAVIEW_HPP

#include "../FESpace/FESpace.hpp"
#include "../util/util.hpp"
#include "../FESpace/macroElement.hpp"

// template<class F>
// class FEMFunction;
//
template<class F>
class FunFEM;

template<class F>
class FEMTimeFunction;

static double paraviewFormat( double x) {
  return (fabs(x) < 1e-20)? 0.: x;
}

/*
 *   Only write P1 solution on the mesh provided
 *
 */
template<class M>
class Paraview {
public :
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename Mesh::Partition Partition;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename Element::Rd Rd;
  typedef FunFEM<Mesh> Fun_h;
  typedef ExpressionVirtual Expression;


  // Fun& u;
  const FESpace & Vh;
  int nbDataFile = 0;
  std::string outFile_;
  int ntCut = -1;
  Fun_h* levelSet = nullptr;


  Paraview(const FESpace&, Fun_h* ls, std::string name);
  // Paraview(const FESpace&, Fun_h &, std::string name = "my_output.vtk");


  void writeFileMesh();
  void writeFileCell();

  void add(Fun_h&, std::string, int, int);
  void add(const ExpressionVirtual& fh, std::string);
  void writeFileScalarData(const ExpressionVirtual&, std::string);
  void writeFileVectorData(Fun_h&, int, std::string);

  void add(const ExpressionVirtual& fh, std::string, MacroElement& macro);
  void add(Fun_h&, std::string, int, int, MacroElement& macro);
  void writeFileVectorData(Fun_h&, int, std::string, MacroElement& macro);
  void writeFileScalarData(const ExpressionVirtual&, std::string, MacroElement& macro);


  int get_ntCut();
  // int get_ntCut(Interface2& interface, Interface2& marker);
  virtual int numCell() =0;
  virtual int nbOfNode() =0;
  virtual int idxVTK(int i) = 0;

};


template<class M>
Paraview<M>::Paraview(const GFESpace<M> & vh, Fun_h* ls, std::string name) :  Vh(vh), levelSet(ls)  {
  outFile_ = name;
}
// template<class M>
// Paraview<M>::Paraview(const GFESpace<M> & vh, std::string name) :  Vh(vh), levelSet(ls) {
//   outFile_ = name;
//   get_ntCut();
//   cutProblem = true;
// }

class Paraview2 : public Paraview<Mesh2> {

public:
  typedef GFESpace<Mesh2> FESpace;
  typedef Paraview<Mesh2> Paraview;

  const int numCell_ = 5;
  const int nbOfNode_ = 3;
  int idxVTK_[6] = {0,1,2,3,4,5};

  Paraview2(const FESpace & vh, std::string n= "my_output2.vtk") : Paraview(vh,nullptr,n) {
    writeFileMesh();
    writeFileCell();
  }
  Paraview2(const FESpace & vh, Fun_h& ls , std::string n= "my_output2.vtk") : Paraview(vh,&ls,n) {
    writeFileMesh();
    writeFileCell();
  }

  int numCell() {return numCell_;}
  int nbOfNode() {return nbOfNode_;}
  int idxVTK(int i) {return idxVTK_[i];}

};


class Paraview3 : public Paraview<Mesh3> {

public :
typedef GFESpace<Mesh3> FESpace;
typedef Paraview<Mesh3> Paraview;

  const int numCell_ = 10;
  const int nbOfNode_ = 4;
  int idxVTK_[10] = {0,1,2,3,4,7,5,6,8,9};

  Paraview3(const FESpace & vh, std::string n= "my_output3.vtk") : Paraview(vh,nullptr,n) {
    writeFileMesh();
    writeFileCell();
}
Paraview3(const FESpace & vh, Fun_h& ls, std::string n= "my_output3.vtk") : Paraview(vh,&ls,n) {
  writeFileMesh();
  writeFileCell();
}

  int numCell() {return numCell_;}
  int nbOfNode() {return nbOfNode_;}
  int idxVTK(int i) {return idxVTK_[i];}


};

template<class M>
int Paraview<M>::get_ntCut() {
  int nn = 0;
  double loc_ls[Element::nv];
  for(int k=0; k<Vh.NbElement(); ++k) {

    const FElement& FK(Vh[k]);
    levelSet->eval(loc_ls, Vh.Th(FK.T));

    const int domain = FK.whichDomain();
    const Partition& cutK =  Partition(FK.T, loc_ls);
    ElementSignEnum the_part = cutK.what_part(domain);

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
  	it != cutK.element_end(the_part); ++it){
      nn++;
    }
  }
  return nn;
}

// template<class M>
// int Paraview<M>::get_ntCut(Interface2& interface, Interface2& marker) {
//   int nn = 0;
//   double loc_ls[Element::nv];
//   for(int k=0; k<Vh.NbElement(); ++k) {
//
//     const FElement& FK(Vh[k]);
//     const int kb = Vh.Th(FK.T);
//
//     levelSet->eval(loc_ls, Vh.Th(FK.T));
//
//     const int domain = FK.whichDomain();
//     CutData2 cutData(interface.getCutData(kb));     // get the cut data
//     const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
//     ElementSignEnum the_part = cutK.what_part(domain);
//
//
//     for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
//   	it != cutK.element_end(the_part); ++it){
//       nn++;
//     }
//   }
//   return nn;
// }



// Writting the nodes
// ---------------------------------------------------------------------------------------
template<class M>
void Paraview<M>::writeFileMesh() {
  ntCut = (levelSet)? get_ntCut() : Vh.NbElement();
  double loc_ls[Element::nv];

  std::ofstream point(outFile_.c_str(), std::ofstream::out);
  point << "# vtk DataFile Version 1.0" << std::endl
	<< "unstructured Grid" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET UNSTRUCTURED_GRID" << std::endl
	<< "POINTS " << ntCut * Element::nv  << " float " << std::endl;


  for(int k=0;k<Vh.NbElement();++k){

    if(levelSet){
      const FElement& FK(Vh[k]);
      levelSet->eval(loc_ls, Vh.Th(FK.T));

      const int domain = FK.whichDomain();
      const Partition& cutK =  Partition(FK.T, loc_ls);
      ElementSignEnum the_part = cutK.what_part(domain);

      for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
      it != cutK.element_end(the_part); ++it){

        for(int i=0;i<Element::nv;++i) {
          if(Rd::d == 2)point << cutK.get_Vertex(it,i) << "\t" << 0.0 << std::endl;
          else point << cutK.get_Vertex(it,i) << std::endl;
        }
      }
    }
    else{
      if(Rd::d == 3)
         for(int i=0;i<Element::nv;++i) point << (Rd) Vh.Th[k][i] << std::endl;
      else
         for(int i=0;i<Element::nv;++i) point << (Rd) Vh.Th[k][i] << "\t 0.0" <<  std::endl;
    }
  }
  point.close();
}


// Writting the cells
// ---------------------------------------------------------------------------------------
template<class M>
void Paraview<M>::writeFileCell() {

  std::ofstream cell(outFile_.c_str(), std::ofstream::out | std::ofstream::app);
  cell << "CELLS " << ntCut << " " << ntCut*(nbOfNode() + 1) << std::endl;
  int l=0;
  for(int k=0;k<ntCut;++k) {
    const Element& K(Vh.Th[k]);
    cell << nbOfNode();
    for(int i=0;i<Element::nv;++i,++l) cell << "  " << l;//mesh(K[idxVTK(i)]);
    cell << std::endl;
  }

  cell << "CELL_TYPES " << ntCut << std::endl;
  for(int k=0;k<ntCut;++k)  cell << numCell() << std::endl;
  cell.close();
}


template<class M>
void Paraview<M>::add(Fun_h& fh, std::string nameField, int begin_comp, int nb_comp){

  if(nb_comp ==1) {
    ExpressionFunFEM<M> ui(fh, begin_comp, op_id);
    writeFileScalarData(ui, nameField);
  }
  else
    writeFileVectorData(fh, begin_comp, nameField);
}


template<class M>
void Paraview<M>::add(const ExpressionVirtual& fh, std::string nameField){

    writeFileScalarData(fh, nameField);
}


template<class M>
void Paraview<M>::writeFileScalarData(const ExpressionVirtual& fh, std::string name){

  std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);
  if(nbDataFile == 0) data << "POINT_DATA " << ntCut*Element::nv << std::endl;
  data << "SCALARS "+ name +" float" << std::endl;
  data << "LOOKUP_TABLE default" << std::endl;

  double loc_ls[Element::nv];

  if(levelSet){
    for(int k=0; k<Vh.NbElement(); ++k) {
      const FElement& FK(Vh[k]);
      levelSet->eval(loc_ls, Vh.Th(FK.T));

      const int domain = FK.whichDomain();
      const Partition& cutK =  Partition(FK.T, loc_ls);
      ElementSignEnum the_part = cutK.what_part(domain);
      int kback = Vh.idxElementInBackMesh(k);

      if(the_part == NoElement) continue;

      if(cutK.is_cut()) {
        for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
        it != cutK.element_end(the_part); ++it){
          for(int i=0;i<Element::nv;++i) {
            Rd x(cutK.get_Vertex(it,i));
            R1 val = fh.evalOnBackMesh(kback, domain, x);
            data << paraviewFormat(val) << std::endl;;
          }
        }
      }
      else{
        for(int i=0;i<Element::nv;++i) {
          R1 val = fh.evalOnBackMesh(kback, domain, FK.T[i]);
          data << paraviewFormat(val.x) << std::endl;
        }
      }
    }
  }
  else {
    for(int k=0;k<Vh.NbElement();++k) {
      const FElement& FK(Vh[k]);
      int kback = Vh.idxElementInBackMesh(k);
      for(int i=0;i<Element::nv;++i){
        R1 val = fh.evalOnBackMesh(kback, -1, FK.T[i]);
        data << paraviewFormat(val.x) << std::endl;
      }
    }
  }

  data.close();
  nbDataFile++;
}


template<class M>
void Paraview<M>::writeFileVectorData(Fun_h& fh,int c0, std::string name){

  std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);

  if(nbDataFile == 0) data << "POINT_DATA " << ntCut*Element::nv
			   << std::endl;
  data << "VECTORS "+ name +" float" << std::endl;

  if(levelSet) {
    double loc_ls[Element::nv];
    for(int k=0; k<Vh.NbElement(); ++k) {

      const FElement& FK(Vh[k]);
      int kback = Vh.idxElementInBackMesh(k);

      levelSet->eval(loc_ls, Vh.Th(FK.T));
      const int domain = FK.whichDomain();
      const Partition& cutK =  Partition(FK.T, loc_ls);
      ElementSignEnum the_part = cutK.what_part(domain);

      int kf = fh.idxElementFromBackMesh(kback, domain);
      if(the_part == NoElement) continue;

      if(cutK.is_cut()) {
        for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
        it != cutK.element_end(the_part); ++it){
          for(int i=0;i<Element::nv;++i) {
            Rd x(cutK.get_Vertex(it,i));
            for(int  dd=0;dd<Rd::d;++dd){
              R val = fh.eval(kf, x, c0+dd);
              data << paraviewFormat(val) << "\t" ;
            }
            if (Rd::d==2) data << "0.0";
            data << std::endl;
          }
        }
      }
      else
      for(int i=0;i<Element::nv;++i) {
        for(int  dd=0;dd<Rd::d;++dd){
          R val = fh.eval(kf, FK.T[i], c0+dd);
          data << val << "\t" ;
        }
        if (Rd::d==2) data << "0.0";
        data << std::endl;
      }
    }
  }
  else{

    for(int k=0;k<Vh.NbElement();++k) {
      const FElement& FK(Vh[k]);

      int kback = Vh.idxElementInBackMesh(k);
      int kf = fh.idxElementFromBackMesh(kback);

      for(int i=0;i<Element::nv;++i){
        for(int  dd=0;dd<Rd::d;++dd){
          R val = fh.eval(kf, FK.T[i], c0+dd);
          data << paraviewFormat(val) << "\t" ;
        }
        if (Rd::d==2) data << "0.0";
        data << std::endl;
      }
    }
  }
  data.close();
  nbDataFile++;
}



template<class M>
void Paraview<M>::add(Fun_h& fh, std::string nameField, int begin_comp, int nb_comp, MacroElement& macro){

  if(nb_comp ==1) {
    ExpressionFunFEM<M> ui(fh, begin_comp, op_id);
    writeFileScalarData(ui, nameField, macro);
  }
  else
    writeFileVectorData(fh, begin_comp, nameField, macro);
}


template<class M>
void Paraview<M>::add(const ExpressionVirtual& fh, std::string nameField, MacroElement& macro){

    writeFileScalarData(fh, nameField, macro);
}


template<class M>
void Paraview<M>::writeFileVectorData(Fun_h& fh,int c0, std::string name, MacroElement& macro){

  std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);

  if(nbDataFile == 0) data << "POINT_DATA " << ntCut*Element::nv
			   << std::endl;
  data << "VECTORS "+ name +" float" << std::endl;

  if(levelSet) {
    double loc_ls[Element::nv];
    for(int k=0; k<Vh.NbElement(); ++k) {

      const FElement& FK(Vh[k]);
      int kback = Vh.idxElementInBackMesh(k);

      levelSet->eval(loc_ls, Vh.Th(FK.T));
      const int domain = FK.whichDomain();
      const Partition& cutK =  Partition(FK.T, loc_ls);
      ElementSignEnum the_part = cutK.what_part(domain);


      int kf = fh.idxElementFromBackMesh(kback, domain);
      if(!macro.isRootFat(k)) {
        kf = macro.getIndexRootElement(k);
      }

      if(the_part == NoElement) continue;

      if(cutK.is_cut()) {
        for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
        it != cutK.element_end(the_part); ++it){
          for(int i=0;i<Element::nv;++i) {
            Rd x(cutK.get_Vertex(it,i));
            for(int  dd=0;dd<Rd::d;++dd){
              R val = fh.eval(kf, x, c0+dd);
              data << paraviewFormat(val) << "\t" ;
            }
            if (Rd::d==2) data << "0.0";
            data << std::endl;
          }
        }
      }
      else
      for(int i=0;i<Element::nv;++i) {
        for(int  dd=0;dd<Rd::d;++dd){
          R val = fh.eval(kf, FK.T[i], c0+dd);
          data << paraviewFormat(val) << "\t" ;
        }
        if (Rd::d==2) data << "0.0";
        data << std::endl;
      }
    }
  }
  else{

    for(int k=0;k<Vh.NbElement();++k) {
      const FElement& FK(Vh[k]);

      int kback = Vh.idxElementInBackMesh(k);
      int kf = fh.idxElementFromBackMesh(kback);

      for(int i=0;i<Element::nv;++i){
        for(int  dd=0;dd<Rd::d;++dd){
          R val = fh.eval(kf, FK.T[i], c0+dd);
          data << paraviewFormat(val) << "\t" ;
        }
        if (Rd::d==2) data << "0.0";
        data << std::endl;
      }
    }
  }
  data.close();
  nbDataFile++;
}



template<class M>
void Paraview<M>::writeFileScalarData(const ExpressionVirtual& fh, std::string name, MacroElement& macro){

  std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);
  if(nbDataFile == 0) data << "POINT_DATA " << ntCut*Element::nv << std::endl;
  data << "SCALARS "+ name +" float" << std::endl;
  data << "LOOKUP_TABLE default" << std::endl;

  double loc_ls[Element::nv];

  if(levelSet){
    for(int k=0; k<Vh.NbElement(); ++k) {
      const FElement& FK(Vh[k]);
      levelSet->eval(loc_ls, Vh.Th(FK.T));

      const int domain = FK.whichDomain();
      const Partition& cutK =  Partition(FK.T, loc_ls);
      ElementSignEnum the_part = cutK.what_part(domain);
      int kback = Vh.idxElementInBackMesh(k);

      if(!macro.isRootFat(k)) {
        int kk = macro.getIndexRootElement(k);
        kback = Vh.idxElementInBackMesh(kk);
      }

      if(the_part == NoElement) continue;

      if(cutK.is_cut()) {
        for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
        it != cutK.element_end(the_part); ++it){
          for(int i=0;i<Element::nv;++i) {
            Rd x(cutK.get_Vertex(it,i));
            R1 val = fh.evalOnBackMesh(kback, domain, x);
            // if(!macro.isRootFat(k)) val = (2*domain-1)*1e7;
            data << paraviewFormat(val) << std::endl;;
          }
        }
      }
      else{
        for(int i=0;i<Element::nv;++i) {
          R1 val = fh.evalOnBackMesh(kback, domain, FK.T[i]);
          data << paraviewFormat(val.x) << std::endl;
        }
      }
    }
  }
  else {
    for(int k=0;k<Vh.NbElement();++k) {
      const FElement& FK(Vh[k]);
      int kback = Vh.idxElementInBackMesh(k);
      for(int i=0;i<Element::nv;++i){
        R1 val = fh.evalOnBackMesh(kback, -1, FK.T[i]);
        data << paraviewFormat(val.x) << std::endl;
      }
    }
  }

  data.close();
  nbDataFile++;
}


#endif
