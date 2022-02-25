#ifndef PARAVIEW_HPP
#define PARAVIEW_HPP

#include "../common/cut_method.hpp"
#include "FESpace.hpp"
#include "macroElement.hpp"
#include "../util/util.hpp"

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
  // typedef typename Mesh::Partition Partition;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename Element::Rd Rd;
  typedef Partition<Element> PartitionT;

  typedef FunFEM<Mesh> Fun_h;
  typedef ExpressionVirtual Expression;

  const int nv_cut_element = Rd::d +1;

  // Fun& u;
  const FESpace & Vh;
  int nbDataFile = 0;
  std::string outFile_;
  int ntCut = -1;
  int nt_cut, nt_notcut;
  Fun_h* levelSet = nullptr;

  struct ParaviewMesh {

    int ntCut_;
    int ntNotcut_;
    int nv_;

    int numCell_;
    int numCutCell_;

    int nvCutCell_;
    int nvCell_;

    std::map<int,int> element_status;
    std::vector<std::vector<Rd>> mesh_node;
    std::vector<int> idx_in_Vh;
    std::vector<std::pair<int,int>> num_cell; //(nb nodes, num_cell)

    ParaviewMesh() : numCell_(Element::ParaviewNumCell),
                     numCutCell_((Rd::d==2)?Triangle2::ParaviewNumCell : Tet::ParaviewNumCell),
                     nvCutCell_(Rd::d+1),
                     nvCell_(Element::nv)
                      {}


   void build(const FESpace & Vh, Fun_h* levelSet);
   void buildNoCut(const FESpace & Vh, Fun_h* levelSet) {
     ntCut_ = 0;
     ntNotcut_ = 0;
     nv_ = 0;
     mesh_node.resize(2*Vh.NbElement());
     idx_in_Vh.resize(2*Vh.NbElement());
     num_cell.resize(2*Vh.NbElement());
     std::vector<Rd> list_node;
     int kk=0;
     if(levelSet){
       int nn = 0;
       double loc_ls[Element::nv];

       for(int k=0; k<Vh.NbElement(); ++k) {

         const FElement& FK(Vh[k]);
         levelSet->eval(loc_ls, Vh.Th(FK.T));

         const int domain = FK.whichDomain();
         const int kb = Vh.idxElementInBackMesh(k);
         const PartitionT& cutK =  PartitionT(FK.T, loc_ls);
         ElementSignEnum the_part = cutK.what_part(domain);
         int ss = (domain == 0)? 1 : -1;
         if(the_part == AllElement){
           for(int sss=-1;sss<=1;sss+=2){
             cutK.get_list_node(list_node, sss);
             int nv_loc = list_node.size();
             num_cell[kk] = make_pair(nv_loc, 7);
             idx_in_Vh[kk] = k;
             nv_+= nv_loc;
             for(int i=0;i<nv_loc;++i){
               mesh_node[kk].push_back(list_node[i]);
             }
             kk++;
             ntCut_++;
           }
         }
         else{
           cutK.get_list_node(list_node, ss);
           int nv_loc = list_node.size();
           num_cell[kk] = make_pair(nv_loc, 7);
           idx_in_Vh[kk] = k;
           nv_+= nv_loc;
           for(int i=0;i<nv_loc;++i){
             mesh_node[kk].push_back(list_node[i]);
           }
           kk++;
           ntNotcut_++;
         }
       }
     }
     else {

       for(int k=0; k<Vh.NbElement(); ++k) {
         idx_in_Vh[k] = k;
         num_cell[k] = make_pair(nvCell_, numCell_);
         for(int i=0;i<nvCell_;++i) {
           mesh_node[kk].push_back(Vh.Th[k][i]);
         }
         kk++;
         nv_+= nvCell_;
       }
       ntNotcut_ = Vh.NbElement();
       ntCut_ = 0;
     }

     mesh_node.resize(kk+1);
     mesh_node.shrink_to_fit();
     idx_in_Vh.resize(kk+1);
     idx_in_Vh.shrink_to_fit();
     num_cell.resize(kk+1);
     num_cell.shrink_to_fit();
   }
   void buildCut(const FESpace & Vh, Fun_h* levelSet) {
     ntCut_ = 0;
     ntNotcut_ = 0;
     nv_ = 0;
     // std::vector<Rd> list_node;
     int kk=0;
     if(levelSet){
       mesh_node.resize(20*Vh.NbElement());
       idx_in_Vh.resize(20*Vh.NbElement());
       num_cell.resize(20*Vh.NbElement());
       // int nn = 0;
       double loc_ls[Element::nv];

       for(int k=0; k<Vh.NbElement(); ++k) {
         const FElement& FK(Vh[k]);
         levelSet->eval(loc_ls, Vh.Th(FK.T));

         const int domain = FK.whichDomain();
         const int kb = Vh.idxElementInBackMesh(k);
         const PartitionT& cutK =  PartitionT(FK.T, loc_ls);
         ElementSignEnum the_part = cutK.what_part(domain);

         if(the_part == NoElement) continue;

         if(cutK.is_cut()) {
           int iii=0;
           for(typename PartitionT::const_element_iterator it = cutK.element_begin(the_part);
           it != cutK.element_end(the_part); ++it){
             std::cout << iii++ << std::endl;
             std::cout << kk << "\t" << mesh_node.size() << std::endl;
             for(int i=0;i< nvCutCell_;++i) {
               Rd x(cutK.get_Vertex(it,i));
               std::cout << i  << "\t" << x << std::endl;
               mesh_node[kk].push_back(x);
             }
             num_cell[kk] = make_pair(nvCutCell_, numCutCell_);
             idx_in_Vh[kk] = k;
             nv_+= nvCutCell_;
             kk++;
             ntCut_++;
           }
         }

     //     ElementSignEnum the_part = cutK.what_part(domain);
     //     int ss = (domain == 0)? 1 : -1;
     //     if(the_part == AllElement){
     //       for(int sss=-1;sss<=1;sss+=2){
     //         cutK.get_list_node(list_node, sss);
     //         int nv_loc = list_node.size();
     //         num_cell[kk] = make_pair(nv_loc, 7);
     //         idx_in_Vh[kk] = k;
     //         nv_+= nv_loc;
     //         for(int i=0;i<nv_loc;++i){
     //           mesh_node[kk].push_back(list_node[i]);
     //         }
     //         kk++;
     //         ntCut_++;
     //       }
     //     }
     //     else{
     //       cutK.get_list_node(list_node, ss);
     //       int nv_loc = list_node.size();
     //       num_cell[kk] = make_pair(nv_loc, 7);
     //       idx_in_Vh[kk] = k;
     //       nv_+= nv_loc;
     //       for(int i=0;i<nv_loc;++i){
     //         mesh_node[kk].push_back(list_node[i]);
     //       }
     //       kk++;
     //       ntNotcut_++;
     //     }
       }


     //
     // mesh_node.resize(kk+1);
     // mesh_node.shrink_to_fit();
     // idx_in_Vh.resize(kk+1);
     // idx_in_Vh.shrink_to_fit();
     // num_cell.resize(kk+1);
     // num_cell.shrink_to_fit();

     // getchar();
     }
     else {
       // mesh_node.resize(Vh.NbElement());
       // idx_in_Vh.resize(Vh.NbElement());
       // num_cell.resize(Vh.NbElement());
       //
       // for(int k=0; k<Vh.NbElement(); ++k) {
       //   idx_in_Vh[k] = k;
       //   num_cell[k] = make_pair(nvCell_, numCell_);
       //   for(int i=0;i<nvCell_;++i) {
       //     mesh_node[kk].push_back(Vh.Th[k][i]);
       //   }
       //   kk++;
       //   nv_+= nvCell_;
       // }
       // ntNotcut_ = Vh.NbElement();
       // ntCut_ = 0;
     }

   }

   bool isCut(int k) const {
      element_status.find(k);
     return ( element_status.find(k) != element_status.end() );
   }
   int numCell(int k) const{
     if( isCut(k) ) return numCutCell_;
     else return numCell_;
   }
   int nvCell(int k) const {
     if(isCut(k)) return nvCutCell_;
     else return nvCell_;
   }
   int nbElement() const {return ntCut_+ ntNotcut_;}
   int nbNode() const {return nv_;}

   int sizeDataCell() const {return this->nbNode() + this->nbElement() ;}
   Rd node(int k, int i ) const {return mesh_node[k][i];}

  } mesh_data;


  Paraview(const FESpace&, Fun_h* ls, std::string name);
  Paraview(const FESpace&, Fun_h& ls, std::string name);

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


  // int get_ntCut(Interface2& interface, Interface2& marker);
  // virtual int numCutCell() { return (Rd::d==2)?Triangle2::ParaviewNumCell : Tet::ParaviewNumCell; };
  // virtual int numCell() { return Element::ParaviewNumCell; };
  // virtual int nbOfNode() { return Element::nv; };



};


template<class M>
Paraview<M>::Paraview(const GFESpace<M> & vh, Fun_h* ls, std::string name) :  Vh(vh), levelSet(ls)  {
  outFile_ = name;
  mesh_data.build(Vh, levelSet);
  this->writeFileMesh();
  this->writeFileCell();
}
template<class M>
Paraview<M>::Paraview(const GFESpace<M> & vh, Fun_h& ls, std::string name) :  Vh(vh), levelSet(&ls)  {
  outFile_ = name;
  mesh_data.build(Vh, levelSet);
  this->writeFileMesh();
  this->writeFileCell();
}


class Paraview2 : public Paraview<Mesh2> {

public:
  typedef GFESpace<Mesh2> FESpace;
  typedef Paraview<Mesh2> Paraview;

  // const int numCell_ = 5;
  // const int nbOfNode_ = 3;
  // int idxVTK_[6] = {0,1,2,3,4,5};

  Paraview2(const FESpace & vh, std::string n= "my_output2.vtk") : Paraview(vh,nullptr,n) {
    writeFileMesh();
    writeFileCell();
  }
  Paraview2(const FESpace & vh, Fun_h& ls , std::string n= "my_output2.vtk") : Paraview(vh,&ls,n) {
    writeFileMesh();
    writeFileCell();
  }

};
class Paraview3 : public Paraview<Mesh3> {

public :
typedef GFESpace<Mesh3> FESpace;
typedef Paraview<Mesh3> Paraview;

  // const int numCell_ = 10;
  // const int nbOfNode_ = 4;
  // int idxVTK_[10] = {0,1,2,3,4,7,5,6,8,9};

  Paraview3(const FESpace & vh, std::string n= "my_output3.vtk") : Paraview(vh,nullptr,n) {
    writeFileMesh();
    writeFileCell();
}
Paraview3(const FESpace & vh, Fun_h& ls, std::string n= "my_output3.vtk") : Paraview(vh,&ls,n) {
  writeFileMesh();
  writeFileCell();
}

  // int numCell() {return numCell_;}
  // int nbOfNode() {return nbOfNode_;}
  // int idxVTK(int i) {return idxVTK_[i];}


};

// Writting the nodes
// ---------------------------------------------------------------------------------------
template<class M>
void Paraview<M>::writeFileMesh() {

  std::ofstream point(outFile_.c_str(), std::ofstream::out);
  point << "# vtk DataFile Version 1.0" << std::endl
	<< "unstructured Grid" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET UNSTRUCTURED_GRID" << std::endl
  << "POINTS " << mesh_data.nbNode()  << " float " << std::endl;


  for(int k=0;k<mesh_data.nbElement();++k)  {
    for(int i=0;i<mesh_data.mesh_node[k].size();++i){
      if(Rd::d == 3){
        point << mesh_data.node(k,i) << std::endl;
      }
      else {
        point << mesh_data.node(k,i) << "\t 0.0" <<  std::endl;
      }
    }
  }

  point.close();
}


// Writting the cells
// ---------------------------------------------------------------------------------------
template<class M>
void Paraview<M>::writeFileCell() {

  std::ofstream cell(outFile_.c_str(), std::ofstream::out | std::ofstream::app);
  cell << "CELLS " << mesh_data.nbElement() << " " << mesh_data.sizeDataCell()<< std::endl;
  int l=0;
  for(int k=0;k<mesh_data.nbElement();++k)  {
    int nvc = mesh_data.num_cell[k].first;
    cell << nvc;
    for(int i=0;i<nvc ;++i,++l) cell << "  " << l;
    cell << std::endl;
  }

  cell << "CELL_TYPES " << mesh_data.nbElement() << std::endl;
  for(int k=0;k<mesh_data.nbElement();++k)  {
    cell << mesh_data.num_cell[k].second << std::endl;
  }

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
  if(nbDataFile == 0) data << "POINT_DATA " << mesh_data.nbNode() << std::endl;
  data << "SCALARS "+ name +" float" << std::endl;
  data << "LOOKUP_TABLE default" << std::endl;


  for(int k=0;k<mesh_data.nbElement();++k)  {
    int kvh = mesh_data.idx_in_Vh[k];
    int kb = Vh.idxElementInBackMesh(kvh);
    int domain = (levelSet)?Vh[kvh].whichDomain() : -1;

    for(int i=0;i<mesh_data.mesh_node[k].size();++i){

      R1 val = fh.evalOnBackMesh(kb, domain, mesh_data.node(k,i));
      data << paraviewFormat(val) << std::endl;;

    }
  }
  data.close();
  nbDataFile++;
}


template<class M>
void Paraview<M>::writeFileVectorData(Fun_h& fh,int c0, std::string name){

  std::ofstream data(outFile_.c_str(), std::ofstream::out | std::ofstream::app);

  if(nbDataFile == 0) data << "POINT_DATA " << mesh_data.nbNode()
			   << std::endl;
  data << "VECTORS "+ name +" float" << std::endl;


  for(int k=0;k<mesh_data.nbElement();++k)  {
    int kvh = mesh_data.idx_in_Vh[k];
    int kb = Vh.idxElementInBackMesh(kvh);
    int domain = (levelSet)?Vh[kvh].whichDomain() : -1;
    int kf = fh.idxElementFromBackMesh(kb, domain);

    for(int i=0;i<mesh_data.mesh_node[k].size();++i){

      for(int  dd=0;dd<Rd::d;++dd){
        // R1 val = fh.evalOnBackMesh(kb, domain, mesh_data.node(k,i), c0+dd);
        R1 val = fh.eval(kf, mesh_data.node(k,i), c0+dd);
        data << paraviewFormat(val) <<  "\t";
      }
      if (Rd::d==2) data << "0.0";
      data << std::endl;
    }
  }
  // if(levelSet) {
  //   double loc_ls[Element::nv];
  //   for(int k=0; k<Vh.NbElement(); ++k) {
  //
  //     const FElement& FK(Vh[k]);
  //     int kback = Vh.idxElementInBackMesh(k);
  //
  //     levelSet->eval(loc_ls, Vh.Th(FK.T));
  //     const int domain = FK.whichDomain();
  //     const PartitionT& cutK =  PartitionT(FK.T, loc_ls);
  //     ElementSignEnum the_part = cutK.what_part(domain);
  //
  //     int kf = fh.idxElementFromBackMesh(kback, domain);
  //     if(the_part == NoElement) continue;
  //
  //     if(cutK.is_cut()) {
  //       for(typename PartitionT::const_element_iterator it = cutK.element_begin(the_part);
  //       it != cutK.element_end(the_part); ++it){
  //         for(int i=0;i<Element::nv;++i) {
  //           Rd x(cutK.get_Vertex(it,i));
  //           for(int  dd=0;dd<Rd::d;++dd){
  //             R val = fh.eval(kf, x, c0+dd);
  //             data << paraviewFormat(val) << "\t" ;
  //           }
  //           if (Rd::d==2) data << "0.0";
  //           data << std::endl;
  //         }
  //       }
  //     }
  //     else
  //     for(int i=0;i<Element::nv;++i) {
  //       for(int  dd=0;dd<Rd::d;++dd){
  //         R val = fh.eval(kf, FK.T[i], c0+dd);
  //         data << val << "\t" ;
  //       }
  //       if (Rd::d==2) data << "0.0";
  //       data << std::endl;
  //     }
  //   }
  // }
  // else{
  //
  //   for(int k=0;k<Vh.NbElement();++k) {
  //     const FElement& FK(Vh[k]);
  //
  //     int kback = Vh.idxElementInBackMesh(k);
  //     int kf = fh.idxElementFromBackMesh(kback);
  //
  //     for(int i=0;i<Element::nv;++i){
  //       for(int  dd=0;dd<Rd::d;++dd){
  //         R val = fh.eval(kf, FK.T[i], c0+dd);
  //         data << paraviewFormat(val) << "\t" ;
  //       }
  //       if (Rd::d==2) data << "0.0";
  //       data << std::endl;
  //     }
  //   }
  // }
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
      const PartitionT& cutK =  PartitionT(FK.T, loc_ls);
      ElementSignEnum the_part = cutK.what_part(domain);


      int kf = fh.idxElementFromBackMesh(kback, domain);
      if(!macro.isRootFat(k)) {
        kf = macro.getIndexRootElement(k);
      }

      if(the_part == NoElement) continue;

      if(cutK.is_cut()) {
        for(typename PartitionT::const_element_iterator it = cutK.element_begin(the_part);
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
      const PartitionT& cutK =  PartitionT(FK.T, loc_ls);
      ElementSignEnum the_part = cutK.what_part(domain);
      int kback = Vh.idxElementInBackMesh(k);

      if(!macro.isRootFat(k)) {
        int kk = macro.getIndexRootElement(k);
        kback = Vh.idxElementInBackMesh(kk);
      }

      if(the_part == NoElement) continue;

      if(cutK.is_cut()) {
        for(typename PartitionT::const_element_iterator it = cutK.element_begin(the_part);
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
