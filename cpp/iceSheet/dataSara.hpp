#ifndef _DATA_SARA_HPP
#define _DATA_SARA_HPP






class DataSara {


  string path_;


public:
  KN<double> ls_sign;
  KN<int> cut_element;
  KN<R> nodeInterfacex;
  KN<R> nodeInterfacey;


  DataSara (string pp) : path_(pp) {
    load_lsSign();
    load_cutElement();
    load_nodeInterface();
  }
  
  
  void load_lsSign() {
    string filename = path_+"lsSign.dat";
    ifstream f(filename.c_str());
    if(!f) {cerr << "DataSara: cannot open LsSign " << filename<<endl;exit(1);}
    cout << " Read On file \"" <<filename<<"\""<<  endl;
    int nv;
    f >> nv;
    ls_sign.resize(nv);
    for (int i=0;i<nv;i++) {
      f >> ls_sign(i);
      assert(f.good());
    }
  }


  
  void load_cutElement() {
   string filename = path_+"cutElement.dat";
   ifstream f(filename.c_str());
   if(!f) {cerr << "DataSara: cannot open cutElement " << filename<<endl;exit(1);}
    cout << " Read On file \"" <<filename<<"\""<<  endl;
    int nt;
    f >> nt;
    cut_element.resize(nt);
    for (int i=0;i<nt;i++) {
      f >> cut_element(i);
      assert(f.good());
    }
  }

  void load_nodeInterface() {
    string filename = path_+"nodeInterface.dat";
    ifstream f(filename.c_str());
    if(!f) {cerr << "DataSara: cannot open nodeInterface " << filename<<endl;exit(1);}
    cout << " Read On file \"" <<filename<<"\""<<  endl;
    int nvcut;
    f >> nvcut;
    nodeInterfacex.resize(2*nvcut);
    nodeInterfacey.resize(2*nvcut);
    for (int i=0;i<2*nvcut;i++) {
      f >> nodeInterfacex(i) >> nodeInterfacey(i);
      assert(f.good());
    }
  }



  









};






#endif
