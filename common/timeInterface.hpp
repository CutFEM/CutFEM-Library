#ifndef TIME_INTERFACE_HPP
#define TIME_INTERFACE_HPP

// #include "Interface2dn.hpp"
// #include "Interface3dn.hpp"
//
// template<typename M> class FunFEM;
//
// template<typename M>
// class TimeInterface {
// public:
// 	typedef M Mesh;
// 	typedef FunFEM<Mesh> Fun_h;
// 	typedef typename TypeInterface<Mesh::Rd::d>::Interface Interface;
// private:
// 	vector<Interface*> interface;
// 	int n;
//
// public:
//
// 	TimeInterface(int nt) : interface(nt), n(nt) {}
//
// 	void init(int i, const Mesh & Th, const KN<double>& ls) {
// 		assert(0 <= i && i < n);
// 		if(interface[i]) {
// 			delete interface[i];
// 		}
// 		interface[i] = new Interface(Th,ls);
// 	}
//
// 	void init(const Mesh & Th, const vector<Fun_h>& ls) {
// 		for(int i=0;i<ls.size();++i){
// 			assert(0 <= i && i < n);
// 			if(interface[i]) {
// 				delete interface[i];
// 			}
// 			interface[i] = new Interface(Th,ls[i].v);
// 		}
// 	}
//
// 	 Interface* operator[](int i) const{
// 		 assert(0 <= i && i < n);
// 		 return interface[i];
// 	 }
// 	 Interface* operator()(int i) const{
// 		 assert(0 <= i && i < n);
// 		 return interface[i];
// 	 }
//
// 	 int size() const { return interface.size();}
//
// 	 ~TimeInterface(){
// 		 for(int i=0;i<n;++i){
// 			 if(interface[i]) delete interface[i];
// 		 }
// 	 }
//
// private:
// 	 TimeInterface(const TimeInterface&);
// 	 void operator=(const TimeInterface &);
// };
//
// class TimeInterface2 : public TimeInterface<Mesh2> {
// public :
// TimeInterface2(int nt) : TimeInterface<Mesh2>(nt) {}
// private:
// 	TimeInterface2(const TimeInterface2&);
// 	void operator=(const TimeInterface2 &);                    // no copy allowed
//
// };
//
// class TimeInterface3 : public TimeInterface<Mesh3> {
// public :
// TimeInterface3(int nt) : TimeInterface<Mesh3>(nt) {}
// private:
// TimeInterface3(const TimeInterface3&);
// void operator=(const TimeInterface3&);                    // no copy allowed
//
// };
//
//

#endif
