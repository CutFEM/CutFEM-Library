#include "spline.hpp"


Spline::Spline(const vector<double> &x, const vector<double> & y) {
  this->init(x, y);
}
// {
//   X = x;
//   Y = y;
//
//   int n = x.size()-1;
//   vector<double> a;
//   a.insert(a.begin(), y.begin(), y.end());
//   vector<double> b(n);
//   vector<double> d(n);
//   vector<double> h;
//
//   for(int i = 0; i < n; ++i)
//   h.push_back(x[i+1]-x[i]);
//
//   vector<double> alpha;
//   for(int i = 0; i < n; ++i)
//   alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );
//
//   vector<double> c(n+1);
//   vector<double> l(n+1);
//   vector<double> mu(n+1);
//   vector<double> z(n+1);
//   l[0] = 1;
//   mu[0] = 0;
//   z[0] = 0;
//
//   for(int i = 1; i < n; ++i)
//   {
//     l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
//     mu[i] = h[i]/l[i];
//     z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
//   }
//
//   l[n] = 1;
//   z[n] = 0;
//   c[n] = 0;
//
//   for(int j = n-1; j >= 0; --j)
//   {
//     c[j] = z [j] - mu[j] * c[j+1];
//     b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
//     d[j] = (c[j+1]-c[j])/3/h[j];
//   }
//
//   mySplineSet = new SplineSet[n];
//   for(int i = 0; i < n; ++i)
//   {
//     mySplineSet[i].a = a[i];
//     mySplineSet[i].b = b[i];
//     mySplineSet[i].c = c[i];
//     mySplineSet[i].d = d[i];
//     mySplineSet[i].x = x[i];
//   }
//   return;
// }
void Spline::init(const vector<double> &x, const vector<double> & y) {
  X = x;
  Y = y;

  int n = x.size()-1;
  vector<double> a;
  a.insert(a.begin(), y.begin(), y.end());
  vector<double> b(n);
  vector<double> d(n);
  vector<double> h;

  for(int i = 0; i < n; ++i)
  h.push_back(x[i+1]-x[i]);

  vector<double> alpha;
  for(int i = 0; i < n; ++i)
  alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

  vector<double> c(n+1);
  vector<double> l(n+1);
  vector<double> mu(n+1);
  vector<double> z(n+1);
  l[0] = 1;
  mu[0] = 0;
  z[0] = 0;

  for(int i = 1; i < n; ++i)
  {
    l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
    mu[i] = h[i]/l[i];
    z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
  }

  l[n] = 1;
  z[n] = 0;
  c[n] = 0;

  for(int j = n-1; j >= 0; --j)
  {
    c[j] = z [j] - mu[j] * c[j+1];
    b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
    d[j] = (c[j+1]-c[j])/3/h[j];
  }

  mySplineSet = new SplineSet[n];
  for(int i = 0; i < n; ++i)
  {
    mySplineSet[i].a = a[i];
    mySplineSet[i].b = b[i];
    mySplineSet[i].c = c[i];
    mySplineSet[i].d = d[i];
    mySplineSet[i].x = x[i];
  }
  return;
}
double Spline::evaluate(const double x) const {
  int j;
  for ( j = 0; static_cast<unsigned>(j) < X.size()-1; j++ ){
    if( mySplineSet[j].x > x ){
      if( j == 0 ) j++;
      break;
    }
  }
  j--;
  double dx = x - mySplineSet[j].x;
  double y = mySplineSet[j].a + mySplineSet[j].b * dx + mySplineSet[j].c * dx* dx +
  mySplineSet[j].d * dx* dx * dx;
  return y;
}

double Spline::diff(const double x) const {
  int j;
  for ( j = 0; static_cast<unsigned>(j) < X.size()-1; j++ ){
    if( mySplineSet[j].x > x ){
      if( j == 0 ) j++;
      break;
    }
  }
  j--;
  double dx = x - mySplineSet[j].x;
  double y = mySplineSet[j].b + 2* mySplineSet[j].c * dx +
  3*mySplineSet[j].d * dx* dx;
  return y;
}

double Spline::ddiff(const double x) const {
  int j;
  for ( j = 0; static_cast<unsigned>(j) < X.size()-1; j++ ){
    if( mySplineSet[j].x > x ){
      if( j == 0 ) j++;
      break;
    }
  }
  j--;
  double dx = x - mySplineSet[j].x;
  double y = 2* mySplineSet[j].c + 6*mySplineSet[j].d * dx;
  return y;
}

int Spline::gnuplot(int N , string filename){
  ofstream plot, plot1;
  plot1.open("splineT.dat", std::ofstream::out);
  plot.open(filename.c_str());
  if(!plot) return 1;
  for (int j = 0; static_cast<unsigned>(j) < X.size()-1; j++ ){
    double h = (X[j+1] - X[j])/N;
    for(int i=0; i<=N;++i) {
      double x=mySplineSet[j].x+i*h;
      double dx=x - mySplineSet[j].x;
      double y=mySplineSet[j].a + mySplineSet[j].b * dx + mySplineSet[j].c * dx* dx +
      mySplineSet[j].d * dx* dx * dx;
      plot << x << "\t" << y << endl;
    }
    double x = mySplineSet[j].x;
    double y = mySplineSet[j].a;
    plot1 << x << "\t" << y << std::endl;
  }

  plot  << X[X.size()-1] << "\t" << Y[Y.size()-1] << endl;
  plot1 << X[X.size()-2] << "\t" << Y[Y.size()-2] << endl;
  plot1 << X[X.size()-1] << "\t" << Y[Y.size()-1] << endl;

  plot.close();
  return 0;
}


Spline2D::Spline2D( const vector<double> & t, const vector<R2> & P ){
  int n = P.size();
  vector<double> x(n);
  vector<double> y(n);
  for(int i=0;i<n;++i) {
    x[i] = P[i].x;
    y[i] = P[i].y;
    // std::cout << P[i] << std::endl;

  }
  std::cout << " hey" << std::endl;
  splineX.init(t,x);
  splineY.init(t,y);
}

R2 Spline2D::evaluate(double t) const {
  return R2(splineX(t), splineY(t));
}
int Spline2D::gnuplot(int N , string filename){
  ofstream plot, plot1;
  plot1.open("spline2T.dat", std::ofstream::out);
  plot.open(filename.c_str());
  if(!plot) return 1;
  for (int j = 0; static_cast<unsigned>(j) < nb_node - 1; j++ ){
    double h = (splineX.X[j+1] - splineX.X[j])/N;
    for(int i=0; i<=N;++i) {
      double t=splineX.X[j]+i*h;
      // double dx=x - mySplineSet[j].x;
      R2 val = evaluate(t);

      plot << val.x << "\t" << val.y << endl;
      if(i==0) plot1 << val.x << "\t" << val.y << endl;

    }

  }

  // plot  << X[X.size()-1] << "\t" << Y[Y.size()-1] << endl;
  // plot1 << X[X.size()-2] << "\t" << Y[Y.size()-2] << endl;
  // plot1 << X[X.size()-1] << "\t" << Y[Y.size()-1] << endl;

  plot.close();
  return 0;
}

void display_spline2d(const Spline & splineX, const Spline & splineY, string filename, int N){
  ofstream plot;
  plot.open(filename.c_str());
  if(!plot) cerr << "ERROR : " << filename << " cannot be opened" << endl;
  int size = splineX.X.size();
  for (int j = 0; j < size-1; j++ ){
    double h = (splineX.X[j+1] - splineX.X[j])/N;
    assert(h >0);
    for(int i=0; i<N;++i) {
      double x=splineX.X[j]+i*h;
      double dx=x - splineX.X[j];
      double X=splineX.mySplineSet[j].a + splineX.mySplineSet[j].b *dx
      + splineX.mySplineSet[j].c * dx* dx
      + splineX.mySplineSet[j].d * dx* dx * dx;
      double Y=splineY.mySplineSet[j].a + splineY.mySplineSet[j].b *dx
      + splineY.mySplineSet[j].c * dx* dx
      + splineY.mySplineSet[j].d * dx* dx * dx;
      plot << X << "\t" << Y << endl;

    }
  }
  plot.close();
}
