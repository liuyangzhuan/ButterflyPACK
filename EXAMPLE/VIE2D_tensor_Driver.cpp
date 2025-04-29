//------------------------------------------------------------------------------
#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <pthread.h>

#include <cmath>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include <atomic>
#include <mpi.h>
// #include <complex.h>

#include <sstream>
#include <cstring>
#include <getopt.h>
#include <unistd.h>

// #include "math.h"
// #include "stdlib.h"
// #include "stdio.h"
// #include "constants.h"
// #include "functions.h"
// #include "FSM5.h"
// #include "FeikoPQ5.h"
// #include "FeikonalPQ5.h"
// #include "AeikoPQ5.h"
// #include "AeikonalPQ5.h"
// #include "WENO3.h"
// #include "WENO5.h"
// #include "Feiko5_WENO3_LxF.h"
// #include "Feiko5_WENO5_LxF.h"
// #include "Feikonal5_WENO3_Sweeping.h"
// #include "Feikonal5_WENO5_Sweeping.h"
// #include "Extrapolation4.h"
// #include "Extrapolation5.h"
// #include "SmootherVel.h"
// #include "SmootherFunction.h"
// #include "SmootherFunction_x.h"
// #include "SmootherFunction2.h"
// #include "SmootherFunction2_x.h"
// #include "SmootherFunction2_y.h"
// #include "SmootherFunction2_xx.h"
// #include "SmootherFunction2_yy.h"
// #include "time.h"
// #include "FSM5_WENO3_LxF.h"
// #include "FSM5_WENO5_LxF.h"
// #include "FSM5_Feikonal5_WENO3_Sweeping.h"
// #include "FSM5_Feikonal5_WENO5_Sweeping.h"
// #include "SmootherConnection.h"
// #include "SmootherConnection_x.h"
// #include "SmootherConnection_y.h"
// #include "PolyConnection.h"
// #include "PolyConnection_x.h"
// #include "PolyConnection_y.h"
// #include "FSM5_Local_Solver.h"
// #include "FSM5_Sweeping.h"
// #include "FSM5_FeikonalPQ5.h"
// #include "FSM5_AeikonalPQ5.h"
// #include "Amplitude2_Babich.h"
// #include "Amplitude2_WENO3_Babich.h"
// #include "Amplitude2_Babich_Term2.h"
// #include "Amplitude2_WENO3_Babich_Term2.h"
// #include "Amplitude2_Babich_square.h"
// #include "Amplitude2_WENO3_Babich_square.h"
// #include "Amplitude2_WENO5_Babich_square.h"
// #include "Amplitude2_Babich_Term2_square.h"
// #include "Amplitude2_WENO3_Babich_Term2_square.h"
// #include "Gradient.h"
// #include "Laplacian.h"
#include "dC_BPACK_wrapper.h"
// #include "DoCubicInterp2D.h"


#include "G2D/bessel.h"
#include "zC_BPACK_wrapper.h"

#define IVELO9_CONST 1

//------------------------------------------------------------------------------
using namespace std;

const double pi = 4.0*atan(1.0);
const _Complex double Im={0.0,1.0};
double timer=0.0;

extern "C" {
      ///////////////////////////////////////////////
      ////// BLACS //////////////////////////////////
      ///////////////////////////////////////////////
      void Cblacs_exit(int);

      extern void dgemv_(char *, int *, int *, double *, double *a, int *,
                  double *, int *, double *, double *, int *, int);
      extern void dgemm_(const char*, const char*, const int*, const int*, const int*,
                  const double*, const double*, const int*, const double*,
                  const int*, const double*, double*, const int*, int, int);

}


namespace my {
std::string to_string( double d ) {

    std::ostringstream stm ;
    stm << std::setprecision(std::numeric_limits<double>::digits10) << d ;
    return stm.str() ;
}
}


int product(int arr[], int n) {
    int out = 1;  // Initialize product to 1
    for (int i = 0; i < n; i++) {
        out *= arr[i];
    }
    return out;
}


void crossProduct3D(double v_A[], double v_B[], double c_P[]) {
   c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

double l2_norm(double const* u, int n) {
    double accum = 0.;
    for (int i = 0; i < n; ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}


int ceil_safe(double x){
  int x_round=round(x);
  if(fabs(x-x_round)<1e-13){
    return x_round;
  }else{
    return ceil(x);
  }
}
int floor_safe(double x){
  int x_round=round(x);
  if(fabs(x-x_round)<1e-13){
    return x_round;
  }else{
    return floor(x);
  }
}

int subdomain_detection(double x, double y, double center[], int shape, double radius_max, double L, double H){
  double point[2]={x-center[0],y-center[1]};
  if(shape==1){
    if(l2_norm(point, 2)<radius_max){
      return 0; // inner domain
    }else{
      return 1; // outer domain
    }
  }else if(shape==4){
    if((fabs(point[0])<L/2 && fabs(point[1])<H/2)){
      return 0; // inner domain
    }else{
      return 1; // outer domain
    }
  }else{
    std::cout<<"not implemented yet"<<std::endl;
  }
}

double slowness(double x,double y, double slow_x0, double slow_y0,int ivelo)
{
  double g1, g2, s0;
  g1=0; g2=-1.0/2.0;
  s0=2;  // This is at the domain reference point (slow_x0,slow_y0)

  if(ivelo==9){
#ifdef IVELO9_CONST
    s0 = 1.0;
#else
    s0 =1/((1.0/s0+g1*(x-slow_x0)+g2*(y-slow_y0)));
#endif
  }else if(ivelo==1){

  }else{
    cout<<"not supported ivelo: "<<ivelo<<endl;
    exit(0);
  }
  return s0;
}

//symmetric about x=0
double TukeyWindow(double x, double constwidth, double taperwidth){
    double tmp,out;
    if(x<constwidth){
      out = 1.0;
    }else if(x>constwidth+taperwidth){
      out = 0.0;
    }else{
      tmp = (x-constwidth)/taperwidth*pi;
      out = (1+cos(tmp))/2.0;
    }
    return out;
}

_Complex double source_function(double x,double y, double xs0, double ys0, int nth_rhs, double h, double omega)
{
  _Complex double out = {0.0,0.0};
  double window;
  // the first rhs is point source
  if(nth_rhs==0){
  if((fabs(x-xs0)<1e-10 && fabs(y-ys0)<1e-10)){
    out = {1.0,0.0};
  }
  }

  // the second rhs is gaussian pulse
  if(nth_rhs==1){
    double sigma=0.15;
    out = exp(-(pow(x-xs0,2)+pow(y-ys0,2))/pow(sigma,2)/2)*pow(h,2);
    double taperwidth=0.05;
    double constwidth=0.1;
    double dist = sqrt(pow(x-xs0,2)+pow(y-ys0,2));
    double window =TukeyWindow(dist,constwidth,taperwidth);
    out = out*window;
  }

  // the third rhs is all-1 inside a small circle
  if(nth_rhs==2){
    double dist = sqrt(pow(x-xs0,2)+pow(y-ys0,2));
    double taperwidth=0.1;
    double constwidth=0.25;
    out = TukeyWindow(dist,constwidth,taperwidth);
    out *= pow(h,2);
  }

  // the fourth rhs is gaussian wave packet
  if(nth_rhs==3){
    double omega0=omega*0.9;
    double dist = sqrt(pow(x-xs0,2)+pow(y-ys0,2));
    double d[2];
    d[0]=1/sqrt(2.0);
    d[1]=1/sqrt(2.0);

    double sigma=0.15;
    double phase = omega0*(x*d[0]+y*d[1]);
    out = exp(-(pow(x-xs0,2)+pow(y-ys0,2))/pow(sigma,2)/2)*(cos(phase)+Im*sin(phase))*pow(h,2);
    double taperwidth=0.05; //0.1; //0.01;
    double constwidth=0.1; //0.3; //0.05;
    double window =TukeyWindow(dist,constwidth,taperwidth);
    out = out*window;


    // double phase = omega0*(x*d[0]+y*d[1]);
    // out = exp(-4*omega0*pow(dist,2))*(cos(phase)+Im*sin(phase))*pow(h,2);
    // double tmp = 1/sqrt(omega0*8)*4;
    // double taperwidth=tmp/2;
    // double constwidth=tmp;
    // double window =TukeyWindow(dist,constwidth,taperwidth);
    // out = out*window;
  }

  return out;
}







_Complex double myhankel(double v, double z){
  double eps1 = 1e-13;
  _Complex double out =  {0.0,0.0};
  if(v<0){
      out = myhankel(-v,z);
      out = out * (cos(-v*pi)+Im*sin(-v*pi));
  }else if(fabs(round(v)-v)<eps1){   //v is integer
      out = bessj(round(v),z) +Im*bessy(round(v),z);
  }else if(fabs(v-0.5)<eps1){ // v=0.5
      out = -Im*sqrt(2.0/(pi*z))*(cos(z)+Im*sin(z));
  }else{
      printf("wrong order in myhankel \n");
      exit(0);
  }
  return out;
}


_Complex double Babich(int m, double w, double v0, double v1, double tau){
  _Complex double out;
  double q = -(m-2.0)/2.0;
  _Complex double t0 = Im*sqrt(pi)/2.0*(cos(q*pi)+Im*sin(q*pi))*pow(2.0*tau/w,q)*(myhankel(q,w*tau));
  q = 1.0-(m-2.0)/2.0;
  _Complex double t1 = Im*sqrt(pi)/2.0*(cos(q*pi)+Im*sin(q*pi))*pow(2.0*tau/w,q)*(myhankel(q,w*tau));
  out =  v0*t0+v1*t1;
  // std::cout << std::fixed;
  // std::cout << std::setprecision(14);
  // cout<<"vt "<<v0<<" "<<v1<<" "<<__real__ t0<<" "<<__imag__ t0<<" "<<t1<<" "<<tau<<" "<<w<<" "<<__real__ myhankel(0.0,w*tau)<<" "<<__imag__ myhankel(0.0,w*tau)<<" "<<w*tau<<" "<<bessj0(w*tau)<<" "<<bessy0(w*tau)<<endl;
  return out;
}


// The object handling kernel parameters and sampling function
class C_QuantApp_BF {
public:
  vector<double> _data;
  vector<double> _data_m;
  vector<double> _panelnodes;
  vector<double> _panelnorms;
  vector<int> _v_sub2glo;
  int _d = 0;   // data dimension 2 or 3
  int _n = 0;   // max dimension: max(nx,ny)
  int _nx = 0;   // nx
  int _ny = 0;   // ny
  double _w = pi;   //angular frequency

  double _xmin=-0.2/1.0,  _xmax=1.2/1.0;
  double _ymin=-0.2/1.0,  _ymax=1.2/1.0;
  double _x0min=-0.0, _x0max=1.0/1.0;
  double _y0min=-0.0, _y0max=1.0/1.0;
  double _slow_x0 = 0.5;
  double _slow_y0 = 0.5;
  double _h= 0.01;
  double _dl= 0.01;
  double _rmax= 10;
  int _ivelo = 9; // 1 homogenous, 2 constant gradient
  int _TNx, _TNy;
  int _nquad = 4;

  int _verbose=0;
  int _vs=1;


  std::vector<double> _x_cheb;
  std::vector<double> _y_cheb;
  std::vector<double> _u1_square_int_cheb;
  std::vector<double> _D1_int_cheb;
  std::vector<double> _D2_int_cheb;

  std::vector<int> _Hperm;
  std::vector<int> _iHperm;

  std::vector<int> _Hperm_m;
  std::vector<int> _iHperm_m;

  C_QuantApp_BF() = default;

  // constructor for the v2v operator
  C_QuantApp_BF(vector<double> data, int d, int nx, int ny, double w, double xmin, double xmax, double ymin, double ymax, double x0min, double x0max, double y0min, double y0max, double h, double dl, int nquad, int ivelo, int rmax, int verbose, int vs, vector<double> x_cheb, vector<double> y_cheb, vector<double> u1_square_int_cheb, vector<double> D1_int_cheb, vector<double> D2_int_cheb)
    :_data(move(data)),_d(d), _nx(nx), _ny(ny), _w(w),_xmin(xmin),_xmax(xmax),_ymin(ymin),_ymax(ymax),_x0min(x0min),_x0max(x0max),_y0min(y0min),_y0max(y0max),_h(h),_dl(dl),_nquad(nquad),_ivelo(ivelo),_rmax(rmax),_verbose(verbose),_vs(vs),_x_cheb(move(x_cheb)),_y_cheb(move(y_cheb)),_u1_square_int_cheb(move(u1_square_int_cheb)),_D1_int_cheb(move(D1_int_cheb)),_D2_int_cheb(move(D2_int_cheb)){
      _TNx = _x_cheb.size();
      _TNy = _y_cheb.size();
      _slow_x0 = (_x0min+_x0max)/2.0;
      _slow_y0 = (_y0min+_y0max)/2.0;
      _n = max(_nx,_ny);
	}


  // inline void Sample_noself(double x1, double y1, double x2, double y2, _Complex double* val){

  //   double fr[3];
  //   lagrange_interp2D_vec(_x_cheb.data(), _y_cheb.data(), x1,y1,x2,y2,_u1_square_int_cheb.data(),_D1_int_cheb.data(),_D2_int_cheb.data(), _TNx,_TNy,fr);
  //   double tau_square=fr[0];
  //   double D1 =fr[1];
  //   double D2 =fr[2];

  //   // double tau_square1 = (pow(x1-x2,2)+pow(y1-y2,2))*4.0;
  //   // if(fabs(tau_square-tau_square1)/tau_square1>1e-10){
  //   //   cout<<tau_square<<" "<<fabs(tau_square-tau_square1)/tau_square1<<" "<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<endl;
  //   //   exit(0);
  //   // }


  //   double s0 = slowness(x1,y1, _slow_x0, _slow_y0,_ivelo);
  //   double s1 = slowness(x2,y2, _slow_x0, _slow_y0,_ivelo);
  //   double tau = sqrt(pow(s0,2)* (pow(x1-x2,2) + pow(y1-y2,2)));
  //   // if(tau_square>0){
  //   // if(tau>2*pi/(_w)/50){
  //   if(tau>2*pi/(_w)/20){
  //     tau = sqrt(max(tau_square,1e-30));
  //   }

  //   // tau= acosh(((pow(x1-x2,2) + pow(y1-y2,2)))*pow(0.25,2)*s0*s1*0.5+1)/0.25;

  //   // if(tau_square<0){
  //   //   cout<<"be careful, tau_square: "<<tau_square<<endl;
  //   // }
  //   // double D1 = lagrange_interp2D(_x_cheb.data(), _y_cheb.data(), x1,y1,x2,y2,_D1_int_cheb.data(), _TNx,_TNy);
  //   // double D2 = lagrange_interp2D(_x_cheb.data(), _y_cheb.data(), x1,y1,x2,y2,_D2_int_cheb.data(), _TNx,_TNy);

  //   //  cout<<m<<" "<<n<<" "<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<D1<<" "<<D2<<" "<<tau<<" gagag"<<endl;

  //   _Complex double out =Babich(_d, _w, D1, D2, tau);

  //   // if(out != out)  // this seems to have no effect and causes some runtime error sometimes
  //   //   cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<D1<<" "<<D2<<" "<<tau<<" "<<tau_square<<" "<<__real__ out<<" "<<__imag__ out<<endl;
  //   *val = out;
  //   // exit(0);

  // }



  inline void SampleSelf(double x1, double y1, double x2, double y2, _Complex double* val){

    int self = sqrt(pow(x1-x2,2)+pow(y1-y2,2))<1e-20? 1:0;
    if(_vs==1){
      // int closeby = fabs(x1-x2)<_h*2+1e-20 && fabs(y1-y2)<_h*2+1e-20?1:0;
      if(self==1){
        double s0 = slowness(x1,y1, _slow_x0, _slow_y0,_ivelo);
        double gamma = 1.781072418;

        // 7-point Legendre-Gauss Quadrature on [0,pi/4] for the integral due to Jianliang Qian
        double nodes[7] = {0.765412887308718, 0.683897697334573, 0.552074099956508, 0.392699081698724, 0.233324063440941, 0.101500466062876, 0.019985276088730};
        double weights[7] = {0.050848627308305, 0.109840050384021, 0.149944310198338, 0.164132187616120, 0.149944310198338, 0.109840050384021, 0.050848627308305};
        *val = 0;
        double tt = _h*_w*s0;
        for (int i=0;i<7;i++){
          double x = nodes[i];
          *val += weights[i]*(tt/(2*cos(x))*(bessj(1,tt/(2*cos(x))) +Im*bessy(1,tt/(2*cos(x)))) + Im*2/pi);  // here -Im instead of Im is needed due to the use of IEEE imaginary unit definition
        }
        *val*=8/pow(tt,2)/4.0*Im;

        // std::cout<<__real__(*val)<<" "<<__imag__(*val) << std::endl;

      }else{
        std::cout<<"should not arrive here for _vs==1 "<< std::endl;
        // Sample_noself(x1, y1, x2, y2,val);
      }
    }else{
      // int closeby = sqrt(pow(x1-x2,2)+pow(y1-y2,2))<1e-10+_dl? 1:0;
      if(self==1){
        double s0 = slowness(x1,y1, _slow_x0, _slow_y0,_ivelo);
        double gamma = 1.781072418;
        // *val = 0;
        *val = Im*_dl/4.0*(1+Im*2/pi*(log(gamma*_w*s0*_dl/4.0)-1));
      }else{
        // Sample_noself(x1, y1, x2, y2,val);
        // *val = *val*_dl;
        std::cout<<"should not arrive here for _vs==0 "<< std::endl;
      }
    }
  }
};



// Assemble a block of matrix entries from interpolated D1, D2, tau
void assemble_fromD1D2Tau(double x1,double x2,double y1,double y2, _Complex double* output, C_QuantApp_BF* Q){

    int self = sqrt(pow(x1-x2,2)+pow(y1-y2,2))<1e-20? 1:0;
    if(self==1){
      Q->SampleSelf(x1, y1, x2, y2, output);
    }else{
      // double tau_square=fr[idxr+idxc*nr];
      double D1 =1.0/2.0/sqrt(pi); //fr[nr*nc + idxr+idxc*nr];
      double D2 =0;// fr[nr*nc*2 + idxr+idxc*nr];
      double s0 = slowness(x1,y1, Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
      // double s1 = slowness(x2[idxc],y2[idxc], Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
      // cout<<s0<<" "<<s1<<" "<<Q->_w<<endl;
      double tau = sqrt(pow(s0,2)* (pow(x1-x2,2) + pow(y1-y2,2)));
      // if(tau_square>0){
      // if(tau>2*pi/(Q->_w)/50){
      // if(tau>2*pi/(Q->_w)/20){
      //   tau = sqrt(max(tau_square,1e-30));
      // }
      // tau= acosh(((pow(x1[idxr]-x2[idxc],2) + pow(y1[idxr]-y2[idxc],2)))*pow(0.25,2)*s0*s1*0.5+1)/0.25;
      *output =Babich(Q->_d, Q->_w, D1, D2, tau);
    }
}


// Assemble a block of matrix entries from interpolated D1, D2, tau
void assemble_fromD1D2Tau_s2s(double x1,double x2,double y1,double y2, _Complex double* output, C_QuantApp_BF* Q){

    double s1 = slowness(x2,y2, Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
    double s0=2;
    double k0 = s0*Q->_w;
    double coef = pow(k0,2.0)*(pow(s1/s0,2.0)-1);
    int self = sqrt(pow(x1-x2,2)+pow(y1-y2,2))<1e-20? 1:0;
    if(self==1){
      Q->SampleSelf(x1, y1, x2, y2, output);
      *output = -*output*coef + 1.0/pow(Q->_h,2.0);
    }else{
      // double tau_square=fr[idxr+idxc*nr];
      double D1 =1.0/2.0/sqrt(pi); //fr[nr*nc + idxr+idxc*nr];
      double D2 =0;// fr[nr*nc*2 + idxr+idxc*nr];
      // double s0 = slowness(x1,y1, Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
      // double s1 = slowness(x2[idxc],y2[idxc], Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
      // cout<<s0<<" "<<s1<<" "<<Q->_w<<endl;
      double tau = sqrt(pow(s0,2)* (pow(x1-x2,2) + pow(y1-y2,2)));
      // if(tau_square>0){
      // if(tau>2*pi/(Q->_w)/50){
      // if(tau>2*pi/(Q->_w)/20){
      //   tau = sqrt(max(tau_square,1e-30));
      // }
      // tau= acosh(((pow(x1[idxr]-x2[idxc],2) + pow(y1[idxr]-y2[idxc],2)))*pow(0.25,2)*s0*s1*0.5+1)/0.25;
      *output =-coef*Babich(Q->_d, Q->_w, D1, D2, tau);
    }
}



// Assemble a block of matrix entries from interpolated D1, D2, tau
void assemble_fromD1D2Tau_block(int nth, int nr, int nc, double* x1,double* x2,double* y1,double* y2, _Complex double* alldat_loc, int64_t* idx_val_map, double *fr, C_QuantApp_BF* Q){
    #pragma omp parallel for
    for (int idxrc=0;idxrc<nr*nc;idxrc++){
      int idxr = idxrc%nr;
      int idxc = idxrc/nr;
      int self = sqrt(pow(x1[idxr]-x2[idxc],2)+pow(y1[idxr]-y2[idxc],2))<1e-20? 1:0;
      // int closeby;
      // if(Q->_vs==1){
      //   closeby = fabs(x1[idxr]-x2[idxc])<Q->_h*2+1e-20 && fabs(y1[idxr]-y2[idxc])<Q->_h*2+1e-20?1:0;
      // }{
      //   closeby = sqrt(pow(x1[idxr]-x2[idxc],2)+pow(y1[idxr]-y2[idxc],2))<1e-10+Q->_dl? 1:0;
      // }
      if(self==1){
        _Complex double valtmp;
        Q->SampleSelf(x1[idxr], y1[idxr], x2[idxc], y2[idxc], &valtmp);
        alldat_loc[idx_val_map[nth]+idxr+idxc*nr]=valtmp;
      }else{
        // double tau_square=fr[idxr+idxc*nr];
        double D1 =1.0/2.0/sqrt(pi); //fr[nr*nc + idxr+idxc*nr];
        double D2 =0;// fr[nr*nc*2 + idxr+idxc*nr];



        double s0 = slowness(x1[idxr],y1[idxr], Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
        // double s1 = slowness(x2[idxc],y2[idxc], Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
        // cout<<s0<<" "<<s1<<" "<<Q->_w<<endl;
        double tau = sqrt(pow(s0,2)* (pow(x1[idxr]-x2[idxc],2) + pow(y1[idxr]-y2[idxc],2)));
        // if(tau_square>0){
        // if(tau>2*pi/(Q->_w)/50){
        // if(tau>2*pi/(Q->_w)/20){
        //   tau = sqrt(max(tau_square,1e-30));
        // }
        // tau= acosh(((pow(x1[idxr]-x2[idxc],2) + pow(y1[idxr]-y2[idxc],2)))*pow(0.25,2)*s0*s1*0.5+1)/0.25;
        if(Q->_vs==1){
          alldat_loc[idx_val_map[nth]+idxr+idxc*nr] =Babich(Q->_d, Q->_w, D1, D2, tau);
        }else{
          alldat_loc[idx_val_map[nth]+idxr+idxc*nr] =Q->_dl*Babich(Q->_d, Q->_w, D1, D2, tau);
        }
      }
    }
}





// Assemble a block of matrix entries from interpolated D1, D2, tau
void assemble_fromD1D2Tau_block_s2s(int nth, int nr, int nc, double* x1,double* x2,double* y1,double* y2, _Complex double* alldat_loc, int64_t* idx_val_map, double *fr, C_QuantApp_BF* Q){
    #pragma omp parallel for
    for (int idxrc=0;idxrc<nr*nc;idxrc++){
      int idxr = idxrc%nr;
      int idxc = idxrc/nr;
      int self = sqrt(pow(x1[idxr]-x2[idxc],2)+pow(y1[idxr]-y2[idxc],2))<1e-20? 1:0;
      // int closeby;
      // if(Q->_vs==1){
      //   closeby = fabs(x1[idxr]-x2[idxc])<Q->_h*2+1e-20 && fabs(y1[idxr]-y2[idxc])<Q->_h*2+1e-20?1:0;
      // }{
      //   closeby = sqrt(pow(x1[idxr]-x2[idxc],2)+pow(y1[idxr]-y2[idxc],2))<1e-10+Q->_dl? 1:0;
      // }

      double s1 = slowness(x2[idxc],y2[idxc], Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
      double s0=2;
      double k0 = s0*Q->_w;
      double coef = pow(k0,2.0)*(pow(s1/s0,2.0)-1);


      if(self==1){
        _Complex double valtmp;
        Q->SampleSelf(x1[idxr], y1[idxr], x2[idxc], y2[idxc], &valtmp);
        alldat_loc[idx_val_map[nth]+idxr+idxc*nr]=-valtmp*coef + 1.0/pow(Q->_h,2.0);
      }else{
        // double tau_square=fr[idxr+idxc*nr];
        double D1 =1.0/2.0/sqrt(pi); //fr[nr*nc + idxr+idxc*nr];
        double D2 =0;// fr[nr*nc*2 + idxr+idxc*nr];

        // double s0 = slowness(x1[idxr],y1[idxr], Q->_slow_x0, Q->_slow_y0,Q->_ivelo);

        double tau = sqrt(pow(s0,2)* (pow(x1[idxr]-x2[idxc],2) + pow(y1[idxr]-y2[idxc],2)));
        // if(tau_square>0){
        // if(tau>2*pi/(Q->_w)/50){
        // if(tau>2*pi/(Q->_w)/20){
        //   tau = sqrt(max(tau_square,1e-30));
        // }
        // tau= acosh(((pow(x1[idxr]-x2[idxc],2) + pow(y1[idxr]-y2[idxc],2)))*pow(0.25,2)*s0*s1*0.5+1)/0.25;
        alldat_loc[idx_val_map[nth]+idxr+idxc*nr] =-coef*Babich(Q->_d, Q->_w, D1, D2, tau);
      }
    }
}


// The distance function wrapper required by the Fortran HODLR code
inline void C_FuncDistmn_dummy(int *m, int *n, double *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}

// The compressibility function wrapper required by the Fortran HODLR code
inline void C_FuncNearFar_dummy(int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBZmn(int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  std::cout<<"C_FuncBZmn is no more used, use C_FuncBZmnBlock instead"<<std::endl;

}



// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn_BF(int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  // double x1 = Q->_data[(*m-1) * Q->_d];
  // double y1 = Q->_data[(*m-1) * Q->_d+1];
  // double x2 = Q->_data[(*n-1) * Q->_d];
  // double y2 = Q->_data[(*n-1) * Q->_d+1];
  // Q->Sample(x1, y1, x2, y2,val);

  std::cout<<"C_FuncZmn_BF is no more used, use C_FuncZmnBlock_BF instead"<<std::endl;

}


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn_BF_V2V_MD(int* Ndim, int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  double x1 = Q->_data[(m[0]-1) * Q->_d];
  double y1 = Q->_data[(m[1]-1) * Q->_d+1];
  double x2 = Q->_data[(n[0]-1) * Q->_d];
  double y2 = Q->_data[(n[1]-1) * Q->_d+1];
  assemble_fromD1D2Tau(x1,x2,y1,y2,val, Q);
}



// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn_BF_S2S_MD(int* Ndim, int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  double x1 = Q->_data[(m[0]-1) * Q->_d];
  double y1 = Q->_data[(m[1]-1) * Q->_d+1];
  double x2 = Q->_data[(n[0]-1) * Q->_d];
  double y2 = Q->_data[(n[1]-1) * Q->_d+1];
  assemble_fromD1D2Tau_s2s(x1,x2,y1,y2,val, Q);
}



// The distance function wrapper required by the Fortran HODLR code
inline void C_FuncDistmn_BF(int *m, int *n, double *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}

// The compressibility function wrapper required by the Fortran HODLR code
inline void C_FuncNearFar_BF(int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}


// The compressibility function (tensor) wrapper required by the Fortran HODLR code
inline void C_FuncNearFar_BF_MD(int* Ndim, int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}



// The entry evaluation routine for a list of subtensors for the background operator
inline void C_FuncZmnBlock_BF_V2V_MD(int* Ndim, int* Ninter, int* Nrow_max, int* Ncol_max, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

    int myrank, size;                     // Store values of processor rank and total no of procs requested
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    int64_t idx_row_1=0;
    int64_t idx_row_2=0;
    int64_t idx_col_1=0;
    int64_t idx_col_2=0;
    int64_t idx_val=0;
    int nrmax=0;
    int ncmax=0;
    int nvalmax=0;
    vector<int64_t> idx_row_map_1(*Ninter,0);
    vector<int64_t> idx_row_map_2(*Ninter,0);
    vector<int64_t> idx_col_map_1(*Ninter,0);
    vector<int64_t> idx_col_map_2(*Ninter,0);
    vector<int64_t> idx_val_map(*Ninter,0);

    time_t start, end;
    time(&start);
    for (int nn=0;nn<*Ninter;nn++){
      int pp = pgidx[nn];
      int nprow = pmaps[pp];
      int npcol = pmaps[*Npmap+pp];
      int pid = pmaps[(*Npmap)*2+pp];
      int nr_1 = rowidx[*Ndim*nn];
      int nr_2 = rowidx[*Ndim*nn+1];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];

      if(nprow*npcol==1){
        idx_row_map_1[nn]=idx_row_1;
        idx_row_map_2[nn]=idx_row_2;
        idx_col_map_1[nn]=idx_col_1;
        idx_col_map_2[nn]=idx_col_2;
        idx_val_map[nn]=idx_val;
        idx_val+=nr_1*nr_2*nc_1*nc_2;
        idx_row_1+=nr_1;
        idx_row_2+=nr_2;
        idx_col_1+=nc_1;
        idx_col_2+=nc_2;
        nrmax = max(nr_1*nr_2,nrmax);
        ncmax = max(nc_1*nc_2,ncmax);
        nvalmax = max(nr_1*nr_2*nc_1*nc_2,nvalmax);
      }else{
        std::cout<<"nprow*npcol>1 in C_FuncZmnBlock_BF_V2V_MD"<<std::endl;
        exit(0);
      }
    }

    double *x1,*y1;
    double *x2,*y2;
    double *fr;
    // cout<<" "<<myrank<<" "<<NinterNew<<endl;
    // #pragma omp parallel private(x1,y1,x2,y2,fr)
    {
    x1= (double*)malloc(nrmax*sizeof(double));
    y1= (double*)malloc(nrmax*sizeof(double));
    x2= (double*)malloc(ncmax*sizeof(double));
    y2= (double*)malloc(ncmax*sizeof(double));
    fr= (double*)malloc(nvalmax*3*sizeof(double));
    // #pragma omp for
    for (int nn=0;nn<*Ninter;nn++){
      int nr_1 = rowidx[*Ndim*nn];
      int nr_2 = rowidx[*Ndim*nn+1];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];

      for (int idxr_1=0;idxr_1<nr_1;idxr_1++){
      for (int idxr_2=0;idxr_2<nr_2;idxr_2++){
        int idxr = idxr_1 + nr_1*idxr_2;
        int m_1=allrows[idx_row_map_1[nn]+idxr_1];
        int m_2=allrows[*Nrow_max + idx_row_map_2[nn]+idxr_2];
        x1[idxr] = Q->_data[(m_1-1) * Q->_d];
        y1[idxr] = Q->_data[(m_2-1) * Q->_d+1];
      }
      }

      for (int idxc_1=0;idxc_1<nc_1;idxc_1++){
      for (int idxc_2=0;idxc_2<nc_2;idxc_2++){
        int idxc = idxc_1 + nc_1*idxc_2;
        int n_1=allcols[idx_col_map_1[nn]+idxc_1];
        int n_2=allcols[*Ncol_max + idx_col_map_2[nn]+idxc_2];
        x2[idxc] = Q->_data[(n_1-1) * Q->_d];
        y2[idxc] = Q->_data[(n_2-1) * Q->_d+1];
      }
      }

      // lagrange_interp2D_vec_block(Q->_x_cheb.data(), Q->_y_cheb.data(), nr, nc, x1,y1,x2,y2,Q->_u1_square_int_cheb.data(),Q->_D1_int_cheb.data(),Q->_D2_int_cheb.data(), Q->_TNx,Q->_TNy, fr);

      assemble_fromD1D2Tau_block(nn,nr_1*nr_2, nc_1*nc_2, x1,x2,y1,y2, alldat_loc, idx_val_map.data(),fr, Q);

    }

    free(x1);
    free(y1);
    free(x2);
    free(y2);
    free(fr);
    }
time(&end);
timer += difftime(end,start);

}



// The entry evaluation routine for a list of subtensors for the VIE operator
inline void C_FuncZmnBlock_BF_S2S_MD(int* Ndim, int* Ninter, int* Nrow_max, int* Ncol_max, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

    int myrank, size;                     // Store values of processor rank and total no of procs requested
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    int64_t idx_row_1=0;
    int64_t idx_row_2=0;
    int64_t idx_col_1=0;
    int64_t idx_col_2=0;
    int64_t idx_val=0;
    int nrmax=0;
    int ncmax=0;
    int nvalmax=0;
    vector<int64_t> idx_row_map_1(*Ninter,0);
    vector<int64_t> idx_row_map_2(*Ninter,0);
    vector<int64_t> idx_col_map_1(*Ninter,0);
    vector<int64_t> idx_col_map_2(*Ninter,0);
    vector<int64_t> idx_val_map(*Ninter,0);

    time_t start, end;
    time(&start);
    for (int nn=0;nn<*Ninter;nn++){
      int pp = pgidx[nn];
      int nprow = pmaps[pp];
      int npcol = pmaps[*Npmap+pp];
      int pid = pmaps[(*Npmap)*2+pp];
      int nr_1 = rowidx[*Ndim*nn];
      int nr_2 = rowidx[*Ndim*nn+1];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];

      if(nprow*npcol==1){
        idx_row_map_1[nn]=idx_row_1;
        idx_row_map_2[nn]=idx_row_2;
        idx_col_map_1[nn]=idx_col_1;
        idx_col_map_2[nn]=idx_col_2;
        idx_val_map[nn]=idx_val;
        idx_val+=nr_1*nr_2*nc_1*nc_2;
        idx_row_1+=nr_1;
        idx_row_2+=nr_2;
        idx_col_1+=nc_1;
        idx_col_2+=nc_2;
        nrmax = max(nr_1*nr_2,nrmax);
        ncmax = max(nc_1*nc_2,ncmax);
        nvalmax = max(nr_1*nr_2*nc_1*nc_2,nvalmax);
      }else{
        std::cout<<"nprow*npcol>1 in C_FuncZmnBlock_BF_S2S_MD"<<std::endl;
        exit(0);
      }
    }

    double *x1,*y1;
    double *x2,*y2;
    double *fr;
    // cout<<" "<<myrank<<" "<<NinterNew<<endl;
    // #pragma omp parallel private(x1,y1,x2,y2,fr)
    {
    x1= (double*)malloc(nrmax*sizeof(double));
    y1= (double*)malloc(nrmax*sizeof(double));
    x2= (double*)malloc(ncmax*sizeof(double));
    y2= (double*)malloc(ncmax*sizeof(double));
    fr= (double*)malloc(nvalmax*3*sizeof(double));
    // #pragma omp for
    for (int nn=0;nn<*Ninter;nn++){
      int nr_1 = rowidx[*Ndim*nn];
      int nr_2 = rowidx[*Ndim*nn+1];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];

      for (int idxr_1=0;idxr_1<nr_1;idxr_1++){
      for (int idxr_2=0;idxr_2<nr_2;idxr_2++){
        int idxr = idxr_1 + nr_1*idxr_2;
        int m_1=allrows[idx_row_map_1[nn]+idxr_1];
        int m_2=allrows[*Nrow_max + idx_row_map_2[nn]+idxr_2];
        x1[idxr] = Q->_data[(m_1-1) * Q->_d];
        y1[idxr] = Q->_data[(m_2-1) * Q->_d+1];
      }
      }

      for (int idxc_1=0;idxc_1<nc_1;idxc_1++){
      for (int idxc_2=0;idxc_2<nc_2;idxc_2++){
        int idxc = idxc_1 + nc_1*idxc_2;
        int n_1=allcols[idx_col_map_1[nn]+idxc_1];
        int n_2=allcols[*Ncol_max + idx_col_map_2[nn]+idxc_2];
        x2[idxc] = Q->_data[(n_1-1) * Q->_d];
        y2[idxc] = Q->_data[(n_2-1) * Q->_d+1];
      }
      }

      // lagrange_interp2D_vec_block(Q->_x_cheb.data(), Q->_y_cheb.data(), nr, nc, x1,y1,x2,y2,Q->_u1_square_int_cheb.data(),Q->_D1_int_cheb.data(),Q->_D2_int_cheb.data(), Q->_TNx,Q->_TNy, fr);

      assemble_fromD1D2Tau_block_s2s(nn,nr_1*nr_2, nc_1*nc_2, x1,x2,y1,y2, alldat_loc, idx_val_map.data(),fr, Q);

    }

    free(x1);
    free(y1);
    free(x2);
    free(y2);
    free(fr);
    }
time(&end);
timer += difftime(end,start);

}



// The command line parser for the example related parameters
void set_option_from_command_line(int argc, const char* const* cargv,F2Cptr option0) {
    double opt_d;
    int opt_i;
    std::vector<std::unique_ptr<char[]>> argv_data(argc);
    std::vector<char*> argv(argc);
    for (int i=0; i<argc; i++) {
      argv_data[i].reset(new char[strlen(cargv[i])+1]);
      argv[i] = argv_data[i].get();
      strcpy(argv[i], cargv[i]);
    }
    option long_options[] =
      {{"nmin_leaf",                     required_argument, 0, 1},
       {"tol_comp",                   required_argument, 0, 2},
       {"tol_rand",                   required_argument, 0, 3},
       {"tol_Rdetect",             required_argument, 0, 4},
       {"tol_itersol",             required_argument, 0, 5},
       {"n_iter",          required_argument, 0, 6},
       {"level_check",         required_argument, 0, 7},
       {"precon",                  required_argument, 0, 8},
       {"xyzsort",      required_argument, 0, 9},
       {"lrlevel",     required_argument, 0, 10},
       {"errfillfull",       required_argument, 0, 11},
       {"baca_batch",      required_argument, 0, 12},
       {"reclr_leaf",      required_argument, 0, 13},
       {"nogeo",     required_argument, 0, 14},
       {"less_adapt",            required_argument, 0, 15},
       {"errsol",           required_argument, 0, 16},
       {"lr_blk_num",                  required_argument, 0, 17},
       {"rank0",  required_argument, 0, 18},
       {"rankrate", required_argument, 0, 19},
       {"itermax",               required_argument, 0, 20},
       {"powiter",  required_argument, 0, 21},
       {"ilu", required_argument, 0, 22},
       {"nbundle",     required_argument, 0, 23},
       {"near_para",  required_argument, 0, 24},
       {"format",  required_argument, 0, 25},
       {"verbosity", required_argument, 0, 26},
       {"rmax", required_argument, 0, 27},
       {"sample_para", required_argument, 0, 28},
       {"pat_comp",    required_argument, 0, 29},
       {"knn",         required_argument, 0, 30},
       {"knn_near_para",         required_argument, 0, 31},
       {"forwardN15flag",         required_argument, 0, 32},
       {"sample_para_outer",         required_argument, 0, 33},
       {"elem_extract",         required_argument, 0, 34},
       {"fastsample_tensor",         required_argument, 0, 35},
       {"use_zfp",         required_argument, 0, 36},
       {"use_qtt",         required_argument, 0, 37},
       {NULL, 0, NULL, 0}
      };
    int c, option_index = 0;
    // bool unrecognized_options = false;
    opterr = optind = 0;
    while ((c = getopt_long_only
            (argc, argv.data(), "",
             long_options, &option_index)) != -1) {
      switch (c) {
      case 1: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "Nmin_leaf", opt_i);
      } break;
      case 2: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_comp", opt_d);
        z_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
        z_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d*0.1);
      } break;
      case 3: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
      } break;
      case 4: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d);
      } break;
      case 5: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_itersol", opt_d);
      } break;
      case 6: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "n_iter", opt_i);
      } break;
      case 7: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "level_check", opt_i);
      } break;
      case 8: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "precon", opt_i);
      } break;
      case 9: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "xyzsort", opt_i);
      } break;
      case 10: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "LRlevel", opt_i);
      } break;
      case 11: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "ErrFillFull", opt_i);
      } break;
      case 12: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "BACA_Batch", opt_i);
      } break;
      case 13: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "RecLR_leaf", opt_i);
      } break;
      case 14: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "nogeo", opt_i);
      } break;
      case 15: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "less_adapt", opt_i);
      } break;
      case 16: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "ErrSol", opt_i);
      } break;
      case 17: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "LR_BLK_NUM", opt_i);
      } break;
      case 18: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "rank0", opt_i);
      } break;
      case 19: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "rankrate", opt_d);
      } break;
      case 20: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "itermax", opt_i);
      } break;
      case 21: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "powiter", opt_i);
      } break;
      case 22: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "ILU", opt_i);
      } break;
      case 23: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "Nbundle", opt_i);
      } break;
      case 24: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "near_para", opt_d);
      } break;
      case 25: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "format", opt_i);
      } break;
      case 26: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "verbosity", opt_i);
      } break;
      case 27: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "rmax", opt_i);
      } break;
      case 28: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "sample_para", opt_d);
      } break;
      case 29: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "pat_comp", opt_i);
      } break;
      case 30: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "knn", opt_i);
      } break;
      case 31: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "knn_near_para", opt_d);
      } break;
      case 32: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "forwardN15flag", opt_i);
      } break;
      case 33: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "sample_para_outer", opt_d);
      } break;
      case 34: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "elem_extract", opt_i);
      } break;
      case 35: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "fastsample_tensor", opt_i);
      } break;
      case 36: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "use_zfp", opt_i);
      } break;
      case 37: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "use_qtt", opt_i);
      } break;
      default: break;
      }
    }
  }



// This example uses Chebshev interpolation from WENO to compute the phase and amplitudes, and then use entry evaluation to compute the hodbf
////////////////////////////////////////////////////////////////////////////////
// --------------------------- Main Code Starts Here ------------------------ //

int main(int argc, char* argv[])
{

    int myrank, size;                     // Store values of processor rank and total no of procs requestedss
    int master_rank = 0;
	  MPI_Init(&argc, &argv); 	                            // Initialize MPI, called only once
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    MPI_Op op;


	int Ndim=1; //data dimension
	double starttime, endtime;
	double* dat_ptr_m, *dat_ptr_n, *dat_ptr_k, *dat_ptr_l;
	int* nns_ptr_m, *nns_ptr_n, *nns_ptr_k, *nns_ptr_l;
	int nogeo=0;  // 1: no geometrical information passed to hodlr, dat_ptr and Ndim are dummy

	int Nmin=4; //finest leafsize
	double tol=1e-4; //compression tolerance
	double sample_para=2.0; //oversampling factor in entry evaluation
	double sample_para_outer=2.0; //oversampling factor in entry evaluation
	int com_opt=5; //1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton 3:ACA_parallel
	int sort_opt=1; //0:natural order 1:kd-tree 2:cobble-like ordering 3:gram distance-based cobble-like ordering
	int checkerr = 0; //1: check compression quality
	int batch = 16; //batch size for BACA
	int bnum = 1; //sqrt of #of subblocks in H-BACA
	int knn=0; //k nearest neighbours stored per point
	int v_major,v_minor,v_bugfix; //version numbers
  int ker=3;
  int vs=1;   //volume or surface tests
  int shape=-1;
  int readtable=0;

  // int tst = 2;
	int lrlevel=100;
  int verbose=0;
  int elem_extract=0;

if(myrank==master_rank){
	z_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}

	/*****************************************************************/
  int Nx, Ny; // number of grid points per x or y direction
  int Nmax; // max between Nx and Ny

  // define the inner domain
  double x0min=-0.0, x0max=2.0/1.0;
  double y0min=-0.0, y0max=2.0/1.0;

  // define the outer domain
  double xouter=0.2, youter=0.2;
  double xmin=x0min-xouter,  xmax=x0max+xouter;
  double ymin=y0min-youter,  ymax=y0max+youter;

  double h=0.02;
  double h0=0.01;
  int ivelo = 1;
  int ivelo_o = 1; // this is assumed fixed
  int rmax = 20;
  double w = 5*pi;
  int TNx =10, TNy=10;

  double slow_x0;
  double slow_y0;



  FILE *fout1;

  //getting the example configurations from command line
  std::vector<std::unique_ptr<char[]>> argv_data(argc);
  std::vector<char*> argv1(argc);
  for (int i=0; i<argc; i++) {
    argv_data[i].reset(new char[strlen(argv[i])+1]);
    argv1[i] = argv_data[i].get();
    strcpy(argv1[i], argv[i]);
  }


  option long_options[] =
    {{"tst",                     required_argument, 0, 1},
      {"vs",                   required_argument, 0, 2},
      // {"N",                   required_argument, 0, 3},
      // {"K",             required_argument, 0, 4},
      // {"L",             required_argument, 0, 5},
      {"TNx",          required_argument, 0, 6},
      {"xmin",          required_argument, 0, 7},
      {"xmax",          required_argument, 0, 8},
      {"ymin",          required_argument, 0, 9},
      {"ymax",          required_argument, 0, 10},
      {"x0min",          required_argument, 0, 11},
      {"x0max",          required_argument, 0, 12},
      {"y0min",          required_argument, 0, 13},
      {"y0max",          required_argument, 0, 14},
      {"h",          required_argument, 0, 15},
      {"ivelo",          required_argument, 0, 16},
      {"omega",          required_argument, 0, 17},
      {"h0",          required_argument, 0, 18},
      {"TNy",          required_argument, 0, 19},
      {"shape",        required_argument, 0, 20},
      {"readtable",        required_argument, 0, 21},
      {NULL, 0, NULL, 0}
    };
  int c, option_index = 0;
  opterr = optind = 0;
  while ((c = getopt_long_only
          (argc, argv1.data(), "",
            long_options, &option_index)) != -1) {

    switch (c) {
    case 1: {
      // std::istringstream iss(optarg);
      // iss >> tst;
    } break;
    case 2: {
      std::istringstream iss(optarg);
      iss >> vs;
    } break;
    case 3: {
      // std::istringstream iss(optarg);
      // iss >> N;
    } break;
    case 4: {
      // std::istringstream iss(optarg);
      // iss >> K;
    } break;
    case 5: {
      // std::istringstream iss(optarg);
      // iss >> L;
    } break;
    case 6: {
      std::istringstream iss(optarg);
      iss >> TNx;
    } break;
    case 7: {
      std::istringstream iss(optarg);
      iss >> xmin;
    } break;
    case 8: {
      std::istringstream iss(optarg);
      iss >> xmax;
    } break;
    case 9: {
      std::istringstream iss(optarg);
      iss >> ymin;
    } break;
    case 10: {
      std::istringstream iss(optarg);
      iss >> ymax;
    } break;
    case 11: {
      std::istringstream iss(optarg);
      iss >> x0min;
    } break;
    case 12: {
      std::istringstream iss(optarg);
      iss >> x0max;
    } break;
    case 13: {
      std::istringstream iss(optarg);
      iss >> y0min;
    } break;
    case 14: {
      std::istringstream iss(optarg);
      iss >> y0max;
    } break;
    case 15: {
      std::istringstream iss(optarg);
      iss >> h;
    } break;
    case 16: {
      std::istringstream iss(optarg);
      iss >> ivelo;
    } break;
    case 17: {
      std::istringstream iss(optarg);
      iss >> w;
    } break;
    case 18: {
      std::istringstream iss(optarg);
      iss >> h0;
    } break;
    case 19: {
      std::istringstream iss(optarg);
      iss >> TNy;
    } break;
    case 20: {
      std::istringstream iss(optarg);
      iss >> shape;
    } break;
    case 21: {
      std::istringstream iss(optarg);
      iss >> readtable;
    } break;
    default: break;
    }
  }

	Ndim=2; //data dimension
  int Ns[Ndim];
  Nx = round((x0max-x0min)/h+1);
  Ny = round((y0max-y0min)/h+1);
  Ns[0]=Nx;
  Ns[1]=Ny;
  Nmax = max(Nx,Ny);


  // if(tst ==1){
	// 	vector<double> matU(M*rank_rand);
	// 	for (int i=0; i<M*rank_rand; i++)
	// 	matU[i] = (double)rand() / RAND_MAX;
	// 	MPI_Bcast(matU.data(), M*rank_rand, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// 	vector<double> matV(N*rank_rand);
	// 	for (int i=0; i<N*rank_rand; i++)
	// 	matV[i] = (double)rand() / RAND_MAX;
	// 	MPI_Bcast(matV.data(), N*rank_rand, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// 	quant_ptr_s2v=new C_QuantApp(M, N, rank_rand, 0, matU, matV);
	// }else if(tst ==2){
    // define the reference point for slowness
    slow_x0 = round((x0min+x0max)/2/h)*h;
    slow_y0 = round((y0min+y0max)/2/h)*h;


    // generate the chebyshev nodes
    double txmin=0, txmax=pi;
    double tymin=0, tymax=pi;
    double thx=(txmax-txmin)/(TNx-1);
    double thy=(tymax-tymin)/(TNy-1);
    double *tx=(double *)malloc(TNx*sizeof(double));
    double *ty=(double *)malloc(TNy*sizeof(double));

    vector<double> x_cheb(TNx);
    vector<double> y_cheb(TNy);
    for(int i=0;i<TNx;i++){
      tx[i]=txmin+i*thx;
      x_cheb[i]=cos(tx[i]);
    }
    for(int j=0;j<TNy;j++){
      ty[j]=tymin+j*thy;
      y_cheb[j]=cos(ty[j]);
    }
    //
    /**********************************************
    [-1 1]^2 to [x0min x0max]x[y0min y0max]
    ***********************************************/
    for(int i=0;i<TNx;i++){
      x_cheb[i]=(x0min+x0max)/2.0+(x0max-x0min)/2.0*x_cheb[i];
      // cout<<x_cheb[i]<<"\t";
    }
    // exit(0);
    for(int j=0;j<TNy;j++){
      y_cheb[j]=(y0min+y0max)/2.0+(y0max-y0min)/2.0*y_cheb[j];
    }

    vector<double> u1_square_int_cheb(TNx*TNy*TNx*TNy,0);
    vector<double> D1_int_cheb(TNx*TNy*TNx*TNy,0);
    vector<double> D2_int_cheb(TNx*TNy*TNx*TNy,0);


    string filename,filename1, str;
    filename ="ivelo_"+to_string(ivelo)+"_TNx_"+to_string(TNx)+"_TNy_"+to_string(TNy);
    std::ostringstream streamObj;
    streamObj << h0;
    str=streamObj.str();
    // str=to_string(h); // this only has 6-digit precision
    str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
    filename += "_h0_"+str;

	/*****************************************************************/
	int* groups = new int[size];
	int i_opt;
	double d_opt;
	int cpp=1; //1: use user-defined cpp/c functions for construction

	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);

	for (int i = 0; i < size; i++)groups[i]=i;
	// create hodlr data structures

	// double* dat_ptr;
	int* nns_ptr;
  int nquad=4;

	Nmin=100; //finest leafsize
	tol=1e-4; //compression tolerance
  double tol_rand=1e-2; //factorization tolerence
	com_opt=5; //1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton
	sort_opt=1; //0:natural order 1:kd-tree 2:cobble-like ordering 3:gram distance-based cobble-like ordering
	checkerr = 0; //1: check compression quality
	batch = 100; //batch size for BACA
	bnum = 1; //sqrt of #of subblocks in H-BACA
	knn=0; //k nearest neighbours stored per point
	double eps1 =1e-6;

	//quantities for the first holdr
	F2Cptr bmat_bf;  //hierarchical matrix returned by Fortran code
	F2Cptr option_bf;     //option structure returned by Fortran code
	F2Cptr stats_bf;      //statistics structure returned by Fortran code
	F2Cptr msh_bf;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant_bf;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree_bf;      //process tree returned by Fortran code


	z_c_bpack_createoption(&option_bf);

	// set hodlr options
	// z_c_bpack_set_D_option(&option_bf, "tol_comp", tol*1e-1);      // bf_a, bf_b, bf_c uses this tolerance
	z_c_bpack_set_D_option(&option_bf, "tol_comp", tol);      // bf_a, bf_b, bf_c uses this tolerance
	z_c_bpack_set_D_option(&option_bf, "tol_rand", tol_rand);           // bf_mv uses this tolerance
	z_c_bpack_set_D_option(&option_bf, "tol_Rdetect", tol_rand*3e-1);   // bf_mv uses this tolerance
	z_c_bpack_set_D_option(&option_bf, "sample_para", sample_para);
	z_c_bpack_set_D_option(&option_bf, "sample_para_outer", sample_para_outer);
	z_c_bpack_set_I_option(&option_bf, "nogeo", nogeo);
	z_c_bpack_set_I_option(&option_bf, "Nmin_leaf", Nmin);
	z_c_bpack_set_I_option(&option_bf, "RecLR_leaf", com_opt);
	z_c_bpack_set_I_option(&option_bf, "xyzsort", sort_opt);
	z_c_bpack_set_I_option(&option_bf, "ErrFillFull", checkerr);
	z_c_bpack_set_I_option(&option_bf, "BACA_Batch", batch);
	z_c_bpack_set_I_option(&option_bf, "LR_BLK_NUM", bnum);
	z_c_bpack_set_I_option(&option_bf, "cpp", cpp);
	z_c_bpack_set_I_option(&option_bf, "LRlevel", lrlevel);
	z_c_bpack_set_I_option(&option_bf, "knn", knn);
	z_c_bpack_set_I_option(&option_bf, "verbosity", verbose);
	z_c_bpack_set_I_option(&option_bf, "less_adapt", 1);
	z_c_bpack_set_I_option(&option_bf, "itermax", 10);
	z_c_bpack_set_I_option(&option_bf, "rmax", rmax);
	z_c_bpack_set_I_option(&option_bf, "ErrSol", 1);
	z_c_bpack_set_I_option(&option_bf, "elem_extract", elem_extract);
  set_option_from_command_line(argc, argv, option_bf);


  const int I=round((xmax-xmin)/h+1),  J=round((ymax-ymin)/h+1);
  const int Iint=round((x0max-x0min)/h+1),  Jint=round((y0max-y0min)/h+1);


  double smax=0;
  for(int i=0;i<Iint;i++){
    for(int j=0;j<Jint;j++){
      int ij = j*Iint+i;
      smax = max(smax,slowness(i*h+x0min,j*h+y0min,slow_x0, slow_y0,ivelo));
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,&smax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);



	/*****************************************************************/
	/* Test the Taylor+WENO-based 2D Green function */
	vector<double> data_geo;
	vector<double> panelnodes;
  double dl=0.01;
  int knn_pre;
  int layer=10;
  knn_pre= (2*layer+1)*(2*layer+1);

  double center[2];
  // center[0]=(x0min+x0max)/2.0;
  // center[1]=(y0min+y0max)/2.0;
  center[0]=0.4;
  center[1]=0.4;
  double radius_max=0.3;
  double L = 0.6;
  double H = 0.6;

    // radius_max=0.5;
    // center[0]=1.0;
    // center[1]=1.0;


  if(vs==1){
    data_geo.resize(Ndim*Nmax);
    for(int ii=0;ii<Nx;ii++){
      data_geo[(ii) * Ndim] = ii*h+x0min;
    }
    for(int ii=0;ii<Ny;ii++){
      data_geo[(ii) * Ndim+1] = ii*h+y0min;
    }
    if(myrank==0){
      cout<<"smax: "<<smax<<" PPW: "<<2*pi/(w*smax)/h<<" From: "<< Nx <<" x "<< Ny <<" To: "<< Nx <<" x "<< Ny <<endl;
    }
  }else{

  }

	/*****************************************************************/
	int myseg[Ndim];     // local number of unknowns
	int* perms_bf = new int[Nmax*Ndim]; //permutation vector returned by HODLR
	int nlevel = 0; // 0: tree level, nonzero if a tree is provided
	int* tree_bf = new int[(int)pow(2,nlevel)*2]; //user provided array containing size of each leaf node of each dimension, not used if nlevel=0
	tree_bf[0] = Nx;
	tree_bf[1] = Ny;

  if(vs==1){

  	C_QuantApp_BF *quant_ptr_bf;
    // create hodlr data structures
    z_c_bpack_createptree(&size, groups, &Fcomm, &ptree_bf);
    z_c_bpack_createstats(&stats_bf);
    z_c_bpack_set_I_option(&option_bf, "cpp", cpp);

    quant_ptr_bf=new C_QuantApp_BF(data_geo, Ndim, Nx, Ny, w, xmin, xmax, ymin, ymax, x0min, x0max, y0min, y0max, h, dl, nquad, 1,rmax,verbose,vs, x_cheb,y_cheb,u1_square_int_cheb,D1_int_cheb,D2_int_cheb);
      // construct hodlr with geometrical points
    z_c_bpack_md_construct_init(Ns,&Nmax, &Ndim, data_geo.data(), perms_bf, myseg, &bmat_bf, &option_bf, &stats_bf, &msh_bf, &kerquant_bf, &ptree_bf, &C_FuncNearFar_BF_MD, quant_ptr_bf);

    quant_ptr_bf->_Hperm.resize(Nmax*Ndim);
    std::copy(perms_bf, perms_bf + Nmax*Ndim, quant_ptr_bf->_Hperm.begin());

	  z_c_bpack_printoption(&option_bf,&ptree_bf);
    z_c_bpack_md_construct_element_compute(&Ndim,&bmat_bf, &option_bf, &stats_bf, &msh_bf, &kerquant_bf, &ptree_bf, &C_FuncZmn_BF_V2V_MD, &C_FuncZmnBlock_BF_V2V_MD, quant_ptr_bf);



    if(myrank==master_rank)std::cout<<"\n\nGenerating the incident fields: "<<std::endl;
    int nvec=2; // the 4th rhs makes precon=2 really hard to converge
    vector<_Complex double> b(product(myseg,Ndim)*nvec,{0.0,0.0});
    vector<_Complex double> x(product(myseg,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg,Ndim); i++){
      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_scalar;
      int i_old_md[Ndim];

      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf,i_new_loc_md,i_old_md);
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns,&i_old_scalar,i_old_md);

      double xs = data_geo[(i_old_md[0]-1) * Ndim];
      double ys = data_geo[(i_old_md[1]-1) * Ndim+1];
      double xs0=slow_x0;
      double ys0=slow_y0;
      for (int nth=0; nth<nvec; nth++){
        x.data()[i+nth*product(myseg,Ndim)]=source_function(xs,ys,x0max-0.15,y0max-0.15,nth,h,w);  // generate a source distribution
      }
    }
    z_c_bpack_md_mult(&Ndim,"N",x.data(),b.data(),myseg,myseg,&nvec,&bmat_bf,&option_bf,&stats_bf,&ptree_bf, &msh_bf);
    vector<_Complex double> u_inc_glo(product(Ns,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg,Ndim); i++){
      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_scalar;
      int i_old_md[Ndim];

      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf,i_new_loc_md,i_old_md);
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns,&i_old_scalar,i_old_md);

      for (int nth=0; nth<nvec; nth++){
        u_inc_glo.data()[i_old_scalar-1+nth*product(Ns,Ndim)] = b.data()[i+nth*product(myseg,Ndim)];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,u_inc_glo.data(), product(Ns,Ndim)*nvec, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    if(myrank==master_rank){
      for(int nth=0; nth<nvec; nth++){
        string filename, str;
        filename = "./VIE_F_inc_f_";
        str=to_string(w/2/pi);str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename +=str+"_vs_"+to_string(vs)+"_ivelo_"+to_string(ivelo);
        std::ostringstream streamObj;
        streamObj << h;
        str=streamObj.str();
        // str=to_string(h); // this only has 6-digit precision
        str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename += "_h_"+str;
        double opt_d;
        z_c_bpack_getoption(&option_bf, "tol_comp", &opt_d);

        str=my::to_string(opt_d);//str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename += "_tol_"+str+"_nth_"+to_string(nth)+"_tensor.bin";
        fout1=fopen(filename.c_str(),"wb");

        int nx = round((x0max-x0min)/h+1);
        int ny = round((y0max-y0min)/h+1);
        fwrite(&nx,sizeof(int),1,fout1);
        fwrite(&ny,sizeof(int),1,fout1);
        fwrite(&h,sizeof(double),1,fout1);
        fwrite(&u_inc_glo.data()[nth*product(Ns,Ndim)],sizeof(_Complex double),product(Ns,Ndim),fout1);
        fclose(fout1);
      }
    }


    if(shape!=4){
      cout<<"The tensor VIE cannot handle irregular shape at this moment"<<endl;
      exit(1);
    }
    int idx_off_x = round((-L/2.0 +center[0] - x0min)/h);
    int idx_off_y = round((-H/2.0 +center[1] - y0min)/h);

  int Ns_s[Ndim];
  int Nx_s, Ny_s, Nmax_s;
  Nx_s = round(L/h+1);
  Ny_s = round(H/h+1);
  Ns_s[0]=Nx_s;
  Ns_s[1]=Ny_s;
  Nmax_s = max(Nx_s,Ny_s);

  data_geo.resize(Ndim*Nmax_s);
  for(int ii=0;ii<Nx_s;ii++){
    data_geo[(ii) * Ndim] = ii*h-L/2+center[0];
  }
  for(int ii=0;ii<Ny_s;ii++){
    data_geo[(ii) * Ndim+1] = ii*h-H/2+center[1];
  }

    if(myrank==0){
      cout<<"smax: "<<smax<<" PPW: "<<2*pi/(w*smax)/h<<" From: "<< Nx_s <<" x "<< Ny_s <<" To: "<< Nx_s <<" x "<< Ny_s <<endl;
    }
	/*****************************************************************/
	int* perms_bf_s2s = new int[Nmax_s*Ndim]; //permutation vector returned by HODLR

  	C_QuantApp_BF *quant_ptr_bf_s2s;


    F2Cptr bmat_bf_s2s;  //hierarchical matrix returned by Fortran code
    F2Cptr stats_bf_s2s;      //statistics structure returned by Fortran code
    F2Cptr msh_bf_s2s;		   //d_mesh structure returned by Fortran code
    F2Cptr kerquant_bf_s2s;   //kernel quantities structure returned by Fortran code
    int myseg_s2s[Ndim];

    // create hodlr data structures
    z_c_bpack_createstats(&stats_bf_s2s);
    quant_ptr_bf_s2s=new C_QuantApp_BF(data_geo, Ndim, Nx_s, Ny_s, w, xmin, xmax, ymin, ymax, x0min, x0max, y0min, y0max, h, dl, nquad, ivelo,rmax,verbose,vs, x_cheb,y_cheb,u1_square_int_cheb,D1_int_cheb,D2_int_cheb);


    z_c_bpack_md_construct_init(Ns_s,&Nmax_s, &Ndim, data_geo.data(), perms_bf_s2s, myseg_s2s, &bmat_bf_s2s, &option_bf, &stats_bf_s2s, &msh_bf_s2s, &kerquant_bf_s2s, &ptree_bf, &C_FuncNearFar_BF_MD, quant_ptr_bf_s2s);

    quant_ptr_bf_s2s->_Hperm.resize(Nmax_s*Ndim);
    std::copy(perms_bf_s2s, perms_bf_s2s + Nmax_s*Ndim, quant_ptr_bf_s2s->_Hperm.begin());

	  z_c_bpack_printoption(&option_bf,&ptree_bf);
    z_c_bpack_md_construct_element_compute(&Ndim,&bmat_bf_s2s, &option_bf, &stats_bf_s2s, &msh_bf_s2s, &kerquant_bf_s2s, &ptree_bf, &C_FuncZmn_BF_S2S_MD, C_FuncZmnBlock_BF_S2S_MD, quant_ptr_bf_s2s);


#if 0
    if(myrank==master_rank)std::cout<<"\n\nFactoring the scatterer-scatterer operator: "<<std::endl;
    // factor hodlr
    z_c_bpack_factor(&bmat_bf_s2s,&option_bf,&stats_bf_s2s,&ptree_bf,&msh_bf_s2s);
#endif

    if(myrank==master_rank)std::cout<<"\n\nSolving the volume IE: "<<std::endl;
    vector<_Complex double> b_s(product(myseg_s2s,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg_s2s,Ndim); i++){

      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_scalar;
      int i_old_md[Ndim];

      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg_s2s,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf_s2s,i_new_loc_md,i_old_md);
      i_old_md[0] += idx_off_x;
      i_old_md[1] += idx_off_y;
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns,&i_old_scalar,i_old_md);
      for (int nth=0; nth<nvec; nth++){
        b_s[i+nth*product(myseg_s2s,Ndim)]=u_inc_glo[i_old_scalar-1+nth*product(Ns,Ndim)];
      }
    }
    int ErrSol=0;
    z_c_bpack_set_I_option(&option_bf, "ErrSol", ErrSol);
    vector<_Complex double> x_s(product(myseg_s2s,Ndim)*nvec,{0.0,0.0});
    z_c_bpack_md_solve(&Ndim, x_s.data(),b_s.data(),myseg_s2s,&nvec,&bmat_bf_s2s,&option_bf,&stats_bf_s2s,&ptree_bf, &msh_bf_s2s);


    vector<_Complex double> x_v_glo(product(Ns,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg_s2s,Ndim); i++){
      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_scalar;
      int i_old_md[Ndim];

      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg_s2s,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf_s2s,i_new_loc_md,i_old_md);
      double xs = data_geo[(i_old_md[0]-1) * Ndim];
      double ys = data_geo[(i_old_md[1]-1) * Ndim+1];
      i_old_md[0] += idx_off_x;
      i_old_md[1] += idx_off_y;
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns,&i_old_scalar,i_old_md);

      for (int nth=0; nth<nvec; nth++){
        double ss = slowness(xs,ys,slow_x0, slow_y0,ivelo);
        double s0=2;
        double k0 = s0*w;
        double coef = pow(k0,2.0)*(pow(ss/s0,2.0)-1);
        x_v_glo[i_old_scalar-1+nth*product(Ns,Ndim)]=x_s[i+nth*product(myseg_s2s,Ndim)]*coef;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,x_v_glo.data(), product(Ns,Ndim)*nvec, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);


    vector<_Complex double> x_v(product(myseg,Ndim)*nvec,{0.0,0.0}),b_v(product(myseg,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg,Ndim); i++){

      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_scalar;
      int i_old_md[Ndim];

      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf,i_new_loc_md,i_old_md);
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns,&i_old_scalar,i_old_md);
      for (int nth=0; nth<nvec; nth++){
        x_v[i+nth*product(myseg,Ndim)] = x_v_glo[i_old_scalar-1+nth*product(Ns,Ndim)];
      }
    }

    z_c_bpack_md_mult(&Ndim,"N",x_v.data(),b_v.data(),myseg,myseg,&nvec,&bmat_bf,&option_bf,&stats_bf,&ptree_bf,&msh_bf);
    vector<_Complex double> u_sca_glo(product(Ns,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg,Ndim); i++){
      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_scalar;
      int i_old_md[Ndim];

      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf,i_new_loc_md,i_old_md);
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns,&i_old_scalar,i_old_md);
      for (int nth=0; nth<nvec; nth++){
        u_sca_glo.data()[i_old_scalar-1+nth*product(Ns,Ndim)] = b_v.data()[i+nth*product(myseg,Ndim)];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,u_sca_glo.data(), product(Ns,Ndim)*nvec, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);


    if(myrank==master_rank){
      for(int nth=0; nth<nvec; nth++){
        string filename, str;
        filename = "./VIE_F_sca_f_";
        str=to_string(w/2/pi);str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename +=str+"_vs_"+to_string((int)0)+"_ivelo_"+to_string(ivelo);
        if(shape>0)
          filename +="_shape_"+to_string(shape);
        std::ostringstream streamObj;
        streamObj << h;
        str=streamObj.str();
        // str=to_string(h); // this only has 6-digit precision
        str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename += "_h_"+str;
        double opt_d;
        z_c_bpack_getoption(&option_bf, "tol_comp", &opt_d);

        str=my::to_string(opt_d);//str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename += "_tol_"+str+"_nth_"+to_string(nth)+"_tensor.bin";
#ifdef IVELO9_CONST
        filename +="_ivelo9_const";
#endif
        fout1=fopen(filename.c_str(),"wb");

        int nx = round((x0max-x0min)/h+1);
        int ny = round((y0max-y0min)/h+1);
        fwrite(&nx,sizeof(int),1,fout1);
        fwrite(&ny,sizeof(int),1,fout1);
        fwrite(&h,sizeof(double),1,fout1);
        fwrite(&u_sca_glo.data()[nth*product(Ns,Ndim)],sizeof(_Complex double),product(Ns,Ndim),fout1);
        fclose(fout1);
      }
    }

    if(myrank==master_rank)std::cout<<"\n\nPrinting stats of the volume-volume operator: "<<std::endl;
    z_c_bpack_printstats(&stats_bf,&ptree_bf);

    if(myrank==master_rank)std::cout<<"\n\nPrinting stats of the scatterer-scatterer operator: "<<std::endl;
    z_c_bpack_printstats(&stats_bf_s2s,&ptree_bf);

    delete quant_ptr_bf_s2s;
    delete quant_ptr_bf;
    z_c_bpack_deletestats(&stats_bf_s2s);
    z_c_bpack_md_deletemesh(&Ndim,&msh_bf_s2s);
    z_c_bpack_deletekernelquant(&kerquant_bf_s2s);
    z_c_bpack_delete(&bmat_bf_s2s);
    delete[] perms_bf_s2s;



  }else{

    cout<<"shouldn't reach here"<<endl;

  }

	// // solve the system
	// int nrhs=1;
	// _Complex double* b = new _Complex double[nrhs*myseg];
	// _Complex double* x = new _Complex double[nrhs*myseg];

	// for (int i = 0; i < nrhs*myseg; i++){
	// 	b[i]=1;
	// }
	// z_c_bpack_solve(x,b,&myseg,&nrhs,&bmat,&option,&stats,&ptree);




	z_c_bpack_deletestats(&stats_bf);
	z_c_bpack_deleteproctree(&ptree_bf);
	z_c_bpack_md_deletemesh(&Ndim,&msh_bf);
	z_c_bpack_deletekernelquant(&kerquant_bf);
	z_c_bpack_delete(&bmat_bf);
	z_c_bpack_deleteoption(&option_bf);



	delete[] perms_bf;
	delete[] tree_bf;
	delete[] groups;
	free(tx);
	free(ty);


	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
// // //------------------------------------------------------------------------------











// int main()
// {
//   // define the sloweness type
//   int ivelo = 9;
//   int ker=3;
//   int verbose=1;

//   // define the inner domain
//   double x0min=-0.0, x0max=1.0/1.0;
//   double y0min=-0.0, y0max=1.0/1.0;

//   // define the outer domain
//   double xouter=0.2, youter=0.2;
//   double xmin=x0min-xouter,  xmax=x0max+xouter;
//   double ymin=y0min-youter,  ymax=y0max+youter;

//   // define the reference point for slowness
//   double slow_x0 = (x0min+x0max)/2.0;
//   double slow_y0 = (y0min+y0max)/2.0;

//   // define the grid size
//   int I=141,  J=I;
//   double h=(xmax-xmin)/(I-1);

//   // generate the chebyshev nodes
//   int TN=10;
//   int TNx=TN;
//   int TNy=TN;
//   double txmin=0, txmax=PI;
//   double tymin=0, tymax=PI;
//   double thx=(txmax-txmin)/(TNx-1);
//   double thy=(tymax-tymin)/(TNy-1);
//   double *tx=(double *)malloc(TNx*sizeof(double));
//   double *ty=(double *)malloc(TNy*sizeof(double));
//   double *x=(double *)malloc(TNx*sizeof(double));
//   double *y=(double *)malloc(TNy*sizeof(double));
//   for(int i=0;i<TNx;i++){
//      tx[i]=txmin+i*thx;
//      x[i]=cos(tx[i]);
//   }
//   for(int j=0;j<TNy;j++){
//      ty[j]=tymin+j*thy;
//      y[j]=cos(ty[j]);
//   }
//   //
//   /**********************************************
//   [-1 1]^3 to [x0min x0max]x[y0min y0max]x[z0min z0max]
//   ***********************************************/
//   for(int i=0;i<TNx;i++){
//      x[i]=(x0min+x0max)/2.0+(x0max-x0min)/2.0*x[i];
//   }
//   for(int j=0;j<TNy;j++){
//      y[j]=(y0min+y0max)/2.0+(y0max-y0min)/2.0*y[j];
//   }

// double u1_square_int_cheb[TNx*TNy*TNx*TNy];
// double D1_int_cheb[TNx*TNy*TNx*TNy];
// double D2_int_cheb[TNx*TNy*TNx*TNy];

// for(int si=0;si<TNx;si++){
//   for(int sj=0;sj<TNy;sj++){
//   // for(int si=0;si<1;si++){
//   //   for(int sj=0;sj<1;sj++){
//       int sij = sj*TNx+si;
//       double sx = x[si];
//       double sy = y[sj];

//       // define shifted inner domain
//       double x0min1=sx - ceil_safe((sx-x0min)/h)*h;
//       double x0max1=sx + ceil_safe((x0max-sx)/h)*h;
//       double y0min1=sy - ceil_safe((sy-y0min)/h)*h;
//       double y0max1=sy + ceil_safe((y0max-sy)/h)*h;

//       // define the outer domain
//       double xmin1=x0min1-xouter,  xmax1=x0max1+xouter;
//       double ymin1=y0min1-youter,  ymax1=y0max1+youter;

//       // define the grid size
//       const int I1=(xmax1-xmin1)/h+1,  J1=(ymax1-ymin1)/h+1;
//       const int Iint1=(x0max1-x0min1)/h+1,  Jint1=(y0max1-y0min1)/h+1;


//       // deine the coordinates of the inner domain grids
//       double *x1=(double *)malloc(Iint1*sizeof(double));
//       double *y1=(double *)malloc(Jint1*sizeof(double));

//       for(int i=0;i<Iint1;i++){
//           x1[i] = i*h+x0min1;
//       }
//       for(int j=0;j<Jint1;j++){
//           y1[j] = j*h+y0min1;
//       }

//       int six,siy;

//       six = (sx-xmin1)/h;
//       siy = (sy-ymin1)/h;

//       double u1_square_int[Jint1*Iint1];
//       double D1_int[Jint1*Iint1];
//       double D2_int[Jint1*Iint1];
//       cout<<six<<" "<<siy<<" "<<I1<<" "<<J1<<endl;
//       compute_one_col(xmin1,xmax1,ymin1,ymax1,x0min1,x0max1,y0min1,y0max1,slow_x0,slow_y0,I1,J1,six,siy,ivelo,u1_square_int,D1_int,D2_int,Iint1,Jint1,ker,verbose);

//       double* u1_square_int_cheb1 = u1_square_int_cheb+TNx*TNy*sij;
//       CubicInterp2D(x1, y1, u1_square_int, Iint1, Jint1, x, y, u1_square_int_cheb1, TNx, TNy);

//       double* D1_int_cheb1 = D1_int_cheb+TNx*TNy*sij;
//       CubicInterp2D(x1, y1, D1_int, Iint1, Jint1, x, y, D1_int_cheb1, TNx, TNy);

//       double* D2_int_cheb1 = D2_int_cheb+TNx*TNy*sij;
//       CubicInterp2D(x1, y1, D2_int, Iint1, Jint1, x, y, D2_int_cheb1, TNx, TNy);


//       // for(int i=0; i<TNx*TNy; i++)
//       // {
//       //         cout<<u1_square_int_cheb1[i]<<"\t"<<x[i%TNx]<<"\t"<<y[i/TNx]<<endl;
//       // }
//       // cout<<endl;

//       delete[] x1;
//       delete[] y1;
//     }
//   }

//   //test the interpolation
//   I=401;
//   J=I;
//   h=(xmax-xmin)/(I-1);
//   int Iint=(x0max-x0min)/h+1;
//   int Jint=(y0max-y0min)/h+1;
//   int idx = (int)((double) rand() / (RAND_MAX)*Iint*Jint); // generate a random source point

//   int idx_x = idx%Iint;
//   int idx_y = idx/Iint;
//   int six,siy;
//   six = (idx_x*h+x0min-xmin)/h;
//   siy = (idx_y*h+y0min-ymin)/h;
//   double* u1_square_int_ref = new double[Jint*Iint];
//   double* D1_int_ref = new double[Jint*Iint];
//   double* D2_int_ref = new double[Jint*Iint];
//   compute_one_col(xmin,xmax,ymin,ymax,x0min,x0max,y0min,y0max,slow_x0, slow_y0, I,J,six,siy,ivelo,u1_square_int_ref,D1_int_ref,D2_int_ref,Iint,Jint,ker, verbose);

//   {
//   double* u1_square_int_1 = new double[Jint*Iint];
//   for(int i=0;i<Iint;i++){
//     for(int j=0;j<Jint;j++){
//       int ij = j*Iint+i;
//       u1_square_int_1[ij] = lagrange_interp2D(x, y, i*h+x0min,j*h+y0min,idx_x*h+x0min,idx_y*h+y0min,u1_square_int_cheb, TNx,TNy);
//     }
//   }
//   double norm1=0;
//   double norm2=0;
//   for(int i=0;i<Jint*Iint;i++){
//     norm1 += u1_square_int_ref[i]*u1_square_int_ref[i];
//     norm2 += (u1_square_int_ref[i]-u1_square_int_1[i])*(u1_square_int_ref[i]-u1_square_int_1[i]);
//   }
//   cout<<"interpolation error (u1_square): "<<norm2/norm1<<endl;
//   delete[] u1_square_int_1;
//   }

//   {
//   double* D1_int_1 = new double[Jint*Iint];
//   for(int i=0;i<Iint;i++){
//     for(int j=0;j<Jint;j++){
//       int ij = j*Iint+i;
//       D1_int_1[ij] = lagrange_interp2D(x, y, i*h+x0min,j*h+y0min,idx_x*h+x0min,idx_y*h+y0min,D1_int_cheb, TNx,TNy);
//     }
//   }
//   double norm1=0;
//   double norm2=0;
//   for(int i=0;i<Jint*Iint;i++){
//     norm1 += D1_int_ref[i]*D1_int_ref[i];
//     norm2 += (D1_int_ref[i]-D1_int_1[i])*(D1_int_ref[i]-D1_int_1[i]);
//   }
//   cout<<"interpolation error (D1): "<<norm2/norm1<<endl;
//   delete[] D1_int_1;
//   }


//   {
//   double* D2_int_1 = new double[Jint*Iint];
//   for(int i=0;i<Iint;i++){
//     for(int j=0;j<Jint;j++){
//       int ij = j*Iint+i;
//       D2_int_1[ij] = lagrange_interp2D(x, y, i*h+x0min,j*h+y0min,idx_x*h+x0min,idx_y*h+y0min,D2_int_cheb, TNx,TNy);
//     }
//   }
//   double norm1=0;
//   double norm2=0;
//   for(int i=0;i<Jint*Iint;i++){
//     norm1 += D2_int_ref[i]*D2_int_ref[i];
//     norm2 += (D2_int_ref[i]-D2_int_1[i])*(D2_int_ref[i]-D2_int_1[i]);
//   }
//   cout<<"interpolation error (D2): "<<norm2/norm1<<endl;
//   delete[] D2_int_1;
//   }



//   delete[] u1_square_int_ref;
//   delete[] D1_int_ref;
//   delete[] D2_int_ref;

// }