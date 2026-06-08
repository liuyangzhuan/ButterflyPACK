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
#include <unordered_map>

#include <cmath>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include <atomic>
#include <limits>
#include <mpi.h>
// #include <complex.h>

#include <sstream>
#include <cstring>
#include <getopt.h>
#include <unistd.h>

#include "G2D/bessel.h"
#include "zBPACK_wrapper.h"

// #define IVELO9_CONST 1

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

int product(const int arr[], int n) {
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

// This pads n to something like 2^k*m, where m is an approximate leaf size with leafsize<=m<leafsize*2.
int pad_to_leaf_power2(int n, int leafsize){
  if(leafsize <= 1) return n;
  long long best = std::numeric_limits<long long>::max();
  for(long long pow2=1; pow2<=std::numeric_limits<int>::max(); pow2*=2){
    long long m = (n + pow2 - 1)/pow2;
    m = max<long long>(m, leafsize);
    if(m < 2LL*leafsize){
      long long padded = m*pow2;
      if(padded >= n) best = min(best,padded);
    }
    if(best != std::numeric_limits<long long>::max() && pow2 > n) break;
    if(pow2 > std::numeric_limits<int>::max()/2) break;
  }
  if(best > std::numeric_limits<int>::max()){
    cout<<"pad_to_leaf_power2 overflow for n="<<n<<" leafsize="<<leafsize<<endl;
    exit(1);
  }
  return static_cast<int>(best);
}

int effective_leaf_size_for_padded_dim(int n, int leafsize){
  if(leafsize <= 1) return n;
  int pow2 = 1;
  while(n%(2*pow2)==0 && n/(2*pow2) >= leafsize){
    if(pow2 > std::numeric_limits<int>::max()/2) break;
    pow2 *= 2;
  }
  int m = n/pow2;
  if(m < leafsize || m >= 2*leafsize){
    cout<<"invalid padded dimension "<<n<<" for leafsize "<<leafsize<<endl;
    exit(1);
  }
  return m;
}

int tree_levels_for_padded_dim(int n, int leafsize){
  if(leafsize <= 1) return 0;
  int pow2 = 1;
  int levels = 0;
  while(n%(2*pow2)==0 && n/(2*pow2) >= leafsize){
    if(pow2 > std::numeric_limits<int>::max()/2) break;
    pow2 *= 2;
    levels++;
  }
  return levels;
}

bool in_physical_tensor_domain(int const md[], int nx, int ny, int nz){
  return md[0] >= 1 && md[0] <= nx &&
         md[1] >= 1 && md[1] <= ny &&
         md[2] >= 1 && md[2] <= nz;
}

int physical_tensor_scalar(int const md[], int nx, int ny){
  return (md[2]-1)*nx*ny + (md[1]-1)*nx + md[0];
}

int plot_slice_index(int n){
  return std::min(n, std::max(1, (n-1)/2 + 1));
}

size_t tensor_slice_size(int nx, int ny, int nz){
  return static_cast<size_t>(ny)*nz + static_cast<size_t>(nx)*nz + static_cast<size_t>(nx)*ny;
}

struct DiagRouteBlock {
  int head;
  int tail;
  int root;
  int gid;
};

void collect_tensor_slices_from_local(int Ndim, const int Ns_phys[], int myseg[],
                                      F2Cptr* msh, const vector<_Complex double>& field_loc,
                                      int nvec, double scale, int ix, int iy, int iz,
                                      vector<_Complex double>& slices_loc){
  int nx = Ns_phys[0];
  int ny = Ns_phys[1];
  int nz = Ns_phys[2];
  size_t one_rhs = tensor_slice_size(nx,ny,nz);
  slices_loc.assign(one_rhs*nvec,{0.0,0.0});

  int nloc = product(myseg,Ndim);
  for(int i=0;i<nloc;i++){
    int i_new_loc_scalar = i+1;
    int i_new_loc_md[3];
    int i_old_md[3];
    z_c_bpack_singleindex_to_multiindex(&Ndim,myseg,&i_new_loc_scalar,i_new_loc_md);
    z_c_bpack_md_new2old(&Ndim,msh,i_new_loc_md,i_old_md);
    if(!in_physical_tensor_domain(i_old_md,nx,ny,nz)) continue;

    int ii = i_old_md[0];
    int jj = i_old_md[1];
    int kk = i_old_md[2];
    for(int nth=0;nth<nvec;nth++){
      _Complex double val = scale*field_loc[i+nth*nloc];
      _Complex double* rhs_slices = slices_loc.data() + nth*one_rhs;
      if(ii == ix){
        size_t off = static_cast<size_t>(jj-1) + static_cast<size_t>(ny)*(kk-1);
        rhs_slices[off] = val;
      }
      if(jj == iy){
        size_t off0 = static_cast<size_t>(ny)*nz;
        size_t off = off0 + static_cast<size_t>(ii-1) + static_cast<size_t>(nx)*(kk-1);
        rhs_slices[off] = val;
      }
      if(kk == iz){
        size_t off0 = static_cast<size_t>(ny)*nz + static_cast<size_t>(nx)*nz;
        size_t off = off0 + static_cast<size_t>(ii-1) + static_cast<size_t>(nx)*(jj-1);
        rhs_slices[off] = val;
      }
    }
  }
}

void write_tensor_slices_file(const string& filename, int nx, int ny, int nz, double h,
                              int ix, int iy, int iz, const _Complex double* slices){
  FILE* fout = fopen(filename.c_str(),"wb");
  if(!fout){
    cout<<"Cannot open slice output file "<<filename<<endl;
    exit(1);
  }
  const int magic = 20260605;
  const int version = 1;
  fwrite(&magic,sizeof(int),1,fout);
  fwrite(&version,sizeof(int),1,fout);
  fwrite(&nx,sizeof(int),1,fout);
  fwrite(&ny,sizeof(int),1,fout);
  fwrite(&nz,sizeof(int),1,fout);
  fwrite(&ix,sizeof(int),1,fout);
  fwrite(&iy,sizeof(int),1,fout);
  fwrite(&iz,sizeof(int),1,fout);
  fwrite(&h,sizeof(double),1,fout);
  fwrite(slices,sizeof(_Complex double),tensor_slice_size(nx,ny,nz),fout);
  fclose(fout);
}



double slowness(double x,double y, double z, double slow_x0, double slow_y0,double slow_z0, int ivelo, const double* slowness_array, double h, int I, int J, int K)
{
  double g1, g2, g3, s0;
  g1=-0.2; g2=-0.4; g3=-0.35;
  s0=2.0;  // This is at the domain reference point (slow_x0,slow_y0)

  double A = -0.01;
  if(ivelo==9){
  #ifdef IVELO9_CONST
      s0 = 1.0;
  #else
      s0 = 1.68;
      s0 =1/((1.0/s0+g1*(x-slow_x0)+g2*(y-slow_y0)+g3*(z-slow_z0)));
  #endif
  }else if(ivelo==10){
      s0=2.0*sqrt(1+A*z);
  }else if(ivelo==11){
      int i=round(x/h);
      int j=round(y/h);
      int k=round(z/h);
      s0=slowness_array[k*I*J + j*I + i];

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

_Complex double source_function(double x,double y, double z, double xs0, double ys0, double zs0, int nth_rhs, double h, double omega)
{
  _Complex double out = {0.0,0.0};
  double window;
  // the first rhs is point source
  if(nth_rhs==0){
  if((fabs(x-xs0)<1e-10 && fabs(y-ys0)<1e-10 && fabs(z-zs0)<1e-10)){
    out = {1.0,0.0};
  }
  }

  // the second rhs is gaussian pulse
  if(nth_rhs==1){
    double sigma=0.15;
    out = exp(-(pow(x-xs0,2)+pow(y-ys0,2)+pow(z-zs0,2))/pow(sigma,2)/2)*pow(h,3);
    double taperwidth=0.05;
    double constwidth=0.05;
    double dist = sqrt(pow(x-xs0,2)+pow(y-ys0,2)+pow(z-zs0,2));
    double window =TukeyWindow(dist,constwidth,taperwidth);
    out = out*window;
  }

  // the third rhs is all-1 inside a small circle
  if(nth_rhs==2){
    double dist = sqrt(pow(x-xs0,2)+pow(y-ys0,2)+pow(z-zs0,2));
    double taperwidth=0.05;
    double constwidth=0.05;
    out = TukeyWindow(dist,constwidth,taperwidth);
    out *= pow(h,3);
  }

  // the fourth rhs is gaussian wave packet
  if(nth_rhs==3){
    double omega0=omega*0.9;
    double dist = sqrt(pow(x-xs0,2)+pow(y-ys0,2)+pow(z-zs0,2));
    double d[3];
    d[0]=1/sqrt(3.0);
    d[1]=1/sqrt(3.0);
    d[2]=1/sqrt(3.0);

    double sigma=0.15;
    double phase = omega0*(x*d[0]+y*d[1]+z*d[2]);
    out = exp(-(pow(x-xs0,2)+pow(y-ys0,2)+pow(z-zs0,2))/pow(sigma,2)/2)*(cos(phase)+Im*sin(phase))*pow(h,3);
    double taperwidth=0.05;
    double constwidth=0.05;
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


  // // the fifth rhs is all-1 inside a concave kite shaped region
  // if(nth_rhs==4){
  //   double axis_ref[2] = {1, 0};
  //   double origin[2] = {xs0, ys0};
  //   double scale=0.2;
  //   double taperwidth=0.1;

  //   double v1[2]={1,0};
  //   double v2[2]={axis_ref[0],axis_ref[1]};
  //   double ang = atan2(v1[0]*v2[1]-v2[0]*v1[1],v1[0]*v2[0]+v1[1]*v2[1]);
  //   double theta = std::fmod(ang, 2*pi);
  //   double rmax = 2.065670838*scale + taperwidth*1.1;
  //   double v0[2]={x-origin[0],y-origin[1]};
  //   double r = sqrt(pow(v0[0],2)+pow(v0[1],2));

  //   if(r<rmax){
  //     double xm,ym;
  //     int Ncurv=1000;
  //     double angmin=pi;
  //     int imin;
  //     for(int ii=0;ii<Ncurv;ii++){ //loop through all points on the curve, find the one colinear with [x,y]

  //       double t = 2*pi/1000*ii;
  //       double xc = scale*(cos(t)+0.65*cos(2*t)-0.65);
  //       double yc = scale*(1.5*sin(t));
  //       xc = xc*cos(theta)-yc*sin(theta)+origin[0];
  //       yc = xc*sin(theta)+yc*cos(theta)+origin[1];

  //       v1[0] = xc-origin[0];
  //       v1[1] = yc-origin[1];
  //       double ang = fabs(atan2(v1[0]*v0[1]-v0[0]*v1[1],v1[0]*v0[0]+v1[1]*v0[1]));
  //       if(angmin>ang){
  //         imin=ii;
  //         angmin=ang;
  //         xm=xc;
  //         ym=yc;
  //       }
  //     }
  //     double vref[2]={xm-origin[0],ym-origin[1]};
  //     double rref = sqrt(pow(vref[0],2)+pow(vref[1],2));
  //     double constwidth=rref;
  //     out = TukeyWindow(r,constwidth,taperwidth);
  //     out *= pow(h,2);
  //   }else{
  //     out=0;
  //   }
  // }

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
  int _bdiag_precon = 0; // whether to construct block-diagonal preconditioners
  int _d = 0;   // data dimension 2 or 3
  int _n = 0;   // size of the matrix
  int _nx = 0;   // nx
  int _ny = 0;   // ny
  int _nz = 0;   // nz
  int _slow_nx = 0;   // slowness-array x dimension, excluding tensor padding
  int _slow_ny = 0;   // slowness-array y dimension, excluding tensor padding
  int _slow_nz = 0;   // slowness-array z dimension, excluding tensor padding
  int _scaleGreen = 0; // with v=Q*u unknown: whether to apply h^{-3}/Q diagonal in matvec instead of sampled entries
  double _w = pi;   //angular frequency
  int _idx_off_x = 0;
  int _idx_off_y = 0;
  int _idx_off_z = 0;

  double _xmin=-0.2/1.0,  _xmax=1.2/1.0;
  double _ymin=-0.2/1.0,  _ymax=1.2/1.0;
  double _zmin=-0.2/1.0,  _zmax=1.2/1.0;
  double _x0min=-0.0, _x0max=1.0/1.0;
  double _y0min=-0.0, _y0max=1.0/1.0;
  double _z0min=-0.0, _z0max=1.0/1.0;
  double _slow_x0 = 0.5;
  double _slow_y0 = 0.5;
  double _slow_z0 = 0.5;
  double _h= 0.01;
  double _dl= 0.01;
  double _rmax= 10;
  int _ivelo = 9; // 1 homogenous, 2 constant gradient
  int _TNx, _TNy, _TNz;
  int _nquad = 4;

  int _verbose=0;
  int _vs=1;
  F2Cptr* bmat_scaled_green, *option_bf, *stats_bf, *ptree_bf, *msh_bf;
  F2Cptr* bmat_diags, *stats_diags, *msh_diags, *kerquant_diags, *ptree_diags;
  F2Cptr option_diag, ptree_diag;
  C_QuantApp_BF** quant_ptr_diags;



  std::vector<double> _x_cheb;
  std::vector<double> _y_cheb;
  std::vector<double> _z_cheb;
  std::vector<double> _u1_square_int_cheb;
  std::vector<double> _D1_int_cheb;
  std::vector<double> _D2_int_cheb;

  std::vector<int> _Hperm;
  std::vector<int> _iHperm;

  std::vector<int> _Hperm_m;
  std::vector<int> _headarray;
  int _nc_m[3];
  int _diag_block_stride = 0;
  int _nblock_diag = 0;
  int _bdiag_level = -1;
  int _nscatter = 0;
  double _background_s0 = 2.0;
  double _scatterer_tol = 1e-12;
  int _outer_idxs[3] = {1, 1, 1};
  int _outer_idxe_plus1[3] = {1, 1, 1};
  std::vector<int> _diag_block_info;
  std::vector<int> _diag_block_offsets;
  std::vector<int> _diag_tensor_ids;
  std::vector<int> _diag_nloc;
  std::vector<int> _diag_npo;
  std::vector<MPI_Comm> _diag_comms;
  std::vector<MPI_Group> _diag_groups;
  std::vector<int> _iHperm_m;
  std::vector<int> _scatter_to_tensor;
  std::vector<int> _local_tensor_ids;
  std::vector<int> _local_scatter_ids;
  std::vector<int> _local_diag_root;
  std::vector<int> _local_diag_gid;
  std::vector<int> _local_diag_pos;
  const vector<double>* _slowness_array = nullptr;

  C_QuantApp_BF() = default;

  inline bool IsScattererTensorScalar(int tensor_scalar) const {
    return ScatterIndexFromTensorScalar(tensor_scalar) > 0;
  }

  inline int ScatterIndexFromTensorScalar(int tensor_scalar) const {
    auto it = std::lower_bound(_scatter_to_tensor.begin(), _scatter_to_tensor.end(), tensor_scalar);
    if(it != _scatter_to_tensor.end() && *it == tensor_scalar){
      return static_cast<int>(it - _scatter_to_tensor.begin()) + 1;
    }
    return 0;
  }

  inline int LocalScatterIndex(int local_idx) const {
    if(local_idx >= 0 && local_idx < static_cast<int>(_local_scatter_ids.size())){
      return _local_scatter_ids[local_idx];
    }
    return 0;
  }

  // constructor for the tensor operator
  C_QuantApp_BF(vector<double> data, int d, int nx, int ny, int nz, int slow_nx, int slow_ny, int slow_nz, int scaleGreen, double w, double x0min, double x0max, double y0min, double y0max, double z0min, double z0max, double h, double dl, int ivelo, const vector<double>* slowness_array, int rmax, int verbose, int vs, int bdiag_precon, vector<double> x_cheb, vector<double> y_cheb, vector<double> z_cheb, vector<double> u1_square_int_cheb, vector<double> D1_int_cheb, vector<double> D2_int_cheb)
    :_data(move(data)), _d(d), _nx(nx), _ny(ny), _nz(nz), _slow_nx(slow_nx), _slow_ny(slow_ny), _slow_nz(slow_nz), _scaleGreen(scaleGreen), _w(w),_x0min(x0min),_x0max(x0max),_y0min(y0min),_y0max(y0max),_z0min(z0min),_z0max(z0max),_h(h),_dl(dl),_ivelo(ivelo),_slowness_array(slowness_array), _rmax(rmax),_verbose(verbose),_vs(vs),_bdiag_precon(bdiag_precon),_x_cheb(move(x_cheb)),_y_cheb(move(y_cheb)),_z_cheb(move(z_cheb)),_u1_square_int_cheb(move(u1_square_int_cheb)),_D1_int_cheb(move(D1_int_cheb)),_D2_int_cheb(move(D2_int_cheb)){
      _TNx = _x_cheb.size();
      _TNy = _y_cheb.size();
      _TNz = _z_cheb.size();
      _slow_x0 = (_x0min+_x0max)/2.0;
      _slow_y0 = (_y0min+_y0max)/2.0;
      _slow_z0 = (_z0min+_z0max)/2.0;
      _u1_square_int_cheb.insert(_u1_square_int_cheb.end(),_D1_int_cheb.begin(),_D1_int_cheb.end());   // concatenate D1 and D2 into u1_square for later interpolation convenience
      _u1_square_int_cheb.insert(_u1_square_int_cheb.end(),_D2_int_cheb.begin(),_D2_int_cheb.end());
	}


  // constructor for the matrix operator
  C_QuantApp_BF(vector<double> data, int d, int scaleGreen, double w, double x0min, double x0max, double y0min, double y0max, double z0min, double z0max, double h, double dl, int ivelo, const vector<double>* slowness_array, int rmax, int verbose, int vs, vector<double> x_cheb, vector<double> y_cheb, vector<double> z_cheb, vector<double> u1_square_int_cheb, vector<double> D1_int_cheb, vector<double> D2_int_cheb)
    :_data(move(data)), _d(d), _n(_data.size() / _d), _scaleGreen(scaleGreen), _w(w),_x0min(x0min),_x0max(x0max),_y0min(y0min),_y0max(y0max),_z0min(z0min),_z0max(z0max),_h(h),_dl(dl),_ivelo(ivelo),_slowness_array(slowness_array), _rmax(rmax),_verbose(verbose),_vs(vs),_x_cheb(move(x_cheb)),_y_cheb(move(y_cheb)),_z_cheb(move(z_cheb)),_u1_square_int_cheb(move(u1_square_int_cheb)),_D1_int_cheb(move(D1_int_cheb)),_D2_int_cheb(move(D2_int_cheb)){
      _TNx = _x_cheb.size();
      _TNy = _y_cheb.size();
      _TNz = _z_cheb.size();
      _slow_x0 = (_x0min+_x0max)/2.0;
      _slow_y0 = (_y0min+_y0max)/2.0;
      _slow_z0 = (_z0min+_z0max)/2.0;
      _nx = round((_x0max-_x0min)/_h);
      _ny = round((_y0max-_y0min)/_h);
      _nz = round((_z0max-_z0min)/_h);
      _slow_nx = _nx;
      _slow_ny = _ny;
      _slow_nz = _nz;
      _u1_square_int_cheb.insert(_u1_square_int_cheb.end(),_D1_int_cheb.begin(),_D1_int_cheb.end());   // concatenate D1 and D2 into u1_square for later interpolation convenience
      _u1_square_int_cheb.insert(_u1_square_int_cheb.end(),_D2_int_cheb.begin(),_D2_int_cheb.end());
	}



  inline void SampleSelf(double x1, double y1, double z1, double x2, double y2,double z2, _Complex double* val){

    int self = sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))<1e-20? 1:0;
    if(_vs==1){
      // int closeby = fabs(x1-x2)<_h*2+1e-20 && fabs(y1-y2)<_h*2+1e-20?1:0;
      if(self==1){
        // double s0 = slowness(x1,y1,z1, _slow_x0, _slow_y0,_slow_z0,_ivelo,_slowness_array.data(),_h, round(_x0max/_h), round(_y0max/_h), round(_z0max/_h));
        double s0 = 2;
        double gamma = 1.781072418;


        // 7-point Legendre-Gauss Quadrature on [0,pi/4] for the integral due to Jianliang Qian
        double nodes_phi[7] = {0.765412887308718, 0.683897697334573, 0.552074099956508, 0.392699081698724, 0.233324063440941, 0.101500466062876, 0.019985276088730};
        double weights_phi[7] = {0.050848627308305, 0.109840050384021, 0.149944310198338, 0.164132187616120, 0.149944310198338, 0.109840050384021, 0.050848627308305};

        // 7-point Legendre-Gauss Quadrature on [0,1] for the integral due to Jianliang Qian
        double nodes_rhop[7] = {0.974553956171379,0.870765592799697,0.702922575688699,0.500000000000000,0.297077424311301,0.129234407200303,0.025446043828621};
        double weights_rhop[7] = {0.064742483084435,0.139852695744638,0.190915025252560,0.208979591836735,0.190915025252560,0.139852695744638,0.064742483084435};

        *val = 0;
        double tt = _h*_w*s0;
        for (int i=0;i<7;i++){
          for (int j=0;j<7;j++){
            double rhop=nodes_rhop[i];
            double phi=nodes_phi[j];
            double rho=rhop/cos(phi);
            double r1=tt/2*sqrt(1+pow(rho,2));
            _Complex double fun=((cos(r1)+Im*sin(r1))*(1-Im*r1)-1)*(rho/pow(1+pow(rho,2),1.5))/cos(phi);
            *val += weights_rhop[i]*weights_phi[j]*fun;
          }
        }
        *val*=-48.0*Im;
        *val*=Im/4.0/pi/pow(tt,2)/_h;

      }else{
        std::cout<<"should not arrive here for _vs==1 "<< std::endl;
        // Sample_noself(x1, y1, x2, y2,val);
      }
    }else{
      std::cout<<"should not arrive here for _vs==0 "<< std::endl;
    }
  }
};



// Assemble a block of matrix entries from interpolated D1, D2, tau
void assemble_fromD1D2Tau(double x1,double x2,double y1,double y2,double z1,double z2, _Complex double* output, C_QuantApp_BF* Q){

    int self = sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))<1e-20? 1:0;
    if(self==1){
      Q->SampleSelf(x1, y1, z1, x2, y2, z2, output);
    }else{
      double s0 = slowness(x1,y1,z1, Q->_slow_x0, Q->_slow_y0,Q->_slow_z0,Q->_ivelo,Q->_slowness_array->data(),Q->_h,Q->_slow_nx, Q->_slow_ny, Q->_slow_nz);
      double D1 =s0/2.0/pi; //fr[nr*nc + idxr+idxc*nr];
      double D2 =0;// fr[nr*nc*2 + idxr+idxc*nr];

      // double s1 = slowness(x2[idxc],y2[idxc], Q->_slow_x0, Q->_slow_y0,Q->_ivelo);
      // cout<<s0<<" "<<s1<<" "<<Q->_w<<endl;
      double tau = sqrt(pow(s0,2)* (pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)));
      // if(tau_square>0){
      // if(tau>2*pi/(Q->_w)/50){
      // if(tau>2*pi/(Q->_w)/20){
      //   tau = sqrt(max(tau_square,1e-30));
      // }
      // tau= acosh(((pow(x1[idxr]-x2[idxc],2) + pow(y1[idxr]-y2[idxc],2)))*pow(0.25,2)*s0*s1*0.5+1)/0.25;
      *output =Babich(Q->_d, Q->_w, D1, D2, tau);
    }
}


// Assemble entries for the scatterer-to-scatterer VIE operator using
// v = Q u as the unknown, where
//   Q(x) = k0^2 * ((s(x)/s0)^2 - 1), s0 = 2, k0 = omega*s0.
//
// Original u-unknown form:
//   (h^{-3} I - G Q) u = b.
// New v=Q*u form:
//   (h^{-3} Q^{-1} - G) v = b.
//
// Therefore the contrast coefficient belongs only to the diagonal/mass term.
// The off-diagonal Green part is always -G.  The scaleGreen option now only
// controls whether the diagonal h^{-3}/Q is inserted here (scaleGreen==0) or
// applied dynamically in C_FuncHMatVec_MD (scaleGreen==1).
void assemble_fromD1D2Tau_s2s(double x1,double x2,double y1,double y2, double z1,double z2, _Complex double* output, C_QuantApp_BF* Q){

    double s0=2;
    int self = sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))<1e-20? 1:0;
    if(self==1){
      Q->SampleSelf(x1, y1, z1, x2, y2, z2, output);
      if(Q->_scaleGreen==0){
        double s1 = slowness(x2,y2,z2, Q->_slow_x0, Q->_slow_y0,Q->_slow_z0,Q->_ivelo,Q->_slowness_array->data(),Q->_h, Q->_slow_nx, Q->_slow_ny, Q->_slow_nz);
        double k0 = s0*Q->_w;
        double coef = pow(k0,2.0)*(pow(s1/s0,2.0)-1); // Q(x2), assumed nonzero
        if(fabs(coef) < 1e-30){
          cout<<"zero contrast in diagonal VIE self term at ("<<x2<<", "<<y2<<", "<<z2<<") with slowness "<<s1<<endl;
          exit(1);
        }
        *output = -*output + 1.0/(pow(Q->_h,3.0)*coef);
      }else{
        *output = -*output ;
      }
    }else{
      double D1 =s0/2.0/pi; //fr[nr*nc + idxr+idxc*nr];
      double D2 =0;// fr[nr*nc*2 + idxr+idxc*nr];
      double tau = sqrt(pow(s0,2)* (pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)));
      *output =-Babich(Q->_d, Q->_w, D1, D2, tau);
    }
}





// Assemble a block of matrix entries from interpolated D1, D2, tau
void assemble_fromD1D2Tau_block(int nth, int nr, int nc, double* x1,double* x2,double* y1,double* y2,double* z1,double* z2, _Complex double* alldat_loc, int64_t* idx_val_map, double *fr, C_QuantApp_BF* Q){
    #pragma omp parallel for
    for (int idxrc=0;idxrc<nr*nc;idxrc++){
      int idxr = idxrc%nr;
      int idxc = idxrc/nr;

      _Complex double valtmp;
      assemble_fromD1D2Tau(x1[idxr],x2[idxc],y1[idxr],y2[idxc],z1[idxr],z2[idxc],&valtmp, Q);
      alldat_loc[idx_val_map[nth]+idxr+idxc*nr]=valtmp;
    }
}


// Assemble a block of matrix entries from interpolated D1, D2, tau
void assemble_fromD1D2Tau_block_s2s(int nth, int nr, int nc, double* x1,double* x2,double* y1,double* y2,double* z1,double* z2, _Complex double* alldat_loc, int64_t* idx_val_map, double *fr, C_QuantApp_BF* Q){
    #pragma omp parallel for
    for (int idxrc=0;idxrc<nr*nc;idxrc++){
      int idxr = idxrc%nr;
      int idxc = idxrc/nr;

      _Complex double valtmp;
      assemble_fromD1D2Tau_s2s(x1[idxr],x2[idxc],y1[idxr],y2[idxc],z1[idxr],z2[idxc],&valtmp, Q);
      alldat_loc[idx_val_map[nth]+idxr+idxc*nr]=valtmp;
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

  std::cout<<"C_FuncZmn_BF is no more used, use C_FuncZmnBlock_BF instead"<<std::endl;

}


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn_BF_V2V_MD(int* Ndim, int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  double x1 = Q->_data[(m[0]-1) * Q->_d];
  double y1 = Q->_data[(m[1]-1) * Q->_d+1];
  double z1 = Q->_data[(m[2]-1) * Q->_d+2];
  double x2 = Q->_data[(n[0]-1) * Q->_d];
  double y2 = Q->_data[(n[1]-1) * Q->_d+1];
  double z2 = Q->_data[(n[2]-1) * Q->_d+2];
  assemble_fromD1D2Tau(x1,x2,y1,y2,z1,z2,val, Q);
}



// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn_BF_S2S_MD(int* Ndim, int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  double x1 = Q->_data[(m[0]-1) * Q->_d];
  double y1 = Q->_data[(m[1]-1) * Q->_d+1];
  double z1 = Q->_data[(m[2]-1) * Q->_d+2];
  double x2 = Q->_data[(n[0]-1) * Q->_d];
  double y2 = Q->_data[(n[1]-1) * Q->_d+1];
  double z2 = Q->_data[(n[2]-1) * Q->_d+2];
  assemble_fromD1D2Tau_s2s(x1,x2,y1,y2,z1,z2, val, Q);
}

// The compressibility function (tensor) wrapper required by the Fortran HODLR code
inline void C_FuncNearFar_BF_MD(int* Ndim, int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}

// The distance function wrapper required by the Fortran HODLR code
inline void C_FuncDistmn_BF(int *m, int *n, double *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}

// The compressibility function wrapper required by the Fortran HODLR code
inline void C_FuncNearFar_BF(int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

}


// The entry evaluation routine for a list of subtensors for the VIE operator
inline void C_FuncZmnBlock_BF_S2S_MD(int* Ndim, int* Ninter, int* Nrow_max, int* Ncol_max, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

    int myrank, size;                     // Store values of processor rank and total no of procs requested
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    int64_t idx_row_1=0;
    int64_t idx_row_2=0;
    int64_t idx_row_3=0;
    int64_t idx_col_1=0;
    int64_t idx_col_2=0;
    int64_t idx_col_3=0;
    int64_t idx_val=0;
    int nrmax=0;
    int ncmax=0;
    int nvalmax=0;
    vector<int64_t> idx_row_map_1(*Ninter,0);
    vector<int64_t> idx_row_map_2(*Ninter,0);
    vector<int64_t> idx_row_map_3(*Ninter,0);
    vector<int64_t> idx_col_map_1(*Ninter,0);
    vector<int64_t> idx_col_map_2(*Ninter,0);
    vector<int64_t> idx_col_map_3(*Ninter,0);
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
      int nr_3 = rowidx[*Ndim*nn+2];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];
      int nc_3 = colidx[*Ndim*nn+2];

      if(nprow*npcol==1){
        idx_row_map_1[nn]=idx_row_1;
        idx_row_map_2[nn]=idx_row_2;
        idx_row_map_3[nn]=idx_row_3;
        idx_col_map_1[nn]=idx_col_1;
        idx_col_map_2[nn]=idx_col_2;
        idx_col_map_3[nn]=idx_col_3;
        idx_val_map[nn]=idx_val;
        idx_val+=nr_1*nr_2*nr_3*nc_1*nc_2*nc_3;
        idx_row_1+=nr_1;
        idx_row_2+=nr_2;
        idx_row_3+=nr_3;
        idx_col_1+=nc_1;
        idx_col_2+=nc_2;
        idx_col_3+=nc_3;
        nrmax = max(nr_1*nr_2*nr_3,nrmax);
        ncmax = max(nc_1*nc_2*nc_3,ncmax);
        nvalmax = max(nr_1*nr_2*nr_3*nc_1*nc_2*nc_3,nvalmax);
      }else{
        std::cout<<"nprow*npcol>1 in C_FuncZmnBlock_BF_S2S_MD"<<std::endl;
        exit(0);
      }
    }

    double *x1,*y1,*z1;
    double *x2,*y2,*z2;
    double *fr;
    // cout<<" "<<myrank<<" "<<NinterNew<<endl;
    #pragma omp parallel private(x1,y1,x2,y2,z1,z2,fr)
    {
    #pragma omp for
    for (int nn=0;nn<*Ninter;nn++){
      x1= (double*)malloc(nrmax*sizeof(double));
      y1= (double*)malloc(nrmax*sizeof(double));
      z1= (double*)malloc(nrmax*sizeof(double));
      x2= (double*)malloc(ncmax*sizeof(double));
      y2= (double*)malloc(ncmax*sizeof(double));
      z2= (double*)malloc(ncmax*sizeof(double));
      fr= (double*)malloc(nvalmax*3*sizeof(double));
      int nr_1 = rowidx[*Ndim*nn];
      int nr_2 = rowidx[*Ndim*nn+1];
      int nr_3 = rowidx[*Ndim*nn+2];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];
      int nc_3 = colidx[*Ndim*nn+2];

      for (int idxr_1=0;idxr_1<nr_1;idxr_1++){
      for (int idxr_2=0;idxr_2<nr_2;idxr_2++){
      for (int idxr_3=0;idxr_3<nr_3;idxr_3++){
        int idxr = idxr_1 + nr_1*idxr_2 + + nr_1*nr_2*idxr_3;
        int m_1=allrows[idx_row_map_1[nn]+idxr_1];
        int m_2=allrows[*Nrow_max + idx_row_map_2[nn]+idxr_2];
        int m_3=allrows[*Nrow_max*2 + idx_row_map_3[nn]+idxr_3];
        x1[idxr] = Q->_data[(m_1-1) * Q->_d];
        y1[idxr] = Q->_data[(m_2-1) * Q->_d+1];
        z1[idxr] = Q->_data[(m_3-1) * Q->_d+2];
      }
      }
      }

      for (int idxc_1=0;idxc_1<nc_1;idxc_1++){
      for (int idxc_2=0;idxc_2<nc_2;idxc_2++){
      for (int idxc_3=0;idxc_3<nc_3;idxc_3++){
        int idxc = idxc_1 + nc_1*idxc_2 + nc_1*nc_2*idxc_3;
        int n_1=allcols[idx_col_map_1[nn]+idxc_1];
        int n_2=allcols[*Ncol_max + idx_col_map_2[nn]+idxc_2];
        int n_3=allcols[*Ncol_max*2 + idx_col_map_3[nn]+idxc_3];
        x2[idxc] = Q->_data[(n_1-1) * Q->_d];
        y2[idxc] = Q->_data[(n_2-1) * Q->_d+1];
        z2[idxc] = Q->_data[(n_3-1) * Q->_d+2];
      }
      }
      }

      // lagrange_interp2D_vec_block(Q->_x_cheb.data(), Q->_y_cheb.data(), nr, nc, x1,y1,x2,y2,Q->_u1_square_int_cheb.data(),Q->_D1_int_cheb.data(),Q->_D2_int_cheb.data(), Q->_TNx,Q->_TNy, fr);

      assemble_fromD1D2Tau_block_s2s(nn,nr_1*nr_2*nr_3, nc_1*nc_2*nc_3, x1,x2,y1,y2,z1,z2, alldat_loc, idx_val_map.data(),fr, Q);
      free(x1);
      free(y1);
      free(z1);
      free(x2);
      free(y2);
      free(z2);
      free(fr);
    }
    }
time(&end);
timer += difftime(end,start);

}



// The entry evaluation routine for a list of subtensors for the background operator
inline void C_FuncZmnBlock_BF_V2V_MD(int* Ndim, int* Ninter, int* Nrow_max, int* Ncol_max, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

    int myrank, size;                     // Store values of processor rank and total no of procs requested
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    int64_t idx_row_1=0;
    int64_t idx_row_2=0;
    int64_t idx_row_3=0;
    int64_t idx_col_1=0;
    int64_t idx_col_2=0;
    int64_t idx_col_3=0;
    int64_t idx_val=0;
    int nrmax=0;
    int ncmax=0;
    int nvalmax=0;
    vector<int64_t> idx_row_map_1(*Ninter,0);
    vector<int64_t> idx_row_map_2(*Ninter,0);
    vector<int64_t> idx_row_map_3(*Ninter,0);
    vector<int64_t> idx_col_map_1(*Ninter,0);
    vector<int64_t> idx_col_map_2(*Ninter,0);
    vector<int64_t> idx_col_map_3(*Ninter,0);
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
      int nr_3 = rowidx[*Ndim*nn+2];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];
      int nc_3 = colidx[*Ndim*nn+2];

      if(nprow*npcol==1){
        idx_row_map_1[nn]=idx_row_1;
        idx_row_map_2[nn]=idx_row_2;
        idx_row_map_3[nn]=idx_row_3;
        idx_col_map_1[nn]=idx_col_1;
        idx_col_map_2[nn]=idx_col_2;
        idx_col_map_3[nn]=idx_col_3;
        idx_val_map[nn]=idx_val;
        idx_val+=nr_1*nr_2*nr_3*nc_1*nc_2*nc_3;
        idx_row_1+=nr_1;
        idx_row_2+=nr_2;
        idx_row_3+=nr_3;
        idx_col_1+=nc_1;
        idx_col_2+=nc_2;
        idx_col_3+=nc_3;
        nrmax = max(nr_1*nr_2*nr_3,nrmax);
        ncmax = max(nc_1*nc_2*nc_3,ncmax);
        nvalmax = max(nr_1*nr_2*nr_3*nc_1*nc_2*nc_3,nvalmax);
      }else{
        std::cout<<"nprow*npcol>1 in C_FuncZmnBlock_BF_V2V_MD"<<std::endl;
        exit(0);
      }
    }

    double *x1,*y1,*z1;
    double *x2,*y2,*z2;
    double *fr;
    // cout<<" "<<myrank<<" "<<NinterNew<<endl;
    #pragma omp parallel private(x1,y1,x2,y2,z1,z2,fr)
    {
    #pragma omp for
    for (int nn=0;nn<*Ninter;nn++){
      x1= (double*)malloc(nrmax*sizeof(double));
      y1= (double*)malloc(nrmax*sizeof(double));
      z1= (double*)malloc(nrmax*sizeof(double));
      x2= (double*)malloc(ncmax*sizeof(double));
      y2= (double*)malloc(ncmax*sizeof(double));
      z2= (double*)malloc(ncmax*sizeof(double));
      fr= (double*)malloc(nvalmax*3*sizeof(double));
      int nr_1 = rowidx[*Ndim*nn];
      int nr_2 = rowidx[*Ndim*nn+1];
      int nr_3 = rowidx[*Ndim*nn+2];
      int nc_1 = colidx[*Ndim*nn];
      int nc_2 = colidx[*Ndim*nn+1];
      int nc_3 = colidx[*Ndim*nn+2];

      for (int idxr_1=0;idxr_1<nr_1;idxr_1++){
      for (int idxr_2=0;idxr_2<nr_2;idxr_2++){
      for (int idxr_3=0;idxr_3<nr_3;idxr_3++){
        int idxr = idxr_1 + nr_1*idxr_2 + + nr_1*nr_2*idxr_3;
        int m_1=allrows[idx_row_map_1[nn]+idxr_1];
        int m_2=allrows[*Nrow_max + idx_row_map_2[nn]+idxr_2];
        int m_3=allrows[*Nrow_max*2 + idx_row_map_3[nn]+idxr_3];
        x1[idxr] = Q->_data[(m_1-1) * Q->_d];
        y1[idxr] = Q->_data[(m_2-1) * Q->_d+1];
        z1[idxr] = Q->_data[(m_3-1) * Q->_d+2];
      }
      }
      }

      for (int idxc_1=0;idxc_1<nc_1;idxc_1++){
      for (int idxc_2=0;idxc_2<nc_2;idxc_2++){
      for (int idxc_3=0;idxc_3<nc_3;idxc_3++){
        int idxc = idxc_1 + nc_1*idxc_2 + nc_1*nc_2*idxc_3;
        int n_1=allcols[idx_col_map_1[nn]+idxc_1];
        int n_2=allcols[*Ncol_max + idx_col_map_2[nn]+idxc_2];
        int n_3=allcols[*Ncol_max*2 + idx_col_map_3[nn]+idxc_3];
        x2[idxc] = Q->_data[(n_1-1) * Q->_d];
        y2[idxc] = Q->_data[(n_2-1) * Q->_d+1];
        z2[idxc] = Q->_data[(n_3-1) * Q->_d+2];
      }
      }
      }

      // lagrange_interp2D_vec_block(Q->_x_cheb.data(), Q->_y_cheb.data(), nr, nc, x1,y1,x2,y2,Q->_u1_square_int_cheb.data(),Q->_D1_int_cheb.data(),Q->_D2_int_cheb.data(), Q->_TNx,Q->_TNy, fr);

      assemble_fromD1D2Tau_block(nn,nr_1*nr_2*nr_3, nc_1*nc_2*nc_3, x1,x2,y1,y2,z1,z2, alldat_loc, idx_val_map.data(),fr, Q);
      free(x1);
      free(y1);
      free(z1);
      free(x2);
      free(y2);
      free(z2);
      free(fr);
    }
    }
time(&end);
timer += difftime(end,start);

}



// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn_BF_S2S(int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  double x1 = Q->_data[(*m-1) * Q->_d];
  double y1 = Q->_data[(*m-1) * Q->_d+1];
  double z1 = Q->_data[(*m-1) * Q->_d+2];
  double x2 = Q->_data[(*n-1) * Q->_d];
  double y2 = Q->_data[(*n-1) * Q->_d+1];
  double z2 = Q->_data[(*n-1) * Q->_d+2];
  assemble_fromD1D2Tau_s2s(x1,x2,y1,y2,z1,z2, val, Q);
}


// The extraction sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmnBlock_BF_S2S(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;
}




inline void C_FuncApplyDiagPrecond(int* Ndim_p, char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;
  int Ndim = *Ndim_p;
  int Npo = product(nout,Ndim);

  if(nout[0]!=nin[0] || nout[1]!=nin[1] || nout[2]!=nin[2]){
    cout<<"nin should equal nout in C_FuncApplyDiagPrecond"<<endl;
    exit(1);
  }

  for(int i0=0;i0<Npo;i0++){
    for (int nth=0; nth<*nvec; nth++){
      xout[i0+nth*Npo] = 0.0;
    }
  }

  if(Q->_bdiag_precon==1){
    int myrank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if(static_cast<int>(Q->_local_diag_root.size()) < Npo ||
       static_cast<int>(Q->_local_diag_gid.size()) < Npo ||
       static_cast<int>(Q->_local_diag_pos.size()) < Npo){
      cout<<"diagonal preconditioner routing has not been initialized"<<endl;
      exit(1);
    }

    vector<vector<int> > send_meta_by_rank(nproc);
    vector<vector<_Complex double> > send_vals_by_rank(nproc);
    for(int i0=0;i0<Npo;i0++){
      int dest = Q->_local_diag_root[i0];
      if(dest >= 0){
        if(dest >= nproc){
          cout<<"invalid diagonal preconditioner destination rank"<<endl;
          exit(1);
        }
        for (int nth=0; nth<*nvec; nth++){
          send_meta_by_rank[dest].push_back(Q->_local_diag_gid[i0]);
          send_meta_by_rank[dest].push_back(Q->_local_diag_pos[i0]);
          send_meta_by_rank[dest].push_back(nth);
          send_meta_by_rank[dest].push_back(i0);
          send_vals_by_rank[dest].push_back(xin[i0+nth*Npo]);
        }
      }
    }

    vector<int> send_counts(nproc,0), recv_counts(nproc,0), send_displs(nproc,0), recv_displs(nproc,0);
    for(int rr=0; rr<nproc; rr++) send_counts[rr] = static_cast<int>(send_vals_by_rank[rr].size());
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int total_send = 0;
    int total_recv = 0;
    for(int rr=0; rr<nproc; rr++){
      send_displs[rr] = total_send;
      recv_displs[rr] = total_recv;
      total_send += send_counts[rr];
      total_recv += recv_counts[rr];
    }

    vector<int> send_meta_counts(nproc,0), recv_meta_counts(nproc,0), send_meta_displs(nproc,0), recv_meta_displs(nproc,0);
    for(int rr=0; rr<nproc; rr++){
      send_meta_counts[rr] = 4*send_counts[rr];
      recv_meta_counts[rr] = 4*recv_counts[rr];
      send_meta_displs[rr] = 4*send_displs[rr];
      recv_meta_displs[rr] = 4*recv_displs[rr];
    }

    vector<int> send_meta;
    vector<_Complex double> send_vals;
    send_meta.reserve(4*total_send);
    send_vals.reserve(total_send);
    for(int rr=0; rr<nproc; rr++){
      send_meta.insert(send_meta.end(), send_meta_by_rank[rr].begin(), send_meta_by_rank[rr].end());
      send_vals.insert(send_vals.end(), send_vals_by_rank[rr].begin(), send_vals_by_rank[rr].end());
    }
    vector<vector<int> >().swap(send_meta_by_rank);
    vector<vector<_Complex double> >().swap(send_vals_by_rank);

    vector<int> recv_meta(4*total_recv);
    vector<_Complex double> recv_vals(total_recv);
    MPI_Alltoallv(send_meta.data(), send_meta_counts.data(), send_meta_displs.data(), MPI_INT,
                  recv_meta.data(), recv_meta_counts.data(), recv_meta_displs.data(), MPI_INT,
                  MPI_COMM_WORLD);
    MPI_Alltoallv(send_vals.data(), send_counts.data(), send_displs.data(), MPI_C_DOUBLE_COMPLEX,
                  recv_vals.data(), recv_counts.data(), recv_displs.data(), MPI_C_DOUBLE_COMPLEX,
                  MPI_COMM_WORLD);

    vector<int> recv_src(total_recv,0);
    for(int rr=0; rr<nproc; rr++){
      for(int pp=recv_displs[rr]; pp<recv_displs[rr]+recv_counts[rr]; pp++){
        recv_src[pp] = rr;
      }
    }
    unordered_map<int, vector<int> > items_by_gid;
    items_by_gid.reserve(total_recv);
    for(int ii=0; ii<total_recv; ii++){
      items_by_gid[recv_meta[4*ii]].push_back(ii);
    }

    vector<vector<int> > ret_meta_by_rank(nproc);
    vector<vector<_Complex double> > ret_vals_by_rank(nproc);

    int nblock = Q->_nblock_diag;
    int info_stride = Q->_diag_block_stride;
    for (int bb=0; bb<nblock; bb++){
      int* info = Q->_diag_block_info.data() + bb*info_stride;
      int head_g = info[0];
      int tail_g = info[1];
      int root_global = info[4];
      int global_block_id = info[7];
      int Npo_from_block = tail_g - head_g;

      int Npo_diag = Q->_diag_npo[bb];
      int Nloc_diag = Q->_diag_nloc[bb];
      if(Npo_diag != Npo_from_block){
        cout<<"inconsistent diagonal block size in C_FuncApplyDiagPrecond"<<endl;
        exit(1);
      }
      int cnt_block = (*nvec)*Npo_diag;
      int cnt_diag = (*nvec)*Nloc_diag;
      vector<_Complex double> block_in(cnt_block,{0.0,0.0});
      vector<_Complex double> block_out(cnt_block,{0.0,0.0});
      vector<_Complex double> xin_diag(cnt_diag,{0.0,0.0});
      vector<_Complex double> xout_diag(cnt_diag,{0.0,0.0});

      if(myrank == root_global){
        auto it = items_by_gid.find(global_block_id);
        if(it != items_by_gid.end()){
          for(size_t qq=0; qq<it->second.size(); qq++){
            int item = it->second[qq];
            int pos = recv_meta[4*item+1];
            int nth = recv_meta[4*item+2];
            if(pos >= 0 && pos < Npo_diag && nth >= 0 && nth < *nvec){
              block_in[pos+nth*Npo_diag] = recv_vals[item];
            }
          }
        }
      }
      MPI_Bcast(block_in.data(), cnt_block, MPI_C_DOUBLE_COMPLEX, 0, Q->_diag_comms[bb]);

      for(int i0=0;i0<Nloc_diag;i0++){
        int i_new_diag_scalar = i0+1;
        int ii;
        z_c_bpack_new2old(Q->quant_ptr_diags[bb]->msh_bf,&i_new_diag_scalar,&ii);
        for (int nth=0; nth<*nvec; nth++){
          xin_diag[i0+nth*Nloc_diag] = block_in[(ii-1)+nth*Npo_diag];
        }
      }

      z_c_bpack_inv_mult("N",xin_diag.data(),xout_diag.data(),&Nloc_diag,&Nloc_diag,nvec,&(Q->bmat_diags[bb]),&(Q->option_diag),&(Q->stats_diags[bb]),&(Q->ptree_diags[bb]));

      for(int i0=0;i0<Nloc_diag;i0++){
        int i_new_diag_scalar = i0+1;
        int ii;
        z_c_bpack_new2old(Q->quant_ptr_diags[bb]->msh_bf,&i_new_diag_scalar,&ii);
        for (int nth=0; nth<*nvec; nth++){
          block_out[(ii-1)+nth*Npo_diag] = xout_diag[i0+nth*Nloc_diag];
        }
      }

      vector<_Complex double> block_out_sum;
      _Complex double* reduce_recv = block_out.data();
      if(myrank == root_global){
        block_out_sum.assign(cnt_block,{0.0,0.0});
        reduce_recv = block_out_sum.data();
      }
      MPI_Reduce(block_out.data(), reduce_recv, cnt_block, MPI_C_DOUBLE_COMPLEX, MPI_SUM, 0, Q->_diag_comms[bb]);

      if(myrank == root_global){
        auto it = items_by_gid.find(global_block_id);
        if(it != items_by_gid.end()){
          for(size_t qq=0; qq<it->second.size(); qq++){
            int item = it->second[qq];
            int src = recv_src[item];
            int pos = recv_meta[4*item+1];
            int nth = recv_meta[4*item+2];
            int local_i0 = recv_meta[4*item+3];
            if(pos >= 0 && pos < Npo_diag && nth >= 0 && nth < *nvec){
              ret_meta_by_rank[src].push_back(local_i0);
              ret_meta_by_rank[src].push_back(nth);
              ret_vals_by_rank[src].push_back(block_out_sum[pos+nth*Npo_diag]);
            }
          }
        }
      }
    }

    vector<int> ret_send_counts(nproc,0), ret_recv_counts(nproc,0), ret_send_displs(nproc,0), ret_recv_displs(nproc,0);
    for(int rr=0; rr<nproc; rr++) ret_send_counts[rr] = static_cast<int>(ret_vals_by_rank[rr].size());
    MPI_Alltoall(ret_send_counts.data(), 1, MPI_INT, ret_recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int ret_total_send = 0;
    int ret_total_recv = 0;
    for(int rr=0; rr<nproc; rr++){
      ret_send_displs[rr] = ret_total_send;
      ret_recv_displs[rr] = ret_total_recv;
      ret_total_send += ret_send_counts[rr];
      ret_total_recv += ret_recv_counts[rr];
    }

    vector<int> ret_send_meta_counts(nproc,0), ret_recv_meta_counts(nproc,0), ret_send_meta_displs(nproc,0), ret_recv_meta_displs(nproc,0);
    for(int rr=0; rr<nproc; rr++){
      ret_send_meta_counts[rr] = 2*ret_send_counts[rr];
      ret_recv_meta_counts[rr] = 2*ret_recv_counts[rr];
      ret_send_meta_displs[rr] = 2*ret_send_displs[rr];
      ret_recv_meta_displs[rr] = 2*ret_recv_displs[rr];
    }

    vector<int> ret_send_meta;
    vector<_Complex double> ret_send_vals;
    ret_send_meta.reserve(2*ret_total_send);
    ret_send_vals.reserve(ret_total_send);
    for(int rr=0; rr<nproc; rr++){
      ret_send_meta.insert(ret_send_meta.end(), ret_meta_by_rank[rr].begin(), ret_meta_by_rank[rr].end());
      ret_send_vals.insert(ret_send_vals.end(), ret_vals_by_rank[rr].begin(), ret_vals_by_rank[rr].end());
    }

    vector<int> ret_recv_meta(2*ret_total_recv);
    vector<_Complex double> ret_recv_vals(ret_total_recv);
    MPI_Alltoallv(ret_send_meta.data(), ret_send_meta_counts.data(), ret_send_meta_displs.data(), MPI_INT,
                  ret_recv_meta.data(), ret_recv_meta_counts.data(), ret_recv_meta_displs.data(), MPI_INT,
                  MPI_COMM_WORLD);
    MPI_Alltoallv(ret_send_vals.data(), ret_send_counts.data(), ret_send_displs.data(), MPI_C_DOUBLE_COMPLEX,
                  ret_recv_vals.data(), ret_recv_counts.data(), ret_recv_displs.data(), MPI_C_DOUBLE_COMPLEX,
                  MPI_COMM_WORLD);

    for(int item=0; item<ret_total_recv; item++){
      int local_i0 = ret_recv_meta[2*item];
      int nth = ret_recv_meta[2*item+1];
      if(local_i0 >= 0 && local_i0 < Npo && nth >= 0 && nth < *nvec){
        xout[local_i0+nth*Npo] = ret_recv_vals[item];
      }
    }
  }else{
    for(int i0=0;i0<Npo;i0++){
      bool is_scatterer = Q->LocalScatterIndex(i0) > 0;
      for (int nth=0; nth<*nvec; nth++){
        xout[i0+nth*Npo] = is_scatterer ? xin[i0+nth*Npo] : 0.0;
      }
    }
  }
}





// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncHMatVec_MD(int* Ndim, char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant) {
  C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;

  int nin_loc = product(nout,*Ndim);
  int64_t cnt = (*nvec)*nin_loc;
  _Complex double* xbuf1 = new _Complex double[cnt];
  _Complex double* xin1 = new _Complex double[cnt];

  if(nout[0]!=nin[0] || nout[1]!=nin[1] || nout[2]!=nin[2]){
    cout<<"nin should equal nout"<<endl;
    exit(1);
  }

  for (int i=0; i<product(nout,*Ndim); i++){

    int i_new_loc_scalar = i+1;
    int i_new_loc_md[*Ndim];
    int i_old_scalar;
    int i_old_md[*Ndim];

    z_c_bpack_singleindex_to_multiindex(Ndim,nout,&i_new_loc_scalar,i_new_loc_md);
    z_c_bpack_md_new2old(Ndim,Q->msh_bf,i_new_loc_md,i_old_md);
    double x1 = Q->_data[(i_old_md[0]-1) * Q->_d];
    double y1 = Q->_data[(i_old_md[1]-1) * Q->_d+1];
    double z1 = Q->_data[(i_old_md[2]-1) * Q->_d+2];

    // double x2 = (i_old_md[0]-1+Q->_idx_off_x) * Q->_h;
    // double y2 = (i_old_md[1]-1+Q->_idx_off_y) * Q->_h;
    // double z2 = (i_old_md[2]-1+Q->_idx_off_z) * Q->_h;

    // if(fabs(x1-x2)>Q->_h/2 || fabs(y1-y2)>Q->_h/2 || fabs(z1-z2)>Q->_h/2){
    //   cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<" "<<z1<<" "<<z2<<endl;
    // }


    bool is_scatterer = Q->LocalScatterIndex(i) > 0;
    double coef = 1.0;
    if(is_scatterer){
      double s1 = slowness(x1,y1,z1,Q->_slow_x0, Q->_slow_y0,Q->_slow_z0,Q->_ivelo,Q->_slowness_array->data(),Q->_h, Q->_slow_nx, Q->_slow_ny, Q->_slow_nz);
      double s0=2;
      double k0 = s0*Q->_w;
      coef = pow(k0,2.0)*(pow(s1/s0,2.0)-1);
    }

    for (int nth=0; nth<*nvec; nth++){
      // v = Q*u is the unknown.  The mass/identity contribution is
      // h^{-3} * Q^{-1} * v.  Keep the sign convention from the sampled
      // matrix: z_c_bpack_md_mult applies the stored -G operator.
      if(is_scatterer){
        xbuf1[i+nth*product(nout,*Ndim)]=xin[i+nth*product(nout,*Ndim)]/(pow(Q->_h,3.0)*coef);
      }else{
        xbuf1[i+nth*product(nout,*Ndim)]=0.0;
      }
    }
    for (int nth=0; nth<*nvec; nth++){
      // The Green operator acts directly on v, not on Q*u.
      xin1[i+nth*product(nout,*Ndim)] = is_scatterer ? xin[i+nth*product(nout,*Ndim)] : 0.0;
    }
  }
  z_c_bpack_md_mult(Ndim,trans,xin1,xout,nin,nout,nvec,Q->bmat_scaled_green,Q->option_bf,Q->stats_bf,Q->ptree_bf,Q->msh_bf);

  for (int i=0; i<product(nout,*Ndim); i++){
    bool is_scatterer = Q->LocalScatterIndex(i) > 0;
    for (int nth=0; nth<*nvec; nth++){
      if(is_scatterer){
        xout[i+nth*product(nout,*Ndim)]=xout[i+nth*product(nout,*Ndim)] + xbuf1[i+nth*product(nout,*Ndim)];
      }else{
        xout[i+nth*product(nout,*Ndim)]=0.0;
      }
    }
  }
  delete[] xbuf1;
  delete[] xin1;
}





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
  double opt_d;

if(myrank==master_rank){
	z_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}

	/*****************************************************************/
  int Nx, Ny, Nz; // number of grid points per x,y,z direction
  int Nmax; // max between Nx, Ny and Nz

  // define the inner domain
  double x0min=-0.0, x0max=1.0/1.0;
  double y0min=-0.0, y0max=1.0/1.0;
  double z0min=-0.0, z0max=1.0/1.0;


  double h=0.02;
  double h0=0.01;
  int ivelo = 1;
  int ivelo_o = 1; // this is assumed fixed
  int rmax = 20;
  double w = 5*pi;
  int TNx =10, TNy=10, TNz=10;

  double slow_x0;
  double slow_y0;
  double slow_z0;

  double L = 0.4;
  double H = 0.4;
  double W = 0.4;

  int scaleGreen=1;
  double smin_ivelo11=1.0;
  double smax_ivelo11=3.0;
  int nshape=200;
  int bdiag_precon = 0 ;
  int bdiag_level = -1 ;
  int pad_s2s = 0;

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
      // {"xmin",          required_argument, 0, 7},
      // {"xmax",          required_argument, 0, 8},
      // {"ymin",          required_argument, 0, 9},
      // {"ymax",          required_argument, 0, 10},
      {"x0min",          required_argument, 0, 11},
      {"x0max",          required_argument, 0, 12},
      {"y0min",          required_argument, 0, 13},
      {"y0max",          required_argument, 0, 14},
      {"h",          required_argument, 0, 15},
      {"ivelo",          required_argument, 0, 16},
      {"omega",          required_argument, 0, 17},
      {"h0",          required_argument, 0, 18},
      {"TNy",          required_argument, 0, 19},
      {"TNz",          required_argument, 0, 20},
      // {"zmin",          required_argument, 0, 21},
      // {"zmax",          required_argument, 0, 22},
      {"z0min",          required_argument, 0, 23},
      {"z0max",          required_argument, 0, 24},
      {"readtable",        required_argument, 0, 25},
      {"shape",        required_argument, 0, 26},
      {"L",        required_argument, 0, 27},
      {"H",        required_argument, 0, 28},
      {"W",        required_argument, 0, 29},
      {"scaleGreen",        required_argument, 0, 30},
      {"nshape",        required_argument, 0, 31},
      {"smin_ivelo11",        required_argument, 0, 32},
      {"smax_ivelo11",        required_argument, 0, 33},
      {"bdiag_precon",        required_argument, 0, 34},
      {"bdiag_level",        required_argument, 0, 35},
      {"pad_s2s",        required_argument, 0, 36},
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
      // std::istringstream iss(optarg);
      // iss >> xmin;
    } break;
    case 8: {
      // std::istringstream iss(optarg);
      // iss >> xmax;
    } break;
    case 9: {
      // std::istringstream iss(optarg);
      // iss >> ymin;
    } break;
    case 10: {
      // std::istringstream iss(optarg);
      // iss >> ymax;
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
      iss >> TNz;
    } break;
    case 21: {
      // std::istringstream iss(optarg);
      // iss >> zmin;
    } break;
    case 22: {
      // std::istringstream iss(optarg);
      // iss >> zmax;
    } break;
    case 23: {
      std::istringstream iss(optarg);
      iss >> z0min;
    } break;
    case 24: {
      std::istringstream iss(optarg);
      iss >> z0max;
    } break;
    case 25: {
      std::istringstream iss(optarg);
      iss >> readtable;
    } break;
    case 26: {
      std::istringstream iss(optarg);
      iss >> shape;
    } break;
    case 27: {
      std::istringstream iss(optarg);
      iss >> L;
    } break;
    case 28: {
      std::istringstream iss(optarg);
      iss >> H;
    } break;
    case 29: {
      std::istringstream iss(optarg);
      iss >> W;
    } break;
    case 30: {
      std::istringstream iss(optarg);
      iss >> scaleGreen;
    } break;
    case 31: {
      std::istringstream iss(optarg);
      iss >> nshape;
    } break;
    case 32: {
      std::istringstream iss(optarg);
      iss >> smin_ivelo11;
    } break;
    case 33: {
      std::istringstream iss(optarg);
      iss >> smax_ivelo11;
    } break;
    case 34: {
      std::istringstream iss(optarg);
      iss >> bdiag_precon;
    } break;
    case 35: {
      std::istringstream iss(optarg);
      iss >> bdiag_level;
    } break;
    case 36: {
      std::istringstream iss(optarg);
      iss >> pad_s2s;
    } break;
    default: break;
    }
  }


	Ndim=3; //data dimension
  int Ns[Ndim];
  Nx = round((x0max-x0min)/h);
  Ny = round((y0max-y0min)/h);
  Nz = round((z0max-z0min)/h);
  Ns[0]=Nx;
  Ns[1]=Ny;
  Ns[2]=Nz;
  Nmax = max(max(Nx,Ny),Nz);

  double center[3]; //geometrical center of the scatterer
  L = x0max - x0min;
  H = y0max - y0min;
  W = z0max - z0min;
  scaleGreen = 1;
  center[0]=(x0min+x0max)/2.0;
  center[1]=(y0min+y0max)/2.0;
  center[2]=(z0min+z0max)/2.0;

  slow_x0 = center[0];
  slow_y0 = center[1];
  slow_z0 = center[2];


  // const int64_t I=round((xmax-xmin)/h),  J=round((ymax-ymin)/h), K=round((zmax-zmin)/h);
  const int64_t Iint=round((x0max-x0min)/h),  Jint=round((y0max-y0min)/h),Kint=round((z0max-z0min)/h);

  int idx_off_x = 0;
  int idx_off_y = 0;
  int idx_off_z = 0;

  int Ns_s[Ndim];
  int Nx_s, Ny_s, Nz_s, Nmax_s;
  Nx_s = Nx;
  Ny_s = Ny;
  Nz_s = Nz;
  Ns_s[0]=Nx_s;
  Ns_s[1]=Ny_s;
  Ns_s[2]=Nz_s;
  Nmax_s = max(max(Nx_s,Ny_s),Nz_s);

  if(shape!=4){
    cout<<"The tensor VIE cannot handle irregular shape at this moment"<<endl;
    exit(1);
  }

    vector<double> slowness_array(Iint*Jint*Kint,2.0);;
    if(ivelo==11){
      string filename_in, str1, str2;
      std::ostringstream streamObj1,streamObj2;
      // cout<<smin_ivelo11 <<" "<<smax_ivelo11<<endl;
      streamObj1 << smin_ivelo11;
      str1=streamObj1.str();
      str1.erase ( str1.find_last_not_of('0') + 1, std::string::npos ); str1.erase ( str1.find_last_not_of('.') + 1, std::string::npos );
      streamObj2 << smax_ivelo11;
      str2=streamObj2.str();
      str2.erase ( str2.find_last_not_of('0') + 1, std::string::npos ); str2.erase ( str2.find_last_not_of('.') + 1, std::string::npos );

      filename_in ="slowness_map_shapeoutter"+to_string(Iint)+"x"+to_string(Jint)+"x"+to_string(Kint)+"_shape"+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+"_off"+to_string(idx_off_x)+"x"+to_string(idx_off_y)+"x"+to_string(idx_off_z)+"_range"+str1+"_"+str2+"_nshape"+to_string(nshape)+".bin";
      // cout<<filename_in<<endl;

      // Open the binary file
      std::ifstream file(filename_in, std::ios::binary);
      if (!file) {
          std::cout << "Unable to open " << filename_in << std::endl;
          exit(1);
      }
      file.read(reinterpret_cast<char*>(slowness_array.data()), Iint*Jint*Kint * sizeof(double));
      file.close();
    }

    // generate the chebyshev nodes
    double txmin=0, txmax=pi;
    double tymin=0, tymax=pi;
    double tzmin=0, tzmax=pi;
    double thx=(txmax-txmin)/(TNx-1);
    double thy=(tymax-tymin)/(TNy-1);
    double thz=(tzmax-tzmin)/(TNz-1);
    double *tx=(double *)malloc(TNx*sizeof(double));
    double *ty=(double *)malloc(TNy*sizeof(double));
    double *tz=(double *)malloc(TNz*sizeof(double));

    vector<double> x_cheb(TNx);
    vector<double> y_cheb(TNy);
    vector<double> z_cheb(TNz);
    for(int i=0;i<TNx;i++){
      tx[i]=txmin+i*thx;
      x_cheb[i]=cos(tx[i]);
    }
    for(int j=0;j<TNy;j++){
      ty[j]=tymin+j*thy;
      y_cheb[j]=cos(ty[j]);
    }
    for(int k=0;k<TNz;k++){
      tz[k]=tzmin+k*thz;
      z_cheb[k]=cos(tz[k]);
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
    for(int k=0;k<TNz;k++){
      z_cheb[k]=(z0min+z0max)/2.0+(z0max-z0min)/2.0*z_cheb[k];
    }
    vector<double> u1_square_int_cheb(TNx*TNy*TNz*TNx*TNy*TNz,0);
    vector<double> D1_int_cheb(TNx*TNy*TNz*TNx*TNy*TNz,0);
    vector<double> D2_int_cheb(TNx*TNy*TNz*TNx*TNy*TNz,0);


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
	int* nns_ptr = nullptr;

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

	//quantities for ButterflyPACK
	F2Cptr option_bf;     //option structure returned by Fortran code
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
  z_c_bpack_set_option_from_command_line(argc, argv, option_bf);

  z_c_bpack_getoption(&option_bf, "Nmin_leaf", &opt_d);
  int tensor_leafsize = max(1, static_cast<int>(round(opt_d)));
  if(pad_s2s){
    Nx_s = pad_to_leaf_power2(Nx, tensor_leafsize);
    Ny_s = pad_to_leaf_power2(Ny, tensor_leafsize);
    Nz_s = pad_to_leaf_power2(Nz, tensor_leafsize);
    tensor_leafsize = effective_leaf_size_for_padded_dim(Nx_s,tensor_leafsize);
    z_c_bpack_set_I_option(&option_bf, "Nmin_leaf", tensor_leafsize);
  }else{
    Nx_s = Nx;
    Ny_s = Ny;
    Nz_s = Nz;
  }
  Ns_s[0]=Nx_s;
  Ns_s[1]=Ny_s;
  Ns_s[2]=Nz_s;
  Nmax_s = max(max(Nx_s,Ny_s),Nz_s);

  double trans_opt = 0.0;
  z_c_bpack_getoption(&option_bf, "trans_invariant", &trans_opt);
  if(myrank==0 && pad_s2s){
    cout<<"Tensor grid padding: physical "<<Nx<<" x "<<Ny<<" x "<<Nz
        <<" padded to "<<Nx_s<<" x "<<Ny_s<<" x "<<Nz_s
        <<" with Nmin_leaf "<<tensor_leafsize
        <<" effective leaves "<<effective_leaf_size_for_padded_dim(Nx_s,tensor_leafsize)
        <<" x "<<effective_leaf_size_for_padded_dim(Ny_s,tensor_leafsize)
        <<" x "<<effective_leaf_size_for_padded_dim(Nz_s,tensor_leafsize)
        <<" and levels "<<tree_levels_for_padded_dim(Nx_s,tensor_leafsize)
        <<" x "<<tree_levels_for_padded_dim(Ny_s,tensor_leafsize)
        <<" x "<<tree_levels_for_padded_dim(Nz_s,tensor_leafsize)<<endl;
  }
  if(myrank==0 && trans_opt >= 0.5 && !pad_s2s){
    cout<<"Warning: trans_invariant is enabled but pad_s2s=0, so the tensor tree may not be uniform/complete."<<endl;
  }



  double smax=2.0;
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      for(int k=0;k<Nz;k++){
        smax = max(smax,slowness(i*h-L/2+center[0],j*h-H/2+center[1],k*h-W/2+center[2],slow_x0, slow_y0,slow_z0,ivelo,slowness_array.data(),h, Iint, Jint, Kint));
      }
    }
  }
  double vtmp=smax;
  MPI_Allreduce(&vtmp,&smax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);



	/*****************************************************************/
	/* Test the free-space 3D Green function */
	vector<double> data_geo;
	vector<double> panelnodes;
  double dl=0.01;
  int knn_pre;
  int layer=3;
  knn_pre= (2*layer+1)*(2*layer+1)*(2*layer+1);



  if(vs==1){
    data_geo.resize(Ndim*Nmax);
    for(int ii=0;ii<Nx;ii++){
      data_geo[(ii) * Ndim] = ii*h+x0min;
    }
    for(int ii=0;ii<Ny;ii++){
      data_geo[(ii) * Ndim+1] = ii*h+y0min;
    }
    for(int ii=0;ii<Nz;ii++){
      data_geo[(ii) * Ndim+2] = ii*h+z0min;
    }
    if(myrank==0){
      cout<<"smax: "<<smax<<" PPW: "<<2*pi/(w*smax)/h<<" From: "<< Nx <<" x "<< Ny <<" x "<< Nz <<" To: "<< Nx <<" x "<< Ny <<" x "<< Nz <<endl;
    }
  }else{
    cout<<"vs=0 is not implemented in 3D code"<<endl;
    exit(0);
  }
  if(vs==1){
    z_c_bpack_createptree(&size, groups, &Fcomm, &ptree_bf);
    z_c_bpack_set_I_option(&option_bf, "cpp", cpp);


    data_geo.resize(Ndim*Nmax_s);
    for(int ii=0;ii<Nx_s;ii++){
      data_geo[(ii) * Ndim] = ii*h-L/2+center[0];
    }
    for(int ii=0;ii<Ny_s;ii++){
      data_geo[(ii) * Ndim+1] = ii*h-H/2+center[1];
    }
    for(int ii=0;ii<Nz_s;ii++){
      data_geo[(ii) * Ndim+2] = ii*h-W/2+center[2];
    }

    int Ntot_s2s = product(Ns_s,Ndim);
    double background_s0 = 2.0;
    double scatterer_tol = 1e-12;
    vector<int> scatter_to_tensor;
    vector<double> data_geo_scatterer;
    for(int kk=0; kk<Nz_s; kk++){
      for(int jj=0; jj<Ny_s; jj++){
        for(int ii=0; ii<Nx_s; ii++){
          double xs = ii*h-L/2+center[0];
          double ys = jj*h-H/2+center[1];
          double zs = kk*h-W/2+center[2];
          int tensor_scalar = kk*Nx_s*Ny_s + jj*Nx_s + ii + 1;
          bool in_physical = ii < Nx && jj < Ny && kk < Nz;
          double s1 = in_physical ? slowness(xs,ys,zs,slow_x0,slow_y0,slow_z0,ivelo,slowness_array.data(),h,Iint,Jint,Kint) : background_s0;
          if(in_physical && fabs(s1-background_s0) > scatterer_tol*max(1.0,fabs(background_s0))){
            scatter_to_tensor.push_back(tensor_scalar);
            data_geo_scatterer.push_back(xs);
            data_geo_scatterer.push_back(ys);
            data_geo_scatterer.push_back(zs);
          }
        }
      }
    }
    int Nscatter = static_cast<int>(scatter_to_tensor.size());

      if(myrank==0){
        cout<<"smax: "<<smax<<" PPW: "<<2*pi/(w*smax)/h<<" From: "<< Nx_s <<" x "<< Ny_s <<" x "<< Nz_s <<" To: "<< Nx_s <<" x "<< Ny_s <<" x "<< Nz_s<<endl;
        cout<<"Scatterer unknowns: "<<Nscatter<<" out of "<<Ntot_s2s<<endl;
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
    quant_ptr_bf_s2s=new C_QuantApp_BF(data_geo, Ndim, Nx_s, Ny_s, Nz_s, static_cast<int>(Iint), static_cast<int>(Jint), static_cast<int>(Kint), 1, w, x0min, x0max, y0min, y0max, z0min, z0max, h, dl, ivelo,&slowness_array,rmax,verbose,vs, bdiag_precon, x_cheb,y_cheb,z_cheb,u1_square_int_cheb,D1_int_cheb,D2_int_cheb);
    quant_ptr_bf_s2s->_bdiag_level = bdiag_level;
    quant_ptr_bf_s2s->_nscatter = Nscatter;
    quant_ptr_bf_s2s->_background_s0 = background_s0;
    quant_ptr_bf_s2s->_scatterer_tol = scatterer_tol;
    quant_ptr_bf_s2s->_scatter_to_tensor = scatter_to_tensor;


    // construct hodlr with geometrical points
    z_c_bpack_md_construct_init(Ns_s,&Nmax_s, &Ndim, data_geo.data(), perms_bf_s2s, myseg_s2s, &bmat_bf_s2s, &option_bf, &stats_bf_s2s, &msh_bf_s2s, &kerquant_bf_s2s, &ptree_bf, &C_FuncNearFar_BF_MD, quant_ptr_bf_s2s);
    quant_ptr_bf_s2s->_Hperm.resize(Nmax_s*Ndim);
    std::copy(perms_bf_s2s, perms_bf_s2s + Nmax_s*Ndim, quant_ptr_bf_s2s->_Hperm.begin());
    quant_ptr_bf_s2s->_local_tensor_ids.resize(product(myseg_s2s,Ndim));
    quant_ptr_bf_s2s->_local_scatter_ids.resize(product(myseg_s2s,Ndim),0);
    for (int i=0; i<product(myseg_s2s,Ndim); i++){
      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_md[Ndim];
      int i_old_scalar;
      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg_s2s,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf_s2s,i_new_loc_md,i_old_md);
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns_s,&i_old_scalar,i_old_md);
      quant_ptr_bf_s2s->_local_tensor_ids[i] = i_old_scalar;
      quant_ptr_bf_s2s->_local_scatter_ids[i] = quant_ptr_bf_s2s->ScatterIndexFromTensorScalar(i_old_scalar);
    }
	  z_c_bpack_printoption(&option_bf,&ptree_bf);
    z_c_bpack_md_construct_element_compute(&Ndim,&bmat_bf_s2s, &option_bf, &stats_bf_s2s, &msh_bf_s2s, &kerquant_bf_s2s, &ptree_bf, &C_FuncZmn_BF_S2S_MD, &C_FuncZmnBlock_BF_S2S_MD, quant_ptr_bf_s2s);

    if(myrank==master_rank)std::cout<<"\n\nGenerating the incident fields: "<<std::endl;
    int nvec=1;
    vector<_Complex double> b(product(myseg_s2s,Ndim)*nvec,{0.0,0.0});
    vector<_Complex double> x(product(myseg_s2s,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg_s2s,Ndim); i++){
      int i_new_loc_scalar = i+1;
      int i_new_loc_md[Ndim];
      int i_old_scalar;
      int i_old_md[Ndim];

      z_c_bpack_singleindex_to_multiindex(&Ndim,myseg_s2s,&i_new_loc_scalar,i_new_loc_md);
      z_c_bpack_md_new2old(&Ndim,&msh_bf_s2s,i_new_loc_md,i_old_md);
      z_c_bpack_multiindex_to_singleindex(&Ndim,Ns_s,&i_old_scalar,i_old_md);

      double xs = data_geo[(i_old_md[0]-1) * Ndim];
      double ys = data_geo[(i_old_md[1]-1) * Ndim+1];
      double zs = data_geo[(i_old_md[2]-1) * Ndim+2];
      bool in_physical = in_physical_tensor_domain(i_old_md,Nx,Ny,Nz);
      // double xs0=x0max-0.1;
      // double ys0=y0max-0.1;
      // double zs0=z0max-0.1;
      double xs0=0.1;
      double ys0=0.1;
      double zs0=0.1;

      for (int nth=0; nth<nvec; nth++){
        x.data()[i+nth*product(myseg_s2s,Ndim)]=in_physical ? source_function(xs,ys,zs, xs0, ys0, zs0, nth,h,w) : 0.0;
      }
    }
    z_c_bpack_md_mult(&Ndim,"N",x.data(),b.data(),myseg_s2s,myseg_s2s,&nvec,&bmat_bf_s2s,&option_bf,&stats_bf_s2s,&ptree_bf, &msh_bf_s2s);
    int out_ix = plot_slice_index(Nx);
    int out_iy = plot_slice_index(Ny);
    int out_iz = plot_slice_index(Nz);
    size_t one_rhs_slice_size = tensor_slice_size(Nx,Ny,Nz);
    int slice_reduce_count = static_cast<int>(one_rhs_slice_size*nvec);
    vector<_Complex double> slices_loc, slices_glo;
    collect_tensor_slices_from_local(Ndim,Ns,myseg_s2s,&msh_bf_s2s,b,nvec,-1.0,out_ix,out_iy,out_iz,slices_loc);
    if(myrank==master_rank) slices_glo.resize(one_rhs_slice_size*nvec);
    MPI_Reduce(slices_loc.data(), myrank==master_rank ? slices_glo.data() : nullptr,
               slice_reduce_count, MPI_C_DOUBLE_COMPLEX, MPI_SUM, master_rank, MPI_COMM_WORLD);

    if(myrank==master_rank){
      for(int nth=0; nth<nvec; nth++){
        string filename, str;
        filename = "./VIE_F_inc_f_";
        str=to_string(w/2/pi);str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename +=str+"_vs_"+to_string(vs)+"_ivelo_"+to_string(ivelo);
        std::ostringstream streamObj;
        streamObj << h;
        str=streamObj.str();
        str.erase ( str.find_last_not_of('0') + 1, std::string::npos ); str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
        filename += "_h_"+str;
        double opt_d;
        z_c_bpack_getoption(&option_bf, "tol_comp", &opt_d);

        str=my::to_string(opt_d);
        filename += "_tol_"+str+"_nth_"+to_string(nth)+"_tensor_slices.bin";
        write_tensor_slices_file(filename,Nx,Ny,Nz,h,out_ix,out_iy,out_iz,
                                 slices_glo.data()+nth*one_rhs_slice_size);
      }
    }

#if 0
    if(myrank==master_rank)std::cout<<"\n\nFactoring the scatterer-scatterer operator: "<<std::endl;
    // factor hodlr

    z_c_bpack_getoption(&option_bf, "precon", &opt_d);
    int precon=round(opt_d);
    if(precon!=2)z_c_bpack_factor(&bmat_bf_s2s,&option_bf,&stats_bf_s2s,&ptree_bf,&msh_bf_s2s);

#else // construct a block-diagonal preconditioner
  MPI_Group world_group;
  if(bdiag_precon==1){
    if(Nscatter <= 0){
      cout<<"No scatterer unknowns found for the block-diagonal preconditioner"<<endl;
      exit(1);
    }
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    z_c_bpack_copyoption(&option_bf,&(quant_ptr_bf_s2s->option_diag));

    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "LRlevel", 0); // lightweight HODLR tree for block discovery
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "format", 1);
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "xyzsort", 1); // KD-tree ordering
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "knn", 0); // avoid global KNN construction over all scatterers
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "precon", 1); // Direct solver
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "elem_extract", 0); // Need to add elem_extract=2 for the matrix solver
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "RecLR_leaf", 5); // sometimes BACA is not accurate enough
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "Nmin_leaf", 64); // leaf size in the matrix solver
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "verbosity", -1); // suppress printing while building the scatterer tree

    F2Cptr bmat_scatterer_blocks;
    F2Cptr stats_scatterer_blocks;
    F2Cptr msh_scatterer_blocks;
    F2Cptr kerquant_scatterer_blocks;
    z_c_bpack_createstats(&stats_scatterer_blocks);
    int myseg_scatterer = 0;
    int* perms_scatterer = new int[Nscatter];
    int nlevel_scatterer = 0;
    int* tree_scatterer = new int[1];
    tree_scatterer[0] = Nscatter;
    z_c_bpack_construct_init(&Nscatter, &Ndim, data_geo_scatterer.data(), nns_ptr, &nlevel_scatterer, tree_scatterer, perms_scatterer, &myseg_scatterer, &bmat_scatterer_blocks, &(quant_ptr_bf_s2s->option_diag), &stats_scatterer_blocks, &msh_scatterer_blocks, &kerquant_scatterer_blocks, &ptree_bf, &C_FuncDistmn_BF, &C_FuncNearFar_BF, quant_ptr_bf_s2s);
    vector<double>().swap(data_geo_scatterer);

    quant_ptr_bf_s2s->_diag_block_stride = 8;
    int mesh_range[2];
    int nblock = 0;
    int dummy_block_info[1];
    z_c_bpack_get_anylevel_blocks(&(quant_ptr_bf_s2s->_bdiag_level), &nblock, dummy_block_info, mesh_range, &bmat_scatterer_blocks, &(quant_ptr_bf_s2s->option_diag), &ptree_bf, &msh_scatterer_blocks);
    quant_ptr_bf_s2s->_nblock_diag = nblock;
    quant_ptr_bf_s2s->_diag_block_info.resize(nblock*quant_ptr_bf_s2s->_diag_block_stride);
    int nblock_capacity = nblock;
    z_c_bpack_get_anylevel_blocks(&(quant_ptr_bf_s2s->_bdiag_level), &nblock_capacity, quant_ptr_bf_s2s->_diag_block_info.data(), mesh_range, &bmat_scatterer_blocks, &(quant_ptr_bf_s2s->option_diag), &ptree_bf, &msh_scatterer_blocks);
    if(nblock_capacity != nblock){
      cout<<"nblock changed between C_BPACK_Get_Anylevel_Blocks calls"<<endl;
      exit(1);
    }

    // The following controls the parameters of each diagonal block. Currently they are hardcoded, but should be made command-line options.   
    // The global matrix tree above is used only to discover diagonal blocks.
    // Enable the high-frequency HODBF options only for the smaller blocks.
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "LRlevel", 2);
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "knn", 20);
    z_c_bpack_set_D_option(&(quant_ptr_bf_s2s->option_diag), "sample_para", 2.01);
    z_c_bpack_set_D_option(&(quant_ptr_bf_s2s->option_diag), "sample_para_outer", 2.01);
    tol_rand = 1e-2;
    z_c_bpack_set_D_option(&(quant_ptr_bf_s2s->option_diag), "tol_rand", tol_rand);
    z_c_bpack_set_D_option(&(quant_ptr_bf_s2s->option_diag), "tol_Rdetect", tol_rand*3e-1);



    quant_ptr_bf_s2s->_outer_idxs[0] = mesh_range[0];
    quant_ptr_bf_s2s->_outer_idxe_plus1[0] = mesh_range[1];

    int local_info_count = nblock*quant_ptr_bf_s2s->_diag_block_stride;
    vector<int> info_counts(size,0), info_displs(size,0);
    MPI_Allgather(&local_info_count, 1, MPI_INT, info_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int total_info_count = 0;
    for(int rr=0; rr<size; rr++){
      info_displs[rr] = total_info_count;
      total_info_count += info_counts[rr];
    }
    vector<int> all_diag_info(total_info_count);
    int* local_info_ptr = local_info_count > 0 ? quant_ptr_bf_s2s->_diag_block_info.data() : nullptr;
    int* all_info_ptr = total_info_count > 0 ? all_diag_info.data() : nullptr;
    MPI_Allgatherv(local_info_ptr, local_info_count, MPI_INT,
                   all_info_ptr, info_counts.data(), info_displs.data(), MPI_INT, MPI_COMM_WORLD);

    vector<DiagRouteBlock> diag_route_blocks;
    diag_route_blocks.reserve(total_info_count/max(1,quant_ptr_bf_s2s->_diag_block_stride));
    unordered_map<int,int> seen_diag_gid;
    for(int off=0; off+quant_ptr_bf_s2s->_diag_block_stride<=total_info_count; off+=quant_ptr_bf_s2s->_diag_block_stride){
      int gid = all_diag_info[off+7];
      if(seen_diag_gid.find(gid) != seen_diag_gid.end()) continue;
      seen_diag_gid[gid] = 1;
      DiagRouteBlock route;
      route.head = all_diag_info[off+0];
      route.tail = all_diag_info[off+1];
      route.root = all_diag_info[off+4];
      route.gid = gid;
      diag_route_blocks.push_back(route);
    }
    std::sort(diag_route_blocks.begin(), diag_route_blocks.end(),
              [](const DiagRouteBlock& a, const DiagRouteBlock& b){return a.head < b.head;});

    int nloc_s2s = product(myseg_s2s,Ndim);
    quant_ptr_bf_s2s->_local_diag_root.assign(nloc_s2s,-1);
    quant_ptr_bf_s2s->_local_diag_gid.assign(nloc_s2s,-1);
    quant_ptr_bf_s2s->_local_diag_pos.assign(nloc_s2s,-1);
    for(int i=0; i<nloc_s2s; i++){
      int scatter_old = quant_ptr_bf_s2s->_local_scatter_ids[i];
      if(scatter_old <= 0) continue;
      int scatter_new = 0;
      z_c_bpack_old2new(&msh_scatterer_blocks, &scatter_old, &scatter_new);

      int lo = 0;
      int hi = static_cast<int>(diag_route_blocks.size());
      while(lo < hi){
        int mid = lo + (hi-lo)/2;
        if(diag_route_blocks[mid].head <= scatter_new){
          lo = mid + 1;
        }else{
          hi = mid;
        }
      }
      int route_idx = lo - 1;
      if(route_idx < 0 || scatter_new >= diag_route_blocks[route_idx].tail){
        cout<<"Rank "<<myrank<<" failed to route scatterer "<<scatter_old
            <<" with new index "<<scatter_new<<" to a diagonal block"<<endl;
        exit(1);
      }
      quant_ptr_bf_s2s->_local_diag_root[i] = diag_route_blocks[route_idx].root;
      quant_ptr_bf_s2s->_local_diag_gid[i] = diag_route_blocks[route_idx].gid;
      quant_ptr_bf_s2s->_local_diag_pos[i] = scatter_new - diag_route_blocks[route_idx].head;
    }

    quant_ptr_bf_s2s->_diag_block_offsets.assign(nblock+1,0);
    for (int bb=0; bb<nblock; bb++){
      int* block_info = quant_ptr_bf_s2s->_diag_block_info.data() + bb*quant_ptr_bf_s2s->_diag_block_stride;
      quant_ptr_bf_s2s->_diag_block_offsets[bb+1] = quant_ptr_bf_s2s->_diag_block_offsets[bb] + block_info[1] - block_info[0];
    }
    quant_ptr_bf_s2s->_diag_tensor_ids.resize(quant_ptr_bf_s2s->_diag_block_offsets[nblock]);

    quant_ptr_bf_s2s->bmat_diags = new F2Cptr[nblock];
    quant_ptr_bf_s2s->kerquant_diags = new F2Cptr[nblock];
    quant_ptr_bf_s2s->stats_diags = new F2Cptr[nblock];
    quant_ptr_bf_s2s->msh_diags = new F2Cptr[nblock];
    quant_ptr_bf_s2s->ptree_diags = new F2Cptr[nblock];
    quant_ptr_bf_s2s->quant_ptr_diags = new C_QuantApp_BF*[nblock];
    quant_ptr_bf_s2s->_diag_nloc.resize(nblock);
    quant_ptr_bf_s2s->_diag_npo.resize(nblock);
    quant_ptr_bf_s2s->_diag_comms.resize(nblock, MPI_COMM_NULL);
    quant_ptr_bf_s2s->_diag_groups.resize(nblock, MPI_GROUP_NULL);


    for (int bb=0; bb<nblock; bb++){
      vector<double> data_geo_diag;
      z_c_bpack_createstats(&(quant_ptr_bf_s2s->stats_diags[bb]));

      int* block_info = quant_ptr_bf_s2s->_diag_block_info.data() + bb*quant_ptr_bf_s2s->_diag_block_stride;
      int head_g = block_info[0];
      int tail_g = block_info[1];
      int Npo = tail_g - head_g;
      int block_offset = quant_ptr_bf_s2s->_diag_block_offsets[bb];
      quant_ptr_bf_s2s->_diag_npo[bb] = Npo;
      data_geo_diag.resize(Ndim*Npo);
      int pgrp_head = block_info[4];
      int pgrp_nproc = block_info[5];
      int global_block_id = block_info[7];

      vector<int> ranks_diag(pgrp_nproc);
      vector<int> groups_diag(pgrp_nproc);
      for(int rr=0; rr<pgrp_nproc; rr++){
        ranks_diag[rr] = pgrp_head + rr;
        groups_diag[rr] = rr;
      }
      MPI_Group_incl(world_group, pgrp_nproc, ranks_diag.data(), &(quant_ptr_bf_s2s->_diag_groups[bb]));
      MPI_Comm_create_group(MPI_COMM_WORLD, quant_ptr_bf_s2s->_diag_groups[bb], global_block_id, &(quant_ptr_bf_s2s->_diag_comms[bb]));
      if(quant_ptr_bf_s2s->_diag_comms[bb] == MPI_COMM_NULL){
        cout<<"Rank "<<myrank<<" failed to create diagonal-block communicator"<<endl;
        exit(1);
      }
      MPI_Fint Fcomm_diag = MPI_Comm_c2f(quant_ptr_bf_s2s->_diag_comms[bb]);
      z_c_bpack_createptree(&pgrp_nproc, groups_diag.data(), &Fcomm_diag, &(quant_ptr_bf_s2s->ptree_diags[bb]));

      if(myrank == pgrp_head){
        std::cout<<"Rank "<< myrank <<": Constructing HODBF for diagonal block: "<< bb+1 << " out of " << nblock <<". Npo: "<< Npo << " MPI ranks: " << pgrp_nproc <<std::endl;
      }

      for(int ii=0;ii<Npo;ii++){
        int scatter_old = perms_scatterer[head_g + ii - 1];
        int tensor_scalar = scatter_to_tensor[scatter_old - 1];
        quant_ptr_bf_s2s->_diag_tensor_ids[block_offset + ii] = tensor_scalar;

        int tensor_scalar_tmp = tensor_scalar;
        int tensor_md[Ndim];
        z_c_bpack_singleindex_to_multiindex(&Ndim,Ns_s,&tensor_scalar_tmp,tensor_md);
        double xs = data_geo[(tensor_md[0]-1) * Ndim];
        double ys = data_geo[(tensor_md[1]-1) * Ndim+1];
        double zs = data_geo[(tensor_md[2]-1) * Ndim+2];

        data_geo_diag[(ii) * Ndim] = xs;
        data_geo_diag[(ii) * Ndim+1] = ys;
        data_geo_diag[(ii) * Ndim+2] = zs;
      }



    // create hodlr data structures
    quant_ptr_bf_s2s->quant_ptr_diags[bb]=new C_QuantApp_BF(data_geo_diag, Ndim, 0, w, x0min, x0max, y0min, y0max, z0min, z0max, h, dl, ivelo,&slowness_array,rmax,verbose,vs, x_cheb,y_cheb,z_cheb,u1_square_int_cheb,D1_int_cheb,D2_int_cheb);


  	int myseg_diag=0;     // local number of unknowns
	  int* perms_diag = new int[Npo]; //permutation vector returned by HODLR
    int nlevel_diag = 0; // 0: tree level, nonzero if a tree is provided
    int* tree_diag = new int[(int)pow(2,nlevel_diag)]; //user provided array containing size of each leaf node, not used if nlevel=0
    tree_diag[0] = Npo;

    // construct hodlr with geometrical points
    z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "verbosity", -1); // suppress prining
    z_c_bpack_construct_init(&Npo, &Ndim, data_geo_diag.data(), nns_ptr, &nlevel_diag, tree_diag, perms_diag, &myseg_diag, &(quant_ptr_bf_s2s->bmat_diags[bb]), &(quant_ptr_bf_s2s->option_diag), &(quant_ptr_bf_s2s->stats_diags[bb]), &(quant_ptr_bf_s2s->msh_diags[bb]), &(quant_ptr_bf_s2s->kerquant_diags[bb]), &(quant_ptr_bf_s2s->ptree_diags[bb]), &C_FuncDistmn_BF, &C_FuncNearFar_BF, quant_ptr_bf_s2s->quant_ptr_diags[bb]);
    quant_ptr_bf_s2s->_diag_nloc[bb] = myseg_diag;

    quant_ptr_bf_s2s->quant_ptr_diags[bb]->_Hperm.resize(Npo);
    std::copy(perms_diag, perms_diag + Npo, quant_ptr_bf_s2s->quant_ptr_diags[bb]->_Hperm.begin());

    if(global_block_id == 1){
      z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "verbosity", 0);
      z_c_bpack_printoption(&(quant_ptr_bf_s2s->option_diag),&(quant_ptr_bf_s2s->ptree_diags[bb]));
      // z_c_bpack_set_I_option(&(quant_ptr_bf_s2s->option_diag), "verbosity", -1);
    }

  	z_c_bpack_construct_element_compute(&(quant_ptr_bf_s2s->bmat_diags[bb]), &(quant_ptr_bf_s2s->option_diag), &(quant_ptr_bf_s2s->stats_diags[bb]), &(quant_ptr_bf_s2s->msh_diags[bb]), &(quant_ptr_bf_s2s->kerquant_diags[bb]), &(quant_ptr_bf_s2s->ptree_diags[bb]), &C_FuncZmn_BF_S2S, &C_FuncZmnBlock_BF_S2S, quant_ptr_bf_s2s->quant_ptr_diags[bb]);
    z_c_bpack_factor(&(quant_ptr_bf_s2s->bmat_diags[bb]),&(quant_ptr_bf_s2s->option_diag),&(quant_ptr_bf_s2s->stats_diags[bb]),&(quant_ptr_bf_s2s->ptree_diags[bb]),&(quant_ptr_bf_s2s->msh_diags[bb]));

    quant_ptr_bf_s2s->quant_ptr_diags[bb]->msh_bf = &(quant_ptr_bf_s2s->msh_diags[bb]);

    delete[] tree_diag;
    delete[] perms_diag;
    }
    z_c_bpack_deletekernelquant(&kerquant_scatterer_blocks);
    z_c_bpack_deletemesh(&msh_scatterer_blocks);
    z_c_bpack_deletestats(&stats_scatterer_blocks);
    z_c_bpack_delete(&bmat_scatterer_blocks);
    delete[] tree_scatterer;
    delete[] perms_scatterer;
  }


#endif





    if(myrank==master_rank)std::cout<<"\n\nSolving the volume IE: "<<std::endl;
    int nin_loc = product(myseg_s2s,Ndim);
    vector<_Complex double> b_s(product(myseg_s2s,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg_s2s,Ndim); i++){
      bool is_scatterer = quant_ptr_bf_s2s->LocalScatterIndex(i) > 0;
      for (int nth=0; nth<nvec; nth++){
        b_s[i+nth*product(myseg_s2s,Ndim)] = is_scatterer ? -b[i+nth*product(myseg_s2s,Ndim)] : 0.0;
      }
    }
    int ErrSol=0;
    z_c_bpack_set_I_option(&option_bf, "ErrSol", ErrSol);
    vector<_Complex double> x_s(product(myseg_s2s,Ndim)*nvec,{0.0,0.0});

    quant_ptr_bf_s2s->bmat_scaled_green = &bmat_bf_s2s;
    quant_ptr_bf_s2s->option_bf = &option_bf;
    quant_ptr_bf_s2s->stats_bf = &stats_bf_s2s;
    quant_ptr_bf_s2s->ptree_bf = &ptree_bf;
    quant_ptr_bf_s2s->msh_bf = &msh_bf_s2s;
    F2Cptr kerquant_s2s;

    z_c_bpack_md_iter_usermatvec_precon(&Ndim, x_s.data(),b_s.data(),myseg_s2s,&nvec,&option_bf, &stats_bf_s2s, &ptree_bf, &kerquant_s2s, &C_FuncHMatVec_MD, &C_FuncApplyDiagPrecond, quant_ptr_bf_s2s);
    z_c_bpack_deletekernelquant(&kerquant_s2s);

    for (int i=0; i<product(myseg_s2s,Ndim); i++){
      if(quant_ptr_bf_s2s->LocalScatterIndex(i) <= 0){
        for (int nth=0; nth<nvec; nth++){
          x_s[i+nth*product(myseg_s2s,Ndim)] = 0.0;
        }
      }
    }

    vector<_Complex double> x_v(product(myseg_s2s,Ndim)*nvec,{0.0,0.0}),b_v(product(myseg_s2s,Ndim)*nvec,{0.0,0.0});
    for (int i=0; i<product(myseg_s2s,Ndim); i++){
      for (int nth=0; nth<nvec; nth++){
        x_v[i+nth*product(myseg_s2s,Ndim)] = x_s[i+nth*product(myseg_s2s,Ndim)];
      }
    }

    z_c_bpack_md_mult(&Ndim,"N",x_v.data(),b_v.data(),myseg_s2s,myseg_s2s,&nvec,&bmat_bf_s2s,&option_bf,&stats_bf_s2s,&ptree_bf,&msh_bf_s2s);
    collect_tensor_slices_from_local(Ndim,Ns,myseg_s2s,&msh_bf_s2s,b_v,nvec,-1.0,out_ix,out_iy,out_iz,slices_loc);
    if(myrank==master_rank) std::fill(slices_glo.begin(), slices_glo.end(), (_Complex double){0.0,0.0});
    MPI_Reduce(slices_loc.data(), myrank==master_rank ? slices_glo.data() : nullptr,
               slice_reduce_count, MPI_C_DOUBLE_COMPLEX, MPI_SUM, master_rank, MPI_COMM_WORLD);

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
        filename += "_tol_"+str+"_nth_"+to_string(nth)+"_tensor_slices.bin";
#ifdef IVELO9_CONST
        filename +="_ivelo9_const";
#endif
        write_tensor_slices_file(filename,Nx,Ny,Nz,h,out_ix,out_iy,out_iz,
                                 slices_glo.data()+nth*one_rhs_slice_size);
      }
    }

    if(myrank==master_rank)std::cout<<"\n\nPrinting stats of the scaled Green/VIE operator: "<<std::endl;
    z_c_bpack_printstats(&stats_bf_s2s,&ptree_bf);


    if(bdiag_precon){
      int nblock = quant_ptr_bf_s2s->_nblock_diag;
      for (int bb=0; bb<nblock; bb++){
        z_c_bpack_deletestats(&(quant_ptr_bf_s2s->stats_diags[bb]));
        z_c_bpack_deletemesh(&(quant_ptr_bf_s2s->msh_diags[bb]));
        z_c_bpack_deletekernelquant(&(quant_ptr_bf_s2s->kerquant_diags[bb]));
        z_c_bpack_delete(&(quant_ptr_bf_s2s->bmat_diags[bb]));
        z_c_bpack_deleteproctree(&(quant_ptr_bf_s2s->ptree_diags[bb]));
        if(quant_ptr_bf_s2s->_diag_comms[bb] != MPI_COMM_NULL) MPI_Comm_free(&(quant_ptr_bf_s2s->_diag_comms[bb]));
        if(quant_ptr_bf_s2s->_diag_groups[bb] != MPI_GROUP_NULL) MPI_Group_free(&(quant_ptr_bf_s2s->_diag_groups[bb]));
        delete quant_ptr_bf_s2s->quant_ptr_diags[bb];
      }
      z_c_bpack_deleteoption(&(quant_ptr_bf_s2s->option_diag));
      delete[] quant_ptr_bf_s2s->bmat_diags;
      delete[] quant_ptr_bf_s2s->kerquant_diags;
      delete[] quant_ptr_bf_s2s->stats_diags;
      delete[] quant_ptr_bf_s2s->msh_diags;
      delete[] quant_ptr_bf_s2s->ptree_diags;
      delete[] quant_ptr_bf_s2s->quant_ptr_diags;
      MPI_Group_free(&world_group);
    }

    delete quant_ptr_bf_s2s;
    z_c_bpack_deletestats(&stats_bf_s2s);
    z_c_bpack_md_deletemesh(&Ndim,&msh_bf_s2s);
    z_c_bpack_deletekernelquant(&kerquant_bf_s2s);
    z_c_bpack_delete(&bmat_bf_s2s);
    delete[] perms_bf_s2s;



  }else{

    cout<<"shouldn't reach here"<<endl;

  }


	z_c_bpack_deleteproctree(&ptree_bf);
	z_c_bpack_deleteoption(&option_bf);

	delete[] groups;
	free(tx);
	free(ty);
	free(tz);


	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
  return 0;
}
// // //------------------------------------------------------------------------------
