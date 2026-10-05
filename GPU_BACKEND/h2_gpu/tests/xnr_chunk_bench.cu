// Does splitting a box's X_NR into chunks cost compression ratio, and how much
// decompression throughput does the extra parallelism buy?
//
// Lossy Bitcomp decompresses a batch in parallel across its blocks, so a wave
// of 37 boxes is a 37-block batch however many bytes it holds.  Chunking each
// box's X_NR into pieces multiplies the block count at the same total bytes and
// the same scratch.  The question is the per-block overhead: each chunk carries
// its own header and its own quantization of the same data.
//
// delta = 2*tol gives max error <= tol (bitcomp rounds delta down to a power of
// two).  The tolerance is per box, from that box's max|.|, as the H2 path does.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>
#include <nvcomp/native/bitcomp.h>

#define CK(x) do{ cudaError_t e=(x); if(e){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__);exit(1);} }while(0)
#define BK(x) do{ bitcompResult_t r=(x); if(r!=BITCOMP_SUCCESS){printf("bitcomp err %d @%d\n",(int)r,__LINE__);exit(1);} }while(0)

static size_t align_up(size_t b, size_t a = 256) { return (b + a - 1) / a * a; }
static double ms_of(cudaEvent_t a, cudaEvent_t b){ float m; cudaEventElapsedTime(&m,a,b); return m; }

struct Box { std::vector<double> a; double maxabs = 0.0; };

static bool load(const char* p, Box& b){
  FILE* f=fopen(p,"rb"); if(!f) return false;
  fseek(f,0,SEEK_END); long sz=ftell(f); fseek(f,0,SEEK_SET);
  if(sz<=0){fclose(f);return false;}
  b.a.resize(sz/sizeof(double));
  size_t got=fread(b.a.data(),sizeof(double),b.a.size(),f); fclose(f);
  if(got!=b.a.size()) return false;
  b.maxabs=0; for(double v:b.a) b.maxabs=fmax(b.maxabs,fabs(v));
  return true;
}

int main(int argc,char**argv){
  const double tolfac = argc>1 ? atof(argv[1]) : 2e-4;  // tol_comp/5 with tol_comp 1e-3
  std::vector<Box> boxes;
  for(int i=0;i<48;i++){
    char p[256];
    snprintf(p,sizeof(p),"/global/cfs/cdirs/m2957/liuyangz/my_research/h2_c1_work/xnr_%02d.f64",i);
    Box b; if(load(p,b)) boxes.push_back(std::move(b));
  }
  if(boxes.empty()){printf("no xnr_*.f64 blocks found\n");return 2;}
  size_t tot_elems=0; for(auto&b:boxes) tot_elems+=b.a.size();
  printf("%zu real X_NR boxes, %.2f MB total, %.2f MB per box, tolfac %.0e\n",
         boxes.size(), tot_elems*8.0/1e6, tot_elems*8.0/1e6/boxes.size(), tolfac);

  cudaStream_t s; CK(cudaStreamCreate(&s));
  cudaEvent_t e0,e1; CK(cudaEventCreate(&e0)); CK(cudaEventCreate(&e1));

  // the boxes, back to back on the device (as the kept block holds them)
  std::vector<double> flat; flat.reserve(tot_elems);
  std::vector<size_t> box_at, box_n;
  for(auto&b:boxes){ box_at.push_back(flat.size()); box_n.push_back(b.a.size());
                     flat.insert(flat.end(), b.a.begin(), b.a.end()); }
  double* d_src; CK(cudaMalloc(&d_src, flat.size()*sizeof(double)));
  CK(cudaMemcpy(d_src, flat.data(), flat.size()*sizeof(double), cudaMemcpyHostToDevice));
  double* d_out; CK(cudaMalloc(&d_out, flat.size()*sizeof(double)));

  printf("\n%10s %8s %9s %8s %10s %10s %9s %12s\n",
         "chunk MB","blocks","stage MB","ratio","comp GB/s","dcmp GB/s","dcmp ms","max_err");

  const double chunks_mb[]={0.0,0.125,0.25,0.5,1.0,2.0,4.0,8.0};  // 0 = whole box, one block
  for(double cmb:chunks_mb){
    // ---- cut the boxes into chunks
    struct Piece { size_t off, elems; double delta; };
    std::vector<Piece> pieces;
    for(size_t i=0;i<boxes.size();++i){
      const double delta = 2.0*boxes[i].maxabs*tolfac;
      const size_t per = cmb<=0.0 ? box_n[i] : std::max<size_t>(1,(size_t)(cmb*1e6/8.0));
      for(size_t o=0;o<box_n[i];o+=per)
        pieces.push_back(Piece{box_at[i]+o, std::min(per, box_n[i]-o), delta});
    }
    const size_t nb=pieces.size();
    std::vector<size_t> nbytes(nb);
    size_t stage=0;
    for(size_t i=0;i<nb;++i){ nbytes[i]=pieces[i].elems*sizeof(double);
                              stage+=align_up(bitcompMaxBuflen(nbytes[i])); }
    char* d_stage; CK(cudaMalloc(&d_stage,stage));
    std::vector<const double*> h_in(nb); std::vector<void*> h_cmp(nb), h_dst(nb);
    std::vector<size_t> h_nb(nb); std::vector<double> h_delta(nb);
    size_t at=0;
    for(size_t i=0;i<nb;++i){
      h_in[i]=(const double*)(d_src+pieces[i].off); h_dst[i]=d_out+pieces[i].off;
      h_cmp[i]=d_stage+at; at+=align_up(bitcompMaxBuflen(nbytes[i]));
      h_nb[i]=nbytes[i]; h_delta[i]=pieces[i].delta;
    }
    const double** d_in; void **d_cmp,**d_dst; size_t *d_nbv,*d_csz; double* d_dlt;
    CK(cudaMalloc(&d_in,nb*sizeof(void*)));  CK(cudaMalloc(&d_cmp,nb*sizeof(void*)));
    CK(cudaMalloc(&d_dst,nb*sizeof(void*))); CK(cudaMalloc(&d_nbv,nb*sizeof(size_t)));
    CK(cudaMalloc(&d_csz,nb*sizeof(size_t))); CK(cudaMalloc(&d_dlt,nb*sizeof(double)));
    CK(cudaMemcpy(d_in,h_in.data(),nb*sizeof(void*),cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_cmp,h_cmp.data(),nb*sizeof(void*),cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_dst,h_dst.data(),nb*sizeof(void*),cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_nbv,h_nb.data(),nb*sizeof(size_t),cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_dlt,h_delta.data(),nb*sizeof(double),cudaMemcpyHostToDevice));

    bitcompHandle_t h;
    BK(bitcompCreateBatchPlan(&h,nb,BITCOMP_FP64_DATA,BITCOMP_LOSSY_FP_TO_SIGNED,BITCOMP_DEFAULT_ALGO));
    BK(bitcompSetStream(h,s));
    BK(bitcompBatchCompressLossy_fp64(h,d_in,d_cmp,d_nbv,d_csz,d_dlt)); CK(cudaStreamSynchronize(s));
    std::vector<size_t> csz(nb); CK(cudaMemcpy(csz.data(),d_csz,nb*sizeof(size_t),cudaMemcpyDeviceToHost));
    size_t ctot=0; for(size_t v:csz) ctot+=v;

    const int reps=20;
    CK(cudaEventRecord(e0,s));
    for(int i=0;i<reps;i++) BK(bitcompBatchCompressLossy_fp64(h,d_in,d_cmp,d_nbv,d_csz,d_dlt));
    CK(cudaEventRecord(e1,s)); CK(cudaStreamSynchronize(s));
    const double tc=ms_of(e0,e1)/reps;
    CK(cudaEventRecord(e0,s));
    for(int i=0;i<reps;i++) BK(bitcompBatchUncompress(h,(const void* const*)d_cmp,d_dst));
    CK(cudaEventRecord(e1,s)); CK(cudaStreamSynchronize(s));
    const double td=ms_of(e0,e1)/reps;

    std::vector<double> got(flat.size());
    CK(cudaMemcpy(got.data(),d_out,flat.size()*sizeof(double),cudaMemcpyDeviceToHost));
    double err=0.0, worst_rel=0.0;
    for(size_t i=0;i<boxes.size();++i){
      const double tol=boxes[i].maxabs*tolfac;
      for(size_t j=0;j<box_n[i];++j){
        const double e=fabs(got[box_at[i]+j]-flat[box_at[i]+j]);
        err=fmax(err,e); worst_rel=fmax(worst_rel,e/tol);
      }
    }
    const double mb=tot_elems*8.0/1e6;
    printf("%10.3f %8zu %9.1f %8.2f %10.1f %10.1f %9.3f %12.3e%s\n",
           cmb, nb, stage/1e6, (double)(tot_elems*8)/ctot, mb/tc, mb/td, td, err,
           worst_rel>1.0 ? "  OVER TOL" : "");
    bitcompDestroyPlan(h);
    cudaFree(d_stage); cudaFree(d_in); cudaFree(d_cmp); cudaFree(d_dst);
    cudaFree(d_nbv); cudaFree(d_csz); cudaFree(d_dlt);
  }
  return 0;
}
