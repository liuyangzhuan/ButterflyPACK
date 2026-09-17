// Rectangular C BP interface regression. Compile with BP_TEST_COMPLEX and/or
// BP_TEST_SINGLE to select z/c/s instead of d. Link against butterflypack.
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
#if defined(BP_TEST_COMPLEX) && defined(BP_TEST_SINGLE)
#include "cBPACK_wrapper.h"
#define API(name) c_##name
typedef __complex__ float Scalar;
#elif defined(BP_TEST_COMPLEX)
#include "zBPACK_wrapper.h"
#define API(name) z_##name
typedef __complex__ double Scalar;
#elif defined(BP_TEST_SINGLE)
#include "sBPACK_wrapper.h"
#define API(name) s_##name
typedef float Scalar;
#else
#include "dBPACK_wrapper.h"
#define API(name) d_##name
typedef double Scalar;
#endif
using Complex = std::complex<double>;
static Scalar scalar(Complex x) {
#ifdef BP_TEST_COMPLEX
    Scalar y; __real__ y = x.real(); __imag__ y = x.imag(); return y;
#else
    return x.real();
#endif
}
static Complex value(Scalar x) {
#ifdef BP_TEST_COMPLEX
    return Complex(__real__ x, __imag__ x);
#else
    return Complex(x, 0);
#endif
}
static Complex entry(int i, int j) {
    return value(scalar(Complex(1.0 + 0.02*i + 0.03*j + 0.001*i*j,
                               0.04*i - 0.015*j + 0.0002*i*j)));
}
static Complex input(int i, int rhs) {
    return value(scalar(Complex(0.2 + 0.01*i + 0.03*rhs, 0.04*rhs - 0.002*i)));
}
struct Context { int m, n, calls, split; };
static void sample(int* m, int* n, Scalar* out, C2Fptr context) {
    Context& ctx = *static_cast<Context*>(context);
    const int row = *m > 0 ? *m : *n;
    const int col = *m < 0 ? -*m : -*n;
    if (row < 1 || row > ctx.m || col < 1 || col > ctx.n) MPI_Abort(MPI_COMM_WORLD, 10);
    ++ctx.calls;
    *out = scalar(entry(row, col));
}
static void distance(int* m, int* n, double* out, C2Fptr) {
    *out = std::abs(std::abs(*m) - std::abs(*n)) + 1.0;
}
static void nearfar(int* m, int* n, int* out, C2Fptr context) {
    // Force some further subdivisions instead of compressing every root block.
    const Context& ctx = *static_cast<Context*>(context);
    *out = ctx.split ? ((*m >= 4 && *n >= 4) && (*m + *n) % 3 != 0) : 1;
}
static void seti(F2Cptr& opt, const char* name, int v) { API(c_bpack_set_I_option)(&opt, name, v); }
static void setd(F2Cptr& opt, const char* name, double v) { API(c_bpack_set_D_option)(&opt, name, v); }
static F2Cptr mesh(int n, F2Cptr& opt, F2Cptr& ptree) {
    F2Cptr stats = nullptr, bmat = nullptr, msh = nullptr, ker = nullptr;
    API(c_bpack_createstats)(&stats);
    int ndim = 1, levels = 0, leaf = n, nloc;
    std::vector<int> perm(n);
    Context ctx{n, n, 0, 0};
    API(c_bpack_construct_init_fortran)(&n, &ndim, nullptr, nullptr, &levels, &leaf,
        perm.data(), &nloc, &bmat, &opt, &stats, &msh, &ker, &ptree, nullptr, nullptr, &ctx);
    API(c_bpack_delete)(&bmat);
    API(c_bpack_deletekernelquant)(&ker);
    API(c_bpack_deletestats)(&stats);
    return msh;
}
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> members(size); for (int i=0; i<size; ++i) members[i]=i;
    MPI_Fint comm = MPI_Comm_c2f(MPI_COMM_WORLD);
    F2Cptr ptree = nullptr, opt = nullptr;
    API(c_bpack_createptree)(&size, members.data(), &comm, &ptree);
    API(c_bpack_createoption)(&opt);
    seti(opt,"cpp",1); seti(opt,"verbosity",-1); seti(opt,"nogeo",1);
    seti(opt,"xyzsort",0); seti(opt,"Nmin_leaf",32); seti(opt,"knn",0);
    seti(opt,"format",1); seti(opt,"LRlevel",0); seti(opt,"elem_extract",0);
#ifdef BP_TEST_SINGLE
    setd(opt,"tol_comp",1e-5); const double tolerance = 2e-3;
#else
    setd(opt,"tol_comp",1e-10); const double tolerance = 1e-7;
#endif
    int m=512, n=768, nrhs=3;
    if (rank==0) std::printf("BP regression: M=%d N=%d Nmin_leaf=32 RHS=%d ranks=%d\n",m,n,nrhs,size);
    F2Cptr mshr=mesh(m,opt,ptree), mshc=mesh(n,opt,ptree);
    int failures=0;
    // HODLR forwardN15flag=2 also exercises its error-check pointer.
    for (int test : {1,2,-2,3,5}) {
        const int format = std::abs(test);
        F2Cptr bp=nullptr, stats=nullptr, msh=nullptr, ker=nullptr;
        API(c_bpack_createstats)(&stats);
        seti(opt,"format",format); seti(opt,"forwardN15flag",format==1 ? 2 : 1);
        seti(opt,"nogeo",format==1 ? 3 : 2); seti(opt,"knn",format==1 ? 1 : 0);
        // Both 2 and -2 use format=2; only -2 disables subdivision.
        Context ctx{m,n,0,test!=-2}; int mloc=0, nloc=0;
        // N>M catches accidental indexing of the row array in the column loop.
        std::vector<int> nnsr(m,1), nnsc(n,1); nnsr[0]=0; nnsc[n-1]=0;
        API(c_bp_construct_init)(&m,&n,&mloc,&nloc,nnsr.data(),nnsc.data(),&mshr,&mshc,
            &bp,&opt,&stats,&msh,&ker,&ptree,distance,nearfar,&ctx);
        API(c_bp_construct_element_compute)(&bp,&opt,&stats,&msh,&ker,&ptree,sample,nullptr,&ctx);
        // A numerical product alone also passes for an entirely dense matrix.
        // Require nonzero compressed storage, a reduced rank, and net storage
        // savings; sum local matrix storage across ranks before comparing.
        double mem_comp=0, mem_direct=0, rank_max=0;
        API(c_bpack_getstats)(&stats,"Mem_Comp_for",&mem_comp);
        API(c_bpack_getstats)(&stats,"Mem_Direct_for",&mem_direct);
        API(c_bpack_getstats)(&stats,"Rank_max_Constr",&rank_max);
        double local_mem[2]={mem_comp,mem_direct}, global_mem[2], global_rank;
        MPI_Allreduce(local_mem,global_mem,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&rank_max,&global_rank,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        const double dense_mb=double(m)*n*sizeof(Scalar)/1024.0/1000.0;
        const double storage_ratio=(global_mem[0]+global_mem[1])/dense_mb;
        const bool compressed=global_mem[0]>0 && global_rank>0 &&
            global_rank<std::min(m,n) && storage_ratio<1;
        if (!compressed) ++failures;
        if (rank==0) std::printf("format=%d split=%d compression=%s max_rank=%.0f "
            "compressed_MB=%.6f direct_MB=%.6f dense_MB=%.6f storage_ratio=%.6f\n",
            format,ctx.split,compressed ? "PASS" : "FAIL",global_rank,
            global_mem[0],global_mem[1],dense_mb,storage_ratio);
        for (char trans : {'N','T','C'}) {
            int ni=trans=='N' ? nloc : mloc, no=trans=='N' ? mloc : nloc;
            F2Cptr inmesh=trans=='N' ? mshc : mshr, outmesh=trans=='N' ? mshr : mshc;
            std::vector<Scalar> x(std::max(1,ni*nrhs)), y(std::max(1,no*nrhs));
            for (int j=0; j<nrhs; ++j) for (int i=0; i<ni; ++i) {
                int iloc=i+1, old; API(c_bpack_new2old)(&inmesh,&iloc,&old);
                x[i+j*ni]=scalar(input(old,j));
            }
            const std::vector<Scalar> before=x;
            API(c_bp_mult)(&trans,x.data(),y.data(),&ni,&no,&nrhs,&bp,&opt,&stats,&ptree);
            double error=0, norm=0;
            for (int j=0; j<nrhs; ++j) for (int i=0; i<no; ++i) {
                int iloc=i+1, old; API(c_bpack_new2old)(&outmesh,&iloc,&old);
                Complex exact=0;
                for (int k=1; k<=(trans=='N' ? n : m); ++k) {
                    Complex a=trans=='N' ? entry(old,k) : entry(k,old);
                    if (trans=='C') a=std::conj(a);
                    exact += a*input(k,j);
                }
                error += std::norm(value(y[i+j*no])-exact); norm += std::norm(exact);
            }
            double sums[2]={error,norm}, all[2]; MPI_Allreduce(sums,all,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            double relative=std::sqrt(all[0]/all[1]);
            if (!(relative<tolerance) || x!=before) ++failures;
            if (rank==0) std::printf("format=%d trans=%c relative_error=%.3e\n",format,trans,relative);
        }
        API(c_bp_delete)(&bp); API(c_bp_delete)(&bp);
        if (bp!=nullptr) ++failures;
        API(c_bpack_deletemesh)(&msh); API(c_bpack_deletekernelquant)(&ker);
        API(c_bpack_deletestats)(&stats);
    }
    API(c_bpack_deletemesh)(&mshr); API(c_bpack_deletemesh)(&mshc);
    API(c_bpack_deleteoption)(&opt); API(c_bpack_deleteproctree)(&ptree);
    int total; MPI_Allreduce(&failures,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (rank==0) std::printf("BP wrapper regression: %s (%d ranks)\n",total ? "FAIL" : "PASS",size);
    MPI_Finalize(); return total ? EXIT_FAILURE : EXIT_SUCCESS;
}
