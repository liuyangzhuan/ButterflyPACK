//
// File: G2D.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Nov-2019 14:31:15
//

// Include Files
#include <cmath>
#include "rt_nonfinite.h"
#include <string.h>
#include <math.h>
#include "rt_defines.h"
#include "G2D.h"

// Variable Definitions
static const double dv0[30] = { 0.6299605249474366, 0.25198420997897464,
  0.15479030041565583, 0.11071306241615901, 0.085730939552739485,
  0.069716131695868433, 0.058608567189371359, 0.050469887353631067,
  0.044260058068915482, 0.039372066154350994, 0.035428319592445537,
  0.032181885750209825, 0.029464624079115768, 0.027158167711293448,
  0.025176827297386177, 0.023457075530607888, 0.021950839013490719,
  0.020621082823564625, 0.019438824089788084, 0.018381063380068317,
  0.017429321323196318, 0.016568583778661234, 0.015786528598791844,
  0.015072950149409559, 0.014419325083995464, 0.013818480573534178,
  0.013264337899427657, 0.012751712197049864, 0.012276154531876277,
  0.011833826239848241 };

static const double dv1[120] = { 1.0, -0.20833333333333334, 0.125,
  0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
  1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
  -11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
  -28.212072558200244, 84.636217674600729, -91.818241543240021,
  42.534998745388457, -7.3687943594796321, 0.22710800170898438,
  212.57013003921713, -765.25246814118168, 1059.9904525279999,
  -699.57962737613252, 218.19051174421159, -26.491430486951554,
  0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
  -13586.550006434138, 11655.393336864534, -5305.646978613403,
  1200.9029132163525, -108.09091978839466, 1.7277275025844574,
  20204.291330966149, -96980.598388637518, 192547.00123253153,
  -203400.17728041555, 122200.46498301746, -41192.65496889755,
  7109.5143024893641, -493.915304773088, 6.074042001273483, -242919.18790055133,
  1.3117636146629772E+6, -2.9980159185381066E+6, 3.7632712976564039E+6,
  -2.8135632265865342E+6, 1.2683652733216248E+6, -331645.17248456361,
  45218.768981362729, -2499.8304818112097, 24.380529699556064,
  3.2844698530720379E+6, -1.9706819118432228E+7, 5.0952602492664643E+7,
  -7.4105148211532652E+7, 6.6344512274729028E+7, -3.7567176660763353E+7,
  1.3288767166421818E+7, -2.7856181280864547E+6, 308186.40461266239,
  -13886.08975371704, 110.01714026924674, -4.932925366450996E+7,
  3.2557307418576574E+8, -9.394623596815784E+8, 1.55359689957058E+9,
  -1.6210805521083372E+9, 1.1068428168230145E+9, -4.9588978427503031E+8,
  1.4206290779753309E+8, -2.447406272573873E+7, 2.2437681779224495E+6,
  -84005.433603024081, 551.33589612202059, 8.1478909611831212E+8,
  -5.8664814920518475E+9, 1.8688207509295826E+10, -3.4632043388158775E+10,
  4.1280185579753975E+10, -3.3026599749800724E+10, 1.79542137311556E+10,
  -6.5632937926192846E+9, 1.5592798648792574E+9, -2.2510566188941526E+8,
  1.7395107553978164E+7, -549842.32757228869, 3038.0905109223841,
  -1.4679261247695616E+10, 1.144982377320258E+11, -3.9909617522446649E+11,
  8.1921866954857727E+11, -1.0983751560812233E+12, 1.0081581068653821E+12,
  -6.4536486924537646E+11, 2.8790064990615057E+11, -8.786707217802327E+10,
  1.7634730606834969E+10, -2.1671649832237949E+9, 1.4315787671888897E+8,
  -3.8718334425726128E+6, 18257.755474293175, 2.86464035717679E+11,
  -2.4062979000285039E+12, 9.1093411852398984E+12, -2.0516899410934438E+13,
  3.056512551993532E+13, -3.166708858478516E+13, 2.334836404458184E+13,
  -1.2320491305598287E+13, 4.6127257808491318E+12, -1.1965528801961816E+12,
  2.0591450323241E+11, -2.1822927757529224E+10, 1.2470092935127103E+9,
  -2.9188388122220814E+7, 118838.42625678325 };

static const double dv2[8] = { 0.57721566490153287, -0.042002635034095237,
  -0.042197734555544333, 0.0072189432466631, -0.00021524167411495098,
  -2.0134854780788239E-5, 1.1330272319816959E-6, 6.1160951044814161E-9 };

// Function Declarations
static double BesselJ0(double x);
static int b_cbknu(const creal_T z, double fnu, int kode, int nin, creal_T y[2],
                   double tol, double elim, double alim);
static int b_ckscl(const creal_T zr, int nin, creal_T y[2], double ascle, double
                   tol, double elim);
static void b_cosh(creal_T *x);
static void b_cunhj(const creal_T z, double fnu, int ipmtr, double tol, creal_T *
                    phi, creal_T *arg, creal_T *zeta1, creal_T *zeta2, creal_T
                    *asum, creal_T *bsum);
static void b_cuni2(const creal_T z, double fnu, int kode, int nin, creal_T y[2],
                    double fnul, double tol, double elim, double alim, int
                    *nlast, int *nz);
static void b_cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double
                    tol, int *init, creal_T cwrk[16], creal_T *phi, creal_T
                    *zeta1, creal_T *zeta2, creal_T *summ);
static int b_cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin,
                   creal_T y[2], double tol, double elim, double alim);
static void b_exp(creal_T *x);
static void b_fix(double *x);
static void b_log(double *x);
static void b_sinh(creal_T *x);
static void b_sqrt(double *x);
static void c_log(creal_T *x);
static void c_sqrt(creal_T *x);
static int cacai(const creal_T z, double fnu, int kode, int mr, creal_T *y,
                 double rl, double tol, double elim, double alim);
static creal_T cairy(const creal_T z, int id, int kode);
static int casyi(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                 double rl, double tol, double elim);
static void cbesj(const creal_T z, double fnu, int kode, creal_T *cy, int *nz,
                  int *ierr);
static int cbinu(const creal_T z, double fnu, int kode, creal_T *cy, double rl,
                 double fnul, double tol, double elim, double alim);
static void cbknu(const creal_T z, double fnu, int kode, double alim, creal_T *y,
                  int *nz);
static void cbuni(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  int nui, double fnul, double tol, double elim, double alim,
                  int *nlast, int *nz);
static int ckscl(const creal_T zr, creal_T *y, double elim);
static int cmlri(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                 double tol);
static void crati(const creal_T z, double fnu, int nin, creal_T *cy, double tol);
static int cseri(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                 double tol, double elim, double alim);
static int cuchk(const creal_T y, double ascle, double tol);
static void cunhj(const creal_T z, double fnu, double tol, creal_T *phi, creal_T
                  *arg, creal_T *zeta1, creal_T *zeta2);
static void cuni1(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double fnul, double tol, double elim, double alim, int *nlast,
                  int *nz);
static void cuni2(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double fnul, double tol, double elim, double alim, int *nlast,
                  int *nz);
static void cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double tol,
                  int init, creal_T *phi, creal_T *zeta1, creal_T *zeta2);
static int cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin,
                 creal_T *y, double tol, double elim, double alim);
static void gammaln(double *x);
static double rt_atan2d_snf(double u0, double u1);
static double rt_hypotd_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);
static double taunum2D(double x, double y, double x0, double b_y0);

// Function Definitions

//
// Arguments    : double x
// Return Type  : double
//
static double BesselJ0(double x)
{
  double J0;
  double ax;
  double z;
  double y;
  ax = std::abs(x);
  if (ax < 8.0) {
    y = x * x;
    J0 = (5.7568490574E+10 + y * (-1.3362690354E+10 + y * (6.516196407E+8 + y *
            (-1.121442418E+7 + y * (77392.33017 + y * -184.9052456))))) /
      (5.7568490411E+10 + y * (1.029532985E+9 + y * (9.494680718E+6 + y *
         (59272.64853 + y * (267.8532712 + y)))));
  } else {
    z = 8.0 / ax;
    y = z * z;
    J0 = std::sqrt(0.636619772 / ax) * (std::cos(ax - 0.785398164) * (1.0 + y *
      (-0.1098628627 + y * (0.2734510407 + y * (-0.2073370639 + y * 0.2093887211))))
      - z * std::sin(ax - 0.785398164) * (-0.1562499995 + y * (0.1430488765 + y *
      (-0.6911147651 + y * (0.7621095161 - y * 0.934935152)))));
  }

  return J0;
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T y[2]
//                double tol
//                double elim
//                double alim
// Return Type  : int
//
static int b_cbknu(const creal_T z, double fnu, int kode, int nin, creal_T y[2],
                   double tol, double elim, double alim)
{
  int nz;
  int n;
  double yy;
  double caz;
  double cssr[3];
  double csrr[3];
  double bry[3];
  int iflag;
  double rz_re;
  double tm;
  double rz_im;
  double g1;
  int inu;
  double fk;
  double dnu;
  boolean_T goto_mw110;
  double dnu2;
  double ak;
  double fhs;
  double s1_re;
  double s1_im;
  double s2_re;
  double s2_im;
  creal_T zd;
  double ck_re;
  double ck_im;
  int inub;
  boolean_T goto_mw225;
  boolean_T goto_mw240;
  boolean_T goto_mw270;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  double fc;
  creal_T p2;
  creal_T cs;
  creal_T fmu;
  double coef_re;
  double coef_im;
  int kflag;
  double t2;
  double a1;
  double t1;
  int i;
  boolean_T exitg1;
  creal_T f;
  creal_T p1;
  int k;
  boolean_T earlyExit;
  int j;
  creal_T cy[2];
  double b_g1;
  boolean_T b_guard1 = false;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  yy = z.im;
  caz = rt_hypotd_snf(z.re, z.im);
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry[0] = 2.2250738585072014E-305 / tol;
  bry[1] = tol / 2.2250738585072014E-305;
  bry[2] = 1.7976931348623157E+308;
  nz = 0;
  iflag = 0;
  if (z.im == 0.0) {
    rz_re = 2.0 / z.re;
    rz_im = 0.0;
  } else if (z.re == 0.0) {
    rz_re = 0.0;
    rz_im = -(2.0 / z.im);
  } else {
    tm = std::abs(z.re);
    g1 = std::abs(z.im);
    if (tm > g1) {
      fk = z.im / z.re;
      g1 = z.re + fk * z.im;
      rz_re = (2.0 + fk * 0.0) / g1;
      rz_im = (0.0 - fk * 2.0) / g1;
    } else if (g1 == tm) {
      if (z.re > 0.0) {
        g1 = 0.5;
      } else {
        g1 = -0.5;
      }

      if (z.im > 0.0) {
        fk = 0.5;
      } else {
        fk = -0.5;
      }

      rz_re = (2.0 * g1 + 0.0 * fk) / tm;
      rz_im = (0.0 * g1 - 2.0 * fk) / tm;
    } else {
      fk = z.re / z.im;
      g1 = z.im + fk * z.re;
      rz_re = fk * 2.0 / g1;
      rz_im = (fk * 0.0 - 2.0) / g1;
    }
  }

  inu = (int)(fnu + 0.5);
  dnu = fnu - (double)inu;
  goto_mw110 = (std::abs(dnu) == 0.5);
  if ((!goto_mw110) && (std::abs(dnu) > tol)) {
    dnu2 = dnu * dnu;
  } else {
    dnu2 = 0.0;
  }

  ak = 1.0;
  fhs = 0.0;
  s1_re = 0.0;
  s1_im = 0.0;
  s2_re = 0.0;
  s2_im = 0.0;
  zd.re = 0.0;
  zd.im = 0.0;
  ck_re = 1.0;
  ck_im = 0.0;
  inub = 1;
  goto_mw225 = false;
  goto_mw240 = false;
  goto_mw270 = false;
  guard1 = false;
  guard2 = false;
  if (goto_mw110 || (caz > 2.0)) {
    goto_mw110 = false;
    cs = z;
    c_sqrt(&cs);
    if (cs.im == 0.0) {
      coef_re = 1.2533141373155001 / cs.re;
      coef_im = 0.0;
    } else if (cs.re == 0.0) {
      coef_re = 0.0;
      coef_im = -(1.2533141373155001 / cs.im);
    } else {
      tm = std::abs(cs.re);
      g1 = std::abs(cs.im);
      if (tm > g1) {
        fk = cs.im / cs.re;
        g1 = cs.re + fk * cs.im;
        coef_re = (1.2533141373155001 + fk * 0.0) / g1;
        coef_im = (0.0 - fk * 1.2533141373155001) / g1;
      } else if (g1 == tm) {
        if (cs.re > 0.0) {
          g1 = 0.5;
        } else {
          g1 = -0.5;
        }

        if (cs.im > 0.0) {
          fk = 0.5;
        } else {
          fk = -0.5;
        }

        coef_re = (1.2533141373155001 * g1 + 0.0 * fk) / tm;
        coef_im = (0.0 * g1 - 1.2533141373155001 * fk) / tm;
      } else {
        fk = cs.re / cs.im;
        g1 = cs.im + fk * cs.re;
        coef_re = fk * 1.2533141373155001 / g1;
        coef_im = (fk * 0.0 - 1.2533141373155001) / g1;
      }
    }

    kflag = 1;
    if (!(kode == 2)) {
      if (z.re > alim) {
        iflag = 1;
      } else {
        a1 = std::exp(-z.re);
        fmu.re = std::cos(z.im) * a1;
        fmu.im = std::sin(z.im) * -a1;
        fc = coef_re;
        coef_re = coef_re * fmu.re - coef_im * fmu.im;
        coef_im = fc * fmu.im + coef_im * fmu.re;
      }
    }

    if (std::abs(dnu) == 0.5) {
      s1_re = coef_re;
      s1_im = coef_im;
      s2_re = coef_re;
      s2_im = coef_im;
      goto_mw110 = true;
    } else {
      ak = std::abs(std::cos(3.1415926535897931 * dnu));
      if (ak == 0.0) {
        s1_re = coef_re;
        s1_im = coef_im;
        s2_re = coef_re;
        s2_im = coef_im;
        goto_mw110 = true;
      } else {
        fhs = std::abs(0.25 - dnu2);
        if (fhs == 0.0) {
          s1_re = coef_re;
          s1_im = coef_im;
          s2_re = coef_re;
          s2_im = coef_im;
          goto_mw110 = true;
        }
      }
    }

    if (!goto_mw110) {
      if (z.re != 0.0) {
        t1 = std::abs(std::atan(z.im / z.re));
      } else {
        t1 = 1.5707963267948966;
      }

      if (28.66666665740641 > caz) {
        g1 = caz;
        b_sqrt(&g1);
        b_sqrt(&g1);
        g1 = 1.8976999933151775 * ak / (tol * g1);
        b_log(&g1);
        ak = (g1 + caz * std::cos(3.0 * t1 / (1.0 + caz)) / (1.0 + 0.008 * caz))
          / std::cos(14.7 * t1 / (28.0 + caz));
        fk = 0.12125 * ak * ak / caz + 1.5;
        guard2 = true;
      } else {
        g1 = ak / (3.1415926535897931 * caz * tol);
        fk = 1.0;
        if (g1 >= 1.0) {
          ak = 2.0;
          fc = (caz + caz) + 2.0;
          a1 = 0.0;
          t2 = 1.0;
          earlyExit = true;
          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i < 31)) {
            tm = t2;
            t2 = fc / (fk + 1.0) * t2 - fhs / ak * a1;
            a1 = tm;
            fc += 2.0;
            ak = ((ak + fk) + fk) + 2.0;
            fhs = (fhs + fk) + fk;
            fk++;
            if (g1 < std::abs(t2) * fk) {
              earlyExit = false;
              exitg1 = true;
            } else {
              i++;
            }
          }

          if (earlyExit) {
            nz = -2;
          } else {
            g1 = 28.66666665740641 / caz;
            b_sqrt(&g1);
            fk += 1.909859317102744 * t1 * g1;
            fhs = std::abs(0.25 - dnu2);
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }
    } else {
      guard1 = true;
    }
  } else {
    fc = 1.0;
    p2.re = rz_re;
    p2.im = rz_im;
    c_log(&p2);
    fmu.re = dnu * p2.re;
    fmu.im = dnu * p2.im;
    if (dnu != 0.0) {
      fc = dnu * 3.1415926535897931;
      fc /= std::sin(fc);
      cs = fmu;
      b_sinh(&cs);
      p2.re = 1.0 / dnu * cs.re;
      p2.im = 1.0 / dnu * cs.im;
    }

    g1 = 1.0 + dnu;
    gammaln(&g1);
    t2 = std::exp(-g1);
    t1 = 1.0 / (t2 * fc);
    if (std::abs(dnu) > 0.1) {
      g1 = (t1 - t2) / (dnu + dnu);
    } else {
      fk = 0.57721566490153287;
      i = 2;
      exitg1 = false;
      while ((!exitg1) && (i < 9)) {
        ak *= dnu2;
        tm = dv2[i - 1] * ak;
        fk += tm;
        if (std::abs(tm) < tol) {
          exitg1 = true;
        } else {
          i++;
        }
      }

      g1 = -fk;
    }

    g1 *= fc;
    cs = fmu;
    b_cosh(&cs);
    f.re = p2.re * (0.5 * (t1 + t2) * fc) + g1 * cs.re;
    f.im = p2.im * (0.5 * (t1 + t2) * fc) + g1 * cs.im;
    b_exp(&fmu);
    p1.re = 0.5 / t2 * fmu.re;
    p1.im = 0.5 / t2 * fmu.im;
    ak = 0.5 / t1;
    if (fmu.im == 0.0) {
      cs.re = ak / fmu.re;
      cs.im = 0.0;
    } else if (fmu.re == 0.0) {
      if (ak == 0.0) {
        cs.re = 0.0 / fmu.im;
        cs.im = 0.0;
      } else {
        cs.re = 0.0;
        cs.im = -(ak / fmu.im);
      }
    } else {
      tm = std::abs(fmu.re);
      g1 = std::abs(fmu.im);
      if (tm > g1) {
        fk = fmu.im / fmu.re;
        g1 = fmu.re + fk * fmu.im;
        cs.re = (ak + fk * 0.0) / g1;
        cs.im = (0.0 - fk * ak) / g1;
      } else if (g1 == tm) {
        if (fmu.re > 0.0) {
          g1 = 0.5;
        } else {
          g1 = -0.5;
        }

        if (fmu.im > 0.0) {
          fk = 0.5;
        } else {
          fk = -0.5;
        }

        cs.re = (ak * g1 + 0.0 * fk) / tm;
        cs.im = (0.0 * g1 - ak * fk) / tm;
      } else {
        fk = fmu.re / fmu.im;
        g1 = fmu.im + fk * fmu.re;
        cs.re = fk * ak / g1;
        cs.im = (fk * 0.0 - ak) / g1;
      }
    }

    s1_re = f.re;
    s1_im = f.im;
    s2_re = p1.re;
    s2_im = p1.im;
    ak = 1.0;
    a1 = 1.0;
    t2 = 1.0 - dnu2;
    if ((inu > 0) || (n > 1)) {
      if (caz >= tol) {
        coef_re = 0.25 * (z.re * z.re - z.im * z.im);
        coef_im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          f.re *= ak;
          f.im *= ak;
          f.re += p1.re;
          f.im += p1.im;
          f.re += cs.re;
          f.im += cs.im;
          f.re *= 1.0 / t2;
          f.im *= 1.0 / t2;
          p1.re *= 1.0 / (ak - dnu);
          p1.im *= 1.0 / (ak - dnu);
          cs.re *= 1.0 / (ak + dnu);
          cs.im *= 1.0 / (ak + dnu);
          g1 = ck_re;
          ck_re = ck_re * coef_re - ck_im * coef_im;
          ck_im = g1 * coef_im + ck_im * coef_re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * f.re - ck_im * f.im;
          s1_im += ck_re * f.im + ck_im * f.re;
          g1 = p1.re - f.re * ak;
          fk = p1.im - f.im * ak;
          s2_re += g1 * ck_re - fk * ck_im;
          s2_im += g1 * ck_im + fk * ck_re;
          a1 = a1 * t1 / ak;
          t2 = ((t2 + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > tol));
      }

      kflag = 1;
      if ((fnu + 1.0) * std::abs(p2.re) > alim) {
        kflag = 2;
      }

      t2 = cssr[kflag] * s2_re;
      s2_im *= cssr[kflag];
      s2_re = t2 * rz_re - s2_im * rz_im;
      s2_im = t2 * rz_im + s2_im * rz_re;
      s1_re *= cssr[kflag];
      s1_im *= cssr[kflag];
      if (kode != 1) {
        f = z;
        b_exp(&f);
        g1 = s1_re;
        s1_re = s1_re * f.re - s1_im * f.im;
        s1_im = g1 * f.im + s1_im * f.re;
        t2 = s2_re;
        s2_re = s2_re * f.re - s2_im * f.im;
        s2_im = t2 * f.im + s2_im * f.re;
      }

      goto_mw110 = true;
      guard1 = true;
    } else {
      if (caz >= tol) {
        coef_re = 0.25 * (z.re * z.re - z.im * z.im);
        coef_im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          f.re *= ak;
          f.im *= ak;
          f.re += p1.re;
          f.im += p1.im;
          f.re += cs.re;
          f.im += cs.im;
          f.re *= 1.0 / t2;
          f.im *= 1.0 / t2;
          p1.re *= 1.0 / (ak - dnu);
          p1.im *= 1.0 / (ak - dnu);
          cs.re *= 1.0 / (ak + dnu);
          cs.im *= 1.0 / (ak + dnu);
          g1 = ck_re;
          ck_re = ck_re * coef_re - ck_im * coef_im;
          ck_im = g1 * coef_im + ck_im * coef_re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * f.re - ck_im * f.im;
          s1_im += ck_re * f.im + ck_im * f.re;
          a1 = a1 * t1 / ak;
          t2 = ((t2 + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > tol));
      }

      y[0].re = s1_re;
      y[0].im = s1_im;
      if (kode != 1) {
        y[0] = z;
        b_exp(&y[0]);
        g1 = y[0].re;
        t2 = y[0].im;
        y[0].re = g1 * s1_re - t2 * s1_im;
        y[0].im = g1 * s1_im + t2 * s1_re;
      }
    }
  }

  if (guard2) {
    b_fix(&fk);
    k = (int)fk;
    fk = k;
    ak = (double)k * (double)k;
    p1.re = 0.0;
    p1.im = 0.0;
    p2.re = tol;
    p2.im = 0.0;
    cs = p2;
    for (i = 1; i <= k; i++) {
      a1 = ak - fk;
      fc = 2.0 / (fk + 1.0);
      fmu = p2;
      t2 = (fk + z.re) * fc;
      g1 = yy * fc;
      fc = p2.re;
      p2.re = p2.re * t2 - p2.im * g1;
      p2.im = fc * g1 + p2.im * t2;
      p2.re -= p1.re;
      p2.im -= p1.im;
      p2.re *= (ak + fk) / (a1 + fhs);
      p2.im *= (ak + fk) / (a1 + fhs);
      p1 = fmu;
      cs.re += p2.re;
      cs.im += p2.im;
      ak = (a1 - fk) + 1.0;
      fk--;
    }

    fmu.re = 1.0 / rt_hypotd_snf(cs.re, cs.im);
    cs.im = -cs.im;
    g1 = cs.re;
    cs.re = cs.re * fmu.re - cs.im * 0.0;
    cs.im = g1 * 0.0 + cs.im * fmu.re;
    t2 = fmu.re * p2.re - 0.0 * p2.im;
    g1 = fmu.re * p2.im + 0.0 * p2.re;
    fc = coef_re * t2 - coef_im * g1;
    coef_im = coef_re * g1 + coef_im * t2;
    s1_re = fc * cs.re - coef_im * cs.im;
    s1_im = fc * cs.im + coef_im * cs.re;
    if ((inu > 0) || (n > 1)) {
      fmu.re = 1.0 / rt_hypotd_snf(p2.re, p2.im);
      g1 = p1.re;
      p1.re = p1.re * fmu.re - p1.im * 0.0;
      p1.im = g1 * 0.0 + p1.im * fmu.re;
      p2.im = -p2.im;
      fc = p2.re;
      p2.re = p2.re * fmu.re - p2.im * 0.0;
      p2.im = fc * 0.0 + p2.im * fmu.re;
      fk = p1.re * p2.im + p1.im * p2.re;
      ak = (dnu + 0.5) - (p1.re * p2.re - p1.im * p2.im);
      t2 = 0.0 - (p1.re * p2.im + p1.im * p2.re);
      if (z.im == 0.0) {
        if (0.0 - fk == 0.0) {
          fc = ak / z.re;
          g1 = 0.0;
        } else if (ak == 0.0) {
          fc = 0.0;
          g1 = (0.0 - fk) / z.re;
        } else {
          fc = ak / z.re;
          g1 = (0.0 - fk) / z.re;
        }
      } else if (z.re == 0.0) {
        if (ak == 0.0) {
          fc = (0.0 - fk) / z.im;
          g1 = 0.0;
        } else if (0.0 - fk == 0.0) {
          fc = 0.0;
          g1 = -(ak / z.im);
        } else {
          fc = (0.0 - fk) / z.im;
          g1 = -(ak / z.im);
        }
      } else {
        tm = std::abs(z.re);
        g1 = std::abs(z.im);
        if (tm > g1) {
          fk = z.im / z.re;
          g1 = z.re + fk * z.im;
          fc = (ak + fk * t2) / g1;
          g1 = (t2 - fk * ak) / g1;
        } else if (g1 == tm) {
          if (z.re > 0.0) {
            g1 = 0.5;
          } else {
            g1 = -0.5;
          }

          if (z.im > 0.0) {
            fk = 0.5;
          } else {
            fk = -0.5;
          }

          fc = (ak * g1 + t2 * fk) / tm;
          g1 = (t2 * g1 - ak * fk) / tm;
        } else {
          fk = z.re / z.im;
          g1 = z.im + fk * z.re;
          fc = (fk * ak + t2) / g1;
          g1 = (fk * t2 - ak) / g1;
        }
      }

      fc++;
      s2_re = s1_re * fc - s1_im * g1;
      s2_im = s1_re * g1 + s1_im * fc;
      goto_mw110 = true;
    } else {
      zd = z;
      if (iflag == 1) {
        goto_mw270 = true;
      } else {
        goto_mw240 = true;
      }
    }

    guard1 = true;
  }

  if (guard1) {
    if (goto_mw240 || goto_mw270) {
    } else if (goto_mw110) {
      ck_re = (dnu + 1.0) * rz_re;
      ck_im = (dnu + 1.0) * rz_im;
      if (n == 1) {
        inu--;
      }

      if (inu > 0) {
        if (iflag == 1) {
          g1 = 0.5 * elim;
          coef_re = std::exp(-elim);
          zd = z;
          fc = z.re;
          k = inu;
          j = 1;
          for (i = 0; i < 2; i++) {
            cy[i].re = 0.0;
            cy[i].im = 0.0;
          }

          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i <= inu)) {
            f.re = s2_re;
            f.im = s2_im;
            t2 = s2_re;
            s2_re = s2_re * ck_re - s2_im * ck_im;
            s2_im = t2 * ck_im + s2_im * ck_re;
            s2_re += s1_re;
            s2_im += s1_im;
            s1_re = f.re;
            s1_im = f.im;
            ck_re += rz_re;
            ck_im += rz_im;
            t2 = rt_hypotd_snf(s2_re, s2_im);
            b_log(&t2);
            b_guard1 = false;
            if (-fc + t2 >= -elim) {
              cs.re = s2_re;
              cs.im = s2_im;
              c_log(&cs);
              p2.re = cs.re + -zd.re;
              p2.im = cs.im + -zd.im;
              p1.re = std::exp(p2.re) / tol * std::cos(p2.im);
              p1.im = std::exp(p2.re) / tol * std::sin(p2.im);
              if (cuchk(p1, bry[0], tol) == 0) {
                j = 1 - j;
                cy[j] = p1;
                if (k == i - 1) {
                  kflag = 0;
                  inub = i + 1;
                  s2_re = cy[j].re;
                  s2_im = cy[j].im;
                  j = 1 - j;
                  s1_re = cy[j].re;
                  s1_im = cy[j].im;
                  if (i + 1 <= inu) {
                    goto_mw225 = true;
                  } else {
                    if (n == 1) {
                      s1_re = s2_re;
                      s1_im = s2_im;
                    }

                    goto_mw240 = true;
                  }

                  exitg1 = true;
                } else {
                  k = i;
                  i++;
                }
              } else {
                b_guard1 = true;
              }
            } else {
              b_guard1 = true;
            }

            if (b_guard1) {
              if (t2 >= g1) {
                fc -= elim;
                s1_re = f.re * coef_re - f.im * 0.0;
                s1_im = f.re * 0.0 + f.im * coef_re;
                t2 = s2_re;
                s2_re = s2_re * coef_re - s2_im * 0.0;
                s2_im = t2 * 0.0 + s2_im * coef_re;
                zd.re = fc;
                zd.im = yy;
              }

              i++;
            }
          }

          if (goto_mw225 || goto_mw240) {
          } else {
            if (n == 1) {
              s1_re = s2_re;
              s1_im = s2_im;
            }

            goto_mw270 = true;
          }
        } else {
          goto_mw225 = true;
        }
      }

      if (goto_mw225 || goto_mw240 || goto_mw270) {
      } else {
        if (n == 1) {
          s1_re = s2_re;
          s1_im = s2_im;
        }

        zd = z;
        if (iflag != 1) {
          goto_mw240 = true;
        }
      }
    } else {
      goto_mw225 = true;
    }

    if (goto_mw225 || goto_mw240) {
      if (goto_mw225) {
        p1.re = csrr[kflag];
        fc = bry[kflag];
        while (inub <= inu) {
          f.re = s2_re;
          f.im = s2_im;
          t2 = s2_re;
          s2_re = s2_re * ck_re - s2_im * ck_im;
          s2_im = t2 * ck_im + s2_im * ck_re;
          s2_re += s1_re;
          s2_im += s1_im;
          s1_re = f.re;
          s1_im = f.im;
          ck_re += rz_re;
          ck_im += rz_im;
          if (kflag + 1 < 3) {
            p2.re = s2_re * p1.re - s2_im * 0.0;
            p2.im = s2_re * 0.0 + s2_im * p1.re;
            g1 = std::abs(p2.re);
            t2 = std::abs(p2.im);
            if ((g1 > t2) || rtIsNaN(t2)) {
              b_g1 = g1;
            } else {
              b_g1 = t2;
            }

            if (b_g1 > fc) {
              kflag++;
              fc = bry[kflag];
              s1_re = cssr[kflag] * (f.re * p1.re - f.im * 0.0);
              s1_im = cssr[kflag] * (f.re * 0.0 + f.im * p1.re);
              s2_re = cssr[kflag] * p2.re;
              s2_im = cssr[kflag] * p2.im;
              p1.re = csrr[kflag];
            }
          }

          inub++;
        }

        if (n == 1) {
          s1_re = s2_re;
          s1_im = s2_im;
        }
      }

      y[0].re = csrr[kflag] * s1_re;
      y[0].im = csrr[kflag] * s1_im;
      if (n == 2) {
        y[1].re = csrr[kflag] * s2_re;
        y[1].im = csrr[kflag] * s2_im;
      }
    } else {
      y[0].re = s1_re;
      y[0].im = s1_im;
      if (n != 1) {
        y[1].re = s2_re;
        y[1].im = s2_im;
      }

      nz = b_ckscl(zd, n, y, bry[0], tol, elim);
      if (n <= nz) {
      } else {
        y[nz].re *= tol;
        y[nz].im *= tol;
        if ((n >= 0) && (nz < n - MAX_int32_T)) {
          k = MAX_int32_T;
        } else if ((n < 0) && (nz > n - MIN_int32_T)) {
          k = MIN_int32_T;
        } else {
          k = n - nz;
        }

        if (k >= 2) {
          y[nz + 1].re *= tol;
          y[nz + 1].im *= tol;
        }
      }
    }
  }

  return nz;
}

//
// Arguments    : const creal_T zr
//                int nin
//                creal_T y[2]
//                double ascle
//                double tol
//                double elim
// Return Type  : int
//
static int b_ckscl(const creal_T zr, int nin, creal_T y[2], double ascle, double
                   tol, double elim)
{
  int nz;
  int n;
  int ic;
  creal_T s1;
  double s1_re;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  ic = 0;
  nz = 0;
  if (n >= 1) {
    s1 = y[0];
    nz = 1;
    y[0].re = 0.0;
    y[0].im = 0.0;
    if (-zr.re + std::log(rt_hypotd_snf(s1.re, s1.im)) >= -elim) {
      c_log(&s1);
      s1.re += -zr.re;
      s1.im += -zr.im;
      s1_re = s1.re;
      s1.re = std::exp(s1.re) / tol * std::cos(s1.im);
      s1.im = std::exp(s1_re) / tol * std::sin(s1.im);
      if (cuchk(s1, ascle, tol) == 0) {
        y[0] = s1;
        nz = 0;
        ic = 1;
      }
    }
  }

  if (n == 2) {
    if (nz == 1) {
      nz = 2;
    } else {
      nz = 1;
    }

    s1 = y[1];
    y[1].re = 0.0;
    y[1].im = 0.0;
    if (-zr.re + std::log(rt_hypotd_snf(s1.re, s1.im)) >= -elim) {
      c_log(&s1);
      s1.re += -zr.re;
      s1.im += -zr.im;
      s1_re = s1.re;
      s1.re = std::exp(s1.re) / tol * std::cos(s1.im);
      s1.im = std::exp(s1_re) / tol * std::sin(s1.im);
      if (cuchk(s1, ascle, tol) == 0) {
        y[1] = s1;
        nz = 0;
        ic = 2;
      }
    }

    if (ic < 2) {
      y[0].re = 0.0;
      y[0].im = 0.0;
      nz = 2;
    }
  }

  return nz;
}

//
// Arguments    : creal_T *x
// Return Type  : void
//
static void b_cosh(creal_T *x)
{
  double x_re;
  double x_im;
  if (x->im == 0.0) {
    x->re = std::cosh(x->re);
    x->im = 0.0;
  } else {
    x_re = x->re;
    x_im = x->im;
    x->re = std::cosh(x->re) * std::cos(x->im);
    x->im = std::sinh(x_re) * std::sin(x_im);
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int ipmtr
//                double tol
//                creal_T *phi
//                creal_T *arg
//                creal_T *zeta1
//                creal_T *zeta2
//                creal_T *asum
//                creal_T *bsum
// Return Type  : void
//
static void b_cunhj(const creal_T z, double fnu, int ipmtr, double tol, creal_T *
                    phi, creal_T *arg, creal_T *zeta1, creal_T *zeta2, creal_T
                    *asum, creal_T *bsum)
{
  double rfnu;
  creal_T up[14];
  double ac;
  creal_T zb;
  double rfnu2;
  double fn23;
  double rfn13_re;
  creal_T w2;
  int k;
  creal_T w;
  creal_T p[30];
  double ap[30];
  double pp;
  creal_T suma;
  int l1;
  boolean_T exitg1;
  double brm;
  creal_T zc;
  double bim;
  double d;
  double zth_re;
  double zth_im;
  static const double BETA[210] = { 0.017998872141355329, 0.0055996491106438812,
    0.0028850140223113277, 0.0018009660676105393, 0.0012475311058919921,
    0.00092287887657293828, 0.00071443042172728737, 0.00057178728178970488,
    0.00046943100760648155, 0.00039323283546291665, 0.00033481888931829768,
    0.00028895214849575154, 0.0002522116155495733, 0.00022228058079888332,
    0.00019754183803306251, 0.00017683685501971802, 0.00015931689966182109,
    0.00014434793019733397, 0.00013144806811996539, 0.00012024544494930288,
    0.00011044914450459939, 0.00010182877074056726, 9.4199822420423752E-5,
    8.7413054575383449E-5, 8.1346626216280142E-5, 7.590022696462193E-5,
    7.0990630063415351E-5, 6.6548287484246817E-5, 6.25146958969275E-5,
    5.8840339442625178E-5, -0.0014928295321342917, -0.00087820470954638936,
    -0.00050291654957203464, -0.000294822138512746, -0.00017546399697078284,
    -0.00010400855046081644, -5.961419530464579E-5, -3.1203892907609836E-5,
    -1.2608973598023005E-5, -2.4289260857573037E-7, 8.059961654142736E-6,
    1.3650700926214739E-5, 1.7396412547292627E-5, 1.9867297884213378E-5,
    2.1446326379082263E-5, 2.2395465923245652E-5, 2.2896778381471263E-5,
    2.3078538981117782E-5, 2.3032197608090914E-5, 2.2823607372034874E-5,
    2.2500588110529241E-5, 2.2098101536199144E-5, 2.164184274481039E-5,
    2.1150764925622083E-5, 2.0638874978217072E-5, 2.0116524199708165E-5,
    1.9591345014117925E-5, 1.9068936791043675E-5, 1.8553371964163667E-5,
    1.8047572225967421E-5, 0.0005522130767212928, 0.00044793258155238465,
    0.00027952065399202059, 0.0001524681561984466, 6.932711056570436E-5,
    1.7625868306999139E-5, -1.3574499634326914E-5, -3.1797241335042717E-5,
    -4.188618616966934E-5, -4.6900488937914104E-5, -4.8766544741378735E-5,
    -4.8701003118673505E-5, -4.7475562089008661E-5, -4.5581305813862843E-5,
    -4.33309644511266E-5, -4.0923019315775034E-5, -3.848226386032213E-5,
    -3.6085716753541052E-5, -3.3779330612336739E-5, -3.1588856077210961E-5,
    -2.9526956175080731E-5, -2.7597891482833575E-5, -2.5800617466688372E-5,
    -2.4130835676128019E-5, -2.2582350951834605E-5, -2.1147965676891298E-5,
    -1.9820063888529493E-5, -1.8590987080106508E-5, -1.7453269984421023E-5,
    -1.63997823854498E-5, -0.0004746177965599598, -0.0004778645671473215,
    -0.00032039022806703763, -0.00016110501611996228, -4.2577810128543523E-5,
    3.4457129429496748E-5, 7.97092684075675E-5, 0.00010313823670827221,
    0.00011246677526220416, 0.00011310364210848139, 0.00010865163484877427,
    0.00010143795159766197, 9.29298396593364E-5, 8.4029313301609E-5,
    7.52727991349134E-5, 6.696325219757309E-5, 5.925645473231947E-5,
    5.2216930882697554E-5, 4.5853948516536063E-5, 4.0144551389148682E-5,
    3.5048173003132809E-5, 3.0515799503434667E-5, 2.6495611995051603E-5,
    2.2936363369099816E-5, 1.9789305666402162E-5, 1.7009198463641262E-5,
    1.45547428261524E-5, 1.2388664099587841E-5, 1.0477587607658323E-5,
    8.7917995497847932E-6, 0.00073646581057257841, 0.000872790805146194,
    0.00062261486257313506, 0.00028599815419430417, 3.8473767287936606E-6,
    -0.00018790600363697156, -0.00029760364659455455, -0.00034599812683265633,
    -0.00035338247091603773, -0.00033571563577504876, -0.00030432112478903981,
    -0.00026672272304761283, -0.00022765421412281953, -0.00018992261185456235,
    -0.00015505891859909386, -0.00012377824076187363, -9.6292614771764412E-5,
    -7.2517832771442527E-5, -5.2207002889563382E-5, -3.5034775051190054E-5,
    -2.0648976103555174E-5, -8.7010609684976711E-6, 1.136986866751003E-6,
    9.1642647412277879E-6, 1.5647778542887261E-5, 2.0822362948246685E-5,
    2.4892338100459516E-5, 2.8034050957414632E-5, 3.0398777462986191E-5,
    3.2115673140670063E-5, -0.0018018219196388571, -0.0024340296293804253,
    -0.001834226635498568, -0.00076220459635400974, 0.00023907947525692722,
    0.00094926611717688109, 0.0013446744970154036, 0.0014845749525944918,
    0.0014473233983061759, 0.0013026826128565718, 0.0011035159737564268,
    0.00088604744041979172, 0.00067307320816566542, 0.00047760387285658237,
    0.00030599192635878935, 0.00016031569459472162, 4.0074955527061327E-5,
    -5.6660746163525162E-5, -0.00013250618677298264, -0.00019029618798961406,
    -0.00023281145037693741, -0.00026262881146466884, -0.00028205046986759866,
    -0.00029308156319286116, -0.00029743596217631662, -0.00029655733423934809,
    -0.00029164736331209088, -0.00028369620383773418, -0.00027351231709567335,
    -0.00026175015580676858, 0.0063858589121205088, 0.00962374215806378,
    0.0076187806120700105, 0.0028321905554562804, -0.0020984135201272008,
    -0.0057382676421662646, -0.0077080424449541465, -0.0082101169226484437,
    -0.0076582452034690543, -0.006472097293910452, -0.0049913241200496648,
    -0.0034561228971313326, -0.0020178558001417079, -0.00075943068678196145,
    0.00028417363152385912, 0.0011089166758633741, 0.0017290149387272878,
    0.0021681259080268472, 0.0024535771049453972, 0.0026128182105833488,
    0.0026714103965627691, 0.0026520307339598045, 0.0025741165287728731,
    0.0024538912623609443, 0.002304600580717955, 0.0021368483768671267,
    0.0019589652847887091, 0.0017773700867945441, 0.0015969028076583906,
    0.0014211197566443854 };

  int l2;
  int ias;
  int ibs;
  int is;
  double rtzta_re;
  double rtzta_im;
  int ks;
  boolean_T exitg2;
  static const double ALFA[180] = { -0.0044444444444444444,
    -0.000922077922077922, -8.8489288489288488E-5, 0.00016592768783244973,
    0.00024669137274179289, 0.00026599558934625478, 0.00026182429706150096,
    0.00024873043734465562, 0.00023272104008323209, 0.00021636248571236508,
    0.00020073885876275234, 0.00018626763663754517, 0.0001730607759178765,
    0.00016109170592901574, 0.00015027477416090814, 0.00014050349739126979,
    0.0001316688165459228, 0.00012366744559825325, 0.00011640527147473791,
    0.00010979829837271337, 0.00010377241042299283, 9.8262607836936344E-5,
    9.321205172495032E-5, 8.857108524787117E-5, 8.4296310571570029E-5,
    8.0349754840779115E-5, 7.6698134535920737E-5, 7.3312215748177779E-5,
    7.0166262516314139E-5, 6.7237563379016026E-5, 0.000693735541354589,
    0.00023224174518292166, -1.419862735566912E-5, -0.00011644493167204864,
    -0.00015080355805304876, -0.00015512192491809622, -0.00014680975664646556,
    -0.00013381550386749137, -0.00011974497568425405, -0.00010618431920797402,
    -9.3769954989119444E-5, -8.2692304558819327E-5, -7.2937434815522126E-5,
    -6.4404235772101633E-5, -5.69611566009369E-5, -5.0473104430356164E-5,
    -4.4813486800888282E-5, -3.9868872771759884E-5, -3.5540053297204251E-5,
    -3.174142566090225E-5, -2.839967939041748E-5, -2.5452272063487058E-5,
    -2.2845929716472455E-5, -2.0535275310648061E-5, -1.848162176276661E-5,
    -1.6651933002139381E-5, -1.5017941298011949E-5, -1.3555403137904052E-5,
    -1.2243474647385812E-5, -1.1064188481130817E-5, -0.00035421197145774384,
    -0.00015616126394515941, 3.0446550359493642E-5, 0.00013019865577324269,
    0.00016747110669971228, 0.00017022258768359256, 0.00015650142760859472,
    0.00013633917097744512, 0.00011488669202982512, 9.4586909303468817E-5,
    7.6449841925089825E-5, 6.0757033496519734E-5, 4.7439429929050881E-5,
    3.6275751200534429E-5, 2.6993971497922491E-5, 1.9321093824793926E-5,
    1.3005667479396321E-5, 7.8262086674449658E-6, 3.5925748581935159E-6,
    1.4404004981425182E-7, -2.6539676969793912E-6, -4.9134686709848593E-6,
    -6.7273929609124832E-6, -8.17269379678658E-6, -9.313047150935612E-6,
    -1.0201141879801643E-5, -1.0880596251059288E-5, -1.1387548150960355E-5,
    -1.1751967567455642E-5, -1.1998736487094414E-5, 0.00037819419920177291,
    0.00020247195276181616, -6.3793850631886236E-5, -0.0002385982306030059,
    -0.00031091625602736159, -0.00031368011524757634, -0.00027895027379132341,
    -0.00022856408261914138, -0.00017524528034084676, -0.00012554406306069035,
    -8.2298287282020835E-5, -4.6286073058811649E-5, -1.7233430236696227E-5,
    5.6069048230460226E-6, 2.313954431482868E-5, 3.6264274585679393E-5,
    4.5800612449018877E-5, 5.2459529495911405E-5, 5.6839620854581527E-5,
    5.9434982039310406E-5, 6.0647852757842175E-5, 6.0802390778843649E-5,
    6.0157789453946036E-5, 5.8919965734469847E-5, 5.72515823777593E-5,
    5.5280437558585257E-5, 5.3106377380288019E-5, 5.080693020123257E-5,
    4.8441864762009484E-5, 4.6056858160747536E-5, -0.00069114139728829421,
    -0.00042997663305887192, 0.000183067735980039, 0.00066008814754201417,
    0.0008759649699511859, 0.00087733523595823551, 0.00074936958537899067,
    0.000563832329756981, 0.00036805931997144317, 0.00018846453551445559,
    3.7066305766490415E-5, -8.28520220232137E-5, -0.000172751952869173,
    -0.00023631487360587297, -0.00027796615069490668, -0.00030207951415545694,
    -0.00031259471264382012, -0.00031287255875806717, -0.0003056780384663244,
    -0.00029322647061455731, -0.0002772556555829348, -0.00025910392846703172,
    -0.00023978401439648034, -0.00022004826004542284, -0.00020044391109497149,
    -0.00018135869221097068, -0.00016305767447865748, -0.00014571267217520584,
    -0.0001294254219839246, -0.00011424569194244596, 0.0019282196424877589,
    0.0013559257630202223, -0.000717858090421303, -0.0025808480257527035,
    -0.0034927113082616847, -0.0034698629934096061, -0.0028228523335131019,
    -0.0018810307640489134, -0.00088953171838394764, 3.8791210263103525E-6,
    0.00072868854011969139, 0.0012656637305345775, 0.0016251815837267443,
    0.0018320315321637317, 0.0019158838899052792, 0.0019058884675554615,
    0.0018279898242182574, 0.0017038950642112153, 0.0015509712717109768,
    0.0013826142185227616, 0.0012088142423006478, 0.0010367653263834496,
    0.00087143791806861914, 0.000716080155297701, 0.00057263700255812935,
    0.00044208981946580229, 0.00032472494850309055, 0.00022034204273024659,
    0.00012841289840135388, 4.8200592455209545E-5 };

  double tfn_re;
  double tfn_im;
  double rzth_re;
  double rzth_im;
  int kp1;
  int l;
  creal_T cr[14];
  creal_T dr[14];
  static const double C[105] = { 1.0, -0.20833333333333334, 0.125,
    0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
    1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
    -11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
    -28.212072558200244, 84.636217674600729, -91.818241543240021,
    42.534998745388457, -7.3687943594796321, 0.22710800170898438,
    212.57013003921713, -765.25246814118168, 1059.9904525279999,
    -699.57962737613252, 218.19051174421159, -26.491430486951554,
    0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
    -13586.550006434138, 11655.393336864534, -5305.646978613403,
    1200.9029132163525, -108.09091978839466, 1.7277275025844574,
    20204.291330966149, -96980.598388637518, 192547.00123253153,
    -203400.17728041555, 122200.46498301746, -41192.65496889755,
    7109.5143024893641, -493.915304773088, 6.074042001273483,
    -242919.18790055133, 1.3117636146629772E+6, -2.9980159185381066E+6,
    3.7632712976564039E+6, -2.8135632265865342E+6, 1.2683652733216248E+6,
    -331645.17248456361, 45218.768981362729, -2499.8304818112097,
    24.380529699556064, 3.2844698530720379E+6, -1.9706819118432228E+7,
    5.0952602492664643E+7, -7.4105148211532652E+7, 6.6344512274729028E+7,
    -3.7567176660763353E+7, 1.3288767166421818E+7, -2.7856181280864547E+6,
    308186.40461266239, -13886.08975371704, 110.01714026924674,
    -4.932925366450996E+7, 3.2557307418576574E+8, -9.394623596815784E+8,
    1.55359689957058E+9, -1.6210805521083372E+9, 1.1068428168230145E+9,
    -4.9588978427503031E+8, 1.4206290779753309E+8, -2.447406272573873E+7,
    2.2437681779224495E+6, -84005.433603024081, 551.33589612202059,
    8.1478909611831212E+8, -5.8664814920518475E+9, 1.8688207509295826E+10,
    -3.4632043388158775E+10, 4.1280185579753975E+10, -3.3026599749800724E+10,
    1.79542137311556E+10, -6.5632937926192846E+9, 1.5592798648792574E+9,
    -2.2510566188941526E+8, 1.7395107553978164E+7, -549842.32757228869,
    3038.0905109223841, -1.4679261247695616E+10, 1.144982377320258E+11,
    -3.9909617522446649E+11, 8.1921866954857727E+11, -1.0983751560812233E+12,
    1.0081581068653821E+12, -6.4536486924537646E+11, 2.8790064990615057E+11,
    -8.786707217802327E+10, 1.7634730606834969E+10, -2.1671649832237949E+9,
    1.4315787671888897E+8, -3.8718334425726128E+6, 18257.755474293175 };

  static const double BR[14] = { 1.0, -0.14583333333333334,
    -0.098741319444444448, -0.14331205391589505, -0.31722720267841353,
    -0.94242914795712029, -3.5112030408263544, -15.727263620368046,
    -82.281439097185938, -492.3553705236705, -3316.2185685479726,
    -24827.674245208589, -204526.58731512979, -1.83844491706821E+6 };

  static const double AR[14] = { 1.0, 0.10416666666666667, 0.083550347222222224,
    0.12822657455632716, 0.29184902646414046, 0.88162726744375763,
    3.3214082818627677, 14.995762986862555, 78.923013011586519,
    474.45153886826432, 3207.490090890662, 24086.549640874004, 198923.1191695098,
    1.7919020077753437E+6 };

  rfnu = 1.0 / fnu;
  memset(&up[0], 0, 14U * sizeof(creal_T));
  ac = fnu * 2.2250738585072014E-305;
  asum->re = 0.0;
  asum->im = 0.0;
  bsum->re = 0.0;
  bsum->im = 0.0;
  if ((std::abs(z.re) <= ac) && (std::abs(z.im) <= ac)) {
    zeta1->re = 1402.9773265065639 + fnu;
    zeta1->im = 0.0;
    zeta2->re = fnu;
    zeta2->im = 0.0;
    phi->re = 1.0;
    phi->im = 0.0;
    arg->re = 1.0;
    arg->im = 0.0;
  } else {
    zb.re = rfnu * z.re;
    zb.im = rfnu * z.im;
    rfnu2 = rfnu * rfnu;
    ac = rt_powd_snf(fnu, 0.33333333333333331);
    fn23 = ac * ac;
    rfn13_re = 1.0 / ac;
    w2.re = 1.0 - (zb.re * zb.re - zb.im * zb.im);
    w2.im = 0.0 - (zb.re * zb.im + zb.im * zb.re);
    ac = rt_hypotd_snf(w2.re, w2.im);
    if (ac > 0.25) {
      w = w2;
      c_sqrt(&w);
      ac = w.re;
      pp = w.im;
      if (w.re < 0.0) {
        ac = 0.0;
      }

      if (w.im < 0.0) {
        pp = 0.0;
      }

      w.re = ac;
      w.im = pp;
      if (zb.im == 0.0) {
        if (pp == 0.0) {
          zc.re = (1.0 + ac) / zb.re;
          zc.im = 0.0;
        } else if (1.0 + ac == 0.0) {
          zc.re = 0.0;
          zc.im = pp / zb.re;
        } else {
          zc.re = (1.0 + ac) / zb.re;
          zc.im = pp / zb.re;
        }
      } else if (zb.re == 0.0) {
        if (1.0 + ac == 0.0) {
          zc.re = pp / zb.im;
          zc.im = 0.0;
        } else if (pp == 0.0) {
          zc.re = 0.0;
          zc.im = -((1.0 + ac) / zb.im);
        } else {
          zc.re = pp / zb.im;
          zc.im = -((1.0 + ac) / zb.im);
        }
      } else {
        brm = std::abs(zb.re);
        bim = std::abs(zb.im);
        if (brm > bim) {
          bim = zb.im / zb.re;
          d = zb.re + bim * zb.im;
          zc.re = ((1.0 + ac) + bim * pp) / d;
          zc.im = (pp - bim * (1.0 + ac)) / d;
        } else if (bim == brm) {
          if (zb.re > 0.0) {
            bim = 0.5;
          } else {
            bim = -0.5;
          }

          if (zb.im > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }

          zc.re = ((1.0 + ac) * bim + pp * d) / brm;
          zc.im = (pp * bim - (1.0 + ac) * d) / brm;
        } else {
          bim = zb.re / zb.im;
          d = zb.im + bim * zb.re;
          zc.re = (bim * (1.0 + ac) + pp) / d;
          zc.im = (bim * pp - (1.0 + ac)) / d;
        }
      }

      c_log(&zc);
      ac = zc.re;
      pp = zc.im;
      if (zc.re < 0.0) {
        ac = 0.0;
      }

      if (zc.im < 0.0) {
        pp = 0.0;
      } else {
        if (zc.im > 1.5707963267948966) {
          pp = 1.5707963267948966;
        }
      }

      zth_re = 1.5 * (ac - w.re);
      zth_im = 1.5 * (pp - w.im);
      zeta1->re = fnu * ac;
      zeta1->im = fnu * pp;
      zeta2->re = fnu * w.re;
      zeta2->im = fnu * w.im;
      if ((zth_re >= 0.0) && (zth_im < 0.0)) {
        ac = 4.71238898038469;
      } else if (zth_re != 0.0) {
        ac = std::atan(zth_im / zth_re);
        if (zth_re < 0.0) {
          ac += 3.1415926535897931;
        }
      } else {
        ac = 1.5707963267948966;
      }

      pp = rt_powd_snf(rt_hypotd_snf(zth_re, zth_im), 0.66666666666666663);
      ac *= 0.66666666666666663;
      zb.re = pp * std::cos(ac);
      zb.im = pp * std::sin(ac);
      if (zb.im < 0.0) {
        zb.im = 0.0;
      }

      arg->re = fn23 * zb.re;
      arg->im = fn23 * zb.im;
      if (zb.im == 0.0) {
        if (zth_im == 0.0) {
          rtzta_re = zth_re / zb.re;
          rtzta_im = 0.0;
        } else if (zth_re == 0.0) {
          rtzta_re = 0.0;
          rtzta_im = zth_im / zb.re;
        } else {
          rtzta_re = zth_re / zb.re;
          rtzta_im = zth_im / zb.re;
        }
      } else if (zb.re == 0.0) {
        if (zth_re == 0.0) {
          rtzta_re = zth_im / zb.im;
          rtzta_im = 0.0;
        } else if (zth_im == 0.0) {
          rtzta_re = 0.0;
          rtzta_im = -(zth_re / zb.im);
        } else {
          rtzta_re = zth_im / zb.im;
          rtzta_im = -(zth_re / zb.im);
        }
      } else {
        brm = std::abs(zb.re);
        bim = zb.im;
        if (brm > bim) {
          bim = zb.im / zb.re;
          d = zb.re + bim * zb.im;
          rtzta_re = (zth_re + bim * zth_im) / d;
          rtzta_im = (zth_im - bim * zth_re) / d;
        } else if (bim == brm) {
          if (zb.re > 0.0) {
            bim = 0.5;
          } else {
            bim = -0.5;
          }

          if (zb.im > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }

          rtzta_re = (zth_re * bim + zth_im * d) / brm;
          rtzta_im = (zth_im * bim - zth_re * d) / brm;
        } else {
          bim = zb.re / zb.im;
          d = zb.im + bim * zb.re;
          rtzta_re = (bim * zth_re + zth_im) / d;
          rtzta_im = (bim * zth_im - zth_re) / d;
        }
      }

      if (w.im == 0.0) {
        if (rtzta_im == 0.0) {
          suma.re = rtzta_re / w.re;
          suma.im = 0.0;
        } else if (rtzta_re == 0.0) {
          suma.re = 0.0;
          suma.im = rtzta_im / w.re;
        } else {
          suma.re = rtzta_re / w.re;
          suma.im = rtzta_im / w.re;
        }
      } else if (w.re == 0.0) {
        if (rtzta_re == 0.0) {
          suma.re = rtzta_im / w.im;
          suma.im = 0.0;
        } else if (rtzta_im == 0.0) {
          suma.re = 0.0;
          suma.im = -(rtzta_re / w.im);
        } else {
          suma.re = rtzta_im / w.im;
          suma.im = -(rtzta_re / w.im);
        }
      } else {
        brm = std::abs(w.re);
        bim = std::abs(w.im);
        if (brm > bim) {
          bim = w.im / w.re;
          d = w.re + bim * w.im;
          suma.re = (rtzta_re + bim * rtzta_im) / d;
          suma.im = (rtzta_im - bim * rtzta_re) / d;
        } else if (bim == brm) {
          if (w.re > 0.0) {
            bim = 0.5;
          } else {
            bim = -0.5;
          }

          if (w.im > 0.0) {
            d = 0.5;
          } else {
            d = -0.5;
          }

          suma.re = (rtzta_re * bim + rtzta_im * d) / brm;
          suma.im = (rtzta_im * bim - rtzta_re * d) / brm;
        } else {
          bim = w.re / w.im;
          d = w.im + bim * w.re;
          suma.re = (bim * rtzta_re + rtzta_im) / d;
          suma.im = (bim * rtzta_im - rtzta_re) / d;
        }
      }

      zb.re = suma.re + suma.re;
      zb.im = suma.im + suma.im;
      c_sqrt(&zb);
      phi->re = zb.re * rfn13_re - zb.im * 0.0;
      phi->im = zb.re * 0.0 + zb.im * rfn13_re;
      if (ipmtr == 1) {
      } else {
        if (w.im == 0.0) {
          tfn_re = rfnu / w.re;
          tfn_im = 0.0;
        } else if (w.re == 0.0) {
          if (rfnu == 0.0) {
            tfn_re = 0.0 / w.im;
            tfn_im = 0.0;
          } else {
            tfn_re = 0.0;
            tfn_im = -(rfnu / w.im);
          }
        } else {
          brm = std::abs(w.re);
          bim = std::abs(w.im);
          if (brm > bim) {
            bim = w.im / w.re;
            d = w.re + bim * w.im;
            tfn_re = (rfnu + bim * 0.0) / d;
            tfn_im = (0.0 - bim * rfnu) / d;
          } else if (bim == brm) {
            if (w.re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (w.im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            tfn_re = (rfnu * bim + 0.0 * d) / brm;
            tfn_im = (0.0 * bim - rfnu * d) / brm;
          } else {
            bim = w.re / w.im;
            d = w.im + bim * w.re;
            tfn_re = bim * rfnu / d;
            tfn_im = (bim * 0.0 - rfnu) / d;
          }
        }

        if (zth_im == 0.0) {
          rzth_re = rfnu / zth_re;
          rzth_im = 0.0;
        } else if (zth_re == 0.0) {
          if (rfnu == 0.0) {
            rzth_re = 0.0 / zth_im;
            rzth_im = 0.0;
          } else {
            rzth_re = 0.0;
            rzth_im = -(rfnu / zth_im);
          }
        } else {
          brm = std::abs(zth_re);
          bim = std::abs(zth_im);
          if (brm > bim) {
            bim = zth_im / zth_re;
            d = zth_re + bim * zth_im;
            rzth_re = (rfnu + bim * 0.0) / d;
            rzth_im = (0.0 - bim * rfnu) / d;
          } else if (bim == brm) {
            if (zth_re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (zth_im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            rzth_re = (rfnu * bim + 0.0 * d) / brm;
            rzth_im = (0.0 * bim - rfnu * d) / brm;
          } else {
            bim = zth_re / zth_im;
            d = zth_im + bim * zth_re;
            rzth_re = bim * rfnu / d;
            rzth_im = (bim * 0.0 - rfnu) / d;
          }
        }

        zc.re = 0.10416666666666667 * rzth_re;
        zc.im = 0.10416666666666667 * rzth_im;
        if (w2.im == 0.0) {
          w.re = 1.0 / w2.re;
          w.im = 0.0;
        } else if (w2.re == 0.0) {
          w.re = 0.0;
          w.im = -(1.0 / w2.im);
        } else {
          brm = std::abs(w2.re);
          bim = std::abs(w2.im);
          if (brm > bim) {
            bim = w2.im / w2.re;
            d = w2.re + bim * w2.im;
            w.re = (1.0 + bim * 0.0) / d;
            w.im = (0.0 - bim) / d;
          } else if (bim == brm) {
            if (w2.re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (w2.im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            w.re = (bim + 0.0 * d) / brm;
            w.im = (0.0 * bim - d) / brm;
          } else {
            bim = w2.re / w2.im;
            d = w2.im + bim * w2.re;
            w.re = bim / d;
            w.im = (bim * 0.0 - 1.0) / d;
          }
        }

        ac = w.re * -0.20833333333333334 + 0.125;
        pp = w.im * -0.20833333333333334;
        up[1].re = ac * tfn_re - pp * tfn_im;
        up[1].im = ac * tfn_im + pp * tfn_re;
        bsum->re = up[1].re + zc.re;
        bsum->im = up[1].im + zc.im;
        if (rfnu >= tol) {
          zth_re = rzth_re;
          zth_im = rzth_im;
          w2.re = tfn_re;
          w2.im = tfn_im;
          up[0].re = 1.0;
          up[0].im = 0.0;
          pp = 1.0;
          fn23 = tol * (std::abs(bsum->re) + std::abs(bsum->im));
          ks = 0;
          kp1 = 2;
          l = 2;
          ias = 0;
          ibs = 0;
          for (l1 = 0; l1 < 14; l1++) {
            cr[l1].re = 0.0;
            cr[l1].im = 0.0;
            dr[l1].re = 0.0;
            dr[l1].im = 0.0;
          }

          is = 2;
          exitg1 = false;
          while ((!exitg1) && (is < 13)) {
            for (k = is; k <= is + 1; k++) {
              ks++;
              kp1++;
              l++;
              suma.re = C[l];
              suma.im = 0.0;
              for (l1 = 2; l1 <= kp1; l1++) {
                l++;
                ac = suma.re;
                suma.re = (suma.re * w.re - suma.im * w.im) + C[l];
                suma.im = ac * w.im + suma.im * w.re;
              }

              ac = w2.re;
              w2.re = w2.re * tfn_re - w2.im * tfn_im;
              w2.im = ac * tfn_im + w2.im * tfn_re;
              up[kp1 - 1].re = w2.re * suma.re - w2.im * suma.im;
              up[kp1 - 1].im = w2.re * suma.im + w2.im * suma.re;
              if (ks > 2147483646) {
                l1 = MAX_int32_T;
              } else {
                l1 = ks + 1;
              }

              cr[ks - 1].re = BR[l1 - 1] * zth_re;
              if (ks > 2147483646) {
                l1 = MAX_int32_T;
              } else {
                l1 = ks + 1;
              }

              cr[ks - 1].im = BR[l1 - 1] * zth_im;
              ac = zth_re;
              zth_re = zth_re * rzth_re - zth_im * rzth_im;
              zth_im = ac * rzth_im + zth_im * rzth_re;
              if (ks > 2147483645) {
                l1 = MAX_int32_T;
              } else {
                l1 = ks + 2;
              }

              dr[ks - 1].re = AR[l1 - 1] * zth_re;
              if (ks > 2147483645) {
                l1 = MAX_int32_T;
              } else {
                l1 = ks + 2;
              }

              dr[ks - 1].im = AR[l1 - 1] * zth_im;
            }

            pp *= rfnu2;
            if (ias != 1) {
              suma = up[is];
              l1 = is;
              for (l2 = 1; l2 <= is; l2++) {
                l1--;
                suma.re += cr[l2 - 1].re * up[l1].re - cr[l2 - 1].im * up[l1].im;
                suma.im += cr[l2 - 1].re * up[l1].im + cr[l2 - 1].im * up[l1].re;
              }

              asum->re += suma.re;
              asum->im += suma.im;
              if ((pp < tol) && (std::abs(asum->re) + std::abs(asum->im) < tol))
              {
                ias = 1;
              }
            }

            if (ibs != 1) {
              zb.re = up[is + 1].re + (up[is].re * zc.re - up[is].im * zc.im);
              zb.im = up[is + 1].im + (up[is].re * zc.im + up[is].im * zc.re);
              l1 = is;
              for (l2 = 1; l2 <= is; l2++) {
                l1--;
                zb.re += dr[l2 - 1].re * up[l1].re - dr[l2 - 1].im * up[l1].im;
                zb.im += dr[l2 - 1].re * up[l1].im + dr[l2 - 1].im * up[l1].re;
              }

              bsum->re += zb.re;
              bsum->im += zb.im;
              if ((pp < fn23) && (std::abs(bsum->re) + std::abs(bsum->im) < tol))
              {
                ibs = 1;
              }
            }

            if ((ias == 1) && (ibs == 1)) {
              exitg1 = true;
            } else {
              is += 2;
            }
          }
        }

        asum->re++;
        bsum->re = -bsum->re;
        bsum->im = -bsum->im;
        ac = bsum->re;
        pp = bsum->im;
        bsum->re = bsum->re * rfn13_re - bsum->im * 0.0;
        bsum->im = ac * 0.0 + bsum->im * rfn13_re;
        fn23 = ac * rfn13_re - pp * 0.0;
        pp = ac * 0.0 + pp * rfn13_re;
        if (rtzta_im == 0.0) {
          if (bsum->im == 0.0) {
            bsum->re /= rtzta_re;
            bsum->im = 0.0;
          } else if (bsum->re == 0.0) {
            bsum->re = 0.0;
            bsum->im /= rtzta_re;
          } else {
            bsum->re /= rtzta_re;
            bsum->im /= rtzta_re;
          }
        } else if (rtzta_re == 0.0) {
          if (bsum->re == 0.0) {
            bsum->re = bsum->im / rtzta_im;
            bsum->im = 0.0;
          } else if (bsum->im == 0.0) {
            bsum->re = 0.0;
            bsum->im = -(fn23 / rtzta_im);
          } else {
            bsum->re = bsum->im / rtzta_im;
            bsum->im = -(fn23 / rtzta_im);
          }
        } else {
          brm = std::abs(rtzta_re);
          bim = std::abs(rtzta_im);
          if (brm > bim) {
            bim = rtzta_im / rtzta_re;
            d = rtzta_re + bim * rtzta_im;
            bsum->re = (fn23 + bim * pp) / d;
            bsum->im = (pp - bim * fn23) / d;
          } else if (bim == brm) {
            if (rtzta_re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (rtzta_im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            bsum->re = (fn23 * bim + pp * d) / brm;
            bsum->im = (pp * bim - fn23 * d) / brm;
          } else {
            bim = rtzta_re / rtzta_im;
            d = rtzta_im + bim * rtzta_re;
            bsum->re = (bim * fn23 + pp) / d;
            bsum->im = (bim * pp - fn23) / d;
          }
        }
      }
    } else {
      k = 1;
      memset(&p[0], 0, 30U * sizeof(creal_T));
      memset(&ap[0], 0, 30U * sizeof(double));
      p[0].re = 1.0;
      p[0].im = 0.0;
      suma.re = 0.6299605249474366;
      suma.im = 0.0;
      ap[0] = 1.0;
      if (ac >= tol) {
        k = 30;
        l1 = 1;
        exitg1 = false;
        while ((!exitg1) && (l1 + 1 < 31)) {
          pp = p[l1 - 1].re;
          p[l1].re = p[l1 - 1].re * w2.re - p[l1 - 1].im * w2.im;
          p[l1].im = pp * w2.im + p[l1 - 1].im * w2.re;
          suma.re += p[l1].re * dv0[l1];
          suma.im += p[l1].im * dv0[l1];
          ap[l1] = ap[l1 - 1] * ac;
          if (ap[l1] < tol) {
            k = l1 + 1;
            exitg1 = true;
          } else {
            l1++;
          }
        }
      }

      zb.re = w2.re * suma.re - w2.im * suma.im;
      zb.im = w2.re * suma.im + w2.im * suma.re;
      arg->re = fn23 * zb.re;
      arg->im = fn23 * zb.im;
      c_sqrt(&suma);
      c_sqrt(&w2);
      zeta2->re = fnu * w2.re;
      zeta2->im = fnu * w2.im;
      ac = (zb.re * suma.re - zb.im * suma.im) * 0.66666666666666663 + 1.0;
      pp = (zb.re * suma.im + zb.im * suma.re) * 0.66666666666666663;
      zeta1->re = zeta2->re * ac - zeta2->im * pp;
      zeta1->im = zeta2->re * pp + zeta2->im * ac;
      suma.re += suma.re;
      suma.im += suma.im;
      c_sqrt(&suma);
      phi->re = suma.re * rfn13_re - suma.im * 0.0;
      phi->im = suma.re * 0.0 + suma.im * rfn13_re;
      if (ipmtr == 1) {
      } else {
        zb.re = 0.0;
        zb.im = 0.0;
        for (l1 = 0; l1 < k; l1++) {
          zb.re += p[l1].re * BETA[l1];
          zb.im += p[l1].im * BETA[l1];
        }

        *bsum = zb;
        l1 = 0;
        l2 = 30;
        fn23 = tol * rt_hypotd_snf(zb.re, zb.im);
        ac = tol;
        pp = 1.0;
        ias = 0;
        ibs = 0;
        if (rfnu2 >= tol) {
          is = 2;
          exitg1 = false;
          while ((!exitg1) && (is < 8)) {
            ac /= rfnu2;
            pp *= rfnu2;
            if (ias != 1) {
              suma.re = 0.0;
              suma.im = 0.0;
              ks = 0;
              exitg2 = false;
              while ((!exitg2) && (ks + 1 <= k)) {
                suma.re += p[ks].re * ALFA[l1 + ks];
                suma.im += p[ks].im * ALFA[l1 + ks];
                if (ap[ks] < ac) {
                  exitg2 = true;
                } else {
                  ks++;
                }
              }

              asum->re += suma.re * pp;
              asum->im += suma.im * pp;
              if (pp < tol) {
                ias = 1;
              }
            }

            if (ibs != 1) {
              zb.re = 0.0;
              zb.im = 0.0;
              ks = 0;
              exitg2 = false;
              while ((!exitg2) && (ks + 1 <= k)) {
                zb.re += p[ks].re * BETA[l2 + ks];
                zb.im += p[ks].im * BETA[l2 + ks];
                if (ap[ks] < ac) {
                  exitg2 = true;
                } else {
                  ks++;
                }
              }

              bsum->re += zb.re * pp;
              bsum->im += zb.im * pp;
              if (pp < fn23) {
                ibs = 1;
              }
            }

            if ((ias == 1) && (ibs == 1)) {
              exitg1 = true;
            } else {
              l1 += 30;
              l2 += 30;
              is++;
            }
          }
        }

        asum->re++;
        bsum->re *= rfnu * rfn13_re;
        bsum->im *= rfnu * rfn13_re;
      }
    }
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T y[2]
//                double fnul
//                double tol
//                double elim
//                double alim
//                int *nlast
//                int *nz
// Return Type  : void
//
static void b_cuni2(const creal_T z, double fnu, int kode, int nin, creal_T y[2],
                    double fnul, double tol, double elim, double alim, int
                    *nlast, int *nz)
{
  int n;
  double cssr[3];
  double csrr[3];
  double bry1;
  double yy;
  creal_T zn;
  double zb_re;
  double zb_im;
  signed char cid_im;
  double ffnu;
  int inu;
  double c2_re;
  double c2_im;
  double zar_re;
  double zar_im;
  int in;
  static const cint8_T icv1[4] = { { 1,// re
      0                                // im
    }, { 0,                            // re
      1                                // im
    }, { -1,                           // re
      0                                // im
    }, { 0,                            // re
      -1                               // im
    } };

  creal_T dai;
  creal_T s1;
  creal_T zeta1;
  creal_T zeta2;
  creal_T ai;
  creal_T unusedU7;
  double br;
  double bi;
  double rs1;
  double brm;
  int exitg1;
  boolean_T goto_mw120;
  int i;
  double b_br;
  int nn;
  double b_bi;
  boolean_T exitg2;
  double fn;
  creal_T phi;
  creal_T asum;
  creal_T bsum;
  boolean_T guard1 = false;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  *nz = 0;
  *nlast = 0;
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry1 = 2.2250738585072014E-305 / tol;
  yy = z.im;
  zn.re = z.im;
  zn.im = -z.re;
  zb_re = z.re;
  zb_im = z.im;
  cid_im = -1;
  if (fnu < 0.0) {
    ffnu = std::ceil(fnu);
  } else {
    ffnu = std::floor(fnu);
  }

  inu = (int)ffnu - 1;
  ffnu = 1.5707963267948966 * (fnu - ffnu);
  c2_re = std::cos(ffnu);
  c2_im = std::sin(ffnu);
  zar_re = c2_re;
  zar_im = c2_im;
  in = inu + n;
  in -= (in >> 2) << 2;
  in++;
  ffnu = c2_re;
  c2_re = c2_re * (double)icv1[in - 1].re - c2_im * (double)icv1[in - 1].im;
  c2_im = ffnu * (double)icv1[in - 1].im + c2_im * (double)icv1[in - 1].re;
  if (z.im <= 0.0) {
    zn.re = -z.im;
    zn.im = -z.re;
    zb_re = z.re;
    zb_im = -z.im;
    cid_im = 1;
    c2_im = -c2_im;
  }

  if (fnu > 1.0) {
    ffnu = fnu;
  } else {
    ffnu = 1.0;
  }

  b_cunhj(zn, ffnu, 1, tol, &dai, &s1, &zeta1, &zeta2, &ai, &unusedU7);
  if (kode == 1) {
    s1.re = zeta2.re - zeta1.re;
  } else {
    br = zb_re + zeta2.re;
    bi = zb_im + zeta2.im;
    if (bi == 0.0) {
      ffnu = fnu / br;
    } else if (br == 0.0) {
      if (fnu == 0.0) {
        ffnu = 0.0 / bi;
      } else {
        ffnu = 0.0;
      }
    } else {
      brm = std::abs(br);
      ffnu = std::abs(bi);
      if (brm > ffnu) {
        brm = bi / br;
        ffnu = (fnu + brm * 0.0) / (br + brm * bi);
      } else if (ffnu == brm) {
        if (br > 0.0) {
          b_br = 0.5;
        } else {
          b_br = -0.5;
        }

        if (bi > 0.0) {
          b_bi = 0.5;
        } else {
          b_bi = -0.5;
        }

        ffnu = (fnu * b_br + 0.0 * b_bi) / brm;
      } else {
        brm = br / bi;
        ffnu = brm * fnu / (bi + brm * br);
      }
    }

    s1.re = fnu * ffnu - zeta1.re;
  }

  rs1 = s1.re;
  if (std::abs(s1.re) > elim) {
    if (s1.re > 0.0) {
      *nz = -1;
    } else {
      *nz = n;
      for (i = 1; i <= n; i++) {
        y[i - 1].re = 0.0;
        y[i - 1].im = 0.0;
      }
    }
  } else {
    do {
      exitg1 = 0;
      goto_mw120 = false;
      in = -1;
      if (2 < n) {
        nn = 2;
      } else {
        nn = n;
      }

      i = 1;
      exitg2 = false;
      while ((!exitg2) && (i <= nn)) {
        fn = fnu + (double)(n - i);
        b_cunhj(zn, fn, 0, tol, &phi, &unusedU7, &zeta1, &zeta2, &asum, &bsum);
        if (kode == 1) {
          s1.re = zeta2.re - zeta1.re;
          s1.im = zeta2.im - zeta1.im;
        } else {
          br = zb_re + zeta2.re;
          bi = zb_im + zeta2.im;
          if (bi == 0.0) {
            rs1 = fn / br;
            ffnu = 0.0;
          } else if (br == 0.0) {
            if (fn == 0.0) {
              rs1 = 0.0 / bi;
              ffnu = 0.0;
            } else {
              rs1 = 0.0;
              ffnu = -(fn / bi);
            }
          } else {
            brm = std::abs(br);
            ffnu = std::abs(bi);
            if (brm > ffnu) {
              brm = bi / br;
              ffnu = br + brm * bi;
              rs1 = (fn + brm * 0.0) / ffnu;
              ffnu = (0.0 - brm * fn) / ffnu;
            } else if (ffnu == brm) {
              if (br > 0.0) {
                br = 0.5;
              } else {
                br = -0.5;
              }

              if (bi > 0.0) {
                ffnu = 0.5;
              } else {
                ffnu = -0.5;
              }

              rs1 = (fn * br + 0.0 * ffnu) / brm;
              ffnu = (0.0 * br - fn * ffnu) / brm;
            } else {
              brm = br / bi;
              ffnu = bi + brm * br;
              rs1 = brm * fn / ffnu;
              ffnu = (brm * 0.0 - fn) / ffnu;
            }
          }

          s1.re = (fn * rs1 - zeta1.re) + std::abs(yy) * 0.0;
          s1.im = (fn * ffnu - zeta1.im) + std::abs(yy);
        }

        rs1 = s1.re;
        if (std::abs(s1.re) > elim) {
          goto_mw120 = true;
          exitg2 = true;
        } else {
          if (i == 1) {
            in = 1;
          }

          guard1 = false;
          if (std::abs(s1.re) >= alim) {
            rs1 = ((s1.re + std::log(rt_hypotd_snf(phi.re, phi.im))) - 0.25 *
                   std::log(rt_hypotd_snf(unusedU7.re, unusedU7.im))) -
              1.2655121234846454;
            if (std::abs(rs1) > elim) {
              goto_mw120 = true;
              exitg2 = true;
            } else {
              if (i == 1) {
                in = 0;
              }

              if ((rs1 >= 0.0) && (i == 1)) {
                in = 2;
              }

              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            ai = cairy(unusedU7, 0, 2);
            dai = cairy(unusedU7, 1, 2);
            ffnu = (ai.re * asum.re - ai.im * asum.im) + (dai.re * bsum.re -
              dai.im * bsum.im);
            brm = (ai.re * asum.im + ai.im * asum.re) + (dai.re * bsum.im +
              dai.im * bsum.re);
            bi = ffnu * phi.re - brm * phi.im;
            brm = ffnu * phi.im + brm * phi.re;
            ffnu = std::exp(s1.re) * cssr[in] * std::cos(s1.im);
            br = std::exp(s1.re) * cssr[in] * std::sin(s1.im);
            dai.re = bi * ffnu - brm * br;
            dai.im = bi * br + brm * ffnu;
            if ((in + 1 == 1) && (cuchk(dai, bry1, tol) != 0)) {
              goto_mw120 = true;
              exitg2 = true;
            } else {
              if (yy <= 0.0) {
                dai.im = -dai.im;
              }

              y[n - i].re = csrr[in] * (dai.re * c2_re - dai.im * c2_im);
              y[n - i].im = csrr[in] * (dai.re * c2_im + dai.im * c2_re);
              ffnu = c2_re;
              c2_re = c2_re * 0.0 - c2_im * (double)cid_im;
              c2_im = ffnu * (double)cid_im + c2_im * 0.0;
              i++;
            }
          }
        }
      }

      if (!goto_mw120) {
        exitg1 = 1;
      } else if (rs1 > 0.0) {
        *nz = -1;
        exitg1 = 1;
      } else {
        y[n - 1].re = 0.0;
        y[n - 1].im = 0.0;
        (*nz)++;
        n--;
        if (n == 0) {
          exitg1 = 1;
        } else {
          in = b_cuoik(z, fnu, kode, 1, n, y, tol, elim, alim);
          if (in < 0) {
            *nz = -1;
            exitg1 = 1;
          } else {
            n -= in;
            *nz += in;
            if (n == 0) {
              exitg1 = 1;
            } else if ((fnu + (double)n) - 1.0 < fnul) {
              *nlast = n;
              exitg1 = 1;
            } else {
              in = inu + n;
              nn = icv1[in - ((in >> 2) << 2)].re;
              in = icv1[in - ((in >> 2) << 2)].im;
              c2_re = zar_re * (double)nn - zar_im * (double)in;
              c2_im = zar_re * (double)in + zar_im * (double)nn;
              if (yy <= 0.0) {
                c2_im = -c2_im;
              }
            }
          }
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : const creal_T zr
//                double fnu
//                int ikflg
//                int ipmtr
//                double tol
//                int *init
//                creal_T cwrk[16]
//                creal_T *phi
//                creal_T *zeta1
//                creal_T *zeta2
//                creal_T *summ
// Return Type  : void
//
static void b_cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double
                    tol, int *init, creal_T cwrk[16], creal_T *phi, creal_T
                    *zeta1, creal_T *zeta2, creal_T *summ)
{
  double rfn;
  double ac;
  boolean_T guard1 = false;
  creal_T t;
  double s_re;
  double s_im;
  creal_T sr;
  int i;
  double t_re;
  double cfn_im;
  double crfn_im;
  double crfn_re;
  double cfn_re;
  int l;
  int exitg1;
  int j;
  rfn = 1.0 / fnu;
  ac = fnu * 2.2250738585072014E-305;
  guard1 = false;
  if (*init == 0) {
    summ->re = 0.0;
    summ->im = 0.0;
    if ((std::abs(zr.re) > ac) || (std::abs(zr.im) > ac)) {
      t.re = zr.re * rfn - zr.im * 0.0;
      t.im = zr.re * 0.0 + zr.im * rfn;
      s_re = (t.re * t.re - t.im * t.im) + 1.0;
      s_im = t.re * t.im + t.im * t.re;
      sr.re = s_re;
      sr.im = s_im;
      c_sqrt(&sr);
      t_re = t.re;
      ac = 1.0 + sr.re;
      cfn_im = sr.im;
      if (t.im == 0.0) {
        if (cfn_im == 0.0) {
          t.re = ac / t.re;
          t.im = 0.0;
        } else if (ac == 0.0) {
          t.re = 0.0;
          t.im = cfn_im / t_re;
        } else {
          t.re = ac / t.re;
          t.im = cfn_im / t_re;
        }
      } else if (t.re == 0.0) {
        if (ac == 0.0) {
          t.re = cfn_im / t.im;
          t.im = 0.0;
        } else if (cfn_im == 0.0) {
          t.re = 0.0;
          t.im = -(ac / t.im);
        } else {
          t.re = cfn_im / t.im;
          t.im = -(ac / t.im);
        }
      } else {
        crfn_im = std::abs(t.re);
        t_re = std::abs(t.im);
        if (crfn_im > t_re) {
          t_re = t.im / t.re;
          crfn_re = t.re + t_re * t.im;
          t.re = (ac + t_re * cfn_im) / crfn_re;
          t.im = (cfn_im - t_re * ac) / crfn_re;
        } else if (t_re == crfn_im) {
          if (t.re > 0.0) {
            t_re = 0.5;
          } else {
            t_re = -0.5;
          }

          if (t.im > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          t.re = (ac * t_re + cfn_im * crfn_re) / crfn_im;
          t.im = (cfn_im * t_re - ac * crfn_re) / crfn_im;
        } else {
          t_re = t.re / t.im;
          crfn_re = t.im + t_re * t.re;
          t.re = (t_re * ac + cfn_im) / crfn_re;
          t.im = (t_re * cfn_im - ac) / crfn_re;
        }
      }

      c_log(&t);
      zeta1->re = fnu * t.re - 0.0 * t.im;
      zeta1->im = fnu * t.im + 0.0 * t.re;
      zeta2->re = fnu * sr.re - 0.0 * sr.im;
      zeta2->im = fnu * sr.im + 0.0 * sr.re;
      if (sr.im == 0.0) {
        sr.re = rfn / sr.re;
        sr.im = 0.0;
      } else if (sr.re == 0.0) {
        if (rfn == 0.0) {
          sr.re = 0.0 / sr.im;
          sr.im = 0.0;
        } else {
          sr.re = 0.0;
          sr.im = -(rfn / sr.im);
        }
      } else {
        crfn_im = std::abs(sr.re);
        t_re = std::abs(sr.im);
        if (crfn_im > t_re) {
          t_re = sr.im / sr.re;
          crfn_re = sr.re + t_re * sr.im;
          sr.re = (rfn + t_re * 0.0) / crfn_re;
          sr.im = (0.0 - t_re * rfn) / crfn_re;
        } else if (t_re == crfn_im) {
          if (sr.re > 0.0) {
            t_re = 0.5;
          } else {
            t_re = -0.5;
          }

          if (sr.im > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          sr.re = (rfn * t_re + 0.0 * crfn_re) / crfn_im;
          sr.im = (0.0 * t_re - rfn * crfn_re) / crfn_im;
        } else {
          t_re = sr.re / sr.im;
          crfn_re = sr.im + t_re * sr.re;
          sr.re = t_re * rfn / crfn_re;
          sr.im = (t_re * 0.0 - rfn) / crfn_re;
        }
      }

      cwrk[15] = sr;
      c_sqrt(&cwrk[15]);
      phi->re = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].re;
      phi->im = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].im;
      if (ipmtr != 0) {
      } else {
        if (s_im == 0.0) {
          cfn_re = 1.0 / s_re;
          cfn_im = 0.0;
        } else if (s_re == 0.0) {
          cfn_re = 0.0;
          cfn_im = -(1.0 / s_im);
        } else {
          crfn_im = std::abs(s_re);
          t_re = std::abs(s_im);
          if (crfn_im > t_re) {
            t_re = s_im / s_re;
            crfn_re = s_re + t_re * s_im;
            cfn_re = (1.0 + t_re * 0.0) / crfn_re;
            cfn_im = (0.0 - t_re) / crfn_re;
          } else if (t_re == crfn_im) {
            if (s_re > 0.0) {
              t_re = 0.5;
            } else {
              t_re = -0.5;
            }

            if (s_im > 0.0) {
              crfn_re = 0.5;
            } else {
              crfn_re = -0.5;
            }

            cfn_re = (t_re + 0.0 * crfn_re) / crfn_im;
            cfn_im = (0.0 * t_re - crfn_re) / crfn_im;
          } else {
            t_re = s_re / s_im;
            crfn_re = s_im + t_re * s_re;
            cfn_re = t_re / crfn_re;
            cfn_im = (t_re * 0.0 - 1.0) / crfn_re;
          }
        }

        cwrk[0].re = 1.0;
        cwrk[0].im = 0.0;
        crfn_re = 1.0;
        crfn_im = 0.0;
        ac = 1.0;
        l = 0;
        *init = 15;
        i = 1;
        do {
          exitg1 = 0;
          if (i < 15) {
            s_re = 0.0;
            s_im = 0.0;
            for (j = 0; j <= i; j++) {
              l++;
              t_re = s_re;
              s_re = s_re * cfn_re - s_im * cfn_im;
              s_im = t_re * cfn_im + s_im * cfn_re;
              s_re += dv1[l];
            }

            t_re = crfn_re;
            crfn_re = crfn_re * sr.re - crfn_im * sr.im;
            crfn_im = t_re * sr.im + crfn_im * sr.re;
            cwrk[i].re = crfn_re * s_re - crfn_im * s_im;
            cwrk[i].im = crfn_re * s_im + crfn_im * s_re;
            ac *= rfn;
            t = cwrk[i];
            if ((ac < tol) && (std::abs(t.re) + std::abs(t.im) < tol)) {
              *init = i + 1;
              guard1 = true;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            guard1 = true;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }
    } else {
      zeta1->re = 1402.9773265065639 + fnu;
      zeta1->im = 0.0;
      zeta2->re = fnu;
      zeta2->im = 0.0;
      phi->re = 1.0;
      phi->im = 0.0;
    }
  } else {
    zeta1->re = 0.0;
    zeta1->im = 0.0;
    zeta2->re = 0.0;
    zeta2->im = 0.0;
    guard1 = true;
  }

  if (guard1) {
    if (ikflg == 2) {
      s_re = 0.0;
      s_im = 0.0;
      t.re = 1.0;
      t.im = 0.0;
      for (i = 1; i <= *init; i++) {
        s_re += t.re * cwrk[i - 1].re - t.im * cwrk[i - 1].im;
        s_im += t.re * cwrk[i - 1].im + t.im * cwrk[i - 1].re;
        t.re = -t.re;
        t.im = -t.im;
      }

      summ->re = s_re;
      summ->im = s_im;
      phi->re = 1.2533141373155003 * cwrk[15].re;
      phi->im = 1.2533141373155003 * cwrk[15].im;
    } else {
      s_re = 0.0;
      s_im = 0.0;
      for (i = 1; i <= *init; i++) {
        s_re += cwrk[i - 1].re;
        s_im += cwrk[i - 1].im;
      }

      summ->re = s_re;
      summ->im = s_im;
      phi->re = 0.3989422804014327 * cwrk[15].re;
      phi->im = 0.3989422804014327 * cwrk[15].im;
    }
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int ikflg
//                int nin
//                creal_T y[2]
//                double tol
//                double elim
//                double alim
// Return Type  : int
//
static int b_cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin,
                   creal_T y[2], double tol, double elim, double alim)
{
  int nuf;
  int n;
  creal_T zr;
  int iform;
  double gnn;
  double gnu;
  double aarg;
  creal_T arg;
  creal_T an;
  creal_T phi;
  creal_T cz;
  creal_T zeta2;
  boolean_T guard1 = false;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  nuf = 0;
  if (z.re < 0.0) {
    zr.re = -z.re;
    zr.im = -z.im;
  } else {
    zr = z;
  }

  if (std::abs(zr.im) > std::abs(z.re) * 1.7321) {
    iform = 2;
  } else {
    iform = 1;
  }

  if (ikflg == 1) {
    if (fnu < 1.0) {
      gnu = 1.0;
    } else {
      gnu = fnu;
    }
  } else {
    gnn = (fnu + (double)n) - 1.0;
    gnu = n;
    if (gnn > gnu) {
      gnu = gnn;
    }
  }

  aarg = 0.0;
  if (iform == 2) {
    if (zr.im <= 0.0) {
      an.re = -zr.im;
      an.im = -zr.re;
    } else {
      an.re = zr.im;
      an.im = -zr.re;
    }

    cunhj(an, gnu, tol, &phi, &arg, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
    aarg = rt_hypotd_snf(arg.re, arg.im);
  } else {
    arg.re = 0.0;
    arg.im = 0.0;
    cunik(zr, gnu, ikflg, 1, tol, 0, &phi, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
  }

  if (kode == 2) {
    cz.re -= zr.re;
    cz.im -= zr.im;
  }

  if (ikflg == 2) {
    cz.re = -cz.re;
    cz.im = -cz.im;
  }

  gnn = rt_hypotd_snf(phi.re, phi.im);
  if (cz.re >= alim) {
    gnn = cz.re + std::log(gnn);
    if (iform == 2) {
      gnn = (gnn - 0.25 * std::log(aarg)) - 1.2655121234846454;
    }

    if (gnn > elim) {
      nuf = -1;
    }
  } else {
    guard1 = false;
    if (cz.re >= -elim) {
      if (cz.re > -alim) {
      } else {
        gnn = cz.re + std::log(gnn);
        if (iform == 2) {
          gnn = (gnn - 0.25 * std::log(aarg)) - 1.2655121234846454;
        }

        if (gnn > -elim) {
          c_log(&phi);
          cz.im += phi.im;
          if (iform == 2) {
            c_log(&arg);
            cz.im -= 0.25 * arg.im;
          }

          gnn = std::exp(gnn) / tol;
          zr.re = gnn * std::cos(cz.im);
          zr.im = gnn * std::sin(cz.im);
          if (cuchk(zr, 2.2250738585072014E-305 / tol, tol) != 1) {
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      for (iform = 1; iform <= n; iform++) {
        y[iform - 1].re = 0.0;
        y[iform - 1].im = 0.0;
      }

      nuf = n;
    }
  }

  return nuf;
}

//
// Arguments    : creal_T *x
// Return Type  : void
//
static void b_exp(creal_T *x)
{
  double r;
  double x_im;
  if (x->im == 0.0) {
    x->re = std::exp(x->re);
    x->im = 0.0;
  } else if (rtIsInf(x->im) && rtIsInf(x->re) && (x->re < 0.0)) {
    x->re = 0.0;
    x->im = 0.0;
  } else {
    r = std::exp(x->re / 2.0);
    x_im = x->im;
    x->re = r * (r * std::cos(x->im));
    x->im = r * (r * std::sin(x_im));
  }
}

//
// Arguments    : double *x
// Return Type  : void
//
static void b_fix(double *x)
{
  if (*x < 0.0) {
    *x = std::ceil(*x);
  } else {
    *x = std::floor(*x);
  }
}

//
// Arguments    : double *x
// Return Type  : void
//
static void b_log(double *x)
{
  *x = std::log(*x);
}

//
// Arguments    : creal_T *x
// Return Type  : void
//
static void b_sinh(creal_T *x)
{
  double x_re;
  double x_im;
  if (x->im == 0.0) {
    x->re = std::sinh(x->re);
    x->im = 0.0;
  } else {
    x_re = x->re;
    x_im = x->im;
    x->re = std::sinh(x->re) * std::cos(x->im);
    x->im = std::cosh(x_re) * std::sin(x_im);
  }
}

//
// Arguments    : double *x
// Return Type  : void
//
static void b_sqrt(double *x)
{
  *x = std::sqrt(*x);
}

//
// Arguments    : creal_T *x
// Return Type  : void
//
static void c_log(creal_T *x)
{
  double x_im;
  double x_re;
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      x->re = std::log(std::abs(x->re));
      x->im = 3.1415926535897931;
    } else {
      x->re = std::log(std::abs(x->re));
      x->im = 0.0;
    }
  } else if ((std::abs(x->re) > 8.9884656743115785E+307) || (std::abs(x->im) >
              8.9884656743115785E+307)) {
    x_im = x->im;
    x_re = x->re;
    x->re = std::log(rt_hypotd_snf(x->re / 2.0, x->im / 2.0)) +
      0.69314718055994529;
    x->im = rt_atan2d_snf(x_im, x_re);
  } else {
    x_im = x->im;
    x_re = x->re;
    x->re = std::log(rt_hypotd_snf(x->re, x->im));
    x->im = rt_atan2d_snf(x_im, x_re);
  }
}

//
// Arguments    : creal_T *x
// Return Type  : void
//
static void c_sqrt(creal_T *x)
{
  double xr;
  double xi;
  double absxi;
  double absxr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = std::sqrt(-xr);
    } else {
      absxi = std::sqrt(xr);
      xr = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxi = std::sqrt(-xi / 2.0);
      xr = -absxi;
    } else {
      absxi = std::sqrt(xi / 2.0);
      xr = absxi;
    }
  } else if (rtIsNaN(xr)) {
    absxi = xr;
  } else if (rtIsNaN(xi)) {
    absxi = xi;
    xr = xi;
  } else if (rtIsInf(xi)) {
    absxi = std::abs(xi);
    xr = xi;
  } else if (rtIsInf(xr)) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = xi * -xr;
    } else {
      absxi = xr;
      xr = 0.0;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = std::sqrt(absxi) * std::sqrt(1.0 + absxr / absxi);
      } else {
        absxi = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxi = std::sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (xr > 0.0) {
      xr = 0.5 * (xi / absxi);
    } else {
      if (xi < 0.0) {
        xr = -absxi;
      } else {
        xr = absxi;
      }

      absxi = 0.5 * (xi / xr);
    }
  }

  x->re = absxi;
  x->im = xr;
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int mr
//                creal_T *y
//                double rl
//                double tol
//                double elim
//                double alim
// Return Type  : int
//
static int cacai(const creal_T z, double fnu, int kode, int mr, creal_T *y,
                 double rl, double tol, double elim, double alim)
{
  int nz;
  creal_T zn;
  double az;
  boolean_T guard1 = false;
  int nw;
  double crsc_re;
  boolean_T iflag;
  double hz_re;
  double hz_im;
  double cz_re;
  double cz_im;
  double acz;
  creal_T cy[2];
  creal_T s2;
  double ak1_re;
  double ak1_im;
  double s1_re;
  double ascle;
  double s1_im;
  double b_atol;
  double ak;
  double s;
  double rs;
  double b_ak1_re;
  nz = 0;
  zn.re = -z.re;
  zn.im = -z.im;
  az = rt_hypotd_snf(z.re, z.im);
  guard1 = false;
  if ((!(az <= 2.0)) && (az * az * 0.25 > ((fnu + 1.0) - 1.0) + 1.0)) {
    if (az < rl) {
      nw = cmlri(zn, fnu, kode, 1, y, tol);
    } else {
      nw = casyi(zn, fnu, kode, 1, y, rl, tol, elim);
    }

    if (nw < 0) {
      if (nw == -2) {
        nz = -2;
      } else {
        nz = -1;
      }
    } else {
      guard1 = true;
    }
  } else {
    az = rt_hypotd_snf(-z.re, -z.im);
    if (az == 0.0) {
      if (fnu == 0.0) {
        y->re = 1.0;
        y->im = 0.0;
      } else {
        y->re = 0.0;
        y->im = 0.0;
      }
    } else {
      crsc_re = 1.0;
      iflag = false;
      if (az < 2.2250738585072014E-305) {
        if (fnu == 0.0) {
          y->re = 1.0;
          y->im = 0.0;
        } else {
          y->re = 0.0;
          y->im = 0.0;
        }
      } else {
        hz_re = 0.5 * -z.re;
        hz_im = 0.5 * -z.im;
        if (az > 4.7170688552396617E-153) {
          cz_re = hz_re * hz_re - hz_im * hz_im;
          cz_im = hz_re * hz_im + hz_im * hz_re;
          acz = rt_hypotd_snf(cz_re, cz_im);
        } else {
          cz_re = 0.0;
          cz_im = 0.0;
          acz = 0.0;
        }

        s2.re = hz_re;
        s2.im = hz_im;
        c_log(&s2);
        az = ((fnu + 1.0) - 1.0) + 1.0;
        gammaln(&az);
        ak1_re = s2.re * ((fnu + 1.0) - 1.0) - az;
        ak1_im = s2.im * ((fnu + 1.0) - 1.0);
        if (kode == 2) {
          ak1_re -= -z.re;
        }

        if (ak1_re > -elim) {
          ascle = 0.0;
          if (ak1_re <= -alim) {
            iflag = true;
            crsc_re = tol;
            ascle = 2.2250738585072014E-305 / tol;
          }

          az = std::exp(ak1_re);
          if (iflag) {
            az /= tol;
          }

          hz_re = az * std::cos(ak1_im);
          hz_im = az * std::sin(ak1_im);
          b_atol = tol * acz / (((fnu + 1.0) - 1.0) + 1.0);
          s1_re = 1.0;
          s1_im = 0.0;
          if (!(acz < tol * (fnu + 1.0))) {
            ak1_re = 1.0;
            ak1_im = 0.0;
            ak = (fnu + 1.0) + 2.0;
            s = fnu + 1.0;
            az = 2.0;
            do {
              rs = 1.0 / s;
              b_ak1_re = ak1_re;
              ak1_re = ak1_re * cz_re - ak1_im * cz_im;
              ak1_im = b_ak1_re * cz_im + ak1_im * cz_re;
              ak1_re *= rs;
              ak1_im *= rs;
              s1_re += ak1_re;
              s1_im += ak1_im;
              s += ak;
              ak += 2.0;
              az = az * acz * rs;
            } while (!!(az > b_atol));
          }

          s2.re = s1_re * hz_re - s1_im * hz_im;
          s2.im = s1_re * hz_im + s1_im * hz_re;
          if (iflag && (cuchk(s2, ascle, tol) != 0)) {
          } else {
            y->re = s2.re * crsc_re - s2.im * 0.0;
            y->im = s2.re * 0.0 + s2.im * crsc_re;
          }
        } else {
          y->re = 0.0;
          y->im = 0.0;
        }
      }
    }

    guard1 = true;
  }

  if (guard1) {
    for (nw = 0; nw < 2; nw++) {
      cy[nw].re = 0.0;
      cy[nw].im = 0.0;
    }

    nw = b_cbknu(zn, fnu, kode, 1, cy, tol, elim, alim);
    if (nw != 0) {
      if (nw == -2) {
        nz = -2;
      } else {
        nz = -1;
      }
    } else {
      if (mr < 0) {
        hz_im = 3.1415926535897931;
      } else {
        hz_im = -3.1415926535897931;
      }

      ak1_re = 0.0;
      ak1_im = hz_im;
      if (kode != 1) {
        az = std::cos(-(-z.im));
        hz_re = std::sin(-(-z.im));
        ak1_re = 0.0 * az - hz_im * hz_re;
        ak1_im = 0.0 * hz_re + hz_im * az;
      }

      nw = (int)fnu;
      az = (fnu - (double)nw) * hz_im;
      s1_re = std::cos(az);
      s1_im = std::sin(az);
      if ((nw & 1) != 0) {
        s1_re = -s1_re;
        s1_im = -s1_im;
      }

      s2 = cy[0];
      if (kode != 1) {
        ascle = 2.2250738585072014E-305 / tol;
        az = rt_hypotd_snf(cy[0].re, cy[0].im);
        if (az > 0.0) {
          if ((-(-z.re) - (-z.re)) + std::log(az) < -alim) {
            s2.re = 0.0;
            s2.im = 0.0;
            az = 0.0;
          } else {
            s2 = cy[0];
            c_log(&s2);
            s2.re = (s2.re - (-z.re)) - (-z.re);
            s2.im = (s2.im - (-z.im)) - (-z.im);
            if (s2.im == 0.0) {
              s2.re = std::exp(s2.re);
              s2.im = 0.0;
            } else if (rtIsInf(s2.im) && rtIsInf(s2.re) && (s2.re < 0.0)) {
              s2.re = 0.0;
              s2.im = 0.0;
            } else {
              az = std::exp(s2.re / 2.0);
              s2.re = az * (az * std::cos(s2.im));
              s2.im = az * (az * std::sin(s2.im));
            }

            az = rt_hypotd_snf(s2.re, s2.im);
          }
        }

        if ((az > ascle) || (rt_hypotd_snf(y->re, y->im) > ascle)) {
          nz = 0;
        } else {
          s2.re = 0.0;
          s2.im = 0.0;
          y->re = 0.0;
          y->im = 0.0;
          nz = 1;
        }
      }

      az = y->re;
      hz_re = y->im;
      y->re = ak1_re * az - ak1_im * hz_re;
      y->im = ak1_re * hz_re + ak1_im * az;
      y->re += s1_re * s2.re - s1_im * s2.im;
      y->im += s1_re * s2.im + s1_im * s2.re;
    }
  }

  return nz;
}

//
// Arguments    : const creal_T z
//                int id
//                int kode
// Return Type  : creal_T
//
static creal_T cairy(const creal_T z, int id, int kode)
{
  creal_T ai;
  double az;
  double s1_re;
  double az3;
  double s1_im;
  double r;
  creal_T s2;
  creal_T trm2;
  double aa;
  int iflag;
  creal_T trm1;
  double ak;
  double atrm;
  boolean_T guard1 = false;
  double z3_re;
  boolean_T guard2 = false;
  double z3_im;
  boolean_T guard3 = false;
  double d2;
  double ad;
  int i;
  double bk;
  boolean_T exitg1;
  int b_z;
  double b_z3_re;
  double b_z3_im;
  ai.re = 0.0;
  ai.im = 0.0;
  az = rt_hypotd_snf(z.re, z.im);
  if (az > 1.0) {
    az3 = (1.0 + (double)id) / 3.0;
    r = std::log(az);
    trm2 = z;
    c_sqrt(&trm2);
    s2.re = 0.66666666666666663 * (z.re * trm2.re - z.im * trm2.im);
    s2.im = 0.66666666666666663 * (z.re * trm2.im + z.im * trm2.re);
    iflag = 0;
    aa = 1.0;
    ak = s2.im;
    if (z.re < 0.0) {
      s2.re = -std::abs(s2.re);
    }

    if ((z.im != 0.0) || (z.re > 0.0)) {
    } else {
      s2.re = 0.0;
      s2.im = ak;
    }

    guard1 = false;
    guard2 = false;
    guard3 = false;
    if ((s2.re >= 0.0) && (z.re > 0.0)) {
      if ((kode != 2) && (s2.re >= 664.87164553371019)) {
        iflag = 2;
        aa = 4.503599627370496E+15;
        if (-s2.re - 0.25 * r < -700.92179369444591) {
        } else {
          guard3 = true;
        }
      } else {
        guard3 = true;
      }
    } else if ((kode != 2) && (s2.re <= -664.87164553371019)) {
      iflag = 1;
      aa = 2.2204460492503131E-16;
      if (-s2.re + 0.25 * r > 700.92179369444591) {
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard3) {
      cbknu(s2, az3, kode, 664.87164553371019, &trm1, &i);
      guard1 = true;
    }

    if (guard2) {
      trm1.re = 0.0;
      trm1.im = 0.0;
      if (z.im < 0.0) {
        b_z = -1;
      } else {
        b_z = 1;
      }

      i = cacai(s2, az3, kode, b_z, &trm1, 21.784271729432426,
                2.2204460492503131E-16, 700.92179369444591, 664.87164553371019);
      if (i < 0) {
      } else {
        guard1 = true;
      }
    }

    if (guard1) {
      s1_re = 0.18377629847393068 * trm1.re;
      s1_im = 0.18377629847393068 * trm1.im;
      if (iflag != 0) {
        s1_re *= aa;
        s1_im *= aa;
        if (id == 1) {
          s1_re = -s1_re;
          s1_im = -s1_im;
          r = s1_re;
          s1_re = s1_re * z.re - s1_im * z.im;
          s1_im = r * z.im + s1_im * z.re;
        } else {
          r = s1_re;
          s1_re = s1_re * trm2.re - s1_im * trm2.im;
          s1_im = r * trm2.im + s1_im * trm2.re;
        }

        ai.re = 1.0 / aa * s1_re;
        ai.im = 1.0 / aa * s1_im;
      } else if (id == 1) {
        ai.re = -z.re * s1_re - -z.im * s1_im;
        ai.im = -z.re * s1_im + -z.im * s1_re;
      } else {
        ai.re = trm2.re * s1_re - trm2.im * s1_im;
        ai.im = trm2.re * s1_im + trm2.im * s1_re;
      }
    }
  } else {
    s1_re = 1.0;
    s1_im = 0.0;
    s2.re = 1.0;
    s2.im = 0.0;
    if (az < 2.2204460492503131E-16) {
      s1_re = 0.0;
      s1_im = 0.0;
      if (id == 1) {
        if (az > 4.7170688552396617E-153) {
          s1_re = 0.5 * (z.re * z.re - z.im * z.im);
          s1_im = 0.5 * (z.re * z.im + z.im * z.re);
        }

        s1_re *= 0.35502805388781722;
        s1_im *= 0.35502805388781722;
        ai.re = -0.25881940379280682 + s1_re;
        ai.im = s1_im;
      } else {
        if (az > 2.2250738585072014E-305) {
          s1_re = 0.25881940379280682 * z.re;
          s1_im = 0.25881940379280682 * z.im;
        }

        ai.re = 0.35502805388781722 - s1_re;
        ai.im = 0.0 - s1_im;
      }
    } else {
      aa = az * az;
      if (aa >= 2.2204460492503131E-16 / az) {
        trm1.re = 1.0;
        trm1.im = 0.0;
        trm2.re = 1.0;
        trm2.im = 0.0;
        atrm = 1.0;
        r = z.re * z.re - z.im * z.im;
        az3 = z.re * z.im + z.im * z.re;
        z3_re = r * z.re - az3 * z.im;
        z3_im = r * z.im + az3 * z.re;
        az3 = az * aa;
        aa = (2.0 + (double)id) * ((3.0 + (double)id) + (double)id);
        d2 = ((3.0 - (double)id) - (double)id) * (4.0 - (double)id);
        if (aa < d2) {
          ad = aa;
        } else {
          ad = d2;
        }

        ak = 24.0 + 9.0 * (double)id;
        bk = 30.0 - 9.0 * (double)id;
        i = 1;
        exitg1 = false;
        while ((!exitg1) && (i < 26)) {
          if (z3_im == 0.0) {
            b_z3_re = z3_re / aa;
            b_z3_im = 0.0;
          } else if (z3_re == 0.0) {
            b_z3_re = 0.0;
            b_z3_im = z3_im / aa;
          } else {
            b_z3_re = z3_re / aa;
            b_z3_im = z3_im / aa;
          }

          r = trm1.re;
          trm1.re = trm1.re * b_z3_re - trm1.im * b_z3_im;
          trm1.im = r * b_z3_im + trm1.im * b_z3_re;
          s1_re += trm1.re;
          s1_im += trm1.im;
          if (z3_im == 0.0) {
            b_z3_re = z3_re / d2;
            b_z3_im = 0.0;
          } else if (z3_re == 0.0) {
            b_z3_re = 0.0;
            b_z3_im = z3_im / d2;
          } else {
            b_z3_re = z3_re / d2;
            b_z3_im = z3_im / d2;
          }

          r = trm2.re;
          trm2.re = trm2.re * b_z3_re - trm2.im * b_z3_im;
          trm2.im = r * b_z3_im + trm2.im * b_z3_re;
          s2.re += trm2.re;
          s2.im += trm2.im;
          atrm = atrm * az3 / ad;
          aa += ak;
          d2 += bk;
          if ((aa < d2) || rtIsNaN(d2)) {
            ad = aa;
          } else {
            ad = d2;
          }

          if (atrm < 2.2204460492503131E-16 * ad) {
            exitg1 = true;
          } else {
            ak += 18.0;
            bk += 18.0;
            i++;
          }
        }
      }

      if (id == 1) {
        ai.re = -0.25881940379280682 * s2.re;
        ai.im = -0.25881940379280682 * s2.im;
        if (az > 2.2204460492503131E-16) {
          r = z.re * z.re - z.im * z.im;
          az3 = z.re * z.im + z.im * z.re;
          ai.re += (r * s1_re - az3 * s1_im) * 0.17751402694390861;
          ai.im += (r * s1_im + az3 * s1_re) * 0.17751402694390861;
        }

        if (kode == 1) {
        } else {
          s2 = z;
          c_sqrt(&s2);
          r = s2.re * z.im + s2.im * z.re;
          s2.re = 0.66666666666666663 * (s2.re * z.re - s2.im * z.im);
          s2.im = 0.66666666666666663 * r;
          if (s2.im == 0.0) {
            s2.re = std::exp(s2.re);
            s2.im = 0.0;
          } else if (rtIsInf(s2.im) && rtIsInf(s2.re) && (s2.re < 0.0)) {
            s2.re = 0.0;
            s2.im = 0.0;
          } else {
            r = std::exp(s2.re / 2.0);
            s2.re = r * (r * std::cos(s2.im));
            s2.im = r * (r * std::sin(s2.im));
          }

          r = ai.re;
          ai.re = ai.re * s2.re - ai.im * s2.im;
          ai.im = r * s2.im + ai.im * s2.re;
        }
      } else {
        ai.re = s1_re * 0.35502805388781722 - (z.re * s2.re - z.im * s2.im) *
          0.25881940379280682;
        ai.im = s1_im * 0.35502805388781722 - (z.re * s2.im + z.im * s2.re) *
          0.25881940379280682;
        if (kode == 1) {
        } else {
          s2 = z;
          c_sqrt(&s2);
          r = s2.re * z.im + s2.im * z.re;
          s2.re = 0.66666666666666663 * (s2.re * z.re - s2.im * z.im);
          s2.im = 0.66666666666666663 * r;
          if (s2.im == 0.0) {
            s2.re = std::exp(s2.re);
            s2.im = 0.0;
          } else if (rtIsInf(s2.im) && rtIsInf(s2.re) && (s2.re < 0.0)) {
            s2.re = 0.0;
            s2.im = 0.0;
          } else {
            r = std::exp(s2.re / 2.0);
            s2.re = r * (r * std::cos(s2.im));
            s2.im = r * (r * std::sin(s2.im));
          }

          r = ai.re;
          ai.re = ai.re * s2.re - ai.im * s2.im;
          ai.im = r * s2.im + ai.im * s2.re;
        }
      }
    }
  }

  return ai;
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T *y
//                double rl
//                double tol
//                double elim
// Return Type  : int
//
static int casyi(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                 double rl, double tol, double elim)
{
  int nz;
  int n;
  double x;
  creal_T ak1;
  double brm;
  double acz;
  double s;
  double cz_re;
  double cz_im;
  double dnu2;
  double fdn;
  double ez_re;
  double ez_im;
  double aez;
  double b_dnu2;
  int jl;
  double p1_re;
  int inu;
  double bk;
  boolean_T exitg1;
  double sqk;
  double b_atol;
  double sgn;
  double cs1_re;
  double cs1_im;
  double cs2_re;
  double cs2_im;
  double ak;
  double aa;
  double bb;
  double dk_re;
  double dk_im;
  boolean_T errflag;
  boolean_T exitg2;
  double ni;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nz = 0;
  x = z.re;
  if (z.im == 0.0) {
    ak1.re = 0.15915494309189535 / z.re;
    ak1.im = 0.0;
  } else if (z.re == 0.0) {
    ak1.re = 0.0;
    ak1.im = -(0.15915494309189535 / z.im);
  } else {
    brm = std::abs(z.re);
    acz = std::abs(z.im);
    if (brm > acz) {
      s = z.im / z.re;
      acz = z.re + s * z.im;
      ak1.re = (0.15915494309189535 + s * 0.0) / acz;
      ak1.im = (0.0 - s * 0.15915494309189535) / acz;
    } else if (acz == brm) {
      if (z.re > 0.0) {
        acz = 0.5;
      } else {
        acz = -0.5;
      }

      if (z.im > 0.0) {
        dnu2 = 0.5;
      } else {
        dnu2 = -0.5;
      }

      ak1.re = (0.15915494309189535 * acz + 0.0 * dnu2) / brm;
      ak1.im = (0.0 * acz - 0.15915494309189535 * dnu2) / brm;
    } else {
      s = z.re / z.im;
      acz = z.im + s * z.re;
      ak1.re = s * 0.15915494309189535 / acz;
      ak1.im = (s * 0.0 - 0.15915494309189535) / acz;
    }
  }

  c_sqrt(&ak1);
  if (kode == 2) {
    cz_re = 0.0;
    cz_im = z.im;
    acz = 0.0;
  } else {
    cz_re = z.re;
    cz_im = z.im;
    acz = z.re;
  }

  if (std::abs(acz) > elim) {
    nz = -1;
    y->re = rtNaN;
    y->im = 0.0;
  } else {
    dnu2 = fnu + fnu;
    if (cz_im == 0.0) {
      cz_re = std::exp(cz_re);
      cz_im = 0.0;
    } else if (rtIsInf(cz_im) && rtIsInf(cz_re) && (cz_re < 0.0)) {
      cz_re = 0.0;
      cz_im = 0.0;
    } else {
      acz = std::exp(cz_re / 2.0);
      cz_re = acz * (acz * std::cos(cz_im));
      cz_im = acz * (acz * std::sin(cz_im));
    }

    acz = ak1.re;
    ak1.re = ak1.re * cz_re - ak1.im * cz_im;
    ak1.im = acz * cz_im + ak1.im * cz_re;
    fdn = 0.0;
    if (dnu2 > 4.7170688552396617E-153) {
      fdn = dnu2 * dnu2;
    }

    ez_re = 8.0 * z.re;
    ez_im = 8.0 * z.im;
    aez = 8.0 * rt_hypotd_snf(z.re, z.im);
    s = tol / aez;
    dnu2 = rl + rl;
    if (dnu2 < 0.0) {
      b_dnu2 = std::ceil(dnu2);
    } else {
      b_dnu2 = std::floor(dnu2);
    }

    jl = (int)b_dnu2 + 2;
    if (z.im != 0.0) {
      inu = (int)fnu;
      acz = (fnu - (double)inu) * 3.1415926535897931;
      dnu2 = std::sin(acz);
      bk = std::cos(acz);
      if (z.im < 0.0) {
        bk = -bk;
      }

      if ((inu & 1) != 0) {
        p1_re = -(-dnu2);
        bk = -bk;
      } else {
        p1_re = -dnu2;
      }
    } else {
      p1_re = 0.0;
      bk = 0.0;
    }

    inu = 1;
    exitg1 = false;
    while ((!exitg1) && (inu <= n)) {
      sqk = fdn - 1.0;
      b_atol = s * std::abs(fdn - 1.0);
      sgn = 1.0;
      cs1_re = 1.0;
      cs1_im = 0.0;
      cs2_re = 1.0;
      cs2_im = 0.0;
      cz_re = 1.0;
      cz_im = 0.0;
      ak = 0.0;
      aa = 1.0;
      bb = aez;
      dk_re = ez_re;
      dk_im = ez_im;
      errflag = true;
      inu = 1;
      exitg2 = false;
      while ((!exitg2) && (inu <= jl)) {
        cz_re *= sqk;
        cz_im *= sqk;
        acz = cz_re;
        if (dk_im == 0.0) {
          if (cz_im == 0.0) {
            cz_re /= dk_re;
            cz_im = 0.0;
          } else if (cz_re == 0.0) {
            cz_re = 0.0;
            cz_im /= dk_re;
          } else {
            cz_re /= dk_re;
            cz_im /= dk_re;
          }
        } else if (dk_re == 0.0) {
          if (cz_re == 0.0) {
            cz_re = cz_im / dk_im;
            cz_im = 0.0;
          } else if (cz_im == 0.0) {
            cz_re = 0.0;
            cz_im = -(acz / dk_im);
          } else {
            cz_re = cz_im / dk_im;
            cz_im = -(acz / dk_im);
          }
        } else {
          brm = std::abs(dk_re);
          acz = std::abs(dk_im);
          if (brm > acz) {
            dnu2 = dk_im / dk_re;
            acz = dk_re + dnu2 * dk_im;
            ni = cz_im - dnu2 * cz_re;
            cz_re = (cz_re + dnu2 * cz_im) / acz;
            cz_im = ni / acz;
          } else if (acz == brm) {
            if (dk_re > 0.0) {
              acz = 0.5;
            } else {
              acz = -0.5;
            }

            if (dk_im > 0.0) {
              dnu2 = 0.5;
            } else {
              dnu2 = -0.5;
            }

            ni = cz_im * acz - cz_re * dnu2;
            cz_re = (cz_re * acz + cz_im * dnu2) / brm;
            cz_im = ni / brm;
          } else {
            dnu2 = dk_re / dk_im;
            acz = dk_im + dnu2 * dk_re;
            ni = dnu2 * cz_im - cz_re;
            cz_re = (dnu2 * cz_re + cz_im) / acz;
            cz_im = ni / acz;
          }
        }

        cs2_re += cz_re;
        cs2_im += cz_im;
        sgn = -sgn;
        cs1_re += cz_re * sgn;
        cs1_im += cz_im * sgn;
        dk_re += ez_re;
        dk_im += ez_im;
        aa = aa * std::abs(sqk) / bb;
        bb += aez;
        ak += 8.0;
        sqk -= ak;
        if (aa <= b_atol) {
          errflag = false;
          exitg2 = true;
        } else {
          inu++;
        }
      }

      if (errflag) {
        nz = -2;
        exitg1 = true;
      } else {
        if (x + x < elim) {
          cz_re = -2.0 * z.re;
          cz_im = -2.0 * z.im;
          if (cz_im == 0.0) {
            cz_re = std::exp(cz_re);
            cz_im = 0.0;
          } else if (rtIsInf(cz_im) && rtIsInf(cz_re) && (cz_re < 0.0)) {
            cz_re = 0.0;
            cz_im = 0.0;
          } else {
            acz = std::exp(cz_re / 2.0);
            cz_re = acz * (acz * std::cos(cz_im));
            cz_im = acz * (acz * std::sin(cz_im));
          }

          acz = cz_re;
          cz_re = cz_re * cs2_re - cz_im * cs2_im;
          cz_im = acz * cs2_im + cz_im * cs2_re;
          acz = cz_re;
          cz_re = cz_re * p1_re - cz_im * bk;
          cz_im = acz * bk + cz_im * p1_re;
          cs1_re += cz_re;
          cs1_im += cz_im;
        }

        fdn = (fdn + 8.0 * fnu) + 4.0;
        p1_re = -p1_re;
        bk = -bk;
        y->re = cs1_re * ak1.re - cs1_im * ak1.im;
        y->im = cs1_re * ak1.im + cs1_im * ak1.re;
        inu = 2;
      }
    }
  }

  return nz;
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                creal_T *cy
//                int *nz
//                int *ierr
// Return Type  : void
//
static void cbesj(const creal_T z, double fnu, int kode, creal_T *cy, int *nz,
                  int *ierr)
{
  double az;
  int inuh;
  double csgn_re;
  double csgn_im;
  creal_T zn;
  double u1;
  double b_az;
  *ierr = 0;
  az = rt_hypotd_snf(z.re, z.im);
  if ((az > 1.0737418235E+9) || (fnu > 1.0737418235E+9)) {
    *ierr = 4;
  } else {
    if ((az > 32767.999992370605) || (fnu > 32767.999992370605)) {
      *ierr = 3;
    }
  }

  inuh = (int)fnu >> 1;
  az = (fnu - (double)(inuh << 1)) * 1.5707963267948966;
  csgn_re = std::cos(az);
  csgn_im = std::sin(az);
  if ((inuh & 1) == 1) {
    csgn_re = -csgn_re;
    csgn_im = -csgn_im;
  }

  zn.re = -z.re * 0.0 - (-z.im);
  zn.im = -z.re + -z.im * 0.0;
  if (z.im < 0.0) {
    zn.re = -zn.re;
    zn.im = -zn.im;
    csgn_im = -csgn_im;
  }

  cy->re = 0.0;
  cy->im = 0.0;
  *nz = cbinu(zn, fnu, kode, cy, 21.784271729432426, 85.921358647162123,
              2.2204460492503131E-16, 701.61506577445994, 665.56491761372422);
  if (*nz < 0) {
    if (*nz == -2) {
      *nz = 0;
      *ierr = 5;
      cy->re = rtNaN;
      cy->im = rtNaN;
    } else {
      *nz = 0;
      *ierr = 2;
      cy->re = rtInf;
      cy->im = 0.0;
    }
  } else {
    if (*nz != 1) {
      zn = *cy;
      az = std::abs(cy->re);
      u1 = std::abs(cy->im);
      if ((az > u1) || rtIsNaN(u1)) {
        b_az = az;
      } else {
        b_az = u1;
      }

      if (b_az <= 1.0020841800044864E-289) {
        zn.re = 4.503599627370496E+15 * cy->re;
        zn.im = 4.503599627370496E+15 * cy->im;
        az = 2.2204460492503131E-16;
      } else {
        az = 1.0;
      }

      cy->re = az * (zn.re * csgn_re - zn.im * csgn_im);
      cy->im = az * (zn.re * csgn_im + zn.im * csgn_re);
    }
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                creal_T *cy
//                double rl
//                double fnul
//                double tol
//                double elim
//                double alim
// Return Type  : int
//
static int cbinu(const creal_T z, double fnu, int kode, creal_T *cy, double rl,
                 double fnul, double tol, double elim, double alim)
{
  int nz;
  double az;
  int nn;
  double dfnu;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  int nw;
  int i;
  creal_T cw[2];
  int b_dfnu;
  double cinu_re;
  double cinu_im;
  double cscl_re;
  double cy_re;
  double ct_re;
  nz = 0;
  az = rt_hypotd_snf(z.re, z.im);
  nn = 1;
  dfnu = fnu;
  guard1 = false;
  guard2 = false;
  guard3 = false;
  if ((az <= 2.0) || (az * az * 0.25 <= fnu + 1.0)) {
    nw = cseri(z, fnu, kode, 1, cy, tol, elim, alim);
    if (nw <= MIN_int32_T) {
      i = MAX_int32_T;
    } else {
      i = -nw;
    }

    if (nw < 0) {
      nz = i;
    } else {
      nz = nw;
    }

    nn = 1 - nz;
    if ((1 - nz == 0) || (nw >= 0)) {
    } else {
      dfnu = (fnu + (1.0 - (double)nz)) - 1.0;
      guard3 = true;
    }
  } else {
    guard3 = true;
  }

  if (guard3) {
    if (az < rl) {
      if (dfnu <= 1.0) {
        nw = cmlri(z, fnu, kode, nn, cy, tol);
        if (nw < 0) {
          if (nw == -2) {
            nz = -2;
          } else {
            nz = -1;
          }
        } else {
          nz = 0;
        }
      } else {
        guard2 = true;
      }
    } else if ((dfnu <= 1.0) || (!(2.0 * az < dfnu * dfnu))) {
      nw = casyi(z, fnu, kode, nn, cy, rl, tol, elim);
      if (nw < 0) {
        if (nw == -2) {
          nz = -2;
        } else {
          nz = -1;
        }
      } else {
        nz = 0;
      }
    } else {
      guard2 = true;
    }
  }

  if (guard2) {
    nw = cuoik(z, fnu, kode, 1, nn, cy, tol, elim, alim);
    if (nw < 0) {
      if (nw == -2) {
        nz = -2;
      } else {
        nz = -1;
      }
    } else {
      nz += nw;
      nn -= nw;
      if (nn == 0) {
      } else {
        dfnu = (fnu + (double)nn) - 1.0;
        if ((dfnu > fnul) || (az > fnul)) {
          if (dfnu > fnul + 1.0) {
            b_dfnu = 0;
          } else {
            b_dfnu = (int)((fnul - dfnu) + 1.0);
          }

          cbuni(z, fnu, kode, nn, cy, b_dfnu, fnul, tol, elim, alim, &i, &nw);
          if (nw < 0) {
            if (nw == -2) {
              nz = -2;
            } else {
              nz = -1;
            }
          } else {
            nz += nw;
            if (i == 0) {
            } else {
              nn = i;
              guard1 = true;
            }
          }
        } else {
          guard1 = true;
        }
      }
    }
  }

  if (guard1) {
    if (az > rl) {
      for (i = 0; i < 2; i++) {
        cw[i].re = 0.0;
        cw[i].im = 0.0;
      }

      nw = b_cuoik(z, fnu, kode, 2, 2, cw, tol, elim, alim);
      if (nw < 0) {
        nz = nn;
        for (i = 1; i <= nn; i++) {
          cy->re = 0.0;
          cy->im = 0.0;
        }
      } else if (nw > 0) {
        nz = -1;
      } else {
        nw = 0;
        i = b_cbknu(z, fnu, kode, 2, cw, tol, elim, alim);
        if (i != 0) {
          nw = -1;
          if (i == -2) {
            nw = -2;
          }
        } else {
          if (!(nn < 1)) {
            nn = 1;
          }

          crati(z, fnu, nn, cy, tol);
          if (kode == 1) {
            cinu_re = 1.0;
            cinu_im = 0.0;
          } else {
            cinu_re = std::cos(z.im);
            cinu_im = std::sin(z.im);
          }

          az = rt_hypotd_snf(cw[1].re, cw[1].im);
          dfnu = 2.2250738585072014E-305 / tol;
          if (az > dfnu) {
            dfnu = 1.0 / dfnu;
            if (az >= dfnu) {
              cscl_re = tol;
            } else {
              cscl_re = 1.0;
            }
          } else {
            cscl_re = 1.0 / tol;
          }

          az = cw[0].re * cscl_re - cw[0].im * 0.0;
          dfnu = cw[0].re * 0.0 + cw[0].im * cscl_re;
          cy_re = (cy->re * az - cy->im * dfnu) + (cw[1].re * cscl_re - cw[1].im
            * 0.0);
          az = (cy->re * dfnu + cy->im * az) + (cw[1].re * 0.0 + cw[1].im *
            cscl_re);
          ct_re = cy_re * z.re - az * z.im;
          dfnu = cy_re * z.im + az * z.re;
          az = 1.0 / rt_hypotd_snf(ct_re, dfnu);
          dfnu = -dfnu;
          ct_re *= az;
          dfnu *= az;
          cinu_re *= az;
          cinu_im *= az;
          az = cinu_re;
          cinu_re = cinu_re * ct_re - cinu_im * dfnu;
          cinu_im = az * dfnu + cinu_im * ct_re;
          cy->re = cinu_re * cscl_re - cinu_im * 0.0;
          cy->im = cinu_re * 0.0 + cinu_im * cscl_re;
        }

        if (nw < 0) {
          if (nw == -2) {
            nz = -2;
          } else {
            nz = -1;
          }
        } else {
          nz = 0;
        }
      }
    } else {
      nw = cmlri(z, fnu, kode, nn, cy, tol);
      if (nw < 0) {
        if (nw == -2) {
          nz = -2;
        } else {
          nz = -1;
        }
      } else {
        nz = 0;
      }
    }
  }

  return nz;
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                double alim
//                creal_T *y
//                int *nz
// Return Type  : void
//
static void cbknu(const creal_T z, double fnu, int kode, double alim, creal_T *y,
                  int *nz)
{
  double yy;
  double caz;
  int iflag;
  double rz_re;
  double fk;
  double rz_im;
  double g1;
  int inu;
  double dnu;
  boolean_T goto_mw110;
  double dnu2;
  double tm;
  double ak;
  double fhs;
  double s1_re;
  double s1_im;
  double s2_re;
  double s2_im;
  creal_T zd;
  double ck_re;
  double ck_im;
  int inub;
  boolean_T goto_mw225;
  boolean_T goto_mw240;
  boolean_T goto_mw270;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  double fc;
  creal_T smu;
  creal_T coef;
  creal_T fmu;
  int kflag;
  double t2;
  double t1;
  int i;
  boolean_T exitg1;
  creal_T p1;
  double p2_re;
  double p2_im;
  int k;
  double cs_re;
  double a1;
  static const double dv3[3] = { 2.2204460492503131E-16, 1.0,
    4.503599627370496E+15 };

  static const double dv4[3] = { 1.0020841800044864E-289, 9.9792015476736E+288,
    1.7976931348623157E+308 };

  boolean_T earlyExit;
  int j;
  creal_T cy[2];
  static const double dv5[3] = { 4.503599627370496E+15, 1.0,
    2.2204460492503131E-16 };

  double b_g1;
  boolean_T b_guard1 = false;
  y->re = 0.0;
  y->im = 0.0;
  yy = z.im;
  caz = rt_hypotd_snf(z.re, z.im);
  *nz = 0;
  iflag = 0;
  if (z.im == 0.0) {
    rz_re = 2.0 / z.re;
    rz_im = 0.0;
  } else if (z.re == 0.0) {
    rz_re = 0.0;
    rz_im = -(2.0 / z.im);
  } else {
    fk = std::abs(z.re);
    g1 = std::abs(z.im);
    if (fk > g1) {
      fk = z.im / z.re;
      g1 = z.re + fk * z.im;
      rz_re = (2.0 + fk * 0.0) / g1;
      rz_im = (0.0 - fk * 2.0) / g1;
    } else if (g1 == fk) {
      if (z.re > 0.0) {
        g1 = 0.5;
      } else {
        g1 = -0.5;
      }

      if (z.im > 0.0) {
        tm = 0.5;
      } else {
        tm = -0.5;
      }

      rz_re = (2.0 * g1 + 0.0 * tm) / fk;
      rz_im = (0.0 * g1 - 2.0 * tm) / fk;
    } else {
      fk = z.re / z.im;
      g1 = z.im + fk * z.re;
      rz_re = fk * 2.0 / g1;
      rz_im = (fk * 0.0 - 2.0) / g1;
    }
  }

  inu = (int)(fnu + 0.5);
  dnu = fnu - (double)inu;
  goto_mw110 = (std::abs(dnu) == 0.5);
  if ((!goto_mw110) && (std::abs(dnu) > 2.2204460492503131E-16)) {
    dnu2 = dnu * dnu;
  } else {
    dnu2 = 0.0;
  }

  ak = 1.0;
  fhs = 0.0;
  s1_re = 0.0;
  s1_im = 0.0;
  s2_re = 0.0;
  s2_im = 0.0;
  zd.re = 0.0;
  zd.im = 0.0;
  ck_re = 1.0;
  ck_im = 0.0;
  inub = 1;
  goto_mw225 = false;
  goto_mw240 = false;
  goto_mw270 = false;
  guard1 = false;
  guard2 = false;
  if (goto_mw110 || (caz > 2.0)) {
    goto_mw110 = false;
    coef = z;
    c_sqrt(&coef);
    if (coef.im == 0.0) {
      coef.re = 1.2533141373155001 / coef.re;
      coef.im = 0.0;
    } else if (coef.re == 0.0) {
      coef.re = 0.0;
      coef.im = -(1.2533141373155001 / coef.im);
    } else {
      fk = std::abs(coef.re);
      g1 = std::abs(coef.im);
      if (fk > g1) {
        fk = coef.im / coef.re;
        g1 = coef.re + fk * coef.im;
        coef.re = (1.2533141373155001 + fk * 0.0) / g1;
        coef.im = (0.0 - fk * 1.2533141373155001) / g1;
      } else if (g1 == fk) {
        if (coef.re > 0.0) {
          g1 = 0.5;
        } else {
          g1 = -0.5;
        }

        if (coef.im > 0.0) {
          tm = 0.5;
        } else {
          tm = -0.5;
        }

        coef.re = (1.2533141373155001 * g1 + 0.0 * tm) / fk;
        coef.im = (0.0 * g1 - 1.2533141373155001 * tm) / fk;
      } else {
        fk = coef.re / coef.im;
        g1 = coef.im + fk * coef.re;
        coef.re = fk * 1.2533141373155001 / g1;
        coef.im = (fk * 0.0 - 1.2533141373155001) / g1;
      }
    }

    kflag = 1;
    if (!(kode == 2)) {
      if (z.re > alim) {
        iflag = 1;
      } else {
        g1 = std::exp(-z.re);
        fmu.re = std::cos(z.im) * g1;
        fmu.im = std::sin(z.im) * -g1;
        t2 = coef.re;
        coef.re = coef.re * fmu.re - coef.im * fmu.im;
        coef.im = t2 * fmu.im + coef.im * fmu.re;
      }
    }

    if (std::abs(dnu) == 0.5) {
      s1_re = coef.re;
      s1_im = coef.im;
      s2_re = coef.re;
      s2_im = coef.im;
      goto_mw110 = true;
    } else {
      ak = std::abs(std::cos(3.1415926535897931 * dnu));
      if (ak == 0.0) {
        s1_re = coef.re;
        s1_im = coef.im;
        s2_re = coef.re;
        s2_im = coef.im;
        goto_mw110 = true;
      } else {
        fhs = std::abs(0.25 - dnu2);
        if (fhs == 0.0) {
          s1_re = coef.re;
          s1_im = coef.im;
          s2_re = coef.re;
          s2_im = coef.im;
          goto_mw110 = true;
        }
      }
    }

    if (!goto_mw110) {
      if (z.re != 0.0) {
        t1 = std::abs(std::atan(z.im / z.re));
      } else {
        t1 = 1.5707963267948966;
      }

      if (28.66666665740641 > caz) {
        g1 = caz;
        b_sqrt(&g1);
        b_sqrt(&g1);
        g1 = 1.8976999933151775 * ak / (2.2204460492503131E-16 * g1);
        b_log(&g1);
        ak = (g1 + caz * std::cos(3.0 * t1 / (1.0 + caz)) / (1.0 + 0.008 * caz))
          / std::cos(14.7 * t1 / (28.0 + caz));
        fk = 0.12125 * ak * ak / caz + 1.5;
        guard2 = true;
      } else {
        g1 = ak / (3.1415926535897931 * caz * 2.2204460492503131E-16);
        fk = 1.0;
        if (g1 >= 1.0) {
          ak = 2.0;
          t2 = (caz + caz) + 2.0;
          a1 = 0.0;
          fc = 1.0;
          earlyExit = true;
          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i < 31)) {
            tm = fc;
            fc = t2 / (fk + 1.0) * fc - fhs / ak * a1;
            a1 = tm;
            t2 += 2.0;
            ak = ((ak + fk) + fk) + 2.0;
            fhs = (fhs + fk) + fk;
            fk++;
            if (g1 < std::abs(fc) * fk) {
              earlyExit = false;
              exitg1 = true;
            } else {
              i++;
            }
          }

          if (earlyExit) {
            *nz = -2;
          } else {
            g1 = 28.66666665740641 / caz;
            b_sqrt(&g1);
            fk += 1.909859317102744 * t1 * g1;
            fhs = std::abs(0.25 - dnu2);
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }
    } else {
      guard1 = true;
    }
  } else {
    fc = 1.0;
    smu.re = rz_re;
    smu.im = rz_im;
    c_log(&smu);
    fmu.re = dnu * smu.re;
    fmu.im = dnu * smu.im;
    if (dnu != 0.0) {
      fc = dnu * 3.1415926535897931;
      fc /= std::sin(fc);
      coef = fmu;
      b_sinh(&coef);
      smu.re = 1.0 / dnu * coef.re;
      smu.im = 1.0 / dnu * coef.im;
    }

    g1 = 1.0 + dnu;
    gammaln(&g1);
    t2 = std::exp(-g1);
    t1 = 1.0 / (t2 * fc);
    if (std::abs(dnu) > 0.1) {
      g1 = (t1 - t2) / (dnu + dnu);
    } else {
      fk = 0.57721566490153287;
      i = 2;
      exitg1 = false;
      while ((!exitg1) && (i < 9)) {
        ak *= dnu2;
        tm = dv2[i - 1] * ak;
        fk += tm;
        if (std::abs(tm) < 2.2204460492503131E-16) {
          exitg1 = true;
        } else {
          i++;
        }
      }

      g1 = -fk;
    }

    g1 *= fc;
    coef = fmu;
    b_cosh(&coef);
    p1.re = smu.re * (0.5 * (t1 + t2) * fc) + g1 * coef.re;
    p1.im = smu.im * (0.5 * (t1 + t2) * fc) + g1 * coef.im;
    b_exp(&fmu);
    p2_re = 0.5 / t2 * fmu.re;
    p2_im = 0.5 / t2 * fmu.im;
    ak = 0.5 / t1;
    if (fmu.im == 0.0) {
      cs_re = ak / fmu.re;
      tm = 0.0;
    } else if (fmu.re == 0.0) {
      if (ak == 0.0) {
        cs_re = 0.0 / fmu.im;
        tm = 0.0;
      } else {
        cs_re = 0.0;
        tm = -(ak / fmu.im);
      }
    } else {
      fk = std::abs(fmu.re);
      g1 = std::abs(fmu.im);
      if (fk > g1) {
        fk = fmu.im / fmu.re;
        g1 = fmu.re + fk * fmu.im;
        cs_re = (ak + fk * 0.0) / g1;
        tm = (0.0 - fk * ak) / g1;
      } else if (g1 == fk) {
        if (fmu.re > 0.0) {
          g1 = 0.5;
        } else {
          g1 = -0.5;
        }

        if (fmu.im > 0.0) {
          tm = 0.5;
        } else {
          tm = -0.5;
        }

        cs_re = (ak * g1 + 0.0 * tm) / fk;
        tm = (0.0 * g1 - ak * tm) / fk;
      } else {
        fk = fmu.re / fmu.im;
        g1 = fmu.im + fk * fmu.re;
        cs_re = fk * ak / g1;
        tm = (fk * 0.0 - ak) / g1;
      }
    }

    s1_re = p1.re;
    s1_im = p1.im;
    s2_re = p2_re;
    s2_im = p2_im;
    ak = 1.0;
    a1 = 1.0;
    fc = 1.0 - dnu2;
    if (!(inu > 0)) {
      if (caz >= 2.2204460492503131E-16) {
        coef.re = 0.25 * (z.re * z.re - z.im * z.im);
        coef.im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          p1.re *= ak;
          p1.im *= ak;
          p1.re += p2_re;
          p1.im += p2_im;
          p1.re += cs_re;
          p1.im += tm;
          p1.re *= 1.0 / fc;
          p1.im *= 1.0 / fc;
          p2_re *= 1.0 / (ak - dnu);
          p2_im *= 1.0 / (ak - dnu);
          cs_re *= 1.0 / (ak + dnu);
          tm *= 1.0 / (ak + dnu);
          g1 = ck_re;
          ck_re = ck_re * coef.re - ck_im * coef.im;
          ck_im = g1 * coef.im + ck_im * coef.re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * p1.re - ck_im * p1.im;
          s1_im += ck_re * p1.im + ck_im * p1.re;
          a1 = a1 * t1 / ak;
          fc = ((fc + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > 2.2204460492503131E-16));
      }

      y->re = s1_re;
      y->im = s1_im;
      if (kode != 1) {
        coef = z;
        b_exp(&coef);
        y->re = coef.re * s1_re - coef.im * s1_im;
        y->im = coef.re * s1_im + coef.im * s1_re;
      }
    } else {
      if (caz >= 2.2204460492503131E-16) {
        coef.re = 0.25 * (z.re * z.re - z.im * z.im);
        coef.im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          p1.re *= ak;
          p1.im *= ak;
          p1.re += p2_re;
          p1.im += p2_im;
          p1.re += cs_re;
          p1.im += tm;
          p1.re *= 1.0 / fc;
          p1.im *= 1.0 / fc;
          p2_re *= 1.0 / (ak - dnu);
          p2_im *= 1.0 / (ak - dnu);
          cs_re *= 1.0 / (ak + dnu);
          tm *= 1.0 / (ak + dnu);
          g1 = ck_re;
          ck_re = ck_re * coef.re - ck_im * coef.im;
          ck_im = g1 * coef.im + ck_im * coef.re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * p1.re - ck_im * p1.im;
          s1_im += ck_re * p1.im + ck_im * p1.re;
          t2 = p2_re - p1.re * ak;
          g1 = p2_im - p1.im * ak;
          s2_re += t2 * ck_re - g1 * ck_im;
          s2_im += t2 * ck_im + g1 * ck_re;
          a1 = a1 * t1 / ak;
          fc = ((fc + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > 2.2204460492503131E-16));
      }

      kflag = 1;
      if ((fnu + 1.0) * std::abs(smu.re) > alim) {
        kflag = 2;
      }

      g1 = dv5[kflag] * s2_re;
      s2_im *= dv5[kflag];
      s2_re = g1 * rz_re - s2_im * rz_im;
      s2_im = g1 * rz_im + s2_im * rz_re;
      s1_re *= dv5[kflag];
      s1_im *= dv5[kflag];
      if (kode != 1) {
        p1 = z;
        b_exp(&p1);
        g1 = s1_re;
        s1_re = s1_re * p1.re - s1_im * p1.im;
        s1_im = g1 * p1.im + s1_im * p1.re;
        g1 = s2_re;
        s2_re = s2_re * p1.re - s2_im * p1.im;
        s2_im = g1 * p1.im + s2_im * p1.re;
      }

      goto_mw110 = true;
      guard1 = true;
    }
  }

  if (guard2) {
    b_fix(&fk);
    k = (int)fk;
    fk = k;
    ak = (double)k * (double)k;
    p1.re = 0.0;
    p1.im = 0.0;
    p2_re = 2.2204460492503131E-16;
    p2_im = 0.0;
    cs_re = 2.2204460492503131E-16;
    tm = 0.0;
    for (i = 1; i <= k; i++) {
      a1 = ak - fk;
      t2 = 2.0 / (fk + 1.0);
      fmu.re = p2_re;
      fmu.im = p2_im;
      fc = (fk + z.re) * t2;
      g1 = yy * t2;
      t2 = p2_re;
      p2_re = p2_re * fc - p2_im * g1;
      p2_im = t2 * g1 + p2_im * fc;
      p2_re -= p1.re;
      p2_im -= p1.im;
      p2_re *= (ak + fk) / (a1 + fhs);
      p2_im *= (ak + fk) / (a1 + fhs);
      p1 = fmu;
      cs_re += p2_re;
      tm += p2_im;
      ak = (a1 - fk) + 1.0;
      fk--;
    }

    fmu.re = 1.0 / rt_hypotd_snf(cs_re, tm);
    tm = -tm;
    g1 = cs_re;
    cs_re = cs_re * fmu.re - tm * 0.0;
    tm = g1 * 0.0 + tm * fmu.re;
    g1 = fmu.re * p2_re - 0.0 * p2_im;
    fc = fmu.re * p2_im + 0.0 * p2_re;
    t2 = coef.re * g1 - coef.im * fc;
    g1 = coef.re * fc + coef.im * g1;
    s1_re = t2 * cs_re - g1 * tm;
    s1_im = t2 * tm + g1 * cs_re;
    if (!(inu > 0)) {
      zd = z;
      if (iflag == 1) {
        goto_mw270 = true;
      } else {
        goto_mw240 = true;
      }
    } else {
      fmu.re = 1.0 / rt_hypotd_snf(p2_re, p2_im);
      g1 = p1.re;
      p1.re = p1.re * fmu.re - p1.im * 0.0;
      p1.im = g1 * 0.0 + p1.im * fmu.re;
      p2_im = -p2_im;
      t2 = p2_re;
      p2_re = p2_re * fmu.re - p2_im * 0.0;
      p2_im = t2 * 0.0 + p2_im * fmu.re;
      g1 = p1.re * p2_im + p1.im * p2_re;
      ak = (dnu + 0.5) - (p1.re * p2_re - p1.im * p2_im);
      fc = 0.0 - (p1.re * p2_im + p1.im * p2_re);
      if (z.im == 0.0) {
        if (0.0 - g1 == 0.0) {
          t2 = ak / z.re;
          g1 = 0.0;
        } else if (ak == 0.0) {
          t2 = 0.0;
          g1 = (0.0 - g1) / z.re;
        } else {
          t2 = ak / z.re;
          g1 = (0.0 - g1) / z.re;
        }
      } else if (z.re == 0.0) {
        if (ak == 0.0) {
          t2 = (0.0 - g1) / z.im;
          g1 = 0.0;
        } else if (0.0 - g1 == 0.0) {
          t2 = 0.0;
          g1 = -(ak / z.im);
        } else {
          t2 = (0.0 - g1) / z.im;
          g1 = -(ak / z.im);
        }
      } else {
        fk = std::abs(z.re);
        g1 = std::abs(z.im);
        if (fk > g1) {
          fk = z.im / z.re;
          g1 = z.re + fk * z.im;
          t2 = (ak + fk * fc) / g1;
          g1 = (fc - fk * ak) / g1;
        } else if (g1 == fk) {
          if (z.re > 0.0) {
            g1 = 0.5;
          } else {
            g1 = -0.5;
          }

          if (z.im > 0.0) {
            tm = 0.5;
          } else {
            tm = -0.5;
          }

          t2 = (ak * g1 + fc * tm) / fk;
          g1 = (fc * g1 - ak * tm) / fk;
        } else {
          fk = z.re / z.im;
          g1 = z.im + fk * z.re;
          t2 = (fk * ak + fc) / g1;
          g1 = (fk * fc - ak) / g1;
        }
      }

      t2++;
      s2_re = s1_re * t2 - s1_im * g1;
      s2_im = s1_re * g1 + s1_im * t2;
      goto_mw110 = true;
    }

    guard1 = true;
  }

  if (guard1) {
    if (goto_mw240 || goto_mw270) {
    } else if (goto_mw110) {
      ck_re = (dnu + 1.0) * rz_re;
      ck_im = (dnu + 1.0) * rz_im;
      inu--;
      if (inu > 0) {
        if (iflag == 1) {
          zd = z;
          fc = z.re;
          k = inu;
          j = 1;
          for (i = 0; i < 2; i++) {
            cy[i].re = 0.0;
            cy[i].im = 0.0;
          }

          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i <= inu)) {
            cs_re = s2_re;
            tm = s2_im;
            g1 = s2_re;
            s2_re = s2_re * ck_re - s2_im * ck_im;
            s2_im = g1 * ck_im + s2_im * ck_re;
            s2_re += s1_re;
            s2_im += s1_im;
            s1_re = cs_re;
            s1_im = tm;
            ck_re += rz_re;
            ck_im += rz_im;
            g1 = rt_hypotd_snf(s2_re, s2_im);
            b_log(&g1);
            b_guard1 = false;
            if (-fc + g1 >= -700.92179369444591) {
              coef.re = s2_re;
              coef.im = s2_im;
              c_log(&coef);
              p2_re = coef.re + -zd.re;
              p2_im = coef.im + -zd.im;
              p1.re = std::exp(p2_re) / 2.2204460492503131E-16 * std::cos(p2_im);
              p1.im = std::exp(p2_re) / 2.2204460492503131E-16 * std::sin(p2_im);
              if (cuchk(p1, 1.0020841800044864E-289, 2.2204460492503131E-16) ==
                  0) {
                j = 1 - j;
                cy[j] = p1;
                if (k == i - 1) {
                  kflag = 0;
                  inub = i + 1;
                  s2_re = cy[j].re;
                  s2_im = cy[j].im;
                  j = 1 - j;
                  s1_re = cy[j].re;
                  s1_im = cy[j].im;
                  if (i + 1 <= inu) {
                    goto_mw225 = true;
                  } else {
                    s1_re = s2_re;
                    s1_im = s2_im;
                    goto_mw240 = true;
                  }

                  exitg1 = true;
                } else {
                  k = i;
                  i++;
                }
              } else {
                b_guard1 = true;
              }
            } else {
              b_guard1 = true;
            }

            if (b_guard1) {
              if (g1 >= 350.46089684722295) {
                fc -= 700.92179369444591;
                s1_re = cs_re * 3.9222272510438042E-305 - tm * 0.0;
                s1_im = cs_re * 0.0 + tm * 3.9222272510438042E-305;
                g1 = s2_re;
                s2_re = s2_re * 3.9222272510438042E-305 - s2_im * 0.0;
                s2_im = g1 * 0.0 + s2_im * 3.9222272510438042E-305;
                zd.re = fc;
                zd.im = yy;
              }

              i++;
            }
          }

          if (goto_mw225 || goto_mw240) {
          } else {
            s1_re = s2_re;
            s1_im = s2_im;
            goto_mw270 = true;
          }
        } else {
          goto_mw225 = true;
        }
      }

      if (goto_mw225 || goto_mw240 || goto_mw270) {
      } else {
        s1_re = s2_re;
        s1_im = s2_im;
        zd = z;
        if (iflag != 1) {
          goto_mw240 = true;
        }
      }
    } else {
      goto_mw225 = true;
    }

    if (goto_mw225 || goto_mw240) {
      if (goto_mw225) {
        p1.re = dv3[kflag];
        t2 = dv4[kflag];
        while (inub <= inu) {
          cs_re = s2_re;
          tm = s2_im;
          g1 = s2_re;
          s2_re = s2_re * ck_re - s2_im * ck_im;
          s2_im = g1 * ck_im + s2_im * ck_re;
          s2_re += s1_re;
          s2_im += s1_im;
          s1_re = cs_re;
          s1_im = tm;
          ck_re += rz_re;
          ck_im += rz_im;
          if (kflag + 1 < 3) {
            p2_re = s2_re * p1.re - s2_im * 0.0;
            p2_im = s2_re * 0.0 + s2_im * p1.re;
            g1 = std::abs(p2_re);
            fc = std::abs(p2_im);
            if ((g1 > fc) || rtIsNaN(fc)) {
              b_g1 = g1;
            } else {
              b_g1 = fc;
            }

            if (b_g1 > t2) {
              kflag++;
              t2 = dv4[kflag];
              s1_re = dv5[kflag] * (cs_re * p1.re - tm * 0.0);
              s1_im = dv5[kflag] * (cs_re * 0.0 + tm * p1.re);
              s2_re = dv5[kflag] * p2_re;
              s2_im = dv5[kflag] * p2_im;
              p1.re = dv3[kflag];
            }
          }

          inub++;
        }

        s1_re = s2_re;
        s1_im = s2_im;
      }

      y->re = dv3[kflag] * s1_re;
      y->im = dv3[kflag] * s1_im;
    } else {
      y->re = s1_re;
      y->im = s1_im;
      *nz = ckscl(zd, y, 700.92179369444591);
      if (1 <= *nz) {
      } else {
        y->re *= 2.2204460492503131E-16;
        y->im *= 2.2204460492503131E-16;
        if (*nz < -2147483646) {
          k = MAX_int32_T;
        } else {
          k = 1 - *nz;
        }

        if (k >= 2) {
          y->re *= 2.2204460492503131E-16;
          y->im *= 2.2204460492503131E-16;
        }
      }
    }
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T *y
//                int nui
//                double fnul
//                double tol
//                double elim
//                double alim
//                int *nlast
//                int *nz
// Return Type  : void
//
static void cbuni(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  int nui, double fnul, double tol, double elim, double alim,
                  int *nlast, int *nz)
{
  int n;
  int iform;
  double fnui;
  double dfnu;
  int nw;
  double gnu;
  creal_T cy[2];
  int nd;
  double cssr[3];
  double rs1;
  double bry0;
  double csrr[3];
  double bry1;
  int iflag;
  double fn;
  double dscl;
  creal_T cwrk[16];
  creal_T s2;
  creal_T s1;
  creal_T zeta2;
  creal_T rz;
  double dscr;
  double brm;
  int exitg1;
  int nn;
  double b_bry0;
  boolean_T goto_mw110;
  double b_rs1;
  boolean_T exitg2;
  int unusedU3;
  boolean_T guard1 = false;
  double c_rs1;
  double d_rs1;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  *nz = 0;
  iform = 1;
  if (std::abs(z.im) > std::abs(z.re) * 1.7321) {
    iform = 2;
  }

  if (nui == 0) {
    if (iform == 2) {
      cuni2(z, fnu, kode, n, y, fnul, tol, elim, alim, nlast, &nw);
    } else {
      cuni1(z, fnu, kode, n, y, fnul, tol, elim, alim, nlast, &nw);
    }

    if (nw < 0) {
      *nz = -1;
      if (nw == -2) {
        *nz = -2;
      }
    } else {
      *nz = nw;
    }
  } else {
    fnui = nui;
    dfnu = (fnu + (double)n) - 1.0;
    gnu = dfnu + (double)nui;
    if (iform == 2) {
      for (iform = 0; iform < 2; iform++) {
        cy[iform].re = 0.0;
        cy[iform].im = 0.0;
      }

      b_cuni2(z, gnu, kode, 2, cy, fnul, tol, elim, alim, nlast, &nw);
    } else {
      for (iform = 0; iform < 2; iform++) {
        cy[iform].re = 0.0;
        cy[iform].im = 0.0;
      }

      nw = 0;
      nd = 2;
      *nlast = 0;
      cssr[0] = 1.0 / tol;
      cssr[1] = 1.0;
      cssr[2] = tol;
      csrr[0] = tol;
      csrr[1] = 1.0;
      csrr[2] = 1.0 / tol;
      bry1 = 2.2250738585072014E-305 / tol;
      if (gnu > 1.0) {
        fn = gnu;
      } else {
        fn = 1.0;
      }

      iform = 0;
      memset(&cwrk[0], 0, sizeof(creal_T) << 4);
      b_cunik(z, fn, 1, 1, tol, &iform, cwrk, &s2, &s1, &zeta2, &rz);
      if (kode == 1) {
        s1.re = zeta2.re - s1.re;
      } else {
        bry0 = z.re + zeta2.re;
        rs1 = z.im + zeta2.im;
        if (rs1 == 0.0) {
          dscr = fn / bry0;
        } else if (bry0 == 0.0) {
          dscr = 0.0;
        } else {
          brm = std::abs(bry0);
          dscl = std::abs(rs1);
          if (brm > dscl) {
            dscl = rs1 / bry0;
            dscr = (fn + dscl * 0.0) / (bry0 + dscl * rs1);
          } else if (dscl == brm) {
            if (bry0 > 0.0) {
              b_bry0 = 0.5;
            } else {
              b_bry0 = -0.5;
            }

            if (rs1 > 0.0) {
              b_rs1 = 0.5;
            } else {
              b_rs1 = -0.5;
            }

            dscr = (fn * b_bry0 + 0.0 * b_rs1) / brm;
          } else {
            dscl = bry0 / rs1;
            dscr = dscl * fn / (rs1 + dscl * bry0);
          }
        }

        s1.re = fn * dscr - s1.re;
      }

      rs1 = s1.re;
      if (std::abs(s1.re) > elim) {
        if (s1.re > 0.0) {
          nw = -1;
        } else {
          nw = 2;
          for (iform = 0; iform < 2; iform++) {
            cy[iform].re = 0.0;
            cy[iform].im = 0.0;
          }
        }
      } else {
        do {
          exitg1 = 0;
          iflag = -1;
          if (2 < nd) {
            nn = 2;
          } else {
            nn = nd;
          }

          goto_mw110 = false;
          iform = 1;
          exitg2 = false;
          while ((!exitg2) && (iform <= nn)) {
            fn = gnu + (double)(nd - iform);
            unusedU3 = 0;
            b_cunik(z, fn, 1, 0, tol, &unusedU3, cwrk, &s2, &s1, &zeta2, &rz);
            if (kode == 1) {
              s1.re = zeta2.re - s1.re;
              s1.im = zeta2.im - s1.im;
            } else {
              bry0 = z.re + zeta2.re;
              rs1 = z.im + zeta2.im;
              if (rs1 == 0.0) {
                dscr = fn / bry0;
                rs1 = 0.0;
              } else if (bry0 == 0.0) {
                if (fn == 0.0) {
                  dscr = 0.0 / rs1;
                  rs1 = 0.0;
                } else {
                  dscr = 0.0;
                  rs1 = -(fn / rs1);
                }
              } else {
                brm = std::abs(bry0);
                dscl = std::abs(rs1);
                if (brm > dscl) {
                  dscl = rs1 / bry0;
                  bry0 += dscl * rs1;
                  dscr = (fn + dscl * 0.0) / bry0;
                  rs1 = (0.0 - dscl * fn) / bry0;
                } else if (dscl == brm) {
                  if (bry0 > 0.0) {
                    dscl = 0.5;
                  } else {
                    dscl = -0.5;
                  }

                  if (rs1 > 0.0) {
                    bry0 = 0.5;
                  } else {
                    bry0 = -0.5;
                  }

                  dscr = (fn * dscl + 0.0 * bry0) / brm;
                  rs1 = (0.0 * dscl - fn * bry0) / brm;
                } else {
                  dscl = bry0 / rs1;
                  bry0 = rs1 + dscl * bry0;
                  dscr = dscl * fn / bry0;
                  rs1 = (dscl * 0.0 - fn) / bry0;
                }
              }

              s1.re = fn * dscr - s1.re;
              s1.im = fn * rs1 - s1.im;
              s1.im += z.im;
            }

            rs1 = s1.re;
            if (std::abs(s1.re) > elim) {
              goto_mw110 = true;
              exitg2 = true;
            } else {
              if (iform == 1) {
                iflag = 1;
              }

              guard1 = false;
              if (std::abs(s1.re) >= alim) {
                rs1 = s1.re + std::log(rt_hypotd_snf(s2.re, s2.im));
                if (std::abs(rs1) > elim) {
                  goto_mw110 = true;
                  exitg2 = true;
                } else {
                  if (iform == 1) {
                    iflag = 0;
                  }

                  if ((rs1 >= 0.0) && (iform == 1)) {
                    iflag = 2;
                  }

                  guard1 = true;
                }
              } else {
                guard1 = true;
              }

              if (guard1) {
                brm = s2.re * rz.re - s2.im * rz.im;
                bry0 = s2.re * rz.im + s2.im * rz.re;
                dscl = std::exp(s1.re) * cssr[iflag] * std::cos(s1.im);
                dscr = std::exp(s1.re) * cssr[iflag] * std::sin(s1.im);
                s2.re = brm * dscl - bry0 * dscr;
                s2.im = brm * dscr + bry0 * dscl;
                if ((iflag + 1 == 1) && (cuchk(s2, bry1, tol) != 0)) {
                  goto_mw110 = true;
                  exitg2 = true;
                } else {
                  cy[nd - iform].re = csrr[iflag] * s2.re;
                  cy[nd - iform].im = csrr[iflag] * s2.im;
                  iform++;
                }
              }
            }
          }

          if (!goto_mw110) {
            exitg1 = 1;
          } else if (rs1 > 0.0) {
            nw = -1;
            exitg1 = 1;
          } else {
            cy[nd - 1].re = 0.0;
            cy[nd - 1].im = 0.0;
            nw++;
            nd--;
            if (nd == 0) {
              exitg1 = 1;
            } else {
              iform = b_cuoik(z, gnu, kode, 1, nd, cy, tol, elim, alim);
              if (iform < 0) {
                nw = -1;
                exitg1 = 1;
              } else {
                nd -= iform;
                nw += iform;
                if (nd == 0) {
                  exitg1 = 1;
                } else {
                  if (!((gnu + (double)nd) - 1.0 >= fnul)) {
                    *nlast = nd;
                    exitg1 = 1;
                  }
                }
              }
            }
          }
        } while (exitg1 == 0);
      }
    }

    if (nw < 0) {
      *nz = -1;
      if (nw == -2) {
        *nz = -2;
      }
    } else if (nw != 0) {
      *nlast = n;
    } else {
      rs1 = rt_hypotd_snf(cy[0].re, cy[0].im);
      bry0 = 2.2250738585072014E-305 / tol;
      bry1 = 1.0 / bry0;
      cssr[0] = bry0;
      cssr[1] = bry1;
      cssr[2] = bry1;
      iflag = 2;
      fn = 1.0;
      dscl = 1.0;
      if (rs1 > bry0) {
        if (rs1 >= bry1) {
          iflag = 3;
          fn = tol;
          dscl = tol;
        }
      } else {
        iflag = 1;
        bry1 = bry0;
        fn = 1.0 / tol;
        dscl = fn;
      }

      dscr = 1.0 / fn;
      s1.re = dscl * cy[1].re;
      s1.im = dscl * cy[1].im;
      s2.re = dscl * cy[0].re;
      s2.im = dscl * cy[0].im;
      if (z.im == 0.0) {
        rz.re = 2.0 / z.re;
        rz.im = 0.0;
      } else if (z.re == 0.0) {
        rz.re = 0.0;
        rz.im = -(2.0 / z.im);
      } else {
        brm = std::abs(z.re);
        dscl = std::abs(z.im);
        if (brm > dscl) {
          dscl = z.im / z.re;
          bry0 = z.re + dscl * z.im;
          rz.re = (2.0 + dscl * 0.0) / bry0;
          rz.im = (0.0 - dscl * 2.0) / bry0;
        } else if (dscl == brm) {
          if (z.re > 0.0) {
            dscl = 0.5;
          } else {
            dscl = -0.5;
          }

          if (z.im > 0.0) {
            bry0 = 0.5;
          } else {
            bry0 = -0.5;
          }

          rz.re = (2.0 * dscl + 0.0 * bry0) / brm;
          rz.im = (0.0 * dscl - 2.0 * bry0) / brm;
        } else {
          dscl = z.re / z.im;
          bry0 = z.im + dscl * z.re;
          rz.re = dscl * 2.0 / bry0;
          rz.im = (dscl * 0.0 - 2.0) / bry0;
        }
      }

      for (iform = 1; iform <= nui; iform++) {
        *y = s2;
        brm = s2.re;
        s2.re = s2.re * rz.re - s2.im * rz.im;
        s2.im = brm * rz.im + s2.im * rz.re;
        s2.re *= dfnu + fnui;
        s2.im *= dfnu + fnui;
        s2.re += s1.re;
        s2.im += s1.im;
        s1 = *y;
        fnui--;
        if (!(iflag >= 3)) {
          y->re = dscr * s2.re;
          y->im = dscr * s2.im;
          rs1 = std::abs(y->re);
          bry0 = std::abs(y->im);
          if ((rs1 > bry0) || rtIsNaN(bry0)) {
            c_rs1 = rs1;
          } else {
            c_rs1 = bry0;
          }

          if (!(c_rs1 <= bry1)) {
            iflag++;
            bry1 = cssr[iflag - 1];
            s1.re *= dscr;
            s1.im *= dscr;
            fn *= tol;
            dscr = 1.0 / fn;
            s1.re *= fn;
            s1.im *= fn;
            s2.re = fn * y->re;
            s2.im = fn * y->im;
          }
        }
      }

      y->re = dscr * s2.re;
      y->im = dscr * s2.im;
      if (n != 1) {
        fnui = n - 1;
        for (iform = 1; iform < n; iform++) {
          *y = s2;
          brm = s2.re;
          s2.re = s2.re * rz.re - s2.im * rz.im;
          s2.im = brm * rz.im + s2.im * rz.re;
          s2.re *= dfnu + fnui;
          s2.im *= dfnu + fnui;
          s2.re += s1.re;
          s2.im += s1.im;
          s1 = *y;
          y->re = dscr * s2.re;
          y->im = dscr * s2.im;
          fnui--;
          if (!(iflag >= 3)) {
            rs1 = std::abs(y->re);
            bry0 = std::abs(y->im);
            if ((rs1 > bry0) || rtIsNaN(bry0)) {
              d_rs1 = rs1;
            } else {
              d_rs1 = bry0;
            }

            if (!(d_rs1 <= bry1)) {
              iflag++;
              bry1 = cssr[iflag - 1];
              s1.re *= dscr;
              s1.im *= dscr;
              fn *= tol;
              dscr = 1.0 / fn;
              s1.re *= fn;
              s1.im *= fn;
              s2.re = fn * y->re;
              s2.im = fn * y->im;
            }
          }
        }
      }
    }
  }
}

//
// Arguments    : const creal_T zr
//                creal_T *y
//                double elim
// Return Type  : int
//
static int ckscl(const creal_T zr, creal_T *y, double elim)
{
  int nz;
  creal_T s1;
  double s1_re;
  s1 = *y;
  nz = 1;
  y->re = 0.0;
  y->im = 0.0;
  if (-zr.re + std::log(rt_hypotd_snf(s1.re, s1.im)) >= -elim) {
    c_log(&s1);
    s1.re += -zr.re;
    s1.im += -zr.im;
    s1_re = s1.re;
    s1.re = std::exp(s1.re) / 2.2204460492503131E-16 * std::cos(s1.im);
    s1.im = std::exp(s1_re) / 2.2204460492503131E-16 * std::sin(s1.im);
    if (cuchk(s1, 1.0020841800044864E-289, 2.2204460492503131E-16) == 0) {
      *y = s1;
      nz = 0;
    }
  }

  return nz;
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T *y
//                double tol
// Return Type  : int
//
static int cmlri(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                 double tol)
{
  int nz;
  int n;
  double az;
  double flooraz;
  int iaz;
  double fixfnu;
  int ifnu;
  int inu;
  double ck_re;
  double flam;
  double ck_im;
  double rho2;
  creal_T rz;
  double ack;
  double p1_re;
  double p1_im;
  double p2_re;
  double p2_im;
  double tst;
  boolean_T earlyExit;
  int icounter;
  int i;
  boolean_T exitg1;
  int kcounter;
  double pt_re;
  boolean_T guard1 = false;
  double pt_im;
  double fkk;
  int itime;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nz = 0;
  az = rt_hypotd_snf(z.re, z.im);
  flooraz = std::floor(az);
  iaz = (int)flooraz;
  if (fnu < 0.0) {
    fixfnu = std::ceil(fnu);
  } else {
    fixfnu = std::floor(fnu);
  }

  ifnu = (int)fixfnu;
  inu = (ifnu + n) - 1;
  if (z.im == 0.0) {
    ck_re = (flooraz + 1.0) / z.re;
    ck_im = 0.0;
    rz.re = 2.0 / z.re;
    rz.im = 0.0;
  } else {
    if (z.re == 0.0) {
      ck_re = 0.0;
      ck_im = -((flooraz + 1.0) / z.im);
    } else {
      flam = std::abs(z.re);
      rho2 = std::abs(z.im);
      if (flam > rho2) {
        rho2 = z.im / z.re;
        ack = z.re + rho2 * z.im;
        ck_re = ((flooraz + 1.0) + rho2 * 0.0) / ack;
        ck_im = (0.0 - rho2 * (flooraz + 1.0)) / ack;
      } else if (rho2 == flam) {
        if (z.re > 0.0) {
          rho2 = 0.5;
        } else {
          rho2 = -0.5;
        }

        if (z.im > 0.0) {
          ack = 0.5;
        } else {
          ack = -0.5;
        }

        ck_re = ((flooraz + 1.0) * rho2 + 0.0 * ack) / flam;
        ck_im = (0.0 * rho2 - (flooraz + 1.0) * ack) / flam;
      } else {
        rho2 = z.re / z.im;
        ack = z.im + rho2 * z.re;
        ck_re = rho2 * (flooraz + 1.0) / ack;
        ck_im = (rho2 * 0.0 - (flooraz + 1.0)) / ack;
      }
    }

    if (z.re == 0.0) {
      rz.re = 0.0;
      rz.im = -(2.0 / z.im);
    } else {
      flam = std::abs(z.re);
      rho2 = std::abs(z.im);
      if (flam > rho2) {
        rho2 = z.im / z.re;
        ack = z.re + rho2 * z.im;
        rz.re = (2.0 + rho2 * 0.0) / ack;
        rz.im = (0.0 - rho2 * 2.0) / ack;
      } else if (rho2 == flam) {
        if (z.re > 0.0) {
          rho2 = 0.5;
        } else {
          rho2 = -0.5;
        }

        if (z.im > 0.0) {
          ack = 0.5;
        } else {
          ack = -0.5;
        }

        rz.re = (2.0 * rho2 + 0.0 * ack) / flam;
        rz.im = (0.0 * rho2 - 2.0 * ack) / flam;
      } else {
        rho2 = z.re / z.im;
        ack = z.im + rho2 * z.re;
        rz.re = rho2 * 2.0 / ack;
        rz.im = (rho2 * 0.0 - 2.0) / ack;
      }
    }
  }

  p1_re = 0.0;
  p1_im = 0.0;
  p2_re = 1.0;
  p2_im = 0.0;
  ack = ((flooraz + 1.0) + 1.0) / az;
  ack += std::sqrt(ack * ack - 1.0);
  rho2 = ack * ack;
  tst = (rho2 + rho2) / ((rho2 - 1.0) * (ack - 1.0)) / tol;
  flooraz++;
  earlyExit = true;
  icounter = 1;
  i = 1;
  exitg1 = false;
  while ((!exitg1) && (i < 81)) {
    icounter++;
    pt_re = p2_re;
    pt_im = p2_im;
    rho2 = ck_re * p2_im + ck_im * p2_re;
    p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
    p2_im = p1_im - rho2;
    p1_re = pt_re;
    p1_im = pt_im;
    ck_re += rz.re;
    ck_im += rz.im;
    if (rt_hypotd_snf(p2_re, p2_im) > tst * flooraz * flooraz) {
      earlyExit = false;
      exitg1 = true;
    } else {
      flooraz++;
      i++;
    }
  }

  if (earlyExit) {
    nz = -2;
  } else {
    kcounter = 1;
    guard1 = false;
    if (inu >= iaz) {
      p1_re = 0.0;
      p1_im = 0.0;
      p2_re = 1.0;
      p2_im = 0.0;
      if (z.im == 0.0) {
        ck_re = ((double)inu + 1.0) / z.re;
        ck_im = 0.0;
      } else if (z.re == 0.0) {
        if ((double)inu + 1.0 == 0.0) {
          ck_re = 0.0 / z.im;
          ck_im = 0.0;
        } else {
          ck_re = 0.0;
          ck_im = -(((double)inu + 1.0) / z.im);
        }
      } else {
        flam = std::abs(z.re);
        rho2 = std::abs(z.im);
        if (flam > rho2) {
          rho2 = z.im / z.re;
          ack = z.re + rho2 * z.im;
          ck_re = (((double)inu + 1.0) + rho2 * 0.0) / ack;
          ck_im = (0.0 - rho2 * ((double)inu + 1.0)) / ack;
        } else if (rho2 == flam) {
          if (z.re > 0.0) {
            rho2 = 0.5;
          } else {
            rho2 = -0.5;
          }

          if (z.im > 0.0) {
            ack = 0.5;
          } else {
            ack = -0.5;
          }

          ck_re = (((double)inu + 1.0) * rho2 + 0.0 * ack) / flam;
          ck_im = (0.0 * rho2 - ((double)inu + 1.0) * ack) / flam;
        } else {
          rho2 = z.re / z.im;
          ack = z.im + rho2 * z.re;
          ck_re = rho2 * ((double)inu + 1.0) / ack;
          ck_im = (rho2 * 0.0 - ((double)inu + 1.0)) / ack;
        }
      }

      tst = std::sqrt(((double)inu + 1.0) / az / tol);
      itime = 1;
      earlyExit = true;
      i = 1;
      exitg1 = false;
      while ((!exitg1) && (i < 81)) {
        kcounter++;
        pt_re = p2_re;
        pt_im = p2_im;
        rho2 = ck_re * p2_im + ck_im * p2_re;
        p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
        p2_im = p1_im - rho2;
        p1_re = pt_re;
        p1_im = pt_im;
        ck_re += rz.re;
        ck_im += rz.im;
        rho2 = rt_hypotd_snf(p2_re, p2_im);
        if (rho2 >= tst * flooraz * flooraz) {
          if (itime == 2) {
            earlyExit = false;
            exitg1 = true;
          } else {
            ack = rt_hypotd_snf(ck_re, ck_im);
            flam = ack + std::sqrt(ack * ack - 1.0);
            ack = rho2 / rt_hypotd_snf(pt_re, pt_im);
            if ((flam < ack) || rtIsNaN(ack)) {
              ack = flam;
            }

            tst *= std::sqrt(ack / (ack * ack - 1.0));
            itime = 2;
            i++;
          }
        } else {
          i++;
        }
      }

      if (earlyExit) {
        nz = -2;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      icounter += iaz;
      iaz = kcounter + inu;
      if (icounter > iaz) {
        iaz = icounter;
      }

      fkk = iaz;
      p1_re = 0.0;
      p1_im = 0.0;
      p2_re = 2.2250738585072014E-305 / tol;
      p2_im = 0.0;
      tst = fnu - fixfnu;
      az = tst + tst;
      flam = ((double)iaz + az) + 1.0;
      gammaln(&flam);
      rho2 = (double)iaz + 1.0;
      gammaln(&rho2);
      ack = az + 1.0;
      gammaln(&ack);
      flooraz = std::exp((flam - rho2) - ack);
      ck_re = 0.0;
      ck_im = 0.0;
      iaz -= inu;
      for (i = 1; i <= iaz; i++) {
        pt_re = p2_re;
        pt_im = p2_im;
        rho2 = (fkk + tst) * rz.re;
        ack = (fkk + tst) * rz.im;
        flam = rho2 * p2_im + ack * p2_re;
        p2_re = p1_re + (rho2 * p2_re - ack * p2_im);
        p2_im = p1_im + flam;
        p1_re = pt_re;
        p1_im = pt_im;
        ack = flooraz * (1.0 - az / (fkk + az));
        ck_re += (ack + flooraz) * pt_re;
        ck_im += (ack + flooraz) * pt_im;
        flooraz = ack;
        fkk--;
      }

      y->re = p2_re;
      y->im = p2_im;
      if (ifnu > 0) {
        for (i = 1; i <= ifnu; i++) {
          pt_re = p2_re;
          pt_im = p2_im;
          rho2 = (fkk + tst) * rz.re;
          ack = (fkk + tst) * rz.im;
          flam = rho2 * p2_im + ack * p2_re;
          p2_re = p1_re + (rho2 * p2_re - ack * p2_im);
          p2_im = p1_im + flam;
          p1_re = pt_re;
          p1_im = pt_im;
          ack = flooraz * (1.0 - az / (fkk + az));
          ck_re += (ack + flooraz) * pt_re;
          ck_im += (ack + flooraz) * pt_im;
          flooraz = ack;
          fkk--;
        }
      }

      pt_re = z.re;
      pt_im = z.im;
      if (kode == 2) {
        pt_re = z.re - z.re;
        pt_im = z.im;
      }

      c_log(&rz);
      flam = 1.0 + tst;
      gammaln(&flam);
      pt_re = ((-tst * rz.re - -0.0 * rz.im) + pt_re) - flam;
      pt_im += -tst * rz.im + -0.0 * rz.re;
      p2_re += ck_re;
      p2_im += ck_im;
      p1_re = 1.0 / rt_hypotd_snf(p2_re, p2_im);
      if (pt_im == 0.0) {
        ck_re = std::exp(pt_re);
        ck_im = 0.0;
      } else if (rtIsInf(pt_im) && rtIsInf(pt_re) && (pt_re < 0.0)) {
        ck_re = 0.0;
        ck_im = 0.0;
      } else {
        rho2 = std::exp(pt_re / 2.0);
        ck_re = rho2 * (rho2 * std::cos(pt_im));
        ck_im = rho2 * (rho2 * std::sin(pt_im));
      }

      ack = ck_re * p1_re - ck_im * 0.0;
      ck_im = ck_re * 0.0 + ck_im * p1_re;
      rho2 = p2_re * p1_re + p2_im * 0.0;
      p2_im = p2_re * 0.0 - p2_im * p1_re;
      ck_re = ack * rho2 - ck_im * p2_im;
      ck_im = ack * p2_im + ck_im * rho2;
      i = 1;
      while (i <= n) {
        rho2 = y->re;
        ack = y->im;
        y->re = rho2 * ck_re - ack * ck_im;
        y->im = rho2 * ck_im + ack * ck_re;
        i = 2;
      }
    }
  }

  return nz;
}

//
// Arguments    : const creal_T z
//                double fnu
//                int nin
//                creal_T *cy
//                double tol
// Return Type  : void
//
static void crati(const creal_T z, double fnu, int nin, creal_T *cy, double tol)
{
  int magz;
  int idnu;
  int id;
  int k;
  double rz_re;
  double brm;
  double rz_im;
  double flam;
  double test;
  double t1_re;
  double t1_im;
  double p2_re;
  double p2_im;
  double ap2;
  double test1;
  double p1_re;
  double p1_im;
  double ap1;
  double pt_re;
  double pt_im;
  if (nin < 1) {
    magz = nin;
  } else {
    magz = 1;
  }

  if (!(magz < 1)) {
    idnu = (int)fnu;
    magz = (int)rt_hypotd_snf(z.re, z.im);
    if (magz > 2147483646) {
      id = MAX_int32_T;
    } else {
      id = magz + 1;
    }

    if (idnu > id) {
      id = 0;
    } else {
      id -= idnu;
    }

    k = 1;
    if (z.im == 0.0) {
      rz_re = 2.0 / z.re;
      rz_im = 0.0;
    } else if (z.re == 0.0) {
      rz_re = 0.0;
      rz_im = -(2.0 / z.im);
    } else {
      brm = std::abs(z.re);
      flam = std::abs(z.im);
      if (brm > flam) {
        test = z.im / z.re;
        flam = z.re + test * z.im;
        rz_re = (2.0 + test * 0.0) / flam;
        rz_im = (0.0 - test * 2.0) / flam;
      } else if (flam == brm) {
        if (z.re > 0.0) {
          flam = 0.5;
        } else {
          flam = -0.5;
        }

        if (z.im > 0.0) {
          test = 0.5;
        } else {
          test = -0.5;
        }

        rz_re = (2.0 * flam + 0.0 * test) / brm;
        rz_im = (0.0 * flam - 2.0 * test) / brm;
      } else {
        test = z.re / z.im;
        flam = z.im + test * z.re;
        rz_re = test * 2.0 / flam;
        rz_im = (test * 0.0 - 2.0) / flam;
      }
    }

    flam = (double)magz + 1.0;
    test = idnu;
    if (flam > test) {
      test = flam;
    }

    t1_re = test * rz_re;
    flam = (double)magz + 1.0;
    test = idnu;
    if (flam > test) {
      test = flam;
    }

    t1_im = test * rz_im;
    p2_re = -t1_re;
    p2_im = -t1_im;
    t1_re += rz_re;
    t1_im += rz_im;
    ap2 = rt_hypotd_snf(p2_re, p2_im);
    test1 = std::sqrt((ap2 + ap2) / tol);
    test = test1;
    p1_re = 1.0;
    p1_im = 0.0;
    magz = 1;
    while (magz <= 2) {
      k++;
      ap1 = ap2;
      pt_re = p2_re;
      pt_im = p2_im;
      brm = p2_re;
      p2_re = t1_re * p2_re - t1_im * p2_im;
      p2_im = t1_re * p2_im + t1_im * brm;
      p2_re = p1_re - p2_re;
      p2_im = p1_im - p2_im;
      p1_re = pt_re;
      p1_im = pt_im;
      t1_re += rz_re;
      t1_im += rz_im;
      ap2 = rt_hypotd_snf(p2_re, p2_im);
      if (!(ap1 <= test)) {
        if (magz == 1) {
          test = rt_hypotd_snf(t1_re, t1_im) * 0.5;
          flam = test + std::sqrt(test * test - 1.0);
          brm = ap2 / ap1;
          if ((brm < flam) || rtIsNaN(flam)) {
            flam = brm;
          }

          test = test1 * std::sqrt(flam / (flam * flam - 1.0));
        }

        magz++;
      }
    }

    magz = (k + id) + 1;
    t1_re = magz;
    p1_re = 1.0 / ap2;
    p1_im = 0.0;
    p2_re = 0.0;
    p2_im = 0.0;
    for (idnu = 1; idnu <= magz; idnu++) {
      pt_re = p1_re;
      pt_im = p1_im;
      flam = ((fnu + 1.0) - 1.0) + t1_re;
      ap2 = rz_re * flam - rz_im * 0.0;
      test = rz_re * 0.0 + rz_im * flam;
      brm = ap2 * p1_im + test * p1_re;
      p1_re = (ap2 * p1_re - test * p1_im) + p2_re;
      p1_im = brm + p2_im;
      p2_re = pt_re;
      p2_im = pt_im;
      t1_re--;
    }

    if ((p1_re == 0.0) && (p1_im == 0.0)) {
      p1_re = tol;
      p1_im = tol;
    }

    if (p1_im == 0.0) {
      if (p2_im == 0.0) {
        cy->re = p2_re / p1_re;
        cy->im = 0.0;
      } else if (p2_re == 0.0) {
        cy->re = 0.0;
        cy->im = p2_im / p1_re;
      } else {
        cy->re = p2_re / p1_re;
        cy->im = p2_im / p1_re;
      }
    } else if (p1_re == 0.0) {
      if (p2_re == 0.0) {
        cy->re = p2_im / p1_im;
        cy->im = 0.0;
      } else if (p2_im == 0.0) {
        cy->re = 0.0;
        cy->im = -(p2_re / p1_im);
      } else {
        cy->re = p2_im / p1_im;
        cy->im = -(p2_re / p1_im);
      }
    } else {
      brm = std::abs(p1_re);
      flam = std::abs(p1_im);
      if (brm > flam) {
        test = p1_im / p1_re;
        flam = p1_re + test * p1_im;
        cy->re = (p2_re + test * p2_im) / flam;
        cy->im = (p2_im - test * p2_re) / flam;
      } else if (flam == brm) {
        if (p1_re > 0.0) {
          flam = 0.5;
        } else {
          flam = -0.5;
        }

        if (p1_im > 0.0) {
          test = 0.5;
        } else {
          test = -0.5;
        }

        cy->re = (p2_re * flam + p2_im * test) / brm;
        cy->im = (p2_im * flam - p2_re * test) / brm;
      } else {
        test = p1_re / p1_im;
        flam = p1_im + test * p1_re;
        cy->re = (test * p2_re + p2_im) / flam;
        cy->im = (test * p2_im - p2_re) / flam;
      }
    }
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T *y
//                double tol
//                double elim
//                double alim
// Return Type  : int
//
static int cseri(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                 double tol, double elim, double alim)
{
  int nz;
  int n;
  double az;
  double crsc_re;
  boolean_T iflag;
  double hz_re;
  double hz_im;
  double cz_re;
  double cz_im;
  double acz;
  double dfnu;
  creal_T ak1;
  double ascle;
  double aa;
  double coef_re;
  double coef_im;
  double b_atol;
  int i;
  boolean_T exitg1;
  double s1_re;
  double s1_im;
  double ak;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nz = 0;
  az = rt_hypotd_snf(z.re, z.im);
  if (az == 0.0) {
    if (fnu == 0.0) {
      y->re = 1.0;
      y->im = 0.0;
    } else {
      y->re = 0.0;
      y->im = 0.0;
    }
  } else {
    crsc_re = 1.0;
    iflag = false;
    if (az < 2.2250738585072014E-305) {
      nz = n;
      if (fnu == 0.0) {
        nz = n - 1;
        y->re = 1.0;
        y->im = 0.0;
      } else {
        y->re = 0.0;
        y->im = 0.0;
      }
    } else {
      hz_re = 0.5 * z.re;
      hz_im = 0.5 * z.im;
      if (az > 4.7170688552396617E-153) {
        cz_re = hz_re * hz_re - hz_im * hz_im;
        cz_im = hz_re * hz_im + hz_im * hz_re;
        acz = rt_hypotd_snf(cz_re, cz_im);
      } else {
        cz_re = 0.0;
        cz_im = 0.0;
        acz = 0.0;
      }

      dfnu = (fnu + (double)n) - 1.0;
      ak1.re = hz_re;
      ak1.im = hz_im;
      c_log(&ak1);
      az = dfnu + 1.0;
      gammaln(&az);
      ak1.re = ak1.re * dfnu - az;
      ak1.im *= dfnu;
      if (kode == 2) {
        ak1.re -= z.re;
      }

      if (ak1.re > -elim) {
        ascle = 0.0;
        if (ak1.re <= -alim) {
          iflag = true;
          crsc_re = tol;
          ascle = 2.2250738585072014E-305 / tol;
        }

        aa = std::exp(ak1.re);
        if (iflag) {
          aa /= tol;
        }

        coef_re = aa * std::cos(ak1.im);
        coef_im = aa * std::sin(ak1.im);
        b_atol = tol * acz / (dfnu + 1.0);
        i = 1;
        exitg1 = false;
        while ((!exitg1) && (i <= n)) {
          s1_re = 1.0;
          s1_im = 0.0;
          if (!(acz < tol * (fnu + 1.0))) {
            ak1.re = 1.0;
            ak1.im = 0.0;
            ak = (fnu + 1.0) + 2.0;
            az = fnu + 1.0;
            aa = 2.0;
            do {
              hz_re = 1.0 / az;
              hz_im = ak1.re;
              ak1.re = hz_re * (ak1.re * cz_re - ak1.im * cz_im);
              ak1.im = hz_re * (hz_im * cz_im + ak1.im * cz_re);
              s1_re += ak1.re;
              s1_im += ak1.im;
              az += ak;
              ak += 2.0;
              aa = aa * acz * hz_re;
            } while (!!(aa > b_atol));
          }

          ak1.re = s1_re * coef_re - s1_im * coef_im;
          ak1.im = s1_re * coef_im + s1_im * coef_re;
          if (iflag && (cuchk(ak1, ascle, tol) != 0)) {
            exitg1 = true;
          } else {
            y->re = ak1.re * crsc_re - ak1.im * 0.0;
            y->im = ak1.re * 0.0 + ak1.im * crsc_re;
            i = 2;
          }
        }
      } else {
        nz = 1;
        y->re = 0.0;
        y->im = 0.0;
        if (acz > dfnu) {
          nz = -1;
        }
      }
    }
  }

  return nz;
}

//
// Arguments    : const creal_T y
//                double ascle
//                double tol
// Return Type  : int
//
static int cuchk(const creal_T y, double ascle, double tol)
{
  int nz;
  double yr;
  double yi;
  double smallpart;
  yr = std::abs(y.re);
  yi = std::abs(y.im);
  if (yr > yi) {
    smallpart = yi;
    yi = yr;
  } else {
    smallpart = yr;
  }

  if ((smallpart <= ascle) && (yi < smallpart / tol)) {
    nz = 1;
  } else {
    nz = 0;
  }

  return nz;
}

//
// Arguments    : const creal_T z
//                double fnu
//                double tol
//                creal_T *phi
//                creal_T *arg
//                creal_T *zeta1
//                creal_T *zeta2
// Return Type  : void
//
static void cunhj(const creal_T z, double fnu, double tol, creal_T *phi, creal_T
                  *arg, creal_T *zeta1, creal_T *zeta2)
{
  double ac;
  creal_T zb;
  double fn23;
  double rfn13_re;
  creal_T w2;
  creal_T p[30];
  double ap[30];
  double zeta_im;
  double zeta_re;
  int i;
  boolean_T exitg1;
  double brm;
  double s;
  double zth_re;
  double zth_im;
  ac = fnu * 2.2250738585072014E-305;
  if ((std::abs(z.re) <= ac) && (std::abs(z.im) <= ac)) {
    zeta1->re = 1402.9773265065639 + fnu;
    zeta1->im = 0.0;
    zeta2->re = fnu;
    zeta2->im = 0.0;
    phi->re = 1.0;
    phi->im = 0.0;
    arg->re = 1.0;
    arg->im = 0.0;
  } else {
    zb.re = 1.0 / fnu * z.re;
    zb.im = 1.0 / fnu * z.im;
    ac = rt_powd_snf(fnu, 0.33333333333333331);
    fn23 = ac * ac;
    rfn13_re = 1.0 / ac;
    w2.re = 1.0 - (zb.re * zb.re - zb.im * zb.im);
    w2.im = 0.0 - (zb.re * zb.im + zb.im * zb.re);
    ac = rt_hypotd_snf(w2.re, w2.im);
    if (ac > 0.25) {
      c_sqrt(&w2);
      zeta_im = w2.re;
      zeta_re = w2.im;
      if (w2.re < 0.0) {
        zeta_im = 0.0;
      }

      if (w2.im < 0.0) {
        zeta_re = 0.0;
      }

      w2.re = zeta_im;
      w2.im = zeta_re;
      ac = zb.re;
      if (zb.im == 0.0) {
        if (zeta_re == 0.0) {
          zb.re = (1.0 + zeta_im) / zb.re;
          zb.im = 0.0;
        } else if (1.0 + zeta_im == 0.0) {
          zb.re = 0.0;
          zb.im = zeta_re / ac;
        } else {
          zb.re = (1.0 + zeta_im) / zb.re;
          zb.im = zeta_re / ac;
        }
      } else if (zb.re == 0.0) {
        if (1.0 + zeta_im == 0.0) {
          zb.re = zeta_re / zb.im;
          zb.im = 0.0;
        } else if (zeta_re == 0.0) {
          zb.re = 0.0;
          zb.im = -((1.0 + zeta_im) / zb.im);
        } else {
          zb.re = zeta_re / zb.im;
          zb.im = -((1.0 + zeta_im) / zb.im);
        }
      } else {
        brm = std::abs(zb.re);
        ac = std::abs(zb.im);
        if (brm > ac) {
          s = zb.im / zb.re;
          ac = zb.re + s * zb.im;
          zb.re = ((1.0 + zeta_im) + s * zeta_re) / ac;
          zb.im = (zeta_re - s * (1.0 + zeta_im)) / ac;
        } else if (ac == brm) {
          if (zb.re > 0.0) {
            s = 0.5;
          } else {
            s = -0.5;
          }

          if (zb.im > 0.0) {
            ac = 0.5;
          } else {
            ac = -0.5;
          }

          zb.re = ((1.0 + zeta_im) * s + zeta_re * ac) / brm;
          zb.im = (zeta_re * s - (1.0 + zeta_im) * ac) / brm;
        } else {
          s = zb.re / zb.im;
          ac = zb.im + s * zb.re;
          zb.re = (s * (1.0 + zeta_im) + zeta_re) / ac;
          zb.im = (s * zeta_re - (1.0 + zeta_im)) / ac;
        }
      }

      c_log(&zb);
      zeta_im = zb.re;
      zeta_re = zb.im;
      if (zb.re < 0.0) {
        zeta_im = 0.0;
      }

      if (zb.im < 0.0) {
        zeta_re = 0.0;
      } else {
        if (zb.im > 1.5707963267948966) {
          zeta_re = 1.5707963267948966;
        }
      }

      zth_re = 1.5 * (zeta_im - w2.re);
      zth_im = 1.5 * (zeta_re - w2.im);
      zeta1->re = fnu * zeta_im;
      zeta1->im = fnu * zeta_re;
      zeta2->re = fnu * w2.re;
      zeta2->im = fnu * w2.im;
      if ((zth_re >= 0.0) && (zth_im < 0.0)) {
        ac = 4.71238898038469;
      } else if (zth_re != 0.0) {
        ac = std::atan(zth_im / zth_re);
        if (zth_re < 0.0) {
          ac += 3.1415926535897931;
        }
      } else {
        ac = 1.5707963267948966;
      }

      zeta_im = rt_powd_snf(rt_hypotd_snf(zth_re, zth_im), 0.66666666666666663);
      ac *= 0.66666666666666663;
      zeta_re = zeta_im * std::cos(ac);
      zeta_im *= std::sin(ac);
      if (zeta_im < 0.0) {
        zeta_im = 0.0;
      }

      arg->re = fn23 * zeta_re;
      arg->im = fn23 * zeta_im;
      if (zeta_im == 0.0) {
        if (zth_im == 0.0) {
          fn23 = zth_re / zeta_re;
          zth_im = 0.0;
        } else if (zth_re == 0.0) {
          fn23 = 0.0;
          zth_im /= zeta_re;
        } else {
          fn23 = zth_re / zeta_re;
          zth_im /= zeta_re;
        }
      } else if (zeta_re == 0.0) {
        if (zth_re == 0.0) {
          fn23 = zth_im / zeta_im;
          zth_im = 0.0;
        } else if (zth_im == 0.0) {
          fn23 = 0.0;
          zth_im = -(zth_re / zeta_im);
        } else {
          fn23 = zth_im / zeta_im;
          zth_im = -(zth_re / zeta_im);
        }
      } else {
        brm = std::abs(zeta_re);
        if (brm > zeta_im) {
          s = zeta_im / zeta_re;
          ac = zeta_re + s * zeta_im;
          fn23 = (zth_re + s * zth_im) / ac;
          zth_im = (zth_im - s * zth_re) / ac;
        } else if (zeta_im == brm) {
          if (zeta_re > 0.0) {
            s = 0.5;
          } else {
            s = -0.5;
          }

          if (zeta_im > 0.0) {
            ac = 0.5;
          } else {
            ac = -0.5;
          }

          fn23 = (zth_re * s + zth_im * ac) / brm;
          zth_im = (zth_im * s - zth_re * ac) / brm;
        } else {
          s = zeta_re / zeta_im;
          ac = zeta_im + s * zeta_re;
          fn23 = (s * zth_re + zth_im) / ac;
          zth_im = (s * zth_im - zth_re) / ac;
        }
      }

      if (w2.im == 0.0) {
        if (zth_im == 0.0) {
          zb.re = fn23 / w2.re;
          zb.im = 0.0;
        } else if (fn23 == 0.0) {
          zb.re = 0.0;
          zb.im = zth_im / w2.re;
        } else {
          zb.re = fn23 / w2.re;
          zb.im = zth_im / w2.re;
        }
      } else if (w2.re == 0.0) {
        if (fn23 == 0.0) {
          zb.re = zth_im / w2.im;
          zb.im = 0.0;
        } else if (zth_im == 0.0) {
          zb.re = 0.0;
          zb.im = -(fn23 / w2.im);
        } else {
          zb.re = zth_im / w2.im;
          zb.im = -(fn23 / w2.im);
        }
      } else {
        brm = std::abs(w2.re);
        ac = std::abs(w2.im);
        if (brm > ac) {
          s = w2.im / w2.re;
          ac = w2.re + s * w2.im;
          zb.re = (fn23 + s * zth_im) / ac;
          zb.im = (zth_im - s * fn23) / ac;
        } else if (ac == brm) {
          if (w2.re > 0.0) {
            s = 0.5;
          } else {
            s = -0.5;
          }

          if (w2.im > 0.0) {
            ac = 0.5;
          } else {
            ac = -0.5;
          }

          zb.re = (fn23 * s + zth_im * ac) / brm;
          zb.im = (zth_im * s - fn23 * ac) / brm;
        } else {
          s = w2.re / w2.im;
          ac = w2.im + s * w2.re;
          zb.re = (s * fn23 + zth_im) / ac;
          zb.im = (s * zth_im - fn23) / ac;
        }
      }

      zb.re += zb.re;
      zb.im += zb.im;
      c_sqrt(&zb);
      phi->re = zb.re * rfn13_re - zb.im * 0.0;
      phi->im = zb.re * 0.0 + zb.im * rfn13_re;
    } else {
      memset(&p[0], 0, 30U * sizeof(creal_T));
      memset(&ap[0], 0, 30U * sizeof(double));
      p[0].re = 1.0;
      p[0].im = 0.0;
      zb.re = 0.6299605249474366;
      zb.im = 0.0;
      ap[0] = 1.0;
      if (ac >= tol) {
        i = 1;
        exitg1 = false;
        while ((!exitg1) && (i + 1 < 31)) {
          zeta_im = p[i - 1].re;
          p[i].re = p[i - 1].re * w2.re - p[i - 1].im * w2.im;
          p[i].im = zeta_im * w2.im + p[i - 1].im * w2.re;
          zb.re += p[i].re * dv0[i];
          zb.im += p[i].im * dv0[i];
          ap[i] = ap[i - 1] * ac;
          if (ap[i] < tol) {
            exitg1 = true;
          } else {
            i++;
          }
        }
      }

      zeta_re = w2.re * zb.re - w2.im * zb.im;
      zeta_im = w2.re * zb.im + w2.im * zb.re;
      arg->re = fn23 * zeta_re;
      arg->im = fn23 * zeta_im;
      c_sqrt(&zb);
      c_sqrt(&w2);
      zeta2->re = fnu * w2.re;
      zeta2->im = fnu * w2.im;
      ac = (zeta_re * zb.re - zeta_im * zb.im) * 0.66666666666666663 + 1.0;
      zeta_im = (zeta_re * zb.im + zeta_im * zb.re) * 0.66666666666666663;
      zeta1->re = zeta2->re * ac - zeta2->im * zeta_im;
      zeta1->im = zeta2->re * zeta_im + zeta2->im * ac;
      zb.re += zb.re;
      zb.im += zb.im;
      c_sqrt(&zb);
      phi->re = zb.re * rfn13_re - zb.im * 0.0;
      phi->im = zb.re * 0.0 + zb.im * rfn13_re;
    }
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T *y
//                double fnul
//                double tol
//                double elim
//                double alim
//                int *nlast
//                int *nz
// Return Type  : void
//
static void cuni1(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double fnul, double tol, double elim, double alim, int *nlast,
                  int *nz)
{
  int n;
  double cssr[3];
  double csrr[3];
  double bry1;
  double fn;
  int iflag;
  creal_T cwrk[16];
  creal_T s1;
  creal_T zeta1;
  creal_T zeta2;
  creal_T unusedU1;
  double fn_re;
  double bi;
  double rs1;
  double brm;
  int exitg1;
  int i;
  int nn;
  double b_fn_re;
  boolean_T goto_mw110;
  double b_bi;
  boolean_T exitg2;
  int unusedU3;
  creal_T summ;
  boolean_T guard1 = false;
  double sgnbr;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  *nz = 0;
  *nlast = 0;
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry1 = 2.2250738585072014E-305 / tol;
  if (fnu > 1.0) {
    fn = fnu;
  } else {
    fn = 1.0;
  }

  iflag = 0;
  memset(&cwrk[0], 0, sizeof(creal_T) << 4);
  b_cunik(z, fn, 1, 1, tol, &iflag, cwrk, &s1, &zeta1, &zeta2, &unusedU1);
  if (kode == 1) {
    s1.re = zeta2.re - zeta1.re;
  } else {
    fn_re = z.re + zeta2.re;
    bi = z.im + zeta2.im;
    if (bi == 0.0) {
      fn_re = fn / fn_re;
    } else if (fn_re == 0.0) {
      fn_re = 0.0;
    } else {
      brm = std::abs(fn_re);
      rs1 = std::abs(bi);
      if (brm > rs1) {
        brm = bi / fn_re;
        fn_re = (fn + brm * 0.0) / (fn_re + brm * bi);
      } else if (rs1 == brm) {
        if (fn_re > 0.0) {
          b_fn_re = 0.5;
        } else {
          b_fn_re = -0.5;
        }

        if (bi > 0.0) {
          b_bi = 0.5;
        } else {
          b_bi = -0.5;
        }

        fn_re = (fn * b_fn_re + 0.0 * b_bi) / brm;
      } else {
        brm = fn_re / bi;
        fn_re = brm * fn / (bi + brm * fn_re);
      }
    }

    s1.re = fn * fn_re - zeta1.re;
  }

  rs1 = s1.re;
  if (std::abs(s1.re) > elim) {
    if (s1.re > 0.0) {
      *nz = -1;
    } else {
      *nz = n;
      i = 1;
      while (i <= n) {
        y->re = 0.0;
        y->im = 0.0;
        i = 2;
      }
    }
  } else {
    do {
      exitg1 = 0;
      iflag = -1;
      if (2 < n) {
        nn = 2;
      } else {
        nn = n;
      }

      goto_mw110 = false;
      i = 1;
      exitg2 = false;
      while ((!exitg2) && (i <= nn)) {
        fn = fnu + (double)(n - i);
        unusedU3 = 0;
        b_cunik(z, fn, 1, 0, tol, &unusedU3, cwrk, &unusedU1, &zeta1, &zeta2,
                &summ);
        if (kode == 1) {
          s1.re = zeta2.re - zeta1.re;
          s1.im = zeta2.im - zeta1.im;
        } else {
          fn_re = z.re + zeta2.re;
          bi = z.im + zeta2.im;
          if (bi == 0.0) {
            fn_re = fn / fn_re;
            rs1 = 0.0;
          } else if (fn_re == 0.0) {
            if (fn == 0.0) {
              fn_re = 0.0 / bi;
              rs1 = 0.0;
            } else {
              fn_re = 0.0;
              rs1 = -(fn / bi);
            }
          } else {
            brm = std::abs(fn_re);
            rs1 = std::abs(bi);
            if (brm > rs1) {
              brm = bi / fn_re;
              rs1 = fn_re + brm * bi;
              fn_re = (fn + brm * 0.0) / rs1;
              rs1 = (0.0 - brm * fn) / rs1;
            } else if (rs1 == brm) {
              if (fn_re > 0.0) {
                sgnbr = 0.5;
              } else {
                sgnbr = -0.5;
              }

              if (bi > 0.0) {
                rs1 = 0.5;
              } else {
                rs1 = -0.5;
              }

              fn_re = (fn * sgnbr + 0.0 * rs1) / brm;
              rs1 = (0.0 * sgnbr - fn * rs1) / brm;
            } else {
              brm = fn_re / bi;
              rs1 = bi + brm * fn_re;
              fn_re = brm * fn / rs1;
              rs1 = (brm * 0.0 - fn) / rs1;
            }
          }

          s1.re = fn * fn_re - zeta1.re;
          s1.im = (fn * rs1 - zeta1.im) + z.im;
        }

        rs1 = s1.re;
        if (std::abs(s1.re) > elim) {
          goto_mw110 = true;
          exitg2 = true;
        } else {
          if (i == 1) {
            iflag = 1;
          }

          guard1 = false;
          if (std::abs(s1.re) >= alim) {
            rs1 = s1.re + std::log(rt_hypotd_snf(unusedU1.re, unusedU1.im));
            if (std::abs(rs1) > elim) {
              goto_mw110 = true;
              exitg2 = true;
            } else {
              if (i == 1) {
                iflag = 0;
              }

              if ((rs1 >= 0.0) && (i == 1)) {
                iflag = 2;
              }

              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            fn_re = unusedU1.re * summ.re - unusedU1.im * summ.im;
            brm = unusedU1.re * summ.im + unusedU1.im * summ.re;
            sgnbr = std::exp(s1.re) * cssr[iflag] * std::cos(s1.im);
            bi = std::exp(s1.re) * cssr[iflag] * std::sin(s1.im);
            s1.re = fn_re * sgnbr - brm * bi;
            s1.im = fn_re * bi + brm * sgnbr;
            if ((iflag + 1 == 1) && (cuchk(s1, bry1, tol) != 0)) {
              goto_mw110 = true;
              exitg2 = true;
            } else {
              y->re = csrr[iflag] * s1.re;
              y->im = csrr[iflag] * s1.im;
              i++;
            }
          }
        }
      }

      if (!goto_mw110) {
        exitg1 = 1;
      } else if (rs1 > 0.0) {
        *nz = -1;
        exitg1 = 1;
      } else {
        y->re = 0.0;
        y->im = 0.0;
        (*nz)++;
        n--;
        if (n == 0) {
          exitg1 = 1;
        } else {
          iflag = cuoik(z, fnu, kode, 1, n, y, tol, elim, alim);
          if (iflag < 0) {
            *nz = -1;
            exitg1 = 1;
          } else {
            n -= iflag;
            *nz += iflag;
            if (n == 0) {
              exitg1 = 1;
            } else {
              if (!((fnu + (double)n) - 1.0 >= fnul)) {
                *nlast = n;
                exitg1 = 1;
              }
            }
          }
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int nin
//                creal_T *y
//                double fnul
//                double tol
//                double elim
//                double alim
//                int *nlast
//                int *nz
// Return Type  : void
//
static void cuni2(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double fnul, double tol, double elim, double alim, int *nlast,
                  int *nz)
{
  int n;
  double cssr[3];
  double csrr[3];
  double bry1;
  double yy;
  creal_T zn;
  double zb_re;
  double zb_im;
  signed char cid_im;
  double ffnu;
  int inu;
  double c2_re;
  double c2_im;
  double zar_re;
  double zar_im;
  int in;
  static const cint8_T icv0[4] = { { 1,// re
      0                                // im
    }, { 0,                            // re
      1                                // im
    }, { -1,                           // re
      0                                // im
    }, { 0,                            // re
      -1                               // im
    } };

  creal_T dai;
  creal_T s1;
  creal_T zeta1;
  creal_T zeta2;
  creal_T ai;
  creal_T unusedU7;
  double br;
  double bi;
  double rs1;
  double brm;
  int exitg1;
  boolean_T goto_mw120;
  int i;
  double b_br;
  int nn;
  double b_bi;
  boolean_T exitg2;
  double fn;
  creal_T phi;
  creal_T asum;
  creal_T bsum;
  boolean_T guard1 = false;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  *nz = 0;
  *nlast = 0;
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry1 = 2.2250738585072014E-305 / tol;
  yy = z.im;
  zn.re = z.im;
  zn.im = -z.re;
  zb_re = z.re;
  zb_im = z.im;
  cid_im = -1;
  if (fnu < 0.0) {
    ffnu = std::ceil(fnu);
  } else {
    ffnu = std::floor(fnu);
  }

  inu = (int)ffnu - 1;
  ffnu = 1.5707963267948966 * (fnu - ffnu);
  c2_re = std::cos(ffnu);
  c2_im = std::sin(ffnu);
  zar_re = c2_re;
  zar_im = c2_im;
  in = inu + n;
  in -= (in >> 2) << 2;
  in++;
  ffnu = c2_re;
  c2_re = c2_re * (double)icv0[in - 1].re - c2_im * (double)icv0[in - 1].im;
  c2_im = ffnu * (double)icv0[in - 1].im + c2_im * (double)icv0[in - 1].re;
  if (z.im <= 0.0) {
    zn.re = -z.im;
    zn.im = -z.re;
    zb_re = z.re;
    zb_im = -z.im;
    cid_im = 1;
    c2_im = -c2_im;
  }

  if (fnu > 1.0) {
    ffnu = fnu;
  } else {
    ffnu = 1.0;
  }

  b_cunhj(zn, ffnu, 1, tol, &dai, &s1, &zeta1, &zeta2, &ai, &unusedU7);
  if (kode == 1) {
    s1.re = zeta2.re - zeta1.re;
  } else {
    br = zb_re + zeta2.re;
    bi = zb_im + zeta2.im;
    if (bi == 0.0) {
      ffnu = fnu / br;
    } else if (br == 0.0) {
      if (fnu == 0.0) {
        ffnu = 0.0 / bi;
      } else {
        ffnu = 0.0;
      }
    } else {
      brm = std::abs(br);
      ffnu = std::abs(bi);
      if (brm > ffnu) {
        brm = bi / br;
        ffnu = (fnu + brm * 0.0) / (br + brm * bi);
      } else if (ffnu == brm) {
        if (br > 0.0) {
          b_br = 0.5;
        } else {
          b_br = -0.5;
        }

        if (bi > 0.0) {
          b_bi = 0.5;
        } else {
          b_bi = -0.5;
        }

        ffnu = (fnu * b_br + 0.0 * b_bi) / brm;
      } else {
        brm = br / bi;
        ffnu = brm * fnu / (bi + brm * br);
      }
    }

    s1.re = fnu * ffnu - zeta1.re;
  }

  rs1 = s1.re;
  if (std::abs(s1.re) > elim) {
    if (s1.re > 0.0) {
      *nz = -1;
    } else {
      *nz = n;
      i = 1;
      while (i <= n) {
        y->re = 0.0;
        y->im = 0.0;
        i = 2;
      }
    }
  } else {
    do {
      exitg1 = 0;
      goto_mw120 = false;
      in = -1;
      if (2 < n) {
        nn = 2;
      } else {
        nn = n;
      }

      i = 1;
      exitg2 = false;
      while ((!exitg2) && (i <= nn)) {
        fn = fnu + (double)(n - i);
        b_cunhj(zn, fn, 0, tol, &phi, &unusedU7, &zeta1, &zeta2, &asum, &bsum);
        if (kode == 1) {
          s1.re = zeta2.re - zeta1.re;
          s1.im = zeta2.im - zeta1.im;
        } else {
          br = zb_re + zeta2.re;
          bi = zb_im + zeta2.im;
          if (bi == 0.0) {
            rs1 = fn / br;
            ffnu = 0.0;
          } else if (br == 0.0) {
            if (fn == 0.0) {
              rs1 = 0.0 / bi;
              ffnu = 0.0;
            } else {
              rs1 = 0.0;
              ffnu = -(fn / bi);
            }
          } else {
            brm = std::abs(br);
            ffnu = std::abs(bi);
            if (brm > ffnu) {
              brm = bi / br;
              ffnu = br + brm * bi;
              rs1 = (fn + brm * 0.0) / ffnu;
              ffnu = (0.0 - brm * fn) / ffnu;
            } else if (ffnu == brm) {
              if (br > 0.0) {
                br = 0.5;
              } else {
                br = -0.5;
              }

              if (bi > 0.0) {
                ffnu = 0.5;
              } else {
                ffnu = -0.5;
              }

              rs1 = (fn * br + 0.0 * ffnu) / brm;
              ffnu = (0.0 * br - fn * ffnu) / brm;
            } else {
              brm = br / bi;
              ffnu = bi + brm * br;
              rs1 = brm * fn / ffnu;
              ffnu = (brm * 0.0 - fn) / ffnu;
            }
          }

          s1.re = (fn * rs1 - zeta1.re) + std::abs(yy) * 0.0;
          s1.im = (fn * ffnu - zeta1.im) + std::abs(yy);
        }

        rs1 = s1.re;
        if (std::abs(s1.re) > elim) {
          goto_mw120 = true;
          exitg2 = true;
        } else {
          if (i == 1) {
            in = 1;
          }

          guard1 = false;
          if (std::abs(s1.re) >= alim) {
            rs1 = ((s1.re + std::log(rt_hypotd_snf(phi.re, phi.im))) - 0.25 *
                   std::log(rt_hypotd_snf(unusedU7.re, unusedU7.im))) -
              1.2655121234846454;
            if (std::abs(rs1) > elim) {
              goto_mw120 = true;
              exitg2 = true;
            } else {
              if (i == 1) {
                in = 0;
              }

              if ((rs1 >= 0.0) && (i == 1)) {
                in = 2;
              }

              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            ai = cairy(unusedU7, 0, 2);
            dai = cairy(unusedU7, 1, 2);
            ffnu = (ai.re * asum.re - ai.im * asum.im) + (dai.re * bsum.re -
              dai.im * bsum.im);
            brm = (ai.re * asum.im + ai.im * asum.re) + (dai.re * bsum.im +
              dai.im * bsum.re);
            bi = ffnu * phi.re - brm * phi.im;
            brm = ffnu * phi.im + brm * phi.re;
            ffnu = std::exp(s1.re) * cssr[in] * std::cos(s1.im);
            br = std::exp(s1.re) * cssr[in] * std::sin(s1.im);
            dai.re = bi * ffnu - brm * br;
            dai.im = bi * br + brm * ffnu;
            if ((in + 1 == 1) && (cuchk(dai, bry1, tol) != 0)) {
              goto_mw120 = true;
              exitg2 = true;
            } else {
              if (yy <= 0.0) {
                dai.im = -dai.im;
              }

              y->re = csrr[in] * (dai.re * c2_re - dai.im * c2_im);
              y->im = csrr[in] * (dai.re * c2_im + dai.im * c2_re);
              ffnu = c2_re;
              c2_re = c2_re * 0.0 - c2_im * (double)cid_im;
              c2_im = ffnu * (double)cid_im + c2_im * 0.0;
              i++;
            }
          }
        }
      }

      if (!goto_mw120) {
        exitg1 = 1;
      } else if (rs1 > 0.0) {
        *nz = -1;
        exitg1 = 1;
      } else {
        y->re = 0.0;
        y->im = 0.0;
        (*nz)++;
        n--;
        if (n == 0) {
          exitg1 = 1;
        } else {
          in = cuoik(z, fnu, kode, 1, n, y, tol, elim, alim);
          if (in < 0) {
            *nz = -1;
            exitg1 = 1;
          } else {
            n -= in;
            *nz += in;
            if (n == 0) {
              exitg1 = 1;
            } else if ((fnu + (double)n) - 1.0 < fnul) {
              *nlast = n;
              exitg1 = 1;
            } else {
              in = inu + n;
              nn = icv0[in - ((in >> 2) << 2)].re;
              in = icv0[in - ((in >> 2) << 2)].im;
              c2_re = zar_re * (double)nn - zar_im * (double)in;
              c2_im = zar_re * (double)in + zar_im * (double)nn;
              if (yy <= 0.0) {
                c2_im = -c2_im;
              }
            }
          }
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : const creal_T zr
//                double fnu
//                int ikflg
//                int ipmtr
//                double tol
//                int init
//                creal_T *phi
//                creal_T *zeta1
//                creal_T *zeta2
// Return Type  : void
//
static void cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double tol,
                  int init, creal_T *phi, creal_T *zeta1, creal_T *zeta2)
{
  int i;
  double rfn;
  creal_T cwrk[16];
  double ac;
  boolean_T guard1 = false;
  creal_T t;
  double s_re;
  double s_im;
  creal_T sr;
  double t_re;
  double cfn_im;
  double crfn_im;
  double crfn_re;
  double cfn_re;
  int l;
  int exitg1;
  int j;
  for (i = 0; i < 16; i++) {
    cwrk[i].re = 0.0;
    cwrk[i].im = 0.0;
  }

  rfn = 1.0 / fnu;
  ac = fnu * 2.2250738585072014E-305;
  guard1 = false;
  if (init == 0) {
    if ((std::abs(zr.re) > ac) || (std::abs(zr.im) > ac)) {
      t.re = zr.re * rfn - zr.im * 0.0;
      t.im = zr.re * 0.0 + zr.im * rfn;
      s_re = (t.re * t.re - t.im * t.im) + 1.0;
      s_im = t.re * t.im + t.im * t.re;
      sr.re = s_re;
      sr.im = s_im;
      c_sqrt(&sr);
      t_re = t.re;
      ac = 1.0 + sr.re;
      cfn_im = sr.im;
      if (t.im == 0.0) {
        if (cfn_im == 0.0) {
          t.re = ac / t.re;
          t.im = 0.0;
        } else if (ac == 0.0) {
          t.re = 0.0;
          t.im = cfn_im / t_re;
        } else {
          t.re = ac / t.re;
          t.im = cfn_im / t_re;
        }
      } else if (t.re == 0.0) {
        if (ac == 0.0) {
          t.re = cfn_im / t.im;
          t.im = 0.0;
        } else if (cfn_im == 0.0) {
          t.re = 0.0;
          t.im = -(ac / t.im);
        } else {
          t.re = cfn_im / t.im;
          t.im = -(ac / t.im);
        }
      } else {
        crfn_im = std::abs(t.re);
        t_re = std::abs(t.im);
        if (crfn_im > t_re) {
          t_re = t.im / t.re;
          crfn_re = t.re + t_re * t.im;
          t.re = (ac + t_re * cfn_im) / crfn_re;
          t.im = (cfn_im - t_re * ac) / crfn_re;
        } else if (t_re == crfn_im) {
          if (t.re > 0.0) {
            t_re = 0.5;
          } else {
            t_re = -0.5;
          }

          if (t.im > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          t.re = (ac * t_re + cfn_im * crfn_re) / crfn_im;
          t.im = (cfn_im * t_re - ac * crfn_re) / crfn_im;
        } else {
          t_re = t.re / t.im;
          crfn_re = t.im + t_re * t.re;
          t.re = (t_re * ac + cfn_im) / crfn_re;
          t.im = (t_re * cfn_im - ac) / crfn_re;
        }
      }

      c_log(&t);
      zeta1->re = fnu * t.re - 0.0 * t.im;
      zeta1->im = fnu * t.im + 0.0 * t.re;
      zeta2->re = fnu * sr.re - 0.0 * sr.im;
      zeta2->im = fnu * sr.im + 0.0 * sr.re;
      if (sr.im == 0.0) {
        sr.re = rfn / sr.re;
        sr.im = 0.0;
      } else if (sr.re == 0.0) {
        if (rfn == 0.0) {
          sr.re = 0.0 / sr.im;
          sr.im = 0.0;
        } else {
          sr.re = 0.0;
          sr.im = -(rfn / sr.im);
        }
      } else {
        crfn_im = std::abs(sr.re);
        t_re = std::abs(sr.im);
        if (crfn_im > t_re) {
          t_re = sr.im / sr.re;
          crfn_re = sr.re + t_re * sr.im;
          sr.re = (rfn + t_re * 0.0) / crfn_re;
          sr.im = (0.0 - t_re * rfn) / crfn_re;
        } else if (t_re == crfn_im) {
          if (sr.re > 0.0) {
            t_re = 0.5;
          } else {
            t_re = -0.5;
          }

          if (sr.im > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          sr.re = (rfn * t_re + 0.0 * crfn_re) / crfn_im;
          sr.im = (0.0 * t_re - rfn * crfn_re) / crfn_im;
        } else {
          t_re = sr.re / sr.im;
          crfn_re = sr.im + t_re * sr.re;
          sr.re = t_re * rfn / crfn_re;
          sr.im = (t_re * 0.0 - rfn) / crfn_re;
        }
      }

      cwrk[15] = sr;
      c_sqrt(&cwrk[15]);
      phi->re = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].re;
      phi->im = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].im;
      if (ipmtr != 0) {
      } else {
        if (s_im == 0.0) {
          cfn_re = 1.0 / s_re;
          cfn_im = 0.0;
        } else if (s_re == 0.0) {
          cfn_re = 0.0;
          cfn_im = -(1.0 / s_im);
        } else {
          crfn_im = std::abs(s_re);
          t_re = std::abs(s_im);
          if (crfn_im > t_re) {
            t_re = s_im / s_re;
            crfn_re = s_re + t_re * s_im;
            cfn_re = (1.0 + t_re * 0.0) / crfn_re;
            cfn_im = (0.0 - t_re) / crfn_re;
          } else if (t_re == crfn_im) {
            if (s_re > 0.0) {
              t_re = 0.5;
            } else {
              t_re = -0.5;
            }

            if (s_im > 0.0) {
              crfn_re = 0.5;
            } else {
              crfn_re = -0.5;
            }

            cfn_re = (t_re + 0.0 * crfn_re) / crfn_im;
            cfn_im = (0.0 * t_re - crfn_re) / crfn_im;
          } else {
            t_re = s_re / s_im;
            crfn_re = s_im + t_re * s_re;
            cfn_re = t_re / crfn_re;
            cfn_im = (t_re * 0.0 - 1.0) / crfn_re;
          }
        }

        cwrk[0].re = 1.0;
        cwrk[0].im = 0.0;
        crfn_re = 1.0;
        crfn_im = 0.0;
        ac = 1.0;
        l = 0;
        i = 1;
        do {
          exitg1 = 0;
          if (i < 15) {
            s_re = 0.0;
            s_im = 0.0;
            for (j = 0; j <= i; j++) {
              l++;
              t_re = s_re;
              s_re = s_re * cfn_re - s_im * cfn_im;
              s_im = t_re * cfn_im + s_im * cfn_re;
              s_re += dv1[l];
            }

            t_re = crfn_re;
            crfn_re = crfn_re * sr.re - crfn_im * sr.im;
            crfn_im = t_re * sr.im + crfn_im * sr.re;
            cwrk[i].re = crfn_re * s_re - crfn_im * s_im;
            cwrk[i].im = crfn_re * s_im + crfn_im * s_re;
            ac *= rfn;
            if ((ac < tol) && (std::abs(cwrk[i].re) + std::abs(cwrk[i].im) < tol))
            {
              guard1 = true;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            guard1 = true;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }
    } else {
      zeta1->re = 1402.9773265065639 + fnu;
      zeta1->im = 0.0;
      zeta2->re = fnu;
      zeta2->im = 0.0;
      phi->re = 1.0;
      phi->im = 0.0;
    }
  } else {
    zeta1->re = 0.0;
    zeta1->im = 0.0;
    zeta2->re = 0.0;
    zeta2->im = 0.0;
    guard1 = true;
  }

  if (guard1) {
    if (ikflg == 2) {
      phi->re = 1.2533141373155003 * cwrk[15].re;
      phi->im = 1.2533141373155003 * cwrk[15].im;
    } else {
      phi->re = 0.3989422804014327 * cwrk[15].re;
      phi->im = 0.3989422804014327 * cwrk[15].im;
    }
  }
}

//
// Arguments    : const creal_T z
//                double fnu
//                int kode
//                int ikflg
//                int nin
//                creal_T *y
//                double tol
//                double elim
//                double alim
// Return Type  : int
//
static int cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin,
                 creal_T *y, double tol, double elim, double alim)
{
  int nuf;
  int n;
  creal_T zr;
  int iform;
  double gnn;
  double gnu;
  double aarg;
  creal_T arg;
  creal_T an;
  creal_T phi;
  creal_T cz;
  creal_T zeta2;
  boolean_T guard1 = false;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nuf = 0;
  if (z.re < 0.0) {
    zr.re = -z.re;
    zr.im = -z.im;
  } else {
    zr = z;
  }

  if (std::abs(zr.im) > std::abs(z.re) * 1.7321) {
    iform = 2;
  } else {
    iform = 1;
  }

  if (ikflg == 1) {
    if (fnu < 1.0) {
      gnu = 1.0;
    } else {
      gnu = fnu;
    }
  } else {
    gnn = (fnu + (double)n) - 1.0;
    gnu = n;
    if (gnn > gnu) {
      gnu = gnn;
    }
  }

  aarg = 0.0;
  if (iform == 2) {
    if (zr.im <= 0.0) {
      an.re = -zr.im;
      an.im = -zr.re;
    } else {
      an.re = zr.im;
      an.im = -zr.re;
    }

    cunhj(an, gnu, tol, &phi, &arg, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
    aarg = rt_hypotd_snf(arg.re, arg.im);
  } else {
    arg.re = 0.0;
    arg.im = 0.0;
    cunik(zr, gnu, ikflg, 1, tol, 0, &phi, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
  }

  if (kode == 2) {
    cz.re -= zr.re;
    cz.im -= zr.im;
  }

  if (ikflg == 2) {
    cz.re = -cz.re;
    cz.im = -cz.im;
  }

  gnn = rt_hypotd_snf(phi.re, phi.im);
  if (cz.re >= alim) {
    gnn = cz.re + std::log(gnn);
    if (iform == 2) {
      gnn = (gnn - 0.25 * std::log(aarg)) - 1.2655121234846454;
    }

    if (gnn > elim) {
      nuf = -1;
    }
  } else {
    guard1 = false;
    if (cz.re >= -elim) {
      if (cz.re > -alim) {
      } else {
        gnn = cz.re + std::log(gnn);
        if (iform == 2) {
          gnn = (gnn - 0.25 * std::log(aarg)) - 1.2655121234846454;
        }

        if (gnn > -elim) {
          c_log(&phi);
          cz.im += phi.im;
          if (iform == 2) {
            c_log(&arg);
            cz.im -= 0.25 * arg.im;
          }

          gnn = std::exp(gnn) / tol;
          zr.re = gnn * std::cos(cz.im);
          zr.im = gnn * std::sin(cz.im);
          if (cuchk(zr, 2.2250738585072014E-305 / tol, tol) != 1) {
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      iform = 1;
      while (iform <= n) {
        y->re = 0.0;
        y->im = 0.0;
        iform = 2;
      }

      nuf = n;
    }
  }

  return nuf;
}

//
// Arguments    : double *x
// Return Type  : void
//
static void gammaln(double *x)
{
  double t;
  double r;
  static const double table100[100] = { 0.0, 0.0, 0.69314718055994529,
    1.791759469228055, 3.1780538303479458, 4.7874917427820458,
    6.5792512120101012, 8.5251613610654147, 10.604602902745251,
    12.801827480081469, 15.104412573075516, 17.502307845873887,
    19.987214495661885, 22.552163853123425, 25.19122118273868, 27.89927138384089,
    30.671860106080672, 33.505073450136891, 36.395445208033053,
    39.339884187199495, 42.335616460753485, 45.380138898476908,
    48.471181351835227, 51.606675567764377, 54.784729398112319,
    58.003605222980518, 61.261701761002, 64.557538627006338, 67.88974313718154,
    71.257038967168015, 74.658236348830158, 78.0922235533153, 81.557959456115043,
    85.054467017581516, 88.580827542197682, 92.1361756036871, 95.7196945421432,
    99.330612454787428, 102.96819861451381, 106.63176026064346,
    110.32063971475739, 114.03421178146171, 117.77188139974507,
    121.53308151543864, 125.3172711493569, 129.12393363912722,
    132.95257503561632, 136.80272263732635, 140.67392364823425,
    144.5657439463449, 148.47776695177302, 152.40959258449735, 156.3608363030788,
    160.3311282166309, 164.32011226319517, 168.32744544842765,
    172.35279713916279, 176.39584840699735, 180.45629141754378,
    184.53382886144948, 188.6281734236716, 192.7390472878449, 196.86618167289,
    201.00931639928152, 205.1681994826412, 209.34258675253685,
    213.53224149456327, 217.73693411395422, 221.95644181913033,
    226.1905483237276, 230.43904356577696, 234.70172344281826,
    238.97838956183432, 243.26884900298271, 247.57291409618688,
    251.89040220972319, 256.22113555000954, 260.56494097186322,
    264.92164979855278, 269.29109765101981, 273.67312428569369,
    278.06757344036612, 282.4742926876304, 286.893133295427, 291.32395009427029,
    295.76660135076065, 300.22094864701415, 304.68685676566872,
    309.1641935801469, 313.65282994987905, 318.1526396202093, 322.66349912672615,
    327.1852877037752, 331.71788719692847, 336.26118197919845, 340.815058870799,
    345.37940706226686, 349.95411804077025, 354.53908551944079,
    359.1342053695754 };

  int i;
  static const double p1[8] = { 4.9452353592967269, 201.8112620856775,
    2290.8383738313464, 11319.672059033808, 28557.246356716354,
    38484.962284437934, 26377.487876241954, 7225.8139797002877 };

  static const double p2[8] = { 4.974607845568932, 542.4138599891071,
    15506.938649783649, 184793.29044456323, 1.0882047694688288E+6,
    3.33815296798703E+6, 5.1066616789273527E+6, 3.0741090548505397E+6 };

  static const double q1[8] = { 67.482125503037778, 1113.3323938571993,
    7738.7570569353984, 27639.870744033407, 54993.102062261576,
    61611.221800660023, 36351.2759150194, 8785.5363024310136 };

  static const double q2[8] = { 183.03283993705926, 7765.0493214450062,
    133190.38279660742, 1.1367058213219696E+6, 5.2679641174379466E+6,
    1.3467014543111017E+7, 1.7827365303532742E+7, 9.5330955918443538E+6 };

  static const double p4[8] = { 14745.0216605994, 2.4268133694867045E+6,
    1.2147555740450932E+8, 2.6634324496309772E+9, 2.9403789566345539E+10,
    1.7026657377653989E+11, 4.926125793377431E+11, 5.6062518562239514E+11 };

  static const double c[7] = { -0.001910444077728, 0.00084171387781295,
    -0.00059523799130430121, 0.0007936507935003503, -0.0027777777777776816,
    0.083333333333333329, 0.0057083835261 };

  static const double q4[8] = { 2690.5301758708993, 639388.56543000927,
    4.1355999302413881E+7, 1.120872109616148E+9, 1.4886137286788137E+10,
    1.0168035862724382E+11, 3.4174763455073773E+11, 4.4631581874197131E+11 };

  if ((!rtIsNaN(*x)) && (!(*x < 0.0))) {
    if (*x > 2.55E+305) {
      *x = rtInf;
    } else if (*x <= 2.2204460492503131E-16) {
      *x = -std::log(*x);
    } else if (*x <= 0.5) {
      t = 1.0;
      r = 0.0;
      for (i = 0; i < 8; i++) {
        r = r * *x + p1[i];
        t = t * *x + q1[i];
      }

      *x = -std::log(*x) + *x * (-0.57721566490153287 + *x * (r / t));
    } else if (*x <= 0.6796875) {
      t = 1.0;
      r = 0.0;
      for (i = 0; i < 8; i++) {
        r = r * ((*x - 0.5) - 0.5) + p2[i];
        t = t * ((*x - 0.5) - 0.5) + q2[i];
      }

      *x = -std::log(*x) + ((*x - 0.5) - 0.5) * (0.42278433509846713 + ((*x -
        0.5) - 0.5) * (r / t));
    } else if ((*x == std::floor(*x)) && (*x <= 100.0)) {
      *x = table100[(int)*x - 1];
    } else if (*x <= 1.5) {
      t = 1.0;
      r = 0.0;
      for (i = 0; i < 8; i++) {
        r = r * ((*x - 0.5) - 0.5) + p1[i];
        t = t * ((*x - 0.5) - 0.5) + q1[i];
      }

      *x = ((*x - 0.5) - 0.5) * (-0.57721566490153287 + ((*x - 0.5) - 0.5) * (r /
        t));
    } else if (*x <= 4.0) {
      t = 1.0;
      r = 0.0;
      for (i = 0; i < 8; i++) {
        r = r * (*x - 2.0) + p2[i];
        t = t * (*x - 2.0) + q2[i];
      }

      *x = (*x - 2.0) * (0.42278433509846713 + (*x - 2.0) * (r / t));
    } else if (*x <= 12.0) {
      t = -1.0;
      r = 0.0;
      for (i = 0; i < 8; i++) {
        r = r * (*x - 4.0) + p4[i];
        t = t * (*x - 4.0) + q4[i];
      }

      *x = 1.791759469228055 + (*x - 4.0) * (r / t);
    } else {
      if (*x <= 2.25E+76) {
        r = 0.0057083835261;
        t = 1.0 / (*x * *x);
        for (i = 0; i < 6; i++) {
          r = r * t + c[i];
        }

        r /= *x;
      } else {
        r = 0.0;
      }

      t = std::log(*x);
      *x = ((r + 0.91893853320467278) - 0.5 * t) + *x * (t - 1.0);
    }
  }
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2((double)b_u0, (double)b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = 1.0;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

//
// TAUNUM2D
//     TAU = TAUNUM2D(X,Y,X0,Y0)
// Arguments    : double x
//                double y
//                double x0
//                double b_y0
// Return Type  : double
//
static double taunum2D(double x, double y, double x0, double b_y0)
{
  double t2;
  double t3;
  double t4;
  double t5;
  double t11;
  double t8;
  double t10;
  double t14;

  //     This function was generated by the Symbolic Math Toolbox version 8.1.
  //     20-Nov-2019 14:26:33
  t2 = x - x0;
  t3 = y - b_y0;
  t4 = ((((x * x0 * -0.2 - y * b_y0 * 0.6) + x * x * 0.1) + x0 * x0 * 0.1) + y *
        y * 0.3) + b_y0 * b_y0 * 0.3;
  t5 = t2 * t2;
  t2 = t3 * t3;
  t3 = b_y0 * 0.2;
  t11 = y * 0.2;
  t8 = t3 - t11;
  t10 = t5 + t2;
  t14 = 1.0 / (b_y0 * 0.4 - 3.0);
  return std::sqrt((-(t3 - 1.5) * t10 + ((t5 * t2 * 0.01 + t4 * t4 * 0.25) - t8 *
    t8 * (t3 - 1.5) * t10 * t14) / (b_y0 * 0.6 - 4.5)) + (t3 - 1.5) * t10 * t14 *
                   (t3 - t11));
}

//
// G2D
//   2D Green function-based Taylor expansions
// Arguments    : const double z[2]
//                const double z0[2]
//                double w
// Return Type  : creal_T
//
creal_T G2D(const double z[2], const double z0[2], double w)
{
  creal_T out;
  double tau;
  double b_z;
  double c_z;
  double y;
  double Y0;
  creal_T zd;
  int ierr;
  creal_T b_w;
  int unusedU0;
  creal_T c_w;
  double re;
  double im;
  double b_re;

  //  taunum2D(0,0,1e-6,1e-6)
  //  f(-(m-2)/2,w,taunum2D(0,0,1e-6,1e-6))
  //  f(1-(m-2)/2,w,taunum2D(0,0,1e-6,1e-6))
  // V0NUM2D
  //     V0 = V0NUM2D(X,Y,X0,Y0)
  //     This function was generated by the Symbolic Math Toolbox version 8.1.
  //     20-Nov-2019 14:26:33
  tau = taunum2D(z[0], z[1], z0[0], z0[1]);
  b_z = w * tau;

  //  define hankel function using besselj and besselh and bessely are not supported in Matlab Coder 
  // v=0
  if (b_z < 8.0) {
    y = b_z * b_z;
    Y0 = (-2.957821389E+9 + y * (7.062834065E+9 + y * (-5.123598036E+8 + y *
            (1.087988129E+7 + y * (-86327.92757 + y * 228.4622733))))) /
      (4.0076544269E+10 + y * (7.452499648E+8 + y * (7.189466438E+6 + y *
         (47447.2647 + y * (226.1030244 + y))))) + 0.636619772 * BesselJ0(b_z) *
      std::log(b_z);
  } else {
    c_z = 8.0 / b_z;
    y = c_z * c_z;
    Y0 = std::sqrt(0.636619772 / b_z) * (std::sin(b_z - 0.785398164) * (1.0 + y *
      (-0.1098628627 + y * (0.2734510407 + y * (-0.2073370639 + y * 0.2093887211))))
      + c_z * std::cos(b_z - 0.785398164) * (-0.1562499995 + y * (0.1430488765 +
      y * (-0.6911147651 + y * (0.7621095161 - y * 0.934935152)))));
  }

  // V1NUM2D
  //     V1 = V1NUM2D(X,Y,X0,Y0)
  //     This function was generated by the Symbolic Math Toolbox version 8.1.
  //     20-Nov-2019 14:26:33
  y = taunum2D(z[0], z[1], z0[0], z0[1]);
  c_z = w * y;

  //  define hankel function using besselj and besselh and bessely are not supported in Matlab Coder 
  //  v is integer
  //  nonzero even integer v will fail here due to csc(v*pi)!
  //      if(round(v)==v && mod(v,2)==0)
  //          v1 = v + 1e-9;
  //      else
  //          v1 = v;
  //      end
  //      out = i * csc(v1*pi)*(exp(-i*pi*v1)*besselj(v1,z)-besselj(-v1,z)); %
  //      this is correct
  zd.re = c_z;
  zd.im = 0.0;
  ierr = 0;
  if (rtIsNaN(c_z)) {
    b_w.re = rtNaN;
  } else {
    cbesj(zd, 1.0, 1, &b_w, &unusedU0, &ierr);
  }

  if (ierr == 5) {
    b_w.re = rtNaN;
  } else {
    if (ierr == 2) {
      b_w.re = rtInf;
    }
  }

  zd.re = c_z;
  zd.im = 0.0;
  ierr = 0;
  if (rtIsNaN(c_z)) {
    c_w.re = rtNaN;
  } else {
    cbesj(zd, 1.0, 1, &c_w, &unusedU0, &ierr);
    c_w.re = -c_w.re;
  }

  if (ierr == 5) {
    c_w.re = rtNaN;
  } else {
    if (ierr == 2) {
      c_w.re = rtInf;
    }
  }

  //  this is wrong if v is even integer, but matlab coder doesn't support negative real order 
  re = rt_powd_snf(2.0 * tau / w, -0.0) * 0.0;
  im = rt_powd_snf(2.0 * tau / w, -0.0) * 0.88622692545275794;
  c_z = BesselJ0(b_z);
  tau = -b_w.re - c_w.re;
  b_z = -1.2246467991473532E-16 * b_w.re;
  b_re = 0.0 * tau - 8.165619676597685E+15 * b_z;
  b_z = 0.0 * b_z + 8.165619676597685E+15 * tau;
  tau = 2.0 * y / w * -1.08531496757392E-16;
  y = 2.0 * y / w * -0.88622692545275794;
  out.re = ((((((z[0] * z0[0] * -0.037612638903183761 - z[1] * z0[1] *
                 0.037612638903183761) + z[0] * z[0] * 0.01880631945159188) +
               z0[0] * z0[0] * 0.01880631945159188) + z[1] * z[1] *
              0.01880631945159188) + z0[1] * z0[1] * 0.01880631945159188) /
            ((z0[1] * 0.8 - 6.0) * (z0[1] * 2.0 - 15.0)) + 0.2820947917738782) *
    (re * c_z - im * Y0) + 1.0 / rt_powd_snf(z0[1] * 2.0 - 15.0, 3.0) *
    -0.47015798628979688 * (tau * b_re - y * b_z);
  out.im = ((((((z[0] * z0[0] * -0.037612638903183761 - z[1] * z0[1] *
                 0.037612638903183761) + z[0] * z[0] * 0.01880631945159188) +
               z0[0] * z0[0] * 0.01880631945159188) + z[1] * z[1] *
              0.01880631945159188) + z0[1] * z0[1] * 0.01880631945159188) /
            ((z0[1] * 0.8 - 6.0) * (z0[1] * 2.0 - 15.0)) + 0.2820947917738782) *
    (re * Y0 + im * c_z) + 1.0 / rt_powd_snf(z0[1] * 2.0 - 15.0, 3.0) *
    -0.47015798628979688 * (tau * b_z + y * b_re);
  return out;
}

//
// Arguments    : void
// Return Type  : void
//
void G2D_initialize()
{
  rt_InitInfAndNaN(8U);
}

//
// Arguments    : void
// Return Type  : void
//
void G2D_terminate()
{
  // (no terminate code required)
}

//
// File trailer for G2D.cpp
//
// [EOF]
//
