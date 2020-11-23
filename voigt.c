/* voigt.c: computes exact Voigt profiles in terms of the complex
   error function (requires libcerf) */

#include "mex.h"
#include <cerf.h>
#include <math.h>

#define LAMBDAS_ARG    prhs[0]
#define Z_ARG          prhs[1]
#define N_ARG          prhs[2]
#define NUM_LINES_ARG  prhs[3]

#define PROFILE_ARG    plhs[0]

/* number of lines in the Lyman series to consider */
#define NUM_LINES 31

/* note: all units are CGS */

/* physical constants */

static const double c   = 2.99792458e+10;              /* speed of light          cm s⁻¹        */
/* static const double k   = 1.38064852e-16; */        /* Boltzmann constant      erg K⁻¹       */
/* static const double m_p = 1.672621898e-24;*/        /* proton mass             g             */
/* static const double m_e = 9.10938356e-28; */        /* electron mass           g             */
/* e = 1.6021766208e-19 * c / 10; */
/* static const double e   = 4.803204672997660e-10; */ /* elementary charge       statC         */

/* CIV doublet */

static const double transition_wavelengths[] =         /* transition wavelengths  cm            */
  {
    1.54820490e-05,
    1.55077845e-05
  };

static const double oscillator_strengths[] =           /* oscillator strengths    dimensionless */
  {
    0.189900,
    0.094750
  };

static const double Gammas[] =                         /* transition rates        s^-1          */
  {
    2.642e+08,
    2.628e+08
  };

/* assumed constant */
/* static const double T = 1e+04; */                   /* gas temperature         K             */

/* derived constants */

/* b = sqrt(2 * k * T / m_p); */
/* static const double b =
/*     1.28486551932562422e+06; */                       /* Doppler parameter       cm s⁻¹        */

/* sigma = b / M_SQRT2; */
static const double sigma = 9.08537121627923800e+05;   /* Gaussian width          cm s⁻¹        */

//static const double sigma =   1771064895826412132e+6 /* Gaussian width          cm s⁻¹        */



/* leading_constants[i] =
      M_PI * e * e * oscillator_strengths[i] * transition_wavelengths[i] / (m_e * c) ;
*/
static const double leading_constants[] =              /* leading constants       cm²           */
  {     7.802895118381213e-08,
     3.899701297867750e-08
 };

/* gammas[i] = Gammas[i] * transition_wavelengths[i] / (4 * M_PI); */
static const double gammas[] =                         /* Lorentzian widths       cm s⁻¹        */
  {
     3.255002952981575e+02,
     3.243136695286643e+02
  };

/* BOSS spectrograph instrumental broadening */
/* R = 2000; */                                  /* resolving power          dimensionless */

/* width of instrument broadening in pixels */
/* pixel_spacing = 1e-4; */
/* pixel_sigma = 1 / (R * 2 * sqrt(2 * M_LN2) * (pow(10, pixel_spacing) - 1)); */

static const int width = 3;                      /* width of convolution     dimensionless */

/*
  total = 0;
  for (i = -width, j = 0; i <= width; i++, j++) {
    instrument_profile[j] = exp(-0.5 * i * i / (pixel_sigma * pixel_sigma));
    total += instrument_profile[j];
  }

  for (i = 0; i < 2 * width + 1; i++)
    instrument_profile[i] /= total;
*/

static const double instrument_profile[] =
  {
   0.002174609921381,
   0.041162305958045,
   0.240309364651847,
   0.432707438937454,
   0.240309364651847,
   0.041162305958045,
   0.002174609921381
  };

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {

  double *lambdas, *profile, *multipliers, *raw_profile;
  double z, N, velocity, total;
  int num_lines, i, j, k;
  mwSize num_points;

  /* get input */
  lambdas = mxGetPr(LAMBDAS_ARG);                /* wavelengths             Å             */
  z       = mxGetScalar(Z_ARG);                  /* redshift                dimensionless */
  N       = mxGetScalar(N_ARG);                  /* column density          cm⁻²          */

  num_lines = (nrhs > 3) ? (int)(mxGetScalar(NUM_LINES_ARG)) : NUM_LINES;

  num_points = mxGetNumberOfElements(LAMBDAS_ARG);

  /* initialize output */
  PROFILE_ARG = mxCreateDoubleMatrix(num_points - 2 * width, 1, mxREAL);
  profile = mxGetPr(PROFILE_ARG);                /* absorption profile      dimensionless */

  /* to hold the profile before instrumental broadening */
  raw_profile = mxMalloc(num_points * sizeof(double));

  multipliers = mxMalloc(num_lines * sizeof(double));
  for (i = 0; i < num_lines; i++)
    multipliers[i] = c / (transition_wavelengths[i] * (1 + z)) / 1e8;

  /* compute raw Voigt profile */
  for (i = 0; i < num_points; i++) {
    /* apply each absorption line */
    total = 0;
    for (j = 0; j < num_lines; j++) {
      /* velocity relative to transition wavelength */
      velocity = lambdas[i] * multipliers[j] - c;
      total += -leading_constants[j] * voigt(velocity, sigma, gammas[j]);
    }

    raw_profile[i] = exp(N * total);
  }

  num_points = mxGetNumberOfElements(PROFILE_ARG);

  /* instrumental broadening */
  for (i = 0; i < num_points; i++)
    for (j = i, k = 0; j <= i + 2 * width; j++, k++)
      profile[i] += raw_profile[j] * instrument_profile[k];

  mxFree(raw_profile);
  mxFree(multipliers);

}
