#pragma once



const int nparallel=4;



std::map<std::string, int> nbetas = {
  { "0p2000", 9},
  { "0p3000", 8},
  { "0p4000", 8}
};



const int nbins=40;




struct F_Tc {
  V fp;

  F_Tc(const V fp_)
    : fp(fp_)
  {}

  double operator() (const double MB) const {
    return fitter( MB, fp[0], fp[1] );
  }

  double fitter( const double MB, const double b0, const double b1 ) const {
    return b0+b1*MB;
  }
};



struct F_DeltaV{
  V fp;

  F_DeltaV(const V fp_)
    : fp(fp_)
  {}

  double operator() (const double TmTc, const double MB) const {
    return fitter( TmTc, MB, fp[0], fp[1], fp[2], fp[3] );
  }

  double fitter( const double TmTc, const double MB,
                 const double a0,
                 const double b0,
                 const double c0,
                 const double d0
                 ) const {
    return TmTc * ( a0 + b0*MB ) + TmTc*TmTc * ( c0 + d0*MB );
  }

  double dT( const double TmTc, const double MB ) const {
    const double a0 = fp[0];
    const double b0 = fp[1];
    const double c0 = fp[2];
    const double d0 = fp[3];
    return ( a0 + b0*MB ) + 2.0*TmTc * ( c0 + d0*MB );
  }
};


struct F_S3hat{
  V fp;

  F_S3hat(const V fp_)
    : fp(fp_)
  {}

  double operator() (const double TmTc, const double MB) const {
    return fitter( TmTc, MB, fp[0], fp[1], fp[2], fp[3], fp[4] );
  }

  double fitter( const double TmTc, const double MB,
                 const double Mc,
                 const double c0,
                 const double c1,
                 const double c2,
                 const double gam ) const {
    assert(MB>Mc);
    return std::pow(MB-Mc, gam) * ( c0 + c1/TmTc + c2/TmTc/TmTc );
  }

  double dT( const double TmTc, const double MB ) const {
    const double Mc = fp[0];
    const double c0 = fp[1];
    const double c1 = fp[2];
    const double c2 = fp[3];
    const double gam = fp[4];
    assert(MB>Mc);
    return std::pow(MB-Mc, gam) * ( -1.0*c1/TmTc/TmTc -2.0 * c2/std::pow(TmTc,3) );
  }
};


// using Func=std::function<double(const double, const double)>;

// template<class F_S3>
struct F_Gammahat {
  const F_S3hat& f_S3;

  F_Gammahat(const F_S3hat& f_S3_)
    : f_S3(f_S3_)
  {}

  double operator()( const double TmTc, const double MB ) const {
    // if(std::abs(TmTc)<1.0e-13) return 0.0;
    const double S3 = f_S3(TmTc, MB);
    // std::cout << "TmTc,MB = " << TmTc << " " << MB << std::endl;
    // std::cout << "S3 = " << S3 << std::endl;
    assert( S3 > 0 );
    if(S3>1.0e20) return 0.0;
    const double coeff = std::pow( S3/(2.0*M_PI), 1.5 );
    const double efactor = std::exp(-S3);
    return coeff * efactor;
  }

  double dT( const double TmTc, const double MB ) const {
    const double S3 =  f_S3(TmTc, MB);
    if(S3>1.0e20) return 0.0;
    const double DS3 = f_S3.dT(TmTc, MB);

    const double coeff = std::pow( S3/(2.0*M_PI), 1.5 );
    const double dcoeff = 1.5*std::pow( S3/(2.0*M_PI), 0.5 )/(2.0*M_PI);
    const double efactor = std::exp(-S3);
    const double defactor = -std::exp(-S3);

    return (dcoeff*efactor + coeff*defactor)*DS3;
  }
};



double getv(const double alpha){
  assert(alpha>0.0);
  const double numer = std::sqrt(2.0*alpha/3.0 + alpha*alpha) + std::sqrt(1.0/3.0);
  const double denom = 1.0+alpha;
  return numer/denom;
}

// double getv(const double alpha){
//   return 1.0;
// }

// double getv(const double alpha){
//   return 0.2;
// }

// double get_alpha(const double T,
//                  const double Tc,
//                  const double MB,
//                  const F_DeltaV& f_DV,
//                  const F_S3hat& f_S3 ){
//   const double DVPrime = f_DV( T-Tc, MB );
//   const double dDVPrime_dT = f_DV.dT( T-Tc, MB );

//   const double numer = T*dDVPrime_dT;
//   const double denom = rhoRprime(T);

//   return (1.0/3.0) * numer/denom;
// }

double get_alpha(const double T,
                 const double Tc,
                 const double MB,
                 const F_DeltaV& f_DV,
                 const F_S3hat& f_S3 ){
  const double DVPrime = f_DV( T-Tc, MB );
  const double dDVPrime_dT = f_DV.dT( T-Tc, MB );

  const double numer = T*dDVPrime_dT;
  const double denom = T*dDVPrime_dT + 4.0*DVPrime;

  return (1.0/3.0) * numer/denom;
}




double HinGeV(const double TinGeV){
  return 1.0/(std::sqrt(3.0)*MP_GeV) * std::sqrt( M_PI*M_PI/30.0 * gstar(TinGeV) );
}



void get_gw_params_n( double& alpha,
                      double& vJ,
                      double& beta_tilde,
                      double& Tstar,
                      double& Tc,
                      const F_Tc& f_Tc,
                      const F_DeltaV& f_DV,
                      const F_S3hat& f_S3,
                      const double MB,
                      const double sqrt_t0_inv_GeV,
                      const double dT=0.000000001, // in sqrt t0
                      const int imax=100000
                      ){
  Tc = f_Tc(MB);

  const F_Gammahat f_Gammahat = F_Gammahat( f_S3 );

  double integral = 0.0;
  double T, TinGeV, H, Gamma;
  for(int iT=0; iT<imax; iT++){
    T = Tc-(iT+0.5)*dT;
    TinGeV = T * sqrt_t0_inv_GeV;

    // H = HinGeV( TinGeV ) / sqrt_t0_inv_GeV;
    // Gamma = f_Gammahat( T-Tc, MB ) * std::pow( T, 4 );

    H = HinGeV( TinGeV );
    Gamma = f_Gammahat( T-Tc, MB ) * std::pow( TinGeV, 4 );

    assert(Gamma>-1.0e-15);
    integral += dT/T * Gamma / std::pow( H, 4 );
    // std::cout << "debug. integral = " << integral << " " << T << std::endl;
    if(integral>1) break;
  }
  Tstar = T;

  alpha = get_alpha(Tstar, Tc, MB, f_DV, f_S3 );
  vJ = getv(alpha);
  beta_tilde = Tstar * f_S3.dT( Tstar-Tc, MB );
}





void get_freq_amplitude( double& fpeak, double& hsq_Omega_peak,
                         const double alpha,
                         const double vJ,
                         const double beta_tilde,
                         const double Tstar,
                         // const double Tc,
                         const double sqrt_t0_inv_GeV){
  const double kappa_v = std::sqrt(alpha)/(0.135+std::sqrt(0.98+alpha));
  const double Ubarfsq = 3.0/4.0 * alpha/(1.0+alpha) * kappa_v;
  const double tau_sw = 1 - 1.0 / std::sqrt( 1.0 + 2.0 * std::pow(8.0*M_PI, 1.0/3.0) * vJ / ( beta_tilde*std::sqrt(Ubarfsq) ) );
  const double kappa_sw = std::sqrt(tau_sw)*kappa_v;

  fpeak = 1.9e-5 * std::pow(gstar(Tstar)/100.,1.0/6.0) * (Tstar*sqrt_t0_inv_GeV/100.) * (beta_tilde / vJ);
  hsq_Omega_peak = 2.65e-6 * (vJ/beta_tilde) * std::pow(kappa_sw*alpha/(1.0+alpha),2) * std::pow(100./gstar(Tstar),1.0/3.0) * std::pow(Om_HSDM(Tstar),2);
}
