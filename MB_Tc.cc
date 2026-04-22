#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <map>
#include <fstream>
#include <filesystem>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_spline.h>

#include <omp.h>

// #include <highfive/highfive.hpp>
#include <highfive/H5File.hpp>
// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>

// #include "H5Cpp.h"
// #include <highfive/H5DataSpace.hpp>


#include <Eigen/Dense>



using V=std::vector<double>;
using Vk=std::vector<V>;
using Vjk=std::vector<Vk>;
using Vjk_=std::vector<std::vector<std::vector<double>>>;

#include "dof.h"

#include "header1.h"
#include "obs.h"


const double epsilon = 0.8;
const double eps = 0.5;
#define is_eps_trick true


void add_to( V& v, const double c0, const V& v0){
  for(int i=0; i<v.size(); i++) v[i] += c0*v0[i];
}

void rescale( V& v, const double& c0){
  for(int i=0; i<v.size(); i++) v[i] *= c0;
}


V get_mean( const Vk& v ){
  V res( v[0].size(), 0.0 );
  for(int i=0; i<v.size(); i++) {
    // res += v[i];
    add_to( res, 1.0, v[i]);
  }
  // res /= v.size();
  rescale( res, 1.0/v.size() );
  return res;
}


Vk eps_trick( const Vk& v,
              const double eps){
  const V mean = get_mean(v);

  Vk tmp = v;
  for(int i=0; i<v.size(); i++) add_to( tmp[i], -1.0, mean ); // tmp[i] -= mean;
  for(int i=0; i<v.size(); i++) rescale( tmp[i], eps );
  for(int i=0; i<v.size(); i++) add_to( tmp[i], 1.0, mean );
  return tmp;
}




int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);


  // std::cout << gstar_SM(10000000) << std::endl;
  // for(double logT=std::log(TQCD); logT<std::log(500); logT+=0.01){
  //   const double T = std::exp(logT);
  //   std::cout << T << " " << gstar_SM(T) << std::endl;
  // }
  // return 1;

  // --------------------
  // std::ofstream ofs_freq_Omega("freq_Omega.dat");
  // ofs_freq_Omega << std::scientific << std::setprecision(25);
  // std::ofstream ofs_alpha_beta("alpha_beta.dat");
  // ofs_alpha_beta << std::scientific << std::setprecision(25);
  // std::ofstream ofs_freq_Omega("freq_Omega.dat");
  // ofs_freq_Omega << std::scientific << std::setprecision(25);


  V Tcparams;
  Vjk jk_Tcparams_m0p2000, jk_Tcparams_m0p3000, jk_Tcparams_m0p4000;
  {
    const std::string filename = "Tcparams.hdf5";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    Tcparams = f.getDataSet("mean").read<V>();
    jk_Tcparams_m0p2000 = f.getDataSet("jk_m0p2000").read<Vjk>();
    jk_Tcparams_m0p3000 = f.getDataSet("jk_m0p3000").read<Vjk>();
    jk_Tcparams_m0p4000 = f.getDataSet("jk_m0p4000").read<Vjk>();
  }
  V DeltaVparams;
  Vjk jk_DeltaVparams_m0p2000, jk_DeltaVparams_m0p3000, jk_DeltaVparams_m0p4000;
  {
    const std::string filename = "DeltaVparams.hdf5";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    DeltaVparams = f.getDataSet("mean").read<V>();
    jk_DeltaVparams_m0p2000 = f.getDataSet("jk_m0p2000").read<Vjk>();
    jk_DeltaVparams_m0p3000 = f.getDataSet("jk_m0p3000").read<Vjk>();
    jk_DeltaVparams_m0p4000 = f.getDataSet("jk_m0p4000").read<Vjk>();
  }
  V S3params;
  Vjk jk_S3params_m0p2000, jk_S3params_m0p3000, jk_S3params_m0p4000;
  {
    const std::string filename = "S3params.hdf5";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    S3params = f.getDataSet("mean").read<V>();
    jk_S3params_m0p2000 = f.getDataSet("jk_m0p2000").read<Vjk>();
    jk_S3params_m0p3000 = f.getDataSet("jk_m0p3000").read<Vjk>();
    jk_S3params_m0p4000 = f.getDataSet("jk_m0p4000").read<Vjk>();
  }


  // using GWParams=Eigen::Vector2d;
  constexpr int nops=2;
  using GWParams=Eigen::Matrix< double, nops, 1 >; // Eigen::VectorXd;

  const int NTs=2000;
  // const int NTs=200;
  const double threshold = -0.0006;
  const double threshold2 = 0.0006;
  std::vector<std::vector<double>> results;

  for(double MB=2.6; MB<4.2; MB+=0.01){
    // const double sqrt_t0_inv_GeV = MBGeV/MB; // # GeV; free parameter
    // const double MP = MP_GeV / sqrt_t0_inv_GeV;

    // # vary m=0.2
    std::vector<std::vector<GWParams>> jk_fO_m0p2000;
    std::string mass="0p2000";

    for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
      std::vector<GWParams> jk(nbins);
# ifdef is_eps_trick
      const Vk jk_Tcparams = eps_trick( jk_Tcparams_m0p2000[jdrop], eps);
      const Vk jk_DeltaVparams = eps_trick( jk_DeltaVparams_m0p2000[jdrop], eps );
      const Vk jk_S3params = eps_trick( jk_S3params_m0p2000[jdrop], eps );
#else
      const Vk& jk_Tcparams = jk_Tcparams_m0p2000[jdrop];
      const Vk& jk_DeltaVparams = jk_DeltaVparams_m0p2000[jdrop];
      const Vk& jk_S3params = jk_S3params_m0p2000[jdrop];
#endif
      for(int ibin=0; ibin<nbins; ibin++){
        F_Tc f_Tc( jk_Tcparams[ibin] );
        double Tc = f_Tc(MB);
        jk[ibin] << Tc, MB/Tc;
      }
      jk_fO_m0p2000.push_back( jk );
    }

    // # vary m=0.3
    std::vector<std::vector<GWParams>> jk_fO_m0p3000;
    mass="0p3000";

    for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
      std::vector<GWParams> jk(nbins);
      // std::vector<GWParams> jk(nbins, GWParams::Zero(nops));
#ifdef is_eps_trick
      const Vk jk_Tcparams = eps_trick( jk_Tcparams_m0p3000[jdrop], eps );
      const Vk jk_DeltaVparams = eps_trick( jk_DeltaVparams_m0p3000[jdrop], eps );
      const Vk jk_S3params = eps_trick( jk_S3params_m0p3000[jdrop], eps );
#else
      const Vk& jk_Tcparams = jk_Tcparams_m0p3000[jdrop];
      const Vk& jk_DeltaVparams = jk_DeltaVparams_m0p3000[jdrop];
      const Vk& jk_S3params = jk_S3params_m0p3000[jdrop];
#endif

      for(int ibin=0; ibin<nbins; ibin++){
        F_Tc f_Tc( jk_Tcparams[ibin] );
        double Tc = f_Tc(MB);
        jk[ibin] << Tc, MB/Tc;
      }
      jk_fO_m0p3000.push_back( jk );
    }

    // # vary m=0.4
    std::vector<std::vector<GWParams>> jk_fO_m0p4000;
    mass="0p4000";

    for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
      std::vector<GWParams> jk(nbins);
      // std::vector<GWParams> jk(nbins, GWParams::Zero(nops));
      // for(int ibin=0; ibin<nbins; ibin++){
      //   F_Tc f_Tc( jk_Tcparams_m0p4000[jdrop][ibin] );
      //   F_DeltaV f_DV( jk_DeltaVparams_m0p4000[jdrop][ibin] );
      //   F_S3hat f_S3( jk_S3params_m0p4000[jdrop][ibin] );
#ifdef is_eps_trick
      const Vk jk_Tcparams = eps_trick( jk_Tcparams_m0p4000[jdrop], eps );
      const Vk jk_DeltaVparams = eps_trick( jk_DeltaVparams_m0p4000[jdrop], eps );
      const Vk jk_S3params = eps_trick( jk_S3params_m0p4000[jdrop], eps );
#else
      const Vk& jk_Tcparams = jk_Tcparams_m0p4000[jdrop];
      const Vk& jk_DeltaVparams = jk_DeltaVparams_m0p4000[jdrop];
      const Vk& jk_S3params = jk_S3params_m0p4000[jdrop];
#endif

      for(int ibin=0; ibin<nbins; ibin++){
        F_Tc f_Tc( jk_Tcparams[ibin] );
        double Tc = f_Tc(MB);
        jk[ibin] << Tc, MB/Tc;
      }
      jk_fO_m0p4000.push_back( jk );
    }


    GWParams var_tot = GWParams::Zero();

    {
      // # vary m=0.2
      mass="0p2000";
      for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
        JackknifeSimp<GWParams> jk( nbins );
        jk.jack_avg = jk_fO_m0p2000[jdrop];
        jk.finalize();
        var_tot += jk.var;
      }
    }
    {
      // # vary m=0.3
      mass="0p3000";
      for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
        JackknifeSimp<GWParams> jk( nbins );
        jk.jack_avg = jk_fO_m0p3000[jdrop];
        jk.finalize();
        var_tot += jk.var;
      }
    }
    {
      // # vary m=0.4
      mass="0p4000";
      for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
        JackknifeSimp<GWParams> jk( nbins );
        jk.jack_avg = jk_fO_m0p4000[jdrop];
        jk.finalize();
        var_tot += jk.var;
      }
    }

    // #######################MEAN

    F_Tc f_Tc( Tcparams );
    double Tc = f_Tc(MB);
    // jk[ibin] << Tc, MB/Tc;

    std::cout << MB << " "
              << MB/Tc << " " // 0, 1, 2
              << std::sqrt( var_tot[1] )/epsilon/eps << " " // 5, 6
              << std::endl;
  }
  return 0;
}



      // // # vary m=0.3
      // std::vector<GWParams> jk_fO_m0p3000;
      // mass="0p3000";

      // for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
      //   GWParams jk_;
      //   for(int ibin=0; ibin<nbins; ibin++){
      //     // # print(jdrop, ibin)
      //     F_Tc f_Tc( jk_Tcparams_m0p3000[jdrop][ibin] );
      //     F_DeltaV f_DV( jk_DeltaVparams_m0p3000[jdrop][ibin] );
      //     F_S3hat f_S3( jk_S3params_m0p3000[jdrop][ibin] );
      //     // std::cout << "debug. f_Tc = " << f_Tc( 3.2) << std::endl;
      //     // std::cout << "debug. f_DV = " << f_DV(-0.000152, 3.2) << std::endl;
      //     // std::cout << "debug. f_S3 = " << f_S3(-0.000152, 3.2) << std::endl;

      //     double alpha, vJ, beta_tilde, Tstar, Tc;
      //     get_gw_params( alpha, vJ, beta_tilde, Tstar, Tc,
      //                    f_Tc, f_DV, f_S3,
      //                    MB, MP, NTs );
      //     // std::cout << alpha << " "
      //     //           << vJ << " "
      //     //           << beta_tilde << " "
      //     //           << Tstar << " "
      //     //           << Tc << std::endl;
      //     // std::cout << Tstar-Tc << std::endl;
      //     assert(Tstar-Tc > -0.00035);
      //     // meanparams = [alpha, vJ, beta_tilde, Tstar, Tc]

      //     double f_m, om_m;
      //     get_freq_amplitude( f_m, om_m,
      //                         alpha, vJ, beta_tilde, Tstar, // Tc,
      //                         sqrt_t0_inv_GeV);
      //     jk_ << f_m, om_m;
      //     // return 1;
      //   }
      //   jk_fO_m0p3000.push_back( jk_ );
      // }

      // // # vary m=0.4
      // std::vector<GWParams> jk_fO_m0p4000;
      // mass="0p4000";

      // for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
      //   GWParams jk_;
      //   for(int ibin=0; ibin<nbins; ibin++){
      //     // # print(jdrop, ibin)
      //     F_Tc f_Tc( jk_Tcparams_m0p4000[jdrop][ibin] );
      //     F_DeltaV f_DV( jk_DeltaVparams_m0p4000[jdrop][ibin] );
      //     F_S3hat f_S3( jk_S3params_m0p4000[jdrop][ibin] );
      //     // std::cout << "debug. f_Tc = " << f_Tc( 3.2) << std::endl;
      //     // std::cout << "debug. f_DV = " << f_DV(-0.000152, 3.2) << std::endl;
      //     // std::cout << "debug. f_S3 = " << f_S3(-0.000152, 3.2) << std::endl;

      //     double alpha, vJ, beta_tilde, Tstar, Tc;
      //     get_gw_params( alpha, vJ, beta_tilde, Tstar, Tc,
      //                    f_Tc, f_DV, f_S3,
      //                    MB, MP, NTs );
      //     // std::cout << alpha << " "
      //     //           << vJ << " "
      //     //           << beta_tilde << " "
      //     //           << Tstar << " "
      //     //           << Tc << std::endl;
      //     // std::cout << Tstar-Tc << std::endl;
      //     assert(Tstar-Tc > -0.00035);
      //     // meanparams = [alpha, vJ, beta_tilde, Tstar, Tc]

      //     double f_m, om_m;
      //     get_freq_amplitude( f_m, om_m,
      //                         alpha, vJ, beta_tilde, Tstar, // Tc,
      //                         sqrt_t0_inv_GeV);
      //     jk_ << f_m, om_m;
      //     // return 1;
      //   }
      //   jk_fO_m0p4000.push_back( jk_ );
      // }
