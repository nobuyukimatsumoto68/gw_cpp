#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <filesystem>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_spline.h>

#include <omp.h>

#include <highfive/highfive.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

// #include "H5Cpp.h"
// #include <highfive/H5DataSpace.hpp>


#include <Eigen/Dense>




#include "header.h"
#include "obs.h"





int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------

  std::vector<double> Tcparams;
  using Vjk=std::vector<std::vector<std::vector<double>>>;
  Vjk jk_Tcparams_m0p2000, jk_Tcparams_m0p3000, jk_Tcparams_m0p4000;
  {
    const std::string filename = "Tcparams.hdf5";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    Tcparams = f.getDataSet("mean").read<std::vector<double>>();
    jk_Tcparams_m0p2000 = f.getDataSet("jk_m0p2000").read<Vjk>();
    jk_Tcparams_m0p3000 = f.getDataSet("jk_m0p3000").read<Vjk>();
    jk_Tcparams_m0p4000 = f.getDataSet("jk_m0p4000").read<Vjk>();
  }
  std::vector<double> DeltaVparams;
  Vjk jk_DeltaVparams_m0p2000, jk_DeltaVparams_m0p3000, jk_DeltaVparams_m0p4000;
  {
    const std::string filename = "DeltaVparams.hdf5";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    DeltaVparams = f.getDataSet("mean").read<std::vector<double>>();
    jk_DeltaVparams_m0p2000 = f.getDataSet("jk_m0p2000").read<Vjk>();
    jk_DeltaVparams_m0p3000 = f.getDataSet("jk_m0p3000").read<Vjk>();
    jk_DeltaVparams_m0p4000 = f.getDataSet("jk_m0p4000").read<Vjk>();
  }
  std::vector<double> S3params;
  Vjk jk_S3params_m0p2000, jk_S3params_m0p3000, jk_S3params_m0p4000;
  {
    const std::string filename = "S3params.hdf5";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    S3params = f.getDataSet("mean").read<std::vector<double>>();
    jk_S3params_m0p2000 = f.getDataSet("jk_m0p2000").read<Vjk>();
    jk_S3params_m0p3000 = f.getDataSet("jk_m0p3000").read<Vjk>();
    jk_S3params_m0p4000 = f.getDataSet("jk_m0p4000").read<Vjk>();
  }


  using GWParams=Eigen::Vector2d;

  // results_=[];
  // const int NTs=4800;
  // const int NTs=500;
  const int NTs=1000;
  // const double threshold = -0.0005;
  // const double threshold2 = 0.0006;
  const double threshold = -0.0006;
  const double threshold2 = 0.001;
  std::vector<std::vector<double>> results;

  // for(const double MBGeV : std::vector<double>{5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0}){
  for(const double MBGeV : std::vector<double>{10.0, 20.0, 50.0, 100.0, 200.0, 400.0}){
    for(const double MB : std::vector<double>{2.4, 2.8, 3.2, 3.6, 4.0, 4.4}){
      // for(const double MB : std::vector<double>{5.0}){
      const double sqrt_t0_inv_GeV = MBGeV/MB; // # GeV; free parameter
      const double MP = MP_GeV / sqrt_t0_inv_GeV;

      // # vary m=0.2
      std::vector<std::vector<GWParams>> jk_fO_m0p2000;
      std::string mass="0p2000";

      for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
        std::vector<GWParams> jk(nbins);
        for(int ibin=0; ibin<nbins; ibin++){
          F_Tc f_Tc( jk_Tcparams_m0p2000[jdrop][ibin] );
          F_DeltaV f_DV( jk_DeltaVparams_m0p2000[jdrop][ibin] );
          F_S3hat f_S3( jk_S3params_m0p2000[jdrop][ibin] );

          double alpha, vJ, beta_tilde, Tstar, Tc;
          get_gw_params( alpha, vJ, beta_tilde, Tstar, Tc,
                         f_Tc, f_DV, f_S3,
                         MB, MP, NTs, threshold2 );
          assert(Tstar-Tc > threshold);

          double f_m, om_m;
          get_freq_amplitude( f_m, om_m,
                              alpha, vJ, beta_tilde, Tstar, // Tc,
                              sqrt_t0_inv_GeV);
          jk[ibin] << f_m, om_m;
        }
        jk_fO_m0p2000.push_back( jk );
      }

      // # vary m=0.3
      std::vector<std::vector<GWParams>> jk_fO_m0p3000;
      mass="0p3000";

      for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
        std::vector<GWParams> jk(nbins);
        for(int ibin=0; ibin<nbins; ibin++){
          F_Tc f_Tc( jk_Tcparams_m0p3000[jdrop][ibin] );
          F_DeltaV f_DV( jk_DeltaVparams_m0p3000[jdrop][ibin] );
          F_S3hat f_S3( jk_S3params_m0p3000[jdrop][ibin] );

          double alpha, vJ, beta_tilde, Tstar, Tc;
          get_gw_params( alpha, vJ, beta_tilde, Tstar, Tc,
                         f_Tc, f_DV, f_S3,
                         MB, MP, NTs, threshold2 );
          assert(Tstar-Tc > threshold);

          double f_m, om_m;
          get_freq_amplitude( f_m, om_m,
                              alpha, vJ, beta_tilde, Tstar, // Tc,
                              sqrt_t0_inv_GeV);
          jk[ibin] << f_m, om_m;
        }
        jk_fO_m0p3000.push_back( jk );
      }

      // # vary m=0.4
      std::vector<std::vector<GWParams>> jk_fO_m0p4000;
      mass="0p4000";

      for(int jdrop=0; jdrop<nbetas[mass]; jdrop++){
        std::vector<GWParams> jk(nbins);
        for(int ibin=0; ibin<nbins; ibin++){
          F_Tc f_Tc( jk_Tcparams_m0p4000[jdrop][ibin] );
          F_DeltaV f_DV( jk_DeltaVparams_m0p4000[jdrop][ibin] );
          F_S3hat f_S3( jk_S3params_m0p4000[jdrop][ibin] );

          double alpha, vJ, beta_tilde, Tstar, Tc;
          get_gw_params( alpha, vJ, beta_tilde, Tstar, Tc,
                         f_Tc, f_DV, f_S3,
                         MB, MP, NTs, threshold2 );
          assert(Tstar-Tc > threshold);

          double f_m, om_m;
          get_freq_amplitude( f_m, om_m,
                              alpha, vJ, beta_tilde, Tstar, // Tc,
                              sqrt_t0_inv_GeV);
          jk[ibin] << f_m, om_m;
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
      F_DeltaV f_DV( DeltaVparams );
      F_S3hat f_S3( S3params );

      double alpha, vJ, beta_tilde, Tstar, Tc;
      get_gw_params( alpha, vJ, beta_tilde, Tstar, Tc,
                     f_Tc, f_DV, f_S3,
                     MB, MP, NTs );
      assert(Tstar-Tc > threshold);

      double f_m, om_m;
      get_freq_amplitude( f_m, om_m,
                          alpha, vJ, beta_tilde, Tstar, // Tc,
                          sqrt_t0_inv_GeV);

      std::cout << MBGeV << " " << MB << " " << sqrt_t0_inv_GeV << " "
                << f_m << " " << om_m << " "
                << std::sqrt( var_tot[0] ) << " " << std::sqrt( var_tot[1] ) << " "
                << alpha << " " << vJ << " " << " "
                << beta_tilde << " " << Tstar << " " << Tc << std::endl;
    }}
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
