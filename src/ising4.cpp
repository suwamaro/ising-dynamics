/***************************************************************************
* Monte Carlo simulations for Ising models.
*
* Copyright (C) 2021 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*****************************************************************************/

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <random>
#include <filesystem>
#include <mpi.h>
#include "ising_config.h"
#include "cpptoml.h"
#include "square_lattice.h"

/* Calculating the effective local field */
double calc_eff_local_field(std::size_t i, std::vector<int> const& spin, Lattice::Lattice* lattice, double J1=0, double J2=0, double J3=0, double Jp=0){
  double hi = 0;

  auto eff_field = [&](std::size_t d, double J){
    double h_d = 0;    
    if ( std::abs(J) > 1e-12 ) {
      for(std::size_t j=0; j < lattice->n_neighbors(d,i); j++){
	std::size_t sj = lattice->get_neighbor(d,i,j);
	h_d += - J * spin[sj];
      }
    }
    return h_d;
  };
  
  /* Nearest neighbors */
  hi += eff_field(0,J1);

  /* Next nearest neighbors */
  hi += eff_field(1, J2);  

  /* Third nearest neighbors */
  hi += eff_field(2, J3);  

  /* Plaquette neighbors */
  if ( std::abs(Jp) > 1e-12 ) {
    for(std::size_t j=0; j < lattice->n_neighbors_plaquette(i); j++){
      int prod = 1;
      for(std::size_t k=0; k < lattice->n_neighbors_plaquette(i,j); k++){      
	std::size_t sk = lattice->get_neighbor_plaquette(i,j,k);
	prod *= spin[sk];
      }
      hi += - Jp * prod;
    }
  }
  
  return hi;
}

void dump_snapshot(std::size_t nstep, std::size_t dump_period, path const& dump_dir, std::vector<int> const& spin){
  int64_t idx = nstep / dump_period;
  std::stringstream fname;
  fname << "dump" << std::setfill('0') << std::setw(4) << std::to_string(idx) << ".json";
  ofstream dump_file(dump_dir/fname.str());
  dump_file << "{\n" << "\"n_steps\":" << nstep << ",\n";

  /* spin */  
  dump_file << "\"spin\":[";
  for(std::size_t s = 0; s < spin.size(); s++){ 	    
    dump_file << spin[s];
    if (s < spin.size() - 1) dump_file << ",";
  }
  dump_file << "]\n";
  dump_file << "}\n";  
  dump_file.close();
}


int main(int argc, char **argv){
  /* Parameters */
  std::string input_file_name = "parameters.toml";
  if (argc != 2) {
    std::cerr << "Input a path to \"" << input_file_name << "\"." << ".\n";    
    std::exit(EXIT_SUCCESS);
  }
  path base_dir(argv[1]);

  /* Input parameters */
  auto input_file = base_dir / input_file_name;  
  if ( !exists(input_file) ) {
    std::cerr << "Required input file " << input_file << " does not exist.\n";
    std::exit(EXIT_FAILURE);
  }
  auto config = cpptoml::parse_file(input_file.string());
  
  /* Reading a toml file */    
  std::size_t L = config->get_as<int64_t>("L").value_or(8);
  double kT = config->get_as<double>("kT").value_or(2.2);
  double J1 = config->get_as<double>("J1").value_or(1.0);
  double J2 = config->get_as<double>("J2").value_or(0.0);
  double J3 = config->get_as<double>("J3").value_or(0.0);
  double Jp = config->get_as<double>("Jp").value_or(0.0);      
  std::size_t n_mcs = config->get_as<int64_t>("n_mcs").value_or(1024);
  bool dump_config = config->get_as<bool>("dump_config").value_or(false);
  std::size_t n_dump = config->get_as<int64_t>("n_dump").value_or(8);  
  std::size_t dump_period = config->get_as<int64_t>("dump_period").value_or(n_mcs/n_dump);
  if ( dump_period == 0 ) { dump_period = 1; }
  std::size_t t_skip = config->get_as<int64_t>("t_skip").value_or(1);  
  std::uint32_t seed = config->get_as<int64_t>("seed").value_or(0);
  std::size_t n_replicas = config->get_as<int64_t>("n_replicas").value_or(1);  
  std::size_t interaction_distance = 3;
  int precision = 10;
  double beta = 1. / kT;
  std::size_t buffer_size = 1024;
    
  /* MPI parallel */  
  int rank, num_procs;  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  /* Output parameters */
  if ( rank == 0 ) {
    std::cout << "kT = " << kT << std::endl;
    std::cout << "Number of MCSs = " << n_mcs << std::endl;
  }

  /* Lattice */
  auto lattice = Lattice::SquareLattice::mk_square(L, interaction_distance);
  std::size_t N = lattice->n_sites();
  
  /* Tasks */
  std::size_t n_tasks = n_replicas / num_procs;
  std::vector<std::size_t> tasks;
  for(std::size_t i=n_tasks*rank; i < n_tasks*(rank+1); i++){ tasks.push_back(i); }
  std::size_t n_rem = n_replicas % num_procs;
  if ( rank < n_rem ) { tasks.push_back(n_replicas - n_rem + rank); }

  /* For each task */
  for(std::size_t task: tasks){
    /* Creating directories */
    path obs_T("observables"+std::to_string(task));
    path obs_dir = base_dir/"observables"/("observables"+std::to_string(task));
    if ( !is_directory(obs_dir) ) { create_directories(obs_dir); }
    path dump_dir = base_dir/"dump"/("dump"+std::to_string(task));
    if ( dump_config && !is_directory(dump_dir) ) { create_directories(dump_dir); }

    /* Spin configuration */
    std::vector<int> spin(N);
    
    /* random number */
    std::mt19937 rng(seed + task);
    
    /* Setting the initial state random. */    
    std::uniform_int_distribution<> uniform_int(0, 1); /* 0 or 1 */    
    for(std::size_t i=0; i < N; i++){
      spin[i] = 2 * uniform_int(rng) - 1;  // -1 (down) or 1 (up)
    }
  
    /* Output files */
    ofstream e_out(obs_dir/"energy.text");
    ofstream m_out(obs_dir/"mag.text");
    ofstream op_stripe_out(obs_dir/"op_stripe.text");
    ofstream op_4by4_out(obs_dir/"op_4by4.text");
    if ( dump_config ) { dump_snapshot(0, dump_period, dump_dir, spin); }

    /* Setting the precision */
    e_out << std::setprecision( precision );
    m_out << std::setprecision( precision );
    op_stripe_out << std::setprecision( precision );
    op_4by4_out << std::setprecision( precision );        

    /* Result */
    std::vector<double> es(buffer_size);
    std::vector<double> ms(buffer_size);
    std::vector<double> op_stripe(buffer_size);
    std::vector<double> op_4by4(buffer_size);          
  
    /* Monte Carlo update */
    for(std::size_t t=0; t < n_mcs; t++){
      /* Random sequential update */
      std::vector<std::size_t> seq(N);
      std::iota(seq.begin(), seq.end(), 0);    
      std::shuffle(seq.begin(), seq.end(), rng);
      
      for(std::size_t r: seq){
	/* Effective local field */
	double hi = calc_eff_local_field(r, spin, lattice.get(), J1, J2, J3, Jp);

	/* Energy difference */
	double e_diff = 2. * hi * spin[r];

	/* Heat bath */
	bool accept = false;
	if ( beta * e_diff > 30. ) {
	  /* Reject */
	} else if ( beta * e_diff < - 30. ) {
	  /* Accept */
	  accept = true;
	} else {
	  double w_diff = std::exp( - beta * e_diff );
	  std::uniform_real_distribution<double> uniform_01(0, 1); /* [ 0, 1 ) */      
	  if ( uniform_01(rng) * (1. + w_diff ) < w_diff ) {
	    accept = true;
	  }
	}
    
	if ( accept ) {
	  spin[r] *= -1;
	}
      } /* end for i */

      /* Measurement */
      if ( t % t_skip == 0 ) {
	/* Calculating the energy, the magnetization, and structure factors */
	double e = 0;
	double m = 0;
	cx_double m_pi_0 = 0;
	cx_double m_0_pi = 0;
	cx_double m_pi2_pi2 = 0;
	std::vector<double> q_pi_0({M_PI, 0});
	std::vector<double> q_0_pi({0, M_PI});
	std::vector<double> q_pi2_pi2({M_PI/2, M_PI/2});
	for(std::size_t i=0; i < N; i++){
	  double hi = calc_eff_local_field(i, spin, lattice.get(), 0.5*J1, 0.5*J2, 0.5*J3, 0.25*Jp);	  
	  double ei = - hi * spin[i];
	  e += ei;
	  m += (double)spin[i];

	  m_pi_0 += (double)spin[i] * lattice->phase(i, q_pi_0);
	  m_0_pi += (double)spin[i] * lattice->phase(i, q_0_pi);
	  m_pi2_pi2 += (double)spin[i] * lattice->phase(i, q_pi2_pi2);
	}

	/* Storing observables */
	std::size_t t_measure = t / t_skip;
	std::size_t index = t_measure % buffer_size;
	es[index] = e / (double)N;
	ms[index] = m / (double)N;
	m_pi_0 /= (double)N;
	m_0_pi /= (double)N;
	m_pi2_pi2 /= (double)N;
	op_stripe[index] = sqrt(std::norm(m_0_pi) + std::norm(m_pi_0));
	op_4by4[index] = sqrt(4.*std::norm(m_pi2_pi2));  // There are 4 equivalent wave vectors.
	
	auto output_result = [&](std::size_t T1, std::size_t T){
	  for(std::size_t ti=0; ti < T; ti++){
	    std::size_t tp = T1 + ti * t_skip;
	    e_out << tp << std::setw(precision + 10) << es[ti] << std::endl;
	    m_out << tp << std::setw(precision + 10) << ms[ti] << std::endl;
	    op_stripe_out << tp << std::setw(precision + 10) << op_stripe[ti] << std::endl;
	    op_4by4_out << tp << std::setw(precision + 10) << op_4by4[ti] << std::endl;	    	    
	  }
	};

	/* Output */      
	if ( index == buffer_size - 1 ) {
	  output_result( t - (buffer_size - 1) * t_skip, buffer_size );
	} else if ( t == n_mcs - 1 ) {
	  std::size_t rem = t_measure % buffer_size;
	  output_result( t - rem * t_skip, rem + 1 );
	}

	/* Dump the snapshot configuration */
	if ( dump_config && (t + 1) % dump_period == 0 ) {
	  dump_snapshot(t+1, dump_period, dump_dir, spin);
	}
      }
    } /* end for t */

    /* Closing the output files */
    e_out.close();
    m_out.close();
    op_stripe_out.close();
    op_4by4_out.close();        
  } /* end for task */

  /* MPI_Finalize */
  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
