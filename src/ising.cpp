#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <random>
#include <filesystem>
#include <mpi.h>
#include "cpptoml.h"

/* Using std::filesystem */
typedef std::ofstream ofstream;
typedef std::filesystem::path path;

// Square lattice ferromagnetic Ising model
std::tuple<std::size_t, std::size_t> site_to_coordinate(std::size_t i, std::size_t L){
  std::size_t x = i % L;
  std::size_t y = i / L;
  return std::make_tuple(x, y);
}

std::size_t coordinate_to_site(std::size_t L, std::size_t x, std::size_t y){
  return x + y * L;
}

std::vector<std::size_t> neighbors(std::size_t i, std::size_t L){
  std::vector<std::size_t> v;
  std::size_t x, y;
  std::tie(x, y) = site_to_coordinate(i, L);  
  std::size_t x1 = (x + L - 1) % L;
  std::size_t x2 = (x + 1 ) % L;
  std::size_t y1 = (y + L - 1) % L;
  std::size_t y2 = (y + 1 ) % L;
  
  v.push_back( coordinate_to_site(L, x2, y) );
  v.push_back( coordinate_to_site(L, x, y2) );
  v.push_back( coordinate_to_site(L, x1, y) );
  v.push_back( coordinate_to_site(L, x, y1) );
  return v;
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
  double J = config->get_as<double>("J").value_or(1.0);
  std::size_t n_mcs = config->get_as<int64_t>("n_mcs").value_or(1024);
  std::size_t t_skip = config->get_as<int64_t>("t_skip").value_or(1);  
  std::uint32_t seed = config->get_as<int64_t>("seed").value_or(0);
  std::size_t n_replicas = config->get_as<int64_t>("n_replicas").value_or(1);  

  int precision = 10;
  std::size_t N = L * L;
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
  
  /* Tasks */
  std::size_t n_tasks = n_replicas / num_procs;
  std::vector<std::size_t> tasks;
  for(std::size_t i=n_tasks*rank; i < n_tasks*(rank+1); i++){ tasks.push_back(i); }
  std::size_t n_rem = n_replicas % num_procs;
  if ( rank < n_rem ) { tasks.push_back(n_replicas - n_rem + rank); }

  /* For each task */
  for(std::size_t task: tasks){
    /* Creating a directory */
    path obs_T("observables"+std::to_string(task));
    path obs_dir = base_dir/"observables"/("observables"+std::to_string(task));
    if ( !is_directory(obs_dir) ) { create_directories(obs_dir); }
  
    /* random number */
    std::mt19937 rng(seed + task);
    std::uniform_real_distribution<double> uniform_01( 0, 1 ); /* [ 0, 1 ) */
    
    /* Random spin configuration */
    std::vector<int> spin(N);  // -1 (down) or 1 (up)
    for(std::size_t i=0; i < N; i++){
      spin[i] = uniform_01(rng) < 0.5 ? -1 : 1;
    }
  
    /* Output files */
    ofstream e_out(obs_dir/"energy.text");
    ofstream m_out(obs_dir/"mag.text");

    /* Setting the precision */
    e_out << std::setprecision( precision );
    m_out << std::setprecision( precision );

    /* Result */
    std::vector<double> es(buffer_size);
    std::vector<double> ms(buffer_size);  
  
    /* Monte Carlo update */
    for(std::size_t t=0; t < n_mcs; t++){
      /* Random sequential update */
      std::vector<std::size_t> seq(N);
      std::iota(seq.begin(), seq.end(), 0);    
      std::shuffle(seq.begin(), seq.end(), rng);

      for(std::size_t i=0; i < N; i++){
	std::size_t r = seq[i];

	/* Effective local field */
	std::vector<std::size_t> js = neighbors(r, L);
	int hi = 0;
	for(std::size_t j=0; j < js.size(); j++){
	  std::size_t sj = js[j];
	  hi += spin[sj];
	}

	/* Energy difference */
	double e_diff = 2. * J * hi * spin[r];

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
	/* Calculating the energy and the magnetization */
	double e = 0;
	double m = 0;
	for(std::size_t i=0; i < N; i++){
	  std::vector<std::size_t> js = neighbors(i, L);
	  int hi = 0;
	  for(std::size_t j=0; j < js.size(); j++){
	    std::size_t sj = js[j];
	    hi += spin[sj];
	  }
	  double ei = - J * hi * spin[i];
	  e += ei;
	  m += (double)spin[i];
	}
	e *= 0.5;

	/* Storing observables */
	std::size_t t_measure = t / t_skip;
	std::size_t index = t_measure % buffer_size;
	es[index] = e;
	ms[index] = m;

	auto output_result = [&](std::size_t T1, std::size_t T){
	  for(std::size_t ti=0; ti < T; ti++){
	    std::size_t tp = T1 + ti * t_skip;
	    e_out << tp << std::setw(precision + 10) << es[ti] / (double)N << std::endl;
	    m_out << tp << std::setw(precision + 10) << ms[ti] / (double)N << std::endl;
	  }
	};

	/* Output */      
	if ( index == buffer_size - 1 ) {
	  output_result( t - (buffer_size - 1) * t_skip, buffer_size );
	} else if ( t == n_mcs - 1 ) {
	  std::size_t rem = t_measure % buffer_size;
	  output_result( t - rem * t_skip, rem + 1 );
	}
      }
    } /* end for t */

    /* Closing the output files */
    e_out.close();
    m_out.close();
  } /* end for task */

  /* MPI_Finalize */
  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
