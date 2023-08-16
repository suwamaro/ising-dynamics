/***************************************************************************
* Square lattice class.
*
* Copyright (C) 2021 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*****************************************************************************/

#ifndef _SQUARE_LATTICE1_3_
#define _SQUARE_LATTICE1_3_

#include <vector>
#include "lattice1.2.h"

namespace Lattice {
  class SquareLattice : public Lattice {
  public:
    SquareLattice(std::size_t L, std::size_t d, bool plaquette, bool fast_code);
    static std::unique_ptr<SquareLattice> mk_square(std::size_t L, std::size_t d, bool plaquette, bool fast_code);
    void init();
    std::vector<std::size_t> site_to_coordinate(std::size_t i) const override;
    std::vector<double> site_to_real_coordinate(std::size_t i) const override;    
    std::size_t coordinate_to_site(std::vector<std::size_t> r) const override;
    cx_double phase(std::size_t, std::vector<double> const& q) const override;
    void set_neighbors();    
    void set_neighbors_plaquettes();
    std::size_t neighbor_index(std::size_t d, std::size_t i, std::size_t j) const;
    std::size_t neighbor_plaquette_index(std::size_t i, std::size_t j, std::size_t k) const;
    std::size_t n_neighbors(std::size_t d, std::size_t i) const override;    
    std::size_t n_neighbors_plaquette(std::size_t i) const override;
    std::size_t n_neighbors_plaquette(std::size_t i, std::size_t j) const override;    
    std::size_t get_neighbor(std::size_t d, std::size_t i, std::size_t j) const override;        
    std::size_t get_neighbor_plaquette(std::size_t i, std::size_t j, std::size_t k) const override;    
    std::size_t get_neighbor_fast(std::size_t d, std::size_t i, std::size_t j) const;
    std::size_t get_neighbor_plaquette_fast(std::size_t i, std::size_t j, std::size_t k) const;
    std::vector<std::size_t>::const_iterator neighbor_ptr(std::size_t i) const override;
    void make_neighbors_fast();
    bool plaquette_interaction() const { return plaquette_interaction_; }
    bool fast_code() const { return fast_code_; }
    const std::size_t n_neighbors() const { return n_neighbors_; }
    std::size_t total_n_neighbors() const { return total_n_neighbors_; }
    std::size_t total_n_neighbors_plaquette() const { return total_n_neighbors_plaquette_; }
    std::size_t total_n_plaquettes() const { return total_n_plaquettes_; }            
    std::size_t L() const { return L_; }
    std::size_t neighbors_fast(std::size_t i) const { return neighbors_fast_[i]; }
    template<class T> static std::vector<T> make_2by1_state(std::size_t N, T up, T down, int type = 0);
    template<class T> static std::vector<T> make_2by2_state(std::size_t N, T up, T down);    
    template<class T> static std::vector<T> make_4by2_state(std::size_t N, T up, T down, int type = 0);
    template<class T> static std::vector<T> make_4by4_state(std::size_t N, T up, T down, int type = 0);
      
  private:
    bool plaquette_interaction_;
    bool fast_code_;    
    static constexpr std::size_t n_neighbors_ = 4;
    std::size_t total_n_neighbors_;
    std::size_t total_n_neighbors_plaquette_;
    std::size_t total_n_plaquettes_;    
    std::size_t L_;
    std::vector<std::size_t> neighbors_fast_;
  };

  SquareLattice::SquareLattice(std::size_t L, std::size_t d, bool plaquette, bool fast_code):Lattice(L*L,d),plaquette_interaction_(plaquette),fast_code_(fast_code),L_(L){    
    total_n_neighbors_ = int_dist() * n_neighbors();
    init();
  }
  
  std::unique_ptr<SquareLattice> SquareLattice::mk_square(std::size_t L, std::size_t d, bool plaquette, bool fast_code) {
    return std::make_unique<SquareLattice>(L,d, plaquette, fast_code);
  }
  
  void SquareLattice::init(){
    if ( fast_code() ) {
      total_n_neighbors_plaquette_ = 12;
      total_n_plaquettes_ = total_n_neighbors_plaquette_ / (4-1);
      make_neighbors_fast();
    } else {
      total_n_neighbors_plaquette_ = 0;
      total_n_plaquettes_ = 0;      
      set_neighbors();
      if ( plaquette_interaction() ) {
	set_neighbors_plaquettes();
      }
    }
  }

  std::vector<std::size_t> SquareLattice::site_to_coordinate(std::size_t i) const {
    std::size_t x = i % L();
    std::size_t y = i / L();
    std::vector<std::size_t> r;
    r.push_back(x);
    r.push_back(y);  
    return r;
  }

  std::vector<double> SquareLattice::site_to_real_coordinate(std::size_t i) const {
    std::vector<std::size_t> r = site_to_coordinate(i);
    std::vector<double> r2({(double)r[0], (double)r[1]});
    return r2;
  }
  
  cx_double SquareLattice::phase(std::size_t i, std::vector<double> const& q) const {
    std::vector<double> r = site_to_real_coordinate(i);
    double sum = 0;
    for(std::size_t d=0; d < r.size(); d++){
      sum += r[d] * q[d];
    }
    return exp(cx_double(0, sum));
  }
  
  std::size_t SquareLattice::coordinate_to_site(std::vector<std::size_t> r) const {
    return r[0] + r[1] * L();
  }

  void SquareLattice::set_neighbors(){
    neighbors_.clear();
    neighbors_.resize(int_dist());
    for(std::size_t d=0; d < int_dist(); d++){
      std::vector<std::vector<std::size_t>> neighbors_d(n_sites());
      for(std::size_t i=0; i < n_sites(); i++){
	std::vector<std::size_t> r = site_to_coordinate(i);
	std::size_t x = r[0];
	std::size_t y = r[1];    
	std::vector<std::size_t> v;
	if ( d == 0 ) {
	  /* Nearest neighbors */
	  std::size_t x_1 = periodic_shift(x, -1, L());
	  std::size_t x1 = periodic_shift(x,  1, L());
	  std::size_t y_1 = periodic_shift(y, -1, L());
	  std::size_t y1 = periodic_shift(y,  1, L());	  
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_1, y})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y_1})) );
	} else if ( d == 1 ) {
	  /* Next nearest neighbors */
	  std::size_t x_1 = periodic_shift(x, -1, L());
	  std::size_t x1 = periodic_shift(x,  1, L());
	  std::size_t y_1 = periodic_shift(y, -1, L());
	  std::size_t y1 = periodic_shift(y,  1, L());      	  
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_1, y1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_1, y_1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y_1})) );	  
	} else if ( d == 2 ) {
	  /* Third nearest neighbors */
	  std::size_t x_2 = periodic_shift(x, -2, L());
	  std::size_t x2 = periodic_shift(x,  2, L());
	  std::size_t y_2 = periodic_shift(y, -2, L());
	  std::size_t y2 = periodic_shift(y,  2, L());      
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y})) );	  
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_2, y})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y_2})) );
	} else if ( d == 3 ) {
	  std::size_t x_2 = periodic_shift(x, -2, L());
	  std::size_t x_1 = periodic_shift(x, -1, L());
	  std::size_t x1 = periodic_shift(x, 1, L());
	  std::size_t x2 = periodic_shift(x, 2, L());
	  std::size_t y_2 = periodic_shift(y, -2, L());
	  std::size_t y_1 = periodic_shift(y, -1, L());
	  std::size_t y1 = periodic_shift(y, 1, L());
	  std::size_t y2 = periodic_shift(y, 2, L());
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_1, y2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_2, y1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_2, y_1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x_1, y_2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y_2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y_1})) );	  	  	  
	} else {
	  std::cerr << "Interaction distance larger than 4 is not supported." << std::endl;
	  std::exit(EXIT_FAILURE);
	}
  
	neighbors_d[i] = v;
      }

      neighbors_[d] = neighbors_d;
    }
  }
  
  void SquareLattice::set_neighbors_plaquettes(){
    neighbors_plaquette_.clear();
    neighbors_plaquette_.resize(n_sites());
    for(std::size_t i=0; i < n_sites(); i++){      
      std::vector<std::size_t> r = site_to_coordinate(i);
      std::size_t x = r[0];
      std::size_t y = r[1];
      std::size_t x1 = periodic_shift(x, -1, L());
      std::size_t x2 = periodic_shift(x,  1, L());
      std::size_t y1 = periodic_shift(y, -1, L());
      std::size_t y2 = periodic_shift(y,  1, L());      

      std::vector<std::size_t> v1;
      v1.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y})));
      v1.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y2})));
      v1.push_back( coordinate_to_site(std::vector<std::size_t>({x, y2})));

      std::vector<std::size_t> v2;
      v2.push_back( coordinate_to_site(std::vector<std::size_t>({x, y2})));
      v2.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y2})));
      v2.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y})));

      std::vector<std::size_t> v3;
      v3.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y})));
      v3.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y1})));
      v3.push_back( coordinate_to_site(std::vector<std::size_t>({x, y1})));

      std::vector<std::size_t> v4;
      v4.push_back( coordinate_to_site(std::vector<std::size_t>({x, y1})));
      v4.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y1})));
      v4.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y})));                  

      std::vector<std::vector<std::size_t>> v;
      v.push_back(v1);
      v.push_back(v2);
      v.push_back(v3);
      v.push_back(v4);
      
      neighbors_plaquette_[i] = v;
    }
  }

  void SquareLattice::make_neighbors_fast(){
    /* Naive implementation */
    set_neighbors();

    /* Making a one-dimensional array */
    neighbors_fast_.clear();
    neighbors_fast_.resize( (total_n_neighbors() + total_n_neighbors_plaquette()) * n_sites() );
    
    for(std::size_t i=0; i < n_sites(); i++){    
      for(std::size_t d=0; d < int_dist(); d++){
	for(std::size_t j=0; j < n_neighbors(); j++){
	  std::size_t idx = neighbor_index(d, i, j);
	  neighbors_fast_[idx] = Lattice::get_neighbor(d, i, j);
	}
      }
    }

    if ( plaquette_interaction() ) {
      set_neighbors_plaquettes();
      for(std::size_t i=0; i < n_sites(); i++){
	for(std::size_t j=0; j < total_n_plaquettes(); j++){
	  for(std::size_t k=0; k < 3; k++){
	    std::size_t idx = neighbor_plaquette_index(i, j , k);
	    neighbors_fast_[idx] = Lattice::get_neighbor_plaquette(i, j , k);
	  }
	}
      }
    }

    /* Clear naive ones */
    neighbors_.clear();
    neighbors_plaquette_.clear();
  }

  std::size_t SquareLattice::neighbor_index(std::size_t d, std::size_t i, std::size_t j) const {
    return i * (total_n_neighbors() + total_n_neighbors_plaquette()) + d * n_neighbors() + j;
  }
  
  std::size_t SquareLattice::neighbor_plaquette_index(std::size_t i, std::size_t j, std::size_t k) const {
    return i * (total_n_neighbors() + total_n_neighbors_plaquette()) + total_n_neighbors() + j * 3  + k;
  }

  std::size_t SquareLattice::n_neighbors(std::size_t d, std::size_t i) const {
    if ( fast_code() ) {
      return n_neighbors();
    } else {
      return Lattice::n_neighbors(d, i);
    }
  }

  std::size_t SquareLattice::n_neighbors_plaquette(std::size_t i) const {
    if ( fast_code() ) {
      return total_n_plaquettes();
    } else {
      return Lattice::n_neighbors_plaquette(i);
    }    
  }

  std::size_t SquareLattice::n_neighbors_plaquette(std::size_t i, std::size_t j) const {
    if ( fast_code() ) {
      return 3;
    } else {
      return Lattice::n_neighbors_plaquette(i, j);
    }        
  }  
  
  std::size_t SquareLattice::get_neighbor(std::size_t d, std::size_t i, std::size_t j) const {
    if ( fast_code() ) {
      return get_neighbor_fast(d, i, j);
    } else {
      return Lattice::get_neighbor(d, i, j);
    }
  }
  
  std::size_t SquareLattice::get_neighbor_plaquette(std::size_t i, std::size_t j, std::size_t k) const {
    if ( fast_code() ) {
      return get_neighbor_plaquette_fast(i, j, k);
    } else {
      return Lattice::get_neighbor_plaquette(i, j, k);
    }
  }
  
  std::size_t SquareLattice::get_neighbor_fast(std::size_t d, std::size_t i, std::size_t j) const {
    std::size_t idx = neighbor_index(d, i, j);
    return neighbors_fast(idx);
  }
  
  std::size_t SquareLattice::get_neighbor_plaquette_fast(std::size_t i, std::size_t j, std::size_t k) const {
    std::size_t idx = neighbor_plaquette_index(i, j, k);
    return neighbors_fast(idx);
  }

  std::vector<std::size_t>::const_iterator SquareLattice::neighbor_ptr(std::size_t i) const {
    if ( fast_code() ) {
      std::size_t idx = (total_n_neighbors() + total_n_neighbors_plaquette()) * i;
      return neighbors_fast_.begin() + idx;
    } else {
      return Lattice::neighbor_ptr(i);
    }
  }

  template<class T> std::vector<T> SquareLattice::make_2by1_state(std::size_t N, T up, T down, int type){
    std::vector<T> spins(N);
    std::size_t L = (std::size_t)(round(sqrt(N)));
    assert(N == L * L);   // Isotropic system
    
    for(std::size_t s=0; s < N; s++){
      if ( type == 0 ) {
	int x = s % L;      
	spins[s] = x & 1 ? down : up;
      } else if ( type == 1 ) {
	int y = s / L;
	spins[s] = y & 1 ? down : up;
      } else {
	std::cerr << "type must be 0 or 1 in " << __func__ << std::endl;
	std::exit(EXIT_FAILURE);
      }
    }
    
    return spins;
  }
  
  template<class T> std::vector<T> SquareLattice::make_2by2_state(std::size_t N, T up, T down){
    std::vector<T> spins(N);
    std::size_t L = (std::size_t)(round(sqrt(N)));
    assert(N == L * L);   // Isotropic system
    
    for(std::size_t s=0; s < N; s++){
      int x = s % L;
      int y = s / L;
      spins[s] = (x+y) & 1 ? down : up;
    }
    
    return spins;
  }
  
  template<class T> std::vector<T> SquareLattice::make_4by2_state(std::size_t N, T up, T down, int type){
    std::vector<T> spins(N);
    std::size_t L = (std::size_t)(round(sqrt(N)));
    assert(N == L * L);   // Isotropic system

    for(std::size_t s=0; s < N; s++){
      int x = s % L;
      int y = s / L;
      int x2 = 0, y2 = 0, xbi = 0, ybi = 0;
      if ( type == 0 ) {
	x2 = x;
	y2 = y;
	xbi = 1;
	ybi = 0;
      } else if ( type == 1 ) {
	x2 = x + 1;
	y2 = y;
	xbi = 1;
	ybi = 0;
      } else if ( type == 2 ) {
	x2 = x;
	y2 = y;
	xbi = 0;
	ybi = 1;
      } else if ( type == 3 ) {
	x2 = x;
	y2 = y + 1;
	xbi = 0;
	ybi = 1;      
      } else {
	std::cerr << "type must be 0, 1, 2, or 3 in " << __func__ << std::endl;
	std::exit(EXIT_FAILURE);
      }
    
      int xb = get_bit(x, xbi);
      int yb = get_bit(y, ybi);    
      spins[s] = xb == yb ? up : down;
    }

    return spins;
  }
  
  template<class T> std::vector<T> SquareLattice::make_4by4_state(std::size_t N, T up, T down, int type){
    std::vector<T> spins(N);
    std::size_t L = (std::size_t)(round(sqrt(N)));
    assert(N == L * L);   // Isotropic system

    for(std::size_t s=0; s < N; s++){
      int x = s % L;
      int y = s / L;
      int x2 = 0, y2 = 0;
      if ( type == 0 ) {
	x2 = x;
	y2 = y;
      } else if ( type == 1 ) {
	x2 = x + 1;
	y2 = y;
      } else if ( type == 2 ) {
	x2 = x;
	y2 = y + 1;
      } else if ( type == 3 ) {
	x2 = x + 1;
	y2 = y + 1;
      } else {
	std::cerr << "type must be 0, 1, 2, or 3 in " << __func__ << std::endl;
	std::exit(EXIT_FAILURE);
      }

      int xb = get_bit(x2, 1);
      int yb = get_bit(y2, 1);    
      spins[s] = xb == yb ? up : down;
    }
    
    return spins;
  }
}

#endif  // _SQUARE_LATTICE1_3_
