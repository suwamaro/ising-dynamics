/***************************************************************************
* Square lattice class.
*
* Copyright (C) 2021 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*****************************************************************************/

#ifndef _SQUARE_LATTICE_
#define _SQUARE_LATTICE_

#include <vector>
#include "lattice.h"

namespace Lattice {

  class SquareLattice : public Lattice {
  public:
    SquareLattice(std::size_t L, std::size_t d);
    static std::unique_ptr<SquareLattice> mk_square(std::size_t L, std::size_t d);
    void init() override;
    std::vector<std::size_t> site_to_coordinate(std::size_t i) const override;
    std::vector<double> site_to_real_coordinate(std::size_t i) const override;    
    std::size_t coordinate_to_site(std::vector<std::size_t> r) const override;
    cx_double phase(std::size_t, std::vector<double> const& q) const override;
    void set_neighbors();    
    void set_neighbors_plaquettes();
  
    std::size_t L() const { return L_; }
    
  private:
    std::size_t L_;  
  };

  SquareLattice::SquareLattice(std::size_t L, std::size_t d):Lattice(L*L,d),L_(L){
    init();
  }

  void SquareLattice::init(){
    set_neighbors();
    set_neighbors_plaquettes();  
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
	  std::size_t x1 = periodic_shift(x, -1, L());
	  std::size_t x2 = periodic_shift(x,  1, L());
	  std::size_t y1 = periodic_shift(y, -1, L());
	  std::size_t y2 = periodic_shift(y,  1, L());	  
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y1})) );
	} else if ( d == 1 ) {
	  /* Next nearest neighbors */
	  std::size_t x1 = periodic_shift(x, -1, L());
	  std::size_t x2 = periodic_shift(x,  1, L());
	  std::size_t y1 = periodic_shift(y, -1, L());
	  std::size_t y2 = periodic_shift(y,  1, L());      	  
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y2})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y1})) );
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y1})) );	  
	} else if ( d == 2 ) {
	  /* Third nearest neighbors */
	  std::size_t x1 = periodic_shift(x, -2, L());
	  std::size_t x2 = periodic_shift(x,  2, L());
	  std::size_t y1 = periodic_shift(y, -2, L());
	  std::size_t y2 = periodic_shift(y,  2, L());      
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x2, y})));
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y2})));
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x1, y})));
	  v.push_back( coordinate_to_site(std::vector<std::size_t>({x, y1})));	  
	} else {
	  std::cerr << "Interaction distance larger than 2 is not supported." << std::endl;
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

  std::unique_ptr<SquareLattice> SquareLattice::mk_square(std::size_t L, std::size_t d) {
    return std::make_unique<SquareLattice>(L,d);
  }
  
}

#endif  // _SQUARE_LATTICE_
