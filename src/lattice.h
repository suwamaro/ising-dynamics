/***************************************************************************
* Square lattice class.
*
* Copyright (C) 2021 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*****************************************************************************/

#ifndef _LATTICE_
#define _LATTICE_

#include <vector>
#include <cassert>

inline std::size_t periodic_shift(std::size_t i, int d, std::size_t L) {
  assert(i + L + d >= 0);
  return (i + L + d) % L;
}

namespace Lattice {
  
  class Lattice {
  public:
    Lattice(std::size_t N, std::size_t d);
    virtual ~Lattice(){}
    virtual void init() = 0;  
    virtual std::vector<std::size_t> site_to_coordinate(std::size_t i) const = 0;
    virtual std::vector<double> site_to_real_coordinate(std::size_t i) const = 0;    
    virtual std::size_t coordinate_to_site(std::vector<std::size_t> r) const = 0;
    virtual cx_double phase(std::size_t i, std::vector<double> const& q) const = 0;
    std::size_t n_neighbors(std::size_t d, std::size_t i) const;    
    std::size_t n_neighbors_plaquette(std::size_t i) const;
    std::size_t n_neighbors_plaquette(std::size_t i, std::size_t j) const;
    std::size_t get_neighbor(std::size_t d, std::size_t i, std::size_t j) const;        
    std::size_t get_neighbor_plaquette(std::size_t i, std::size_t j, std::size_t k) const;        
    
    std::size_t n_sites() const { return n_sites_; }
    std::size_t int_dist() const { return int_dist_; }      

  protected:
    std::vector<std::vector<std::vector<std::size_t>>> neighbors_;  // Neighbors
    std::vector<std::vector<std::vector<std::size_t>>> neighbors_plaquette_;  // Plaquettes
  
  private:
    std::size_t n_sites_;
    std::size_t int_dist_;
  };

  Lattice::Lattice(std::size_t N, std::size_t d):n_sites_(N),int_dist_(d){}

  std::size_t Lattice::n_neighbors(std::size_t d, std::size_t i) const {
    return neighbors_[d][i].size();
  }
  
  std::size_t Lattice::n_neighbors_plaquette(std::size_t i) const {
    return neighbors_plaquette_[i].size();
  }

  std::size_t Lattice::n_neighbors_plaquette(std::size_t i, std::size_t j) const {
    return neighbors_plaquette_[i][j].size();
  }  

  std::size_t Lattice::get_neighbor(std::size_t d, std::size_t i, std::size_t j) const {
    return neighbors_[d][i][j];
  }

  std::size_t Lattice::get_neighbor_plaquette(std::size_t i, std::size_t j, std::size_t k) const {
    return neighbors_plaquette_[i][j][k];
  }    
  
}

#endif  // _LATTICE_
