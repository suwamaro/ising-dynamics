/***************************************************************************
* Functions related to domains.
*
* Copyright (C) 2021 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*****************************************************************************/

#ifndef _DOMAIN_
#define _DOMAIN_

#include <vector>
#include <stack>
#include "lattice1.2.h"

std::vector<std::size_t> picking_out_largest(std::vector<std::vector<std::size_t>> const& v, std::size_t rank){
  /* Flattening */
  std::vector<std::size_t> comb;
  for(std::size_t type = 0; type < v.size(); type++){
    comb.insert(comb.end(), v[type].begin(), v[type].end());
  }

  /* Sorting */
  std::sort(comb.begin(), comb.end(), std::greater<std::size_t>{});

  /* Picking out the largest ones. */
  std::vector<std::size_t> largest(rank, 0);
  std::size_t copy_range = std::min(comb.size(), rank);      
  std::copy(comb.begin(), comb.begin() + copy_range, largest.begin());
  
  return largest;
}

std::vector<std::size_t> find_largest_domain_sizes(Lattice::Lattice const& lattice, std::vector<int> const& spin, std::vector<int> const& ref, std::size_t range, std::size_t rank){
  int ops[2] = {1, -1};  /* Up or down */
  std::size_t n_sites = spin.size();  
  std::vector<int> spin2(n_sites);
  
  /* Overlap with the reference state. */
  for(std::size_t s = 0; s < n_sites; s++){
    spin2[s] = spin[s] * ref[s];
  }
  
  /* Checking the neighborhood. */
  std::vector<int> local_ops(n_sites, 0);
  for(std::size_t i = 0; i < n_sites; i++){
    for(int op: ops){
      if ( spin2[i] == op ) {
	bool ordered = true;
	for(std::size_t d=0; d < range; d++){
	  for(std::size_t j=0; j < lattice.n_neighbors(d,i); j++){
	    std::size_t sj = lattice.get_neighbor(d,i,j);
	    if ( spin2[sj] != op ) {
	      ordered = false;
	      break;
	    }
	  }
	}

	if ( ordered ) {
	  local_ops[i] = op;
	}
      }
    }
  }

  /* Measuring the domain sizes. */
  std::vector<std::size_t> size_vec;
  for(int op: ops){
    std::vector<bool> is_checked(n_sites, false);
    for(std::size_t i = 0; i < n_sites; i++){
      if ( is_checked[i] ) { continue; }
      is_checked[i] = true;

      /* Starting if site i is ordered. */
      if ( local_ops[i] == op ) {
	std::size_t size = 0;	
	std::stack<std::size_t> stck;
	stck.push(i);

	/* Continue unless stck is empty. */
	while(!stck.empty()){
	  std::size_t t = stck.top();
	  stck.pop();
	  ++size;	  
	  
	  /* Checking the nearest neighbors only. */
	  for(std::size_t j=0; j < lattice.n_neighbors(0,t); j++){
	    std::size_t sj = lattice.get_neighbor(0,t,j);
	    if ( !is_checked[sj] ) {
	      is_checked[sj] = true;

	      /* Adding the neighboring site to the domain if it is ordered. */
	      if ( local_ops[sj] == op ) {
		stck.push(sj);
	      }
	    }	    
	  }
	}

	/* Storing the domain dize. */
	size_vec.push_back(size);
      }
    }
  }

  /* Ordering the domain sizes. */
  std::vector<std::size_t> largest(rank, 0);
  if ( size_vec.size() > 0 ) {
    std::sort(size_vec.begin(), size_vec.end(), std::greater<std::size_t>{});
    std::size_t copy_range = std::min(size_vec.size(), rank);
    std::copy(size_vec.begin(), size_vec.begin() + copy_range, largest.begin());
  }
  
  return largest;
}

#endif  // _DOMAIN_
