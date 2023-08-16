/***************************************************************************
* Monte Carlo simulations for Ising models.
*
* Copyright (C) 2021 Hidemaro Suwa <suwamaro@phys.s.u-tokyo.ac.jp>
*****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <complex>

/* Uncomment if you use Boost Filesystem. */
// #define _USE_BOOST_FILESYSTEM_

#ifdef _USE_BOOST_FILESYSTEM_
/* Using boost::filesystem */
#include <boost/filesystem.hpp>
typedef boost::filesystem::ofstream ofstream;
typedef boost::filesystem::path path;
#else
/* Using std::filesystem */
typedef std::ofstream ofstream;
typedef std::filesystem::path path;
#endif

typedef std::complex<double> cx_double;

constexpr int dump_file_name_digit = 4;
