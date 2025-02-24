// 3-Clause BSD License
//
// Copyright (c) 2023-2025 Renaud AUBIN
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions
//    and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to
//    endorse or promote products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
// WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// This is a C++ implementation of the Frisch-Newton algorithm and the preprocessing phase of the rq
// algorithm as described in Portnoy and Koenker, Statistical Science, (1997) 279-300.
// Original version consists in the Roger Koenker's quantreg R package rq.fit.pfn and rq.fit.fnb
// functions (quantreg.R) and the rqfnb fortran function (rqfnb.f). Quantreg R package repository:
// https://github.com/cran/quantreg.
// Translated and adapted to C++ / Rcpp / RcppArmadillo by Dajun Xu, under the instructions of
// Stephen L. Portnoy, 2015. Repository: https://github.com/xdj701/STAT391-Quantreg.
// Adapted to C++ & Armadillo by Renaud Aubin, 2023.
//
// This file is released under the 3-Clause BSD License, although the original versions were
// not. Permissions was given by Professor Roger Koenker and Dajun Xu to Renaud Aubin on May/June
// 2023 via electronic correspondences.
// Original correspondences can be provided upon request.
//

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "quantreg-cpp"

namespace quantreg {
const double big = 1e20;
const int max_it = 50;

arma::mat cross_chol_inv(const arma::mat& x) {
  arma::mat tmp = arma::trans(x) * x;
  arma::mat rt = arma::inv(arma::chol(tmp));
  return rt;
}

arma::mat times_sqrt(const arma::mat& x, const arma::mat& xxinv, arma::vec p) {
  arma::mat p_arma(p.memptr(), p.n_elem, 1, false);
  arma::mat tmp = arma::pow((x)*xxinv, 2);
  arma::mat rt = arma::sqrt(tmp * p_arma);
  return rt;
}

arma::mat sub_matrix(const arma::mat& mat, const arma::ivec& samp) {
  arma::uword row = samp.n_elem;
  arma::uword col = mat.n_cols;
  arma::mat rt(row, col);
  for (arma::uword i = 0; i < row; i++) {
    rt.row(i) = mat.row(samp(i));
  }
  return rt;
}

arma::vec sub_vector(const arma::vec& vec, const arma::ivec& samp) {
  arma::uword len = samp.n_elem;
  arma::vec rt(len);
  for (arma::uword i = 0; i < len; i++) {
    rt(i) = vec(samp(i));
  }
  return rt;
}

arma::vec mat_times_vec(const arma::mat& x, arma::vec p) {
  arma::mat p_arma(p.memptr(), p.n_elem, 1, false);
  return x * p_arma;
}

arma::ivec trueVector(const arma::uchar_vec& x) {
  std::vector<long long> count;
  for (arma::uword i = 0; i < x.n_elem; i++) {
    if (x[i] == true) count.push_back(i);  // this work because 1 == true is true in C++
  }
  return arma::ivec(count);
}

arma::vec glob_transform(const arma::mat& mat, const arma::uchar_vec& vec) {
  arma::mat sub = sub_matrix(mat, trueVector(vec));
  arma::vec rt(sub.n_cols);

  for (int i = 0; i < rt.n_elem; i++) {
    for (int j = 0; j < sub.n_rows; j++) rt(i) += sub(j, i);
  }
  return rt;
}

std::tuple<arma::vec, arma::vec> rqfnb(const arma::mat& X, const arma::vec& y, double tau,
                                       double beta, double eps);

std::tuple<arma::vec, arma::vec> qr(const arma::mat& X, const arma::vec& y, double tau,
                                    double Mm_factor, int max_bad_fixup, double eps) {
  // Create a copy of x and add a column for intercept if the x matrix doesn't already contain 1s
  // as first column
  arma::mat x = X;
  if (!arma::all(x.col(0) == 1)) {
    x.insert_cols(0, arma::vec(x.n_rows, arma::fill::ones));
  }

  int n = y.n_elem;
  if (x.n_rows != n) {
    throw std::runtime_error("X and y don't match n");
  }
  if (tau < 0 || tau > 1) {
    throw std::runtime_error("tau outside (0,1)");
  }
  int p = x.n_cols;
  int m = round(pow((p + 1.0) * n, 2.0 / 3));
  bool not_optimal = true;
  arma::vec b(p);
  while (not_optimal) {
    if (m >= n) {
      return rqfnb(x, y, tau, 0.99995, eps);
    } else {
      arma::arma_rng::set_seed_random();
      arma::ivec tmp = arma::linspace<arma::ivec>(0, n - 1, n);
      arma::ivec s = arma::shuffle(tmp);
      s.resize(m);

      arma::mat xx = sub_matrix(x, s);
      arma::vec yy = sub_vector(y, s);
      auto z = rqfnb(xx, yy, tau, 0.99995, eps);

      arma::mat xxinv = cross_chol_inv(xx);
      arma::mat band = times_sqrt(x, xxinv, arma::vec(p, arma::fill::ones));
      arma::vec r = y - mat_times_vec(x, std::get<0>(z));
      double M = Mm_factor * m;
      double lo_q = std::max(1.0 / n, tau - M / (2.0 * n));
      double hi_q = std::min(tau + M / (2.0 * n), (n - 1.0) / n);
      arma::vec kappa =
          arma::quantile(r / arma::clamp(band, eps, band.max()), arma::vec{lo_q, hi_q});
      arma::uchar_vec su(r.n_elem, arma::fill::zeros);
      for (arma::uword i = 0; i < su.n_elem; ++i) {
        su(i) = (r(i) > (band(i) * kappa(1)));
      }

      arma::uchar_vec sl(r.n_elem, arma::fill::zeros);
      for (arma::uword i = 0; i < sl.n_elem; ++i) {
        sl(i) = (r(i) < (band(i) * kappa(0)));
      }

      int bad_fixup = 0;
      while (not_optimal && (bad_fixup < max_bad_fixup)) {
        arma::uchar_vec nsu_nsl(su.n_elem);
        for (arma::uword i = 0; i < nsu_nsl.n_elem; ++i) {
          nsu_nsl(i) = !su(i) && !sl(i);
        }
        arma::ivec tv_nsu_nsl = trueVector(nsu_nsl);
        xx = sub_matrix(x, tv_nsu_nsl);
        yy = sub_vector(y, tv_nsu_nsl);

        if (arma::any(su)) {
          arma::vec ghib_x = glob_transform(x, su);
          double ghib_y = arma::sum(sub_vector(y, trueVector(su)));
          xx.insert_rows(xx.n_rows, ghib_x.t());
          yy.resize(yy.n_elem + 1);
          yy(yy.n_elem - 1) = ghib_y;
        }
        if (arma::any(sl)) {
          arma::vec glob_x = glob_transform(x, sl);
          double glob_y = arma::sum(sub_vector(y, trueVector(sl)));
          xx.insert_rows(xx.n_rows, glob_x.t());
          yy.resize(yy.n_elem + 1);
          yy(yy.n_elem - 1) = glob_y;
        }

        z = rqfnb(xx, yy, tau, 0.99995, eps);
        b = std::get<0>(z);
        r = y - mat_times_vec(x, std::get<0>(z));
        arma::uchar_vec su_bad(su.n_elem, arma::fill::zeros);
        for (arma::uword i = 0; i < su_bad.n_elem; ++i) {
          su_bad(i) = ((r(i) < 0) && su(i));
        }
        arma::uchar_vec sl_bad(sl.n_elem, arma::fill::zeros);
        for (arma::uword i = 0; i < sl_bad.n_elem; ++i) {
          sl_bad(i) = ((r(i) > 0) && sl(i));
        }

        if (arma::any(arma::join_cols(su_bad, sl_bad))) {
          arma::uchar_vec or_bad(su_bad.n_elem, arma::fill::zeros);
          for (arma::uword i = 0; i < or_bad.n_elem; ++i) {
            or_bad(i) = (su_bad(i) || sl_bad(i));
          }
          if (arma::sum(or_bad) > 0.1 * M) {
            // warning("Too many fixups: doubling m");
            m = 2 * m;
            break;
          }
          for (arma::uword i = 0; i < su.n_elem; ++i) {
            su(i) = (su(i) && !su_bad(i));
          }
          for (arma::uword i = 0; i < sl.n_elem; ++i) {
            sl(i) = (sl(i) && !sl_bad(i));
          }
          bad_fixup++;
        } else {
          not_optimal = false;
        }
      }
    }
  }
  return {b, y - mat_times_vec(x, b)};
}

void stepy(arma::vec& b, const arma::vec& d, const arma::mat& a, arma::mat& ada) {
  arma::mat tmp = arma::trans(a);
  for (int i = 0; i < (int)tmp.n_rows; i++) {
    tmp.row(i) *= d(i);
  }

  ada = a * tmp;  // dsyr('U',p,d(i),a(1,i),1,ada,p)

  b = arma::solve(ada, b);  // dposv('U',p,1,ada,p,b,p,info)

  ada = arma::chol(ada);
  /* what dposv does; possibly faster
  ada = chol(ada);
  b = arma::solve(trans(trimatu(ada))*trimatu(ada),b);//dposv('U',p,1,ada,p,b,p,info)
  */
}

arma::vec lpfnb(int p, int n, const arma::mat& a, const arma::vec& b, const arma::vec& c,
                arma::vec& x, double beta, double eps) {
  arma::vec nit(3, arma::fill::zeros);
  nit(2) = n;
  arma::vec d(n, arma::fill::ones);  // rep(1,n)
  arma::vec y = a * c;               // dgemv('N',p,n,one,a,p,c,1,zero,y,1)
                                     // actually it's gonna be coefficients
  arma::mat ada(p, p, arma::fill::zeros);

  stepy(y, d, a, ada);

  arma::vec s = -arma::trans(a) * y + c;  // dcopy(n,c,1,s,1) + dgemv('T',p,n,mone,a,p,y,1,one,s,1)

  arma::vec u(n, arma::fill::ones);
  arma::vec z(n, arma::fill::zeros);
  arma::vec w(n, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    if (std::abs(s(i)) < eps) {
      z(i) = std::max(s(i), 0.0) + eps;
      w(i) = std::max(-s(i), 0.0) + eps;
    } else {
      z(i) = std::max(s(i), 0.0);
      w(i) = std::max(-s(i), 0.0);
    }
    s(i) = 1 - x(i);
  }

  double gap = arma::dot(z, x) + arma::dot(w, s);

  arma::vec dx(n, arma::fill::zeros);
  arma::vec ds(n, arma::fill::zeros);
  arma::vec dz(n, arma::fill::zeros);
  arma::vec dy(n, arma::fill::zeros);
  arma::vec rhs(n, arma::fill::zeros);
  arma::vec dw(n, arma::fill::zeros);
  arma::vec dr(n, arma::fill::zeros);

  while ((gap > eps) && (nit(0) < max_it)) {
    nit(0)++;
    for (int i = 0; i < n; i++) {
      d(i) = 1.0 / (z(i) / x(i) + w(i) / s(i));
      ds(i) = z(i) - w(i);
      dz(i) = d(i) * ds(i);
    }

    dy = -a * x + b;   // dcopy(p,b,1,dy,1) + dgemv('N',p,n,mone,a,p,x,1,one,dy,1)
    dy = a * dz + dy;  // dgemv('N',p,n,one,a,p,dz,1,one,dy,1)
    rhs = dy;          // dcopy(p,dy,1,rhs,1)

    stepy(dy, d, a, ada);

    ds = arma::trans(a) * dy - ds;  // call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)

    double deltap = big;
    double deltad = big;
    for (int i = 0; i < n; i++) {
      dx(i) = d(i) * ds(i);
      ds(i) = -dx(i);
      dz(i) = -z(i) * (dx(i) / x(i) + 1.0);
      dw(i) = -w(i) * (ds(i) / s(i) + 1.0);
      if (dx(i) < 0) deltap = std::min(deltap, -x(i) / dx(i));
      if (ds(i) < 0) deltap = std::min(deltap, -s(i) / ds(i));
      if (dz(i) < 0) deltad = std::min(deltad, -z(i) / dz(i));
      if (dw(i) < 0) deltad = std::min(deltad, -w(i) / dw(i));
    }
    deltap = std::min(beta * deltap, 1.0);
    deltad = std::min(beta * deltad, 1.0);

    if (std::min(deltap, deltad) < 1) {
      nit(1)++;
      double mu = arma::dot(z, x) + arma::dot(w, s);
      double g = mu + deltap * arma::dot(dx, z) + deltad * arma::dot(dz, x) +
                 deltap * deltad * arma::dot(dx, dz) + deltap * arma::dot(ds, w) +
                 deltad * arma::dot(dw, s) + deltap * deltad * arma::dot(ds, dw);  // gap is mu
      mu = mu * std::pow(g / mu, 3) / (2.0 * n);
      for (int i = 0; i < n; i++) {
        dr(i) = d(i) * (mu * (1 / s(i) - 1 / x(i)) + dx(i) * dz(i) / x(i) - ds(i) * dw(i) / s(i));
      }
      std::swap(rhs, dy);
      dy = a * dr + dy;                              // dgemv('N',p,n,one,a,p,dr,1,one,dy,1)
      dy = arma::solve(arma::trans(ada) * ada, dy);  // dpotrs('U',p,1,ada,p,dy,p,info)
      u = arma::trans(a) * dy;                       // dgemv('T',p,n,one,a,p,dy,1,zero,u,1)
      deltap = big;
      deltad = big;
      for (int i = 0; i < n; i++) {
        double dxdz = dx(i) * dz(i);
        double dsdw = ds(i) * dw(i);
        dx(i) = d(i) * (u(i) - z(i) + w(i)) - dr(i);
        ds(i) = -dx(i);
        dz(i) = -z(i) + (mu - z(i) * dx(i) - dxdz) / x(i);
        dw(i) = -w(i) + (mu - w(i) * ds(i) - dsdw) / s(i);
        if (dx(i) < 0) deltap = std::min(deltap, -x(i) / dx(i));
        if (ds(i) < 0) deltap = std::min(deltap, -s(i) / ds(i));
        if (dz(i) < 0) deltad = std::min(deltad, -z(i) / dz(i));
        if (dw(i) < 0) deltad = std::min(deltad, -w(i) / dw(i));
      }

      deltap = std::min(beta * deltap, 1.0);
      deltad = std::min(beta * deltad, 1.0);
    }
    x = deltap * dx + x;  // daxpy(n,deltap,dx,1,x,1)
    s = deltap * ds + s;  // daxpy(n,deltap,ds,1,s,1)
    y = deltad * dy + y;  // daxpy(p,deltad,dy,1,y,1)
    z = deltad * dz + z;  // daxpy(n,deltad,dz,1,z,1)
    w = deltad * dw + w;  // daxpy(n,deltad,dw,1,w,1)
    gap = arma::dot(z, x) + arma::dot(w, s);
  }
  /*
      z = -w+z; //daxpy(n,mone,w,1,z,1)
      std::swap(x,z);
  */
  return y;
}

std::tuple<arma::vec, arma::vec> rqfnb(const arma::mat& X, const arma::vec& y, double tau,
                                       double beta, double eps) {
  int n = y.n_elem;
  int p = X.n_cols;

  if (X.n_rows != n) {
    throw std::runtime_error("X and y don't match n");
  }
  if (tau < eps || tau > (1 - eps)) {
    throw std::runtime_error("No parametric Frisch-Newton method. Set tau in (0,1)");
  }
  arma::vec rhs(p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      rhs(j) += X(i, j);
    }
  }
  rhs = rhs * (1 - tau);
  arma::vec wn(n);
  wn.fill(1 - tau);
  arma::mat a = arma::trans(X);
  arma::vec c = -y;
  arma::vec coef = lpfnb(p, n, a, rhs, c, wn, beta, eps);

  arma::vec coefficients = -coef;
  arma::vec residuals = y + X * coef;

  return {coefficients, residuals};
}
}  // namespace quantreg
