// RcppArmadillo.h does not play nicely with R_ext/Applic.h -- hence this file

#include <R_ext/Applic.h>
#include <Rcpp.h>

namespace fhcrc_example {

  struct sbinorm_total { double k, mu1, sd1, mu2, sd2, norm, mean; };
  void binorm_total_norm(double *x, int n, void *ex) {
    sbinorm_total *p = (sbinorm_total *) ex;
    for (int i=0; i<n; ++i) 
      x[i] = R::dnorm(x[i], p->mu1, p->sd1, 0)*R::dnorm(p->k-x[i], p->mu2, p->sd2, 0);
  }
  void binorm_total_mean(double *x, int n, void *ex) {
    sbinorm_total *p = (sbinorm_total *) ex;
    for (int i=0; i<n; ++i) 
      x[i] = x[i]*R::dnorm(x[i], p->mu1, p->sd1, 0)*R::dnorm(p->k-x[i], p->mu2, p->sd2, 0)/p->norm;
  }
  void binorm_total_sd(double *x, int n, void *ex) {
    sbinorm_total *p = (sbinorm_total *) ex;
    for (int i=0; i<n; ++i) 
      x[i] = (x[i]-p->mean)*(x[i]-p->mean)*R::dnorm(x[i], p->mu1, p->sd1, 0)*R::dnorm(p->k-x[i], p->mu2, p->sd2, 0)/p->norm;
  }
  void binorm_total(const double k, const double mu1, const double sd1,
		    const double mu2, const double sd2, double *mean, double *sd) {
    int inf=2, neval=0, ier=0, limit=100, last=0, lenw;
    double epsabs = 1.0e-5, epsrel = 1.0e-5, result = 0.0, 
      abserr = 0.0, bound = 0.0;
    lenw = 4 * limit;
    int iwork[limit];
    double work[lenw];
    sbinorm_total p{k, mu1, sd1, mu2, sd2, 0.0, 0.0};
    Rdqagi(binorm_total_norm, (void *) &p, &bound, &inf,
	   &epsabs, &epsrel,
	   &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last,
	   iwork, work);
    p.norm = result;
    Rdqagi(binorm_total_mean, (void *) &p, &bound, &inf,
	   &epsabs, &epsrel,
	   &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last,
	   iwork, work);
    p.mean = result;
    Rdqagi(binorm_total_sd, (void *) &p, &bound, &inf,
	   &epsabs, &epsrel,
	   &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last,
	   iwork, work);
    *mean = p.mean;
    *sd = result;
  }

}
