module alglib.specialfunctions;
import alglib.ap;
import alglib.internal;

/**
extern(C++,alglib)
{
  double gammafunction(const double x);
  double lngamma(const double x, double &sgngam);
  double errorfunction(const double x);
  double errorfunctionc(const double x);
  double normaldistribution(const double x);
  double inverf(const double e);
  double invnormaldistribution(const double y0);
  double incompletegamma(const double a, const double x);
  double incompletegammac(const double a, const double x);
  double invincompletegammac(const double a, const double y0);
  void airy(const double x, double &ai, double &aip, double &bi, double &bip);
  double besselj0(const double x);
  double besselj1(const double x);
  double besseljn(const ae_int_t n, const double x);
  double bessely0(const double x);
  double bessely1(const double x);
  double besselyn(const ae_int_t n, const double x);
  double besseli0(const double x);
  double besseli1(const double x);
  double besselk0(const double x);
  double besselk1(const double x);
  double besselkn(const ae_int_t nn, const double x);
  double beta(const double a, const double b);
  double incompletebeta(const double a, const double b, const double x);
  double invincompletebeta(const double a, const double b, const double y);
  double binomialdistribution(const ae_int_t k, const ae_int_t n, const double p);
  double binomialcdistribution(const ae_int_t k, const ae_int_t n, const double p);
  double invbinomialdistribution(const ae_int_t k, const ae_int_t n, const double y);
  double chebyshevcalculate(const ae_int_t r, const ae_int_t n, const double x);
  double chebyshevsum(const real_1d_array &c, const ae_int_t r, const ae_int_t n, const double x);
  void chebyshevcoefficients(const ae_int_t n, real_1d_array &c);
  void fromchebyshev(const real_1d_array &a, const ae_int_t n, real_1d_array &b);
  double chisquaredistribution(const double v, const double x);
  double chisquarecdistribution(const double v, const double x);
  double invchisquaredistribution(const double v, const double y);
  double dawsonintegral(const double x);
  double ellipticintegralk(const double m);
  double ellipticintegralkhighprecision(const double m1);
  double incompleteellipticintegralk(const double phi, const double m);
  double ellipticintegrale(const double m);
  double exponentialintegralei(const double x);
  double exponentialintegralen(const double x, const ae_int_t n);
  double fdistribution(const ae_int_t a, const ae_int_t b, const double x);
  double fcdistribution(const ae_int_t a, const ae_int_t b, const double x);
  double invfdistribution(const ae_int_t a, const ae_int_t b, const double y);
  void fresnelintegral(const double x, double &c, double &s);
  double hermitecalculate(const ae_int_t n, const double x);
  double hermitesum(const real_1d_array &c, const ae_int_t n, const double x);
  void hermitecoefficients(const ae_int_t n, real_1d_array &c);
  void jacobianellipticfunctions(const double u, const double m, double &sn, double &cn, double &dn, double &ph);
  double laguerrecalculate(const ae_int_t n, const double x);
  double laguerresum(const real_1d_array &c, const ae_int_t n, const double x);
  void laguerrecoefficients(const ae_int_t n, real_1d_array &c);
  double legendrecalculate(const ae_int_t n, const double x);
  double legendresum(const real_1d_array &c, const ae_int_t n, const double x);
  void legendrecoefficients(const ae_int_t n, real_1d_array &c);
  double poissondistribution(const ae_int_t k, const double m);
  double poissoncdistribution(const ae_int_t k, const double m);
  double invpoissondistribution(const ae_int_t k, const double y);
  double psi(const double x);
  double studenttdistribution(const ae_int_t k, const double t);
  double invstudenttdistribution(const ae_int_t k, const double p);
  void sinecosineintegrals(const double x, double &si, double &ci);
  void hyperbolicsinecosineintegrals(const double x, double &shi, double &chi);
}

namespace alglib_impl
{
  double gammafunction(double x, ae_state *_state);
  double lngamma(double x, double* sgngam, ae_state *_state);
  double errorfunction(double x, ae_state *_state);
  double errorfunctionc(double x, ae_state *_state);
  double normaldistribution(double x, ae_state *_state);
  double inverf(double e, ae_state *_state);
  double invnormaldistribution(double y0, ae_state *_state);
  double incompletegamma(double a, double x, ae_state *_state);
  double incompletegammac(double a, double x, ae_state *_state);
  double invincompletegammac(double a, double y0, ae_state *_state);
  void airy(double x, double* ai, double* aip, double* bi, double* bip, ae_state *_state);
  double besselj0(double x, ae_state *_state);
  double besselj1(double x, ae_state *_state);
  double besseljn(ae_int_t n, double x, ae_state *_state);
  double bessely0(double x, ae_state *_state);
  double bessely1(double x, ae_state *_state);
  double besselyn(ae_int_t n, double x, ae_state *_state);
  double besseli0(double x, ae_state *_state);
  double besseli1(double x, ae_state *_state);
  double besselk0(double x, ae_state *_state);
  double besselk1(double x, ae_state *_state);
  double besselkn(ae_int_t nn, double x, ae_state *_state);
  double beta(double a, double b, ae_state *_state);
  double incompletebeta(double a, double b, double x, ae_state *_state);
  double invincompletebeta(double a, double b, double y, ae_state *_state);
  double binomialdistribution(ae_int_t k, ae_int_t n, double p, ae_state *_state);
  double binomialcdistribution(ae_int_t k, ae_int_t n, double p, ae_state *_state);
  double invbinomialdistribution(ae_int_t k, ae_int_t n, double y, ae_state *_state);
  double chebyshevcalculate(ae_int_t r, ae_int_t n, double x, ae_state *_state);
  double chebyshevsum(ae_vector* c, ae_int_t r, ae_int_t n, double x, ae_state *_state);
  void chebyshevcoefficients(ae_int_t n, ae_vector* c, ae_state *_state);
  void fromchebyshev(ae_vector* a, ae_int_t n,  ae_vector* b, ae_state *_state);
  double chisquaredistribution(double v, double x, ae_state *_state);
  double chisquarecdistribution(double v, double x, ae_state *_state);
  double invchisquaredistribution(double v, double y, ae_state *_state);
  double dawsonintegral(double x, ae_state *_state);
  double ellipticintegralk(double m, ae_state *_state);
  double ellipticintegralkhighprecision(double m1, ae_state *_state);
  double incompleteellipticintegralk(double phi, double m, ae_state *_state);
  double ellipticintegrale(double m, ae_state *_state);
  double incompleteellipticintegrale(double phi, double m, ae_state *_state);
  double exponentialintegralei(double x, ae_state *_state);
  double exponentialintegralen(double x, ae_int_t n, ae_state *_state);
  double fdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
  double fcdistribution(ae_int_t a, ae_int_t b, double x, ae_state *_state);
  double invfdistribution(ae_int_t a, ae_int_t b, double y, ae_state *_state);
  void fresnelintegral(double x, double* c, double* s, ae_state *_state);
  double hermitecalculate(ae_int_t n, double x, ae_state *_state);
  double hermitesum( ae_vector* c, ae_int_t n, double x, ae_state *_state);
  void hermitecoefficients(ae_int_t n,  ae_vector* c, ae_state *_state);
  void jacobianellipticfunctions(double u, double m, double* sn, double* cn, double* dn, double* ph, ae_state *_state);
  double laguerrecalculate(ae_int_t n, double x, ae_state *_state);
  double laguerresum( ae_vector* c, ae_int_t n, double x, ae_state *_state);
  void laguerrecoefficients(ae_int_t n,  ae_vector* c, ae_state *_state);
  double legendrecalculate(ae_int_t n, double x, ae_state *_state);
  double legendresum( ae_vector* c, ae_int_t n, double x, ae_state *_state);
  void legendrecoefficients(ae_int_t n,  ae_vector* c, ae_state *_state);
  double poissondistribution(ae_int_t k, double m, ae_state *_state);
  double poissoncdistribution(ae_int_t k, double m, ae_state *_state);
  double invpoissondistribution(ae_int_t k, double y, ae_state *_state);
  double psi(double x, ae_state *_state);
  double studenttdistribution(ae_int_t k, double t, ae_state *_state);
  double invstudenttdistribution(ae_int_t k, double p, ae_state *_state);
  void sinecosineintegrals(double x, double* si, double* ci, ae_state *_state);
  void hyperbolicsinecosineintegrals(double x, double* shi, double* chi, ae_state *_state);
}

*/