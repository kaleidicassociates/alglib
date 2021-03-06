module alglib.statistics;
import alglib.ap;
import alglib.internal;
import alglib.linalg;
import alglib.specialfunctions;

/**
extern(C++,alglib)
{
  void samplemoments(const real_1d_array &x, const ae_int_t n, double &mean, double &variance, double &skewness, double &kurtosis);
  void samplemoments(const real_1d_array &x, double &mean, double &variance, double &skewness, double &kurtosis);
  double samplemean(const real_1d_array &x, const ae_int_t n);
  double samplemean(const real_1d_array &x);
  double samplevariance(const real_1d_array &x, const ae_int_t n);
  double samplevariance(const real_1d_array &x);
  double sampleskewness(const real_1d_array &x, const ae_int_t n);
  double sampleskewness(const real_1d_array &x);
  double samplekurtosis(const real_1d_array &x, const ae_int_t n);
  double samplekurtosis(const real_1d_array &x);
  void sampleadev(const real_1d_array &x, const ae_int_t n, double &adev);
  void sampleadev(const real_1d_array &x, double &adev);
  void samplemedian(const real_1d_array &x, const ae_int_t n, double &median);
  void samplemedian(const real_1d_array &x, double &median);
  void samplepercentile(const real_1d_array &x, const ae_int_t n, const double p, double &v);
  void samplepercentile(const real_1d_array &x, const double p, double &v);
  double cov2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
  double cov2(const real_1d_array &x, const real_1d_array &y);
  double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
  double pearsoncorr2(const real_1d_array &x, const real_1d_array &y);
  double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
  double spearmancorr2(const real_1d_array &x, const real_1d_array &y);
  void covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
  void smp_covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
  void covm(const real_2d_array &x, real_2d_array &c);
  void smp_covm(const real_2d_array &x, real_2d_array &c);
  void pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
  void smp_pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
  void pearsoncorrm(const real_2d_array &x, real_2d_array &c);
  void smp_pearsoncorrm(const real_2d_array &x, real_2d_array &c);
  void spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
  void smp_spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
  void spearmancorrm(const real_2d_array &x, real_2d_array &c);
  void smp_spearmancorrm(const real_2d_array &x, real_2d_array &c);
  void covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
  void smp_covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
  void covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
  void smp_covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
  void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
  void smp_pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
  void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
  void smp_pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
  void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
  void smp_spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
  void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
  void smp_spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
  void rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
  void smp_rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
  void rankdata(real_2d_array &xy);
  void smp_rankdata(real_2d_array &xy);
  void rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
  void smp_rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
  void rankdatacentered(real_2d_array &xy);
  void smp_rankdatacentered(real_2d_array &xy);
  double pearsoncorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
  double spearmanrankcorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
  void pearsoncorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);
  void spearmanrankcorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);
  void jarqueberatest(const real_1d_array &x, const ae_int_t n, double &p);
  void mannwhitneyutest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
  void onesamplesigntest(const real_1d_array &x, const ae_int_t n, const double median, double &bothtails, double &lefttail, double &righttail);
  void studentttest1(const real_1d_array &x, const ae_int_t n, const double mean, double &bothtails, double &lefttail, double &righttail);
  void studentttest2(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
  void unequalvariancettest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
  void ftest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);
  void onesamplevariancetest(const real_1d_array &x, const ae_int_t n, const double variance, double &bothtails, double &lefttail, double &righttail);
}
*/

extern(C)
{
  void samplemoments(/* Real    */ ae_vector* x, ae_int_t n, double* mean, double* variance, double* skewness, double* kurtosis, ae_state *_state);
  double samplemean(/* Real    */ ae_vector* x, ae_int_t n, ae_state *_state);
  double samplevariance(/* Real    */ ae_vector* x, ae_int_t n, ae_state *_state);
  double sampleskewness(/* Real    */ ae_vector* x, ae_int_t n, ae_state *_state);
  double samplekurtosis(/* Real    */ ae_vector* x, ae_int_t n, ae_state *_state);
  void sampleadev(/* Real    */ ae_vector* x, ae_int_t n, double* adev, ae_state *_state);
  void samplemedian(/* Real    */ ae_vector* x, ae_int_t n, double* median, ae_state *_state);
  void samplepercentile(/* Real    */ ae_vector* x, ae_int_t n, double p, double* v, ae_state *_state);
  double cov2(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_state *_state);
  double pearsoncorr2(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_state *_state);
  double spearmancorr2(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_state *_state);
  void covm(/* Real    */ ae_matrix* x, ae_int_t n, ae_int_t m, /* Real    */ ae_matrix* c, ae_state *_state);
  void _pexec_covm(/* Real    */ ae_matrix* x, ae_int_t n, ae_int_t m, /* Real    */ ae_matrix* c, ae_state *_state);
  void pearsoncorrm(/* Real    */ ae_matrix* x, ae_int_t n, ae_int_t m, /* Real    */ ae_matrix* c, ae_state *_state);
  void _pexec_pearsoncorrm(/* Real    */ ae_matrix* x, ae_int_t n, ae_int_t m, /* Real    */ ae_matrix* c, ae_state *_state);
  void spearmancorrm(/* Real    */ ae_matrix* x, ae_int_t n, ae_int_t m, /* Real    */ ae_matrix* c, ae_state *_state);
  void _pexec_spearmancorrm(/* Real    */ ae_matrix* x, ae_int_t n, ae_int_t m, /* Real    */ ae_matrix* c, ae_state *_state);
  void covm2(/* Real    */ ae_matrix* x, /* Real    */ ae_matrix* y, ae_int_t n, ae_int_t m1, ae_int_t m2, /* Real    */ ae_matrix* c, ae_state *_state);
  void _pexec_covm2(/* Real    */ ae_matrix* x, /* Real    */ ae_matrix* y, ae_int_t n, ae_int_t m1, ae_int_t m2, /* Real    */ ae_matrix* c, ae_state *_state);
  void pearsoncorrm2(/* Real    */ ae_matrix* x, /* Real    */ ae_matrix* y, ae_int_t n, ae_int_t m1, ae_int_t m2, /* Real    */ ae_matrix* c, ae_state *_state);
  void _pexec_pearsoncorrm2(/* Real    */ ae_matrix* x, /* Real    */ ae_matrix* y, ae_int_t n, ae_int_t m1, ae_int_t m2, /* Real    */ ae_matrix* c, ae_state *_state);
  void spearmancorrm2(/* Real    */ ae_matrix* x, /* Real    */ ae_matrix* y, ae_int_t n, ae_int_t m1, ae_int_t m2, /* Real    */ ae_matrix* c, ae_state *_state);
  void _pexec_spearmancorrm2(/* Real    */ ae_matrix* x, /* Real    */ ae_matrix* y, ae_int_t n, ae_int_t m1, ae_int_t m2, /* Real    */ ae_matrix* c, ae_state *_state);
  void rankdata(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
  void _pexec_rankdata(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
  void rankdatacentered(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
  void _pexec_rankdatacentered(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_state *_state);
  double pearsoncorrelation(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_state *_state);
  double spearmanrankcorrelation(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_state *_state);
  void pearsoncorrelationsignificance(double r, ae_int_t n, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void spearmanrankcorrelationsignificance(double r, ae_int_t n, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void jarqueberatest(/* Real    */ ae_vector* x, ae_int_t n, double* p, ae_state *_state);
  void mannwhitneyutest(/* Real    */ ae_vector* x, ae_int_t n, /* Real    */ ae_vector* y, ae_int_t m, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void onesamplesigntest(/* Real    */ ae_vector* x, ae_int_t n, double median, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void studentttest1(/* Real    */ ae_vector* x, ae_int_t n, double mean, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void studentttest2(/* Real    */ ae_vector* x, ae_int_t n, /* Real    */ ae_vector* y, ae_int_t m, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void unequalvariancettest(/* Real    */ ae_vector* x, ae_int_t n, /* Real    */ ae_vector* y, ae_int_t m, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void ftest(/* Real    */ ae_vector* x, ae_int_t n, /* Real    */ ae_vector* y, ae_int_t m, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void onesamplevariancetest(/* Real    */ ae_vector* x, ae_int_t n, double variance, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
  void wilcoxonsignedranktest(/* Real    */ ae_vector* x, ae_int_t n, double e, double* bothtails, double* lefttail, double* righttail, ae_state *_state);
}
