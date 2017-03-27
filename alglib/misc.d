module alglib.misc;
import alglib.ap;
import alglib.internal;

struct hqrndstate
{
    ae_int_t s1;
    ae_int_t s2;
    ae_int_t magicv;
}

struct kdtree
{
    ae_int_t n;
    ae_int_t nx;
    ae_int_t ny;
    ae_int_t normtype;
    ae_matrix xy;
    ae_vector tags;
    ae_vector boxmin;
    ae_vector boxmax;
    ae_vector nodes;
    ae_vector splits;
    ae_vector x;
    ae_int_t kneeded;
    double rneeded;
    ae_bool selfmatch;
    double approxf;
    ae_int_t kcur;
    ae_vector idx;
    ae_vector r;
    ae_vector buf;
    ae_vector curboxmin;
    ae_vector curboxmax;
    double curdist;
    ae_int_t debugcounter;
}

struct xdebugrecord1
{
    ae_int_t i;
    ae_complex c;
    ae_vector a;
} 

/**
extern(C++,alglib)
{
  class _hqrndstate_owner
  {
  public:
      _hqrndstate_owner();
      _hqrndstate_owner(const _hqrndstate_owner &rhs);
      _hqrndstate_owner& operator=(const _hqrndstate_owner &rhs);
      virtual ~_hqrndstate_owner();
      alglib_impl::hqrndstate* c_ptr();
      alglib_impl::hqrndstate* c_ptr() const;
  protected:
      alglib_impl::hqrndstate *p_struct;
  };
  class hqrndstate : public _hqrndstate_owner
  {
  public:
      hqrndstate();
      hqrndstate(const hqrndstate &rhs);
      hqrndstate& operator=(const hqrndstate &rhs);
      virtual ~hqrndstate();

  };

  class _kdtree_owner
  {
  public:
      _kdtree_owner();
      _kdtree_owner(const _kdtree_owner &rhs);
      _kdtree_owner& operator=(const _kdtree_owner &rhs);
      virtual ~_kdtree_owner();
      alglib_impl::kdtree* c_ptr();
      alglib_impl::kdtree* c_ptr() const;
  protected:
      alglib_impl::kdtree *p_struct;
  };
  class kdtree : public _kdtree_owner
  {
  public:
      kdtree();
      kdtree(const kdtree &rhs);
      kdtree& operator=(const kdtree &rhs);
      virtual ~kdtree();

  }

  class _xdebugrecord1_owner
  {
      _xdebugrecord1_owner();
      _xdebugrecord1_owner(const _xdebugrecord1_owner &rhs);
      _xdebugrecord1_owner& operator=(const _xdebugrecord1_owner &rhs);
      virtual ~_xdebugrecord1_owner();
      alglib_impl::xdebugrecord1* c_ptr();
      alglib_impl::xdebugrecord1* c_ptr() const;
    protected:
        alglib_impl::xdebugrecord1 *p_struct;
  }

  class xdebugrecord1 :  _xdebugrecord1_owner
  {
      xdebugrecord1();
      xdebugrecord1(const xdebugrecord1 &rhs);
      xdebugrecord1& operator=(const xdebugrecord1 &rhs);
      virtual ~xdebugrecord1();
      ae_int_t &i;
      alglib::complex &c;
      real_1d_array a;
  }

  void hqrndrandomize(hqrndstate &state);
  void hqrndseed(const ae_int_t s1, const ae_int_t s2, hqrndstate &state);
double hqrnduniformr(const hqrndstate &state);
  ae_int_t hqrnduniformi(const hqrndstate &state, const ae_int_t n);
  double hqrndnormal(const hqrndstate &state);
  void hqrndunit2(const hqrndstate &state, double &x, double &y);
  void hqrndnormal2(const hqrndstate &state, double &x1, double &x2);
  double hqrndexponential(const hqrndstate &state, const double lambdav);
  double hqrnddiscrete(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);
  double hqrndcontinuous(const hqrndstate &state, const real_1d_array &x, const ae_int_t n);
  void kdtreeserialize(kdtree &obj, std::string &s_out);
  void kdtreeunserialize(std::string &s_in, kdtree &obj);
  void kdtreebuild(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
  void kdtreebuild(const real_2d_array &xy, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
  void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t n, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
  void kdtreebuildtagged(const real_2d_array &xy, const integer_1d_array &tags, const ae_int_t nx, const ae_int_t ny, const ae_int_t normtype, kdtree &kdt);
  ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch);
  ae_int_t kdtreequeryknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k);
  ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r, const bool selfmatch);
  ae_int_t kdtreequeryrnn(const kdtree &kdt, const real_1d_array &x, const double r);
  ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const bool selfmatch, const double eps);
  ae_int_t kdtreequeryaknn(const kdtree &kdt, const real_1d_array &x, const ae_int_t k, const double eps);
  void kdtreequeryresultsx(const kdtree &kdt, real_2d_array &x);
  void kdtreequeryresultsxy(const kdtree &kdt, real_2d_array &xy);
  void kdtreequeryresultstags(const kdtree &kdt, integer_1d_array &tags);
  void kdtreequeryresultsdistances(const kdtree &kdt, real_1d_array &r);
  void kdtreequeryresultsxi(const kdtree &kdt, real_2d_array &x);
  void kdtreequeryresultsxyi(const kdtree &kdt, real_2d_array &xy);
  void kdtreequeryresultstagsi(const kdtree &kdt, integer_1d_array &tags);
  void kdtreequeryresultsdistancesi(const kdtree &kdt, real_1d_array &r);
  void xdebuginitrecord1(xdebugrecord1 &rec1);
  ae_int_t xdebugb1count(const boolean_1d_array &a);
  void xdebugb1not(const boolean_1d_array &a);
  void xdebugb1appendcopy(boolean_1d_array &a);
  void xdebugb1outeven(const ae_int_t n, boolean_1d_array &a);
  ae_int_t xdebugi1sum(const integer_1d_array &a);
  void xdebugi1neg(const integer_1d_array &a);
  void xdebugi1appendcopy(integer_1d_array &a);
  void xdebugi1outeven(const ae_int_t n, integer_1d_array &a);
  double xdebugr1sum(const real_1d_array &a);
  void xdebugr1neg(const real_1d_array &a);
  void xdebugr1appendcopy(real_1d_array &a);
  void xdebugr1outeven(const ae_int_t n, real_1d_array &a);
  alglib::complex xdebugc1sum(const complex_1d_array &a);
  void xdebugc1neg(const complex_1d_array &a);
  void xdebugc1appendcopy(complex_1d_array &a);
  void xdebugc1outeven(const ae_int_t n, complex_1d_array &a);
  ae_int_t xdebugb2count(const boolean_2d_array &a);
  void xdebugb2not(const boolean_2d_array &a);
  void xdebugb2transpose(boolean_2d_array &a);
  void xdebugb2outsin(const ae_int_t m, const ae_int_t n, boolean_2d_array &a);
  ae_int_t xdebugi2sum(const integer_2d_array &a);
  void xdebugi2neg(const integer_2d_array &a);
  void xdebugi2transpose(integer_2d_array &a);
  void xdebugi2outsin(const ae_int_t m, const ae_int_t n, integer_2d_array &a);
  double xdebugr2sum(const real_2d_array &a);
  void xdebugr2neg(const real_2d_array &a);
  void xdebugr2transpose(real_2d_array &a);
  void xdebugr2outsin(const ae_int_t m, const ae_int_t n, real_2d_array &a);
  alglib::complex xdebugc2sum(const complex_2d_array &a);
  void xdebugc2neg(const complex_2d_array &a);
  void xdebugc2transpose(complex_2d_array &a);
  void xdebugc2outsincos(const ae_int_t m, const ae_int_t n, complex_2d_array &a);
  double xdebugmaskedbiasedproductsum(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const real_2d_array &b, const boolean_2d_array &c);
}
*/

extern(C)
{
  void hqrndrandomize(hqrndstate* state, ae_state *_state);
  void hqrndseed(ae_int_t s1, ae_int_t s2, hqrndstate* state, ae_state *_state);
  double hqrnduniformr(hqrndstate* state, ae_state *_state);
  ae_int_t hqrnduniformi(hqrndstate* state, ae_int_t n, ae_state *_state);
  double hqrndnormal(hqrndstate* state, ae_state *_state);
  void hqrndunit2(hqrndstate* state, double* x, double* y, ae_state *_state);
  void hqrndnormal2(hqrndstate* state, double* x1, double* x2, ae_state *_state);
  double hqrndexponential(hqrndstate* state, double lambdav, ae_state *_state);
  double hqrnddiscrete(hqrndstate* state, /* Real    */ ae_vector* x, ae_int_t n, ae_state *_state);
  double hqrndcontinuous(hqrndstate* state, /* Real    */ ae_vector* x, ae_int_t n, ae_state *_state);
  void _hqrndstate_init(void* _p, ae_state *_state);
  void _hqrndstate_init_copy(void* _dst, void* _src, ae_state *_state);
  void _hqrndstate_clear(void* _p);
  void _hqrndstate_destroy(void* _p);
  void kdtreebuild(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree* kdt, ae_state *_state);
  void kdtreebuildtagged(/* Real    */ ae_matrix* xy, /* Integer */ ae_vector* tags, ae_int_t n, ae_int_t nx, ae_int_t ny, ae_int_t normtype, kdtree* kdt, ae_state *_state);
  ae_int_t kdtreequeryknn(kdtree* kdt, /* Real    */ ae_vector* x, ae_int_t k, ae_bool selfmatch, ae_state *_state);
  ae_int_t kdtreequeryrnn(kdtree* kdt, /* Real    */ ae_vector* x, double r, ae_bool selfmatch, ae_state *_state);
  ae_int_t kdtreequeryaknn(kdtree* kdt, /* Real    */ ae_vector* x, ae_int_t k, ae_bool selfmatch, double eps, ae_state *_state);
  void kdtreequeryresultsx(kdtree* kdt, /* Real    */ ae_matrix* x, ae_state *_state);
  void kdtreequeryresultsxy(kdtree* kdt, /* Real    */ ae_matrix* xy, ae_state *_state);
  void kdtreequeryresultstags(kdtree* kdt, /* Integer */ ae_vector* tags, ae_state *_state);
  void kdtreequeryresultsdistances(kdtree* kdt, /* Real    */ ae_vector* r, ae_state *_state);
  void kdtreequeryresultsxi(kdtree* kdt, /* Real    */ ae_matrix* x, ae_state *_state);
  void kdtreequeryresultsxyi(kdtree* kdt, /* Real    */ ae_matrix* xy, ae_state *_state);
  void kdtreequeryresultstagsi(kdtree* kdt, /* Integer */ ae_vector* tags, ae_state *_state);
  void kdtreequeryresultsdistancesi(kdtree* kdt, /* Real    */ ae_vector* r, ae_state *_state);
  void kdtreealloc(ae_serializer* s, kdtree* tree, ae_state *_state);
  void kdtreeserialize(ae_serializer* s, kdtree* tree, ae_state *_state);
  void kdtreeunserialize(ae_serializer* s, kdtree* tree, ae_state *_state);
  void _kdtree_init(void* _p, ae_state *_state);
  void _kdtree_init_copy(void* _dst, void* _src, ae_state *_state);
  void _kdtree_clear(void* _p);
  void _kdtree_destroy(void* _p);
  void xdebuginitrecord1(xdebugrecord1* rec1, ae_state *_state);
  ae_int_t xdebugb1count(/* Boolean */ ae_vector* a, ae_state *_state);
  void xdebugb1not(/* Boolean */ ae_vector* a, ae_state *_state);
  void xdebugb1appendcopy(/* Boolean */ ae_vector* a, ae_state *_state);
  void xdebugb1outeven(ae_int_t n, /* Boolean */ ae_vector* a, ae_state *_state);
  ae_int_t xdebugi1sum(/* Integer */ ae_vector* a, ae_state *_state);
  void xdebugi1neg(/* Integer */ ae_vector* a, ae_state *_state);
  void xdebugi1appendcopy(/* Integer */ ae_vector* a, ae_state *_state);
  void xdebugi1outeven(ae_int_t n, /* Integer */ ae_vector* a, ae_state *_state);
  double xdebugr1sum(/* Real    */ ae_vector* a, ae_state *_state);
  void xdebugr1neg(/* Real    */ ae_vector* a, ae_state *_state);
  void xdebugr1appendcopy(/* Real    */ ae_vector* a, ae_state *_state);
  void xdebugr1outeven(ae_int_t n, /* Real    */ ae_vector* a, ae_state *_state);
  ae_complex xdebugc1sum(/* Complex */ ae_vector* a, ae_state *_state);
  void xdebugc1neg(/* Complex */ ae_vector* a, ae_state *_state);
  void xdebugc1appendcopy(/* Complex */ ae_vector* a, ae_state *_state);
  void xdebugc1outeven(ae_int_t n, /* Complex */ ae_vector* a, ae_state *_state);
  ae_int_t xdebugb2count(/* Boolean */ ae_matrix* a, ae_state *_state);
  void xdebugb2not(/* Boolean */ ae_matrix* a, ae_state *_state);
  void xdebugb2transpose(/* Boolean */ ae_matrix* a, ae_state *_state);
  void xdebugb2outsin(ae_int_t m, ae_int_t n, /* Boolean */ ae_matrix* a, ae_state *_state);
  ae_int_t xdebugi2sum(/* Integer */ ae_matrix* a, ae_state *_state);
  void xdebugi2neg(/* Integer */ ae_matrix* a, ae_state *_state);
  void xdebugi2transpose(/* Integer */ ae_matrix* a, ae_state *_state);
  void xdebugi2outsin(ae_int_t m, ae_int_t n, /* Integer */ ae_matrix* a, ae_state *_state);
  double xdebugr2sum(/* Real    */ ae_matrix* a, ae_state *_state);
  void xdebugr2neg(/* Real    */ ae_matrix* a, ae_state *_state);
  void xdebugr2transpose(/* Real    */ ae_matrix* a, ae_state *_state);
  void xdebugr2outsin(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_state *_state);
  ae_complex xdebugc2sum(/* Complex */ ae_matrix* a, ae_state *_state);
  void xdebugc2neg(/* Complex */ ae_matrix* a, ae_state *_state);
  void xdebugc2transpose(/* Complex */ ae_matrix* a, ae_state *_state);
  void xdebugc2outsincos(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_state *_state);
  double xdebugmaskedbiasedproductsum(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, /* Real    */ ae_matrix* b, /* Boolean */ ae_matrix* c, ae_state *_state);
  void _xdebugrecord1_init(void* _p, ae_state *_state);
  void _xdebugrecord1_init_copy(void* _dst, void* _src, ae_state *_state);
  void _xdebugrecord1_clear(void* _p);
  void _xdebugrecord1_destroy(void* _p);
}
