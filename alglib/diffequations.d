module alglib.diffequations;
import alglib.ap;
import alglib.internal;

struct odesolverstate
{
    ae_int_t n;
    ae_int_t m;
    double xscale;
    double h;
    double eps;
    ae_bool fraceps;
    ae_vector yc;
    ae_vector escale;
    ae_vector xg;
    ae_int_t solvertype;
    ae_bool needdy;
    double x;
    ae_vector y;
    ae_vector dy;
    ae_matrix ytbl;
    ae_int_t repterminationtype;
    ae_int_t repnfev;
    ae_vector yn;
    ae_vector yns;
    ae_vector rka;
    ae_vector rkc;
    ae_vector rkcs;
    ae_matrix rkb;
    ae_matrix rkk;
    rcommstate rstate;
}

struct odesolverreport
{
    ae_int_t nfev;
    ae_int_t terminationtype;
}

/**
extern(C++,alglib)
{
  class _odesolverstate_owner
  {
      _odesolverstate_owner();
      _odesolverstate_owner(const _odesolverstate_owner &rhs);
      _odesolverstate_owner& operator=(const _odesolverstate_owner &rhs);
      virtual ~_odesolverstate_owner();
      alglib_impl::odesolverstate* c_ptr();
      alglib_impl::odesolverstate* c_ptr() const;
    protected:
        alglib_impl::odesolverstate *p_struct;
  }

  class odesolverstate : _odesolverstate_owner
  {
      odesolverstate();
      odesolverstate(const odesolverstate &rhs);
      odesolverstate& operator=(const odesolverstate &rhs);
      virtual ~odesolverstate();
      ae_bool &needdy;
      real_1d_array y;
      real_1d_array dy;
      double &x;
  }

  class _odesolverreport_owner
  {
      _odesolverreport_owner();
      _odesolverreport_owner(const _odesolverreport_owner &rhs);
      _odesolverreport_owner& operator=(const _odesolverreport_owner &rhs);
      virtual ~_odesolverreport_owner();
      alglib_impl::odesolverreport* c_ptr();
      alglib_impl::odesolverreport* c_ptr() const;
    protected:
        alglib_impl::odesolverreport *p_struct;
  }

  class odesolverreport : _odesolverreport_owner
  {
      odesolverreport();
      odesolverreport(const odesolverreport &rhs);
      odesolverreport& operator=(const odesolverreport &rhs);
      virtual ~odesolverreport();
      ae_int_t &nfev;
      ae_int_t &terminationtype;
  }
  void odesolverrkck(const real_1d_array &y, const ae_int_t n, const real_1d_array &x, const ae_int_t m, const double eps, const double h, odesolverstate &state);
  void odesolverrkck(const real_1d_array &y, const real_1d_array &x, const double eps, const double h, odesolverstate &state);
  bool odesolveriteration(const odesolverstate &state);
  void odesolversolve(odesolverstate &state, void (*diff)(const real_1d_array &y, double x, real_1d_array &dy, void *ptr), void *ptr = NULL);
  void odesolverresults(const odesolverstate &state, ae_int_t &m, real_1d_array &xtbl, real_2d_array &ytbl, odesolverreport &rep);
}
*/

extern(C)
{
  void odesolverrkck(/* Real    */ ae_vector* y, ae_int_t n, /* Real    */ ae_vector* x, ae_int_t m, double eps, double h, odesolverstate* state, ae_state *_state);
  ae_bool odesolveriteration(odesolverstate* state, ae_state *_state);
  void odesolverresults(odesolverstate* state, ae_int_t* m, /* Real    */ ae_vector* xtbl, /* Real    */ ae_matrix* ytbl, odesolverreport* rep, ae_state *_state);
  void _odesolverstate_init(void* _p, ae_state *_state);
  void _odesolverstate_init_copy(void* _dst, void* _src, ae_state *_state);
  void _odesolverstate_clear(void* _p);
  void _odesolverstate_destroy(void* _p);
  void _odesolverreport_init(void* _p, ae_state *_state);
  void _odesolverreport_init_copy(void* _dst, void* _src, ae_state *_state);
  void _odesolverreport_clear(void* _p);
  void _odesolverreport_destroy(void* _p);
}