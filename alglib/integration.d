module alglib.integration;
import alglib.ap;
import alglib.internal;
import alglib.linalg;
import alglib.specialfunctions;

struct autogkreport
{
    ae_int_t terminationtype;
    ae_int_t nfev;
    ae_int_t nintervals;
}

struct autogkinternalstate
{
    double a;
    double b;
    double eps;
    double xwidth;
    double x;
    double f;
    ae_int_t info;
    double r;
    ae_matrix heap;
    ae_int_t heapsize;
    ae_int_t heapwidth;
    ae_int_t heapused;
    double sumerr;
    double sumabs;
    ae_vector qn;
    ae_vector wg;
    ae_vector wk;
    ae_vector wr;
    ae_int_t n;
    rcommstate rstate;
}

struct autogkstate
{
    double a;
    double b;
    double alpha;
    double beta;
    double xwidth;
    double x;
    double xminusa;
    double bminusx;
    ae_bool needf;
    double f;
    ae_int_t wrappermode;
    autogkinternalstate internalstate;
    rcommstate rstate;
    double v;
    ae_int_t terminationtype;
    ae_int_t nfev;
    ae_int_t nintervals;
}

/**
extern(C++,alglib)
{
    class _autogkreport_owner
    {
        _autogkreport_owner();
        _autogkreport_owner(const _autogkreport_owner &rhs);
        _autogkreport_owner& operator=(const _autogkreport_owner &rhs);
        virtual ~_autogkreport_owner();
        alglib_impl::autogkreport* c_ptr();
        alglib_impl::autogkreport* c_ptr() const;
        protected:
            alglib_impl::autogkreport *p_struct;
    }

    class autogkreport :  _autogkreport_owner
    {
        autogkreport();
        autogkreport(const autogkreport &rhs);
        autogkreport& operator=(const autogkreport &rhs);
        virtual ~autogkreport();
        ae_int_t &terminationtype;
        ae_int_t &nfev;
        ae_int_t &nintervals;
    }

    class _autogkstate_owner
    {
        _autogkstate_owner();
        _autogkstate_owner(const _autogkstate_owner &rhs);
        _autogkstate_owner& operator=(const _autogkstate_owner &rhs);
        virtual ~_autogkstate_owner();
        alglib_impl::autogkstate* c_ptr();
        alglib_impl::autogkstate* c_ptr() const;
        protected:
            alglib_impl::autogkstate *p_struct;
    }

    class autogkstate : _autogkstate_owner
    {
        autogkstate();
        autogkstate(const autogkstate &rhs);
        autogkstate& operator=(const autogkstate &rhs);
        virtual ~autogkstate();
        ae_bool &needf;
        double &x;
        double &xminusa;
        double &bminusx;
        double &f;
    }

    void gqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
    void gqgenerategausslobattorec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const double b, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
    void gqgenerategaussradaurec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const double a, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
    void gqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
    void gqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &w);
    void gqgenerategausslaguerre(const ae_int_t n, const double alpha, ae_int_t &info, real_1d_array &x, real_1d_array &w);
    void gqgenerategausshermite(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &w);
    void gkqgeneraterec(const real_1d_array &alpha, const real_1d_array &beta, const double mu0, const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
    void gkqgenerategausslegendre(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
    void gkqgenerategaussjacobi(const ae_int_t n, const double alpha, const double beta, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
    void gkqlegendrecalc(const ae_int_t n, ae_int_t &info, real_1d_array &x, real_1d_array &wkronrod, real_1d_array &wgauss);
    void autogksmooth(const double a, const double b, autogkstate &state);
    void autogksmoothw(const double a, const double b, const double xwidth, autogkstate &state);
    void autogksingular(const double a, const double b, const double alpha, const double beta, autogkstate &state);
    bool autogkiteration(const autogkstate &state);
    void autogkintegrate(autogkstate &state, void* function(double x, double xminusa, double bminusx, double &y, void *ptr) func, void *ptr = NULL);
    void autogkresults(const autogkstate &state, double &v, autogkreport &rep);
}
*/

extern(C)
{
    void gqgeneraterec(/* Real    */ ae_vector* alpha, /* Real    */ ae_vector* beta, double mu0, ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* w, ae_state *_state);
    void gqgenerategausslobattorec(/* Real    */ ae_vector* alpha, /* Real    */ ae_vector* beta, double mu0, double a, double b, ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* w, ae_state *_state);
    void gqgenerategaussradaurec(/* Real    */ ae_vector* alpha, /* Real    */ ae_vector* beta, double mu0, double a, ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* w, ae_state *_state);
    void gqgenerategausslegendre(ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* w, ae_state *_state);
    void gqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* w, ae_state *_state);
    void gqgenerategausslaguerre(ae_int_t n, double alpha, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* w, ae_state *_state);
    void gqgenerategausshermite(ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* w, ae_state *_state);
    void gkqgeneraterec(/* Real    */ ae_vector* alpha, /* Real    */ ae_vector* beta, double mu0, ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* wkronrod, /* Real    */ ae_vector* wgauss, ae_state *_state);
    void gkqgenerategausslegendre(ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* wkronrod, /* Real    */ ae_vector* wgauss, ae_state *_state);
    void gkqgenerategaussjacobi(ae_int_t n, double alpha, double beta, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* wkronrod, /* Real    */ ae_vector* wgauss, ae_state *_state);
    void gkqlegendrecalc(ae_int_t n, ae_int_t* info, /* Real    */ ae_vector* x, /* Real    */ ae_vector* wkronrod, /* Real    */ ae_vector* wgauss, ae_state *_state);
    void gkqlegendretbl(ae_int_t n, /* Real    */ ae_vector* x, /* Real    */ ae_vector* wkronrod, /* Real    */ ae_vector* wgauss, double* eps, ae_state *_state);
    void autogksmooth(double a, double b, autogkstate* state, ae_state *_state);
    void autogksmoothw(double a, double b, double xwidth, autogkstate* state, ae_state *_state);
    void autogksingular(double a, double b, double alpha, double beta, autogkstate* state, ae_state *_state);
    ae_bool autogkiteration(autogkstate* state, ae_state *_state);
    void autogkresults(autogkstate* state, double* v, autogkreport* rep, ae_state *_state);
    void _autogkreport_init(void* _p, ae_state *_state);
    void _autogkreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _autogkreport_clear(void* _p);
    void _autogkreport_destroy(void* _p);
    void _autogkinternalstate_init(void* _p, ae_state *_state);
    void _autogkinternalstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _autogkinternalstate_clear(void* _p);
    void _autogkinternalstate_destroy(void* _p);
    void _autogkstate_init(void* _p, ae_state *_state);
    void _autogkstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _autogkstate_clear(void* _p);
    void _autogkstate_destroy(void* _p);
}

