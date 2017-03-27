module alglib.solvers;
import alglib.ap;
import alglib.internal;
import alglib.linalg;
import alglib.misc;

struct densesolverreport
{
    double r1;
    double rinf;
}

struct densesolverlsreport
{
    double r2;
    ae_matrix cx;
    ae_int_t n;
    ae_int_t k;
}

struct linlsqrstate
{
    normestimatorstate nes;
    ae_vector rx;
    ae_vector b;
    ae_int_t n;
    ae_int_t m;
    ae_int_t prectype;
    ae_vector ui;
    ae_vector uip1;
    ae_vector vi;
    ae_vector vip1;
    ae_vector omegai;
    ae_vector omegaip1;
    double alphai;
    double alphaip1;
    double betai;
    double betaip1;
    double phibari;
    double phibarip1;
    double phii;
    double rhobari;
    double rhobarip1;
    double rhoi;
    double ci;
    double si;
    double theta;
    double lambdai;
    ae_vector d;
    double anorm;
    double bnorm2;
    double dnorm;
    double r2;
    ae_vector x;
    ae_vector mv;
    ae_vector mtv;
    double epsa;
    double epsb;
    double epsc;
    ae_int_t maxits;
    ae_bool xrep;
    ae_bool xupdated;
    ae_bool needmv;
    ae_bool needmtv;
    ae_bool needmv2;
    ae_bool needvmv;
    ae_bool needprec;
    ae_int_t repiterationscount;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    ae_bool running;
    ae_vector tmpd;
    ae_vector tmpx;
    rcommstate rstate;
}

struct linlsqrreport
{
    ae_int_t iterationscount;
    ae_int_t nmv;
    ae_int_t terminationtype;
}

struct lincgstate
{
    ae_vector rx;
    ae_vector b;
    ae_int_t n;
    ae_int_t prectype;
    ae_vector cx;
    ae_vector cr;
    ae_vector cz;
    ae_vector p;
    ae_vector r;
    ae_vector z;
    double alpha;
    double beta;
    double r2;
    double meritfunction;
    ae_vector x;
    ae_vector mv;
    ae_vector pv;
    double vmv;
    ae_vector startx;
    double epsf;
    ae_int_t maxits;
    ae_int_t itsbeforerestart;
    ae_int_t itsbeforerupdate;
    ae_bool xrep;
    ae_bool xupdated;
    ae_bool needmv;
    ae_bool needmtv;
    ae_bool needmv2;
    ae_bool needvmv;
    ae_bool needprec;
    ae_int_t repiterationscount;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    ae_bool running;
    ae_vector tmpd;
    rcommstate rstate;
}

struct lincgreport
{
    ae_int_t iterationscount;
    ae_int_t nmv;
    ae_int_t terminationtype;
    double r2;
}

struct nleqstate
{
    ae_int_t n;
    ae_int_t m;
    double epsf;
    ae_int_t maxits;
    ae_bool xrep;
    double stpmax;
    ae_vector x;
    double f;
    ae_vector fi;
    ae_matrix j;
    ae_bool needf;
    ae_bool needfij;
    ae_bool xupdated;
    rcommstate rstate;
    ae_int_t repiterationscount;
    ae_int_t repnfunc;
    ae_int_t repnjac;
    ae_int_t repterminationtype;
    ae_vector xbase;
    double fbase;
    double fprev;
    ae_vector candstep;
    ae_vector rightpart;
    ae_vector cgbuf;
}

struct nleqreport
{
    ae_int_t iterationscount;
    ae_int_t nfunc;
    ae_int_t njac;
    ae_int_t terminationtype;
} 

struct polynomialsolverreport
{
    double maxerr;
}

/**
extern(C++,alglib)
{
    class _densesolverreport_owner
    {
        _densesolverreport_owner();
        _densesolverreport_owner(const _densesolverreport_owner &rhs);
        _densesolverreport_owner& operator=(const _densesolverreport_owner &rhs);
        virtual ~_densesolverreport_owner();
        alglib_impl::densesolverreport* c_ptr();
        alglib_impl::densesolverreport* c_ptr() const;
        protected:
            alglib_impl::densesolverreport *p_struct;
    }

    class densesolverreport :  _densesolverreport_owner
    {
        densesolverreport();
        densesolverreport(const densesolverreport &rhs);
        densesolverreport& operator=(const densesolverreport &rhs);
        virtual ~densesolverreport();
        double &r1;
        double &rinf;
    }

    class _densesolverlsreport_owner
    {
        _densesolverlsreport_owner();
        _densesolverlsreport_owner(const _densesolverlsreport_owner &rhs);
        _densesolverlsreport_owner& operator=(const _densesolverlsreport_owner &rhs);
        virtual ~_densesolverlsreport_owner();
        alglib_impl::densesolverlsreport* c_ptr();
        alglib_impl::densesolverlsreport* c_ptr() const;
        protected:
            alglib_impl::densesolverlsreport *p_struct;
    }

    class densesolverlsreport :  _densesolverlsreport_owner
    {
        densesolverlsreport();
        densesolverlsreport(const densesolverlsreport &rhs);
        densesolverlsreport& operator=(const densesolverlsreport &rhs);
        virtual ~densesolverlsreport();
        double &r2;
        real_2d_array cx;
        ae_int_t &n;
        ae_int_t &k;
    }
    
    class _linlsqrstate_owner
    {
        _linlsqrstate_owner();
        _linlsqrstate_owner(const _linlsqrstate_owner &rhs);
        _linlsqrstate_owner& operator=(const _linlsqrstate_owner &rhs);
        virtual ~_linlsqrstate_owner();
        alglib_impl::linlsqrstate* c_ptr();
        alglib_impl::linlsqrstate* c_ptr() const;
        protected:
            alglib_impl::linlsqrstate *p_struct;
    }

    class linlsqrstate :  _linlsqrstate_owner
    {
        linlsqrstate();
        linlsqrstate(const linlsqrstate &rhs);
        linlsqrstate& operator=(const linlsqrstate &rhs);
        virtual ~linlsqrstate();
    }

    class _linlsqrreport_owner
    {
        _linlsqrreport_owner();
        _linlsqrreport_owner(const _linlsqrreport_owner &rhs);
        _linlsqrreport_owner& operator=(const _linlsqrreport_owner &rhs);
        virtual ~_linlsqrreport_owner();
        alglib_impl::linlsqrreport* c_ptr();
        alglib_impl::linlsqrreport* c_ptr() const;
        protected:
            alglib_impl::linlsqrreport *p_struct;
    }

    class linlsqrreport : _linlsqrreport_owner
    {
        linlsqrreport();
        linlsqrreport(const linlsqrreport &rhs);
        linlsqrreport& operator=(const linlsqrreport &rhs);
        virtual ~linlsqrreport();
        ae_int_t &iterationscount;
        ae_int_t &nmv;
        ae_int_t &terminationtype;
    }

    class _lincgstate_owner
    {
        _lincgstate_owner();
        _lincgstate_owner(const _lincgstate_owner &rhs);
        _lincgstate_owner& operator=(const _lincgstate_owner &rhs);
        virtual ~_lincgstate_owner();
        alglib_impl::lincgstate* c_ptr();
        alglib_impl::lincgstate* c_ptr() const;
        protected:
            alglib_impl::lincgstate *p_struct;
    }

    class lincgstate :  _lincgstate_owner
    {
        lincgstate();
        lincgstate(const lincgstate &rhs);
        lincgstate& operator=(const lincgstate &rhs);
        virtual ~lincgstate();
    }

    class _lincgreport_owner
    {
        _lincgreport_owner();
        _lincgreport_owner(const _lincgreport_owner &rhs);
        _lincgreport_owner& operator=(const _lincgreport_owner &rhs);
        virtual ~_lincgreport_owner();
        alglib_impl::lincgreport* c_ptr();
        alglib_impl::lincgreport* c_ptr() const;
        protected:
            alglib_impl::lincgreport *p_struct;
    }

    class lincgreport :  _lincgreport_owner
    {
        lincgreport();
        lincgreport(const lincgreport &rhs);
        lincgreport& operator=(const lincgreport &rhs);
        virtual ~lincgreport();
        ae_int_t &iterationscount;
        ae_int_t &nmv;
        ae_int_t &terminationtype;
        double &r2;
    }

    class _nleqstate_owner
    {
        _nleqstate_owner();
        _nleqstate_owner(const _nleqstate_owner &rhs);
        _nleqstate_owner& operator=(const _nleqstate_owner &rhs);
        virtual ~_nleqstate_owner();
        alglib_impl::nleqstate* c_ptr();
        alglib_impl::nleqstate* c_ptr() const;
        protected:
            alglib_impl::nleqstate *p_struct;
    }

    class nleqstate : _nleqstate_owner
    {
        nleqstate();
        nleqstate(const nleqstate &rhs);
        nleqstate& operator=(const nleqstate &rhs);
        virtual ~nleqstate();
        ae_bool &needf;
        ae_bool &needfij;
        ae_bool &xupdated;
        double &f;
        real_1d_array fi;
        real_2d_array j;
        real_1d_array x;
    }

    class _nleqreport_owner
    {
        _nleqreport_owner();
        _nleqreport_owner(const _nleqreport_owner &rhs);
        _nleqreport_owner& operator=(const _nleqreport_owner &rhs);
        virtual ~_nleqreport_owner();
        alglib_impl::nleqreport* c_ptr();
        alglib_impl::nleqreport* c_ptr() const;
        protected:
            alglib_impl::nleqreport *p_struct;
    }

    class nleqreport :  _nleqreport_owner
    {
        nleqreport();
        nleqreport(const nleqreport &rhs);
        nleqreport& operator=(const nleqreport &rhs);
        virtual ~nleqreport();
        ae_int_t &iterationscount;
        ae_int_t &nfunc;
        ae_int_t &njac;
        ae_int_t &terminationtype;
    }

    class _polynomialsolverreport_owner
    {
        _polynomialsolverreport_owner();
        _polynomialsolverreport_owner(const _polynomialsolverreport_owner &rhs);
        _polynomialsolverreport_owner& operator=(const _polynomialsolverreport_owner &rhs);
        virtual ~_polynomialsolverreport_owner();
        alglib_impl::polynomialsolverreport* c_ptr();
        alglib_impl::polynomialsolverreport* c_ptr() const;
        protected:
            alglib_impl::polynomialsolverreport *p_struct;
    }

    class polynomialsolverreport :  _polynomialsolverreport_owner
    {
        polynomialsolverreport();
        polynomialsolverreport(const polynomialsolverreport &rhs);
        polynomialsolverreport& operator=(const polynomialsolverreport &rhs);
        virtual ~polynomialsolverreport();
        double &maxerr;
    }

    void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
    void smp_rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
    void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
    void smp_rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
    void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void smp_rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
    void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info);
    void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void smp_rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
    void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void smp_cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
    void smp_cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
    void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
    void smp_cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
    void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void smp_cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
    void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info);
    void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
    void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void smp_spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
    void smp_spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
    void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
    void smp_spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
    void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void smp_spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x);
    void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info);
    void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x);
    void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info);
    void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void smp_hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
    void smp_hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
    void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
    void smp_hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
    void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void smp_hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x);
    void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void smp_hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info);
    void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x);
    void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info);
    void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x);
    void smp_rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x);
    void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state);
    void linlsqrsetprecunit(const linlsqrstate &state);
    void linlsqrsetprecdiag(const linlsqrstate &state);
    void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai);
    void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b);
    void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits);
    void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep);
    void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep);
    void lincgcreate(const ae_int_t n, lincgstate &state);
    void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x);
    void lincgsetprecunit(const lincgstate &state);
    void lincgsetprecdiag(const lincgstate &state);
    void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits);
    void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b);
    void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep);
    void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf);
   void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq);
    void lincgsetxrep(const ref lincgstate state, const bool needxrep);
    void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state);
    void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state);
    void nleqsetcond(const ref nleqstate state, const double epsf, const ae_int_t maxits);
    void nleqsetxrep(const ref nleqstate state, const bool needxrep);
    void nleqsetstpmax(const ref nleqstate state, const double stpmax);
    bool nleqiteration(const ref nleqstate state);
    void nleqsolve(ref nleqstate state,
        void* function(const ref real_1d_array x, ref double func, void *ptr) func,
        void* function(const ref real_1d_array x, ref real_1d_array fi, ref real_2d_array jac, void *ptr) jac,
        void* function(const ref real_1d_array x, double func, void *ptr) rep = null,
        void *ptr = null);
    void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep);
    void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep);
    void nleqrestartfrom(const nleqstate &state, const real_1d_array &x);
    void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep);
}
*/
extern(C)
{
    void rmatrixsolve(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void _pexec_rmatrixsolve(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void rmatrixsolvefast(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void _pexec_rmatrixsolvefast(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void rmatrixsolvem(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_bool rfs, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void _pexec_rmatrixsolvem(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_bool rfs, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void rmatrixsolvemfast(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_rmatrixsolvemfast(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void rmatrixlusolve(/* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void rmatrixlusolvefast(/* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void rmatrixlusolvem(/* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void _pexec_rmatrixlusolvem(/* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void rmatrixlusolvemfast(/* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_rmatrixlusolvemfast(/* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void rmatrixmixedsolve(/* Real    */ ae_matrix* a, /* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void rmatrixmixedsolvem(/* Real    */ ae_matrix* a, /* Real    */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void cmatrixsolvem(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_bool rfs, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void _pexec_cmatrixsolvem(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_bool rfs, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void cmatrixsolvemfast(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_cmatrixsolvemfast(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void cmatrixsolve(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_vector* x, ae_state *_state);
    void _pexec_cmatrixsolve(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_vector* x, ae_state *_state);
    void cmatrixsolvefast(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void _pexec_cmatrixsolvefast(/* Complex */ ae_matrix* a, ae_int_t n, /* Complex */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void cmatrixlusolvem(/* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void _pexec_cmatrixlusolvem(/* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void cmatrixlusolvemfast(/* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_cmatrixlusolvemfast(/* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void cmatrixlusolve(/* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_vector* x, ae_state *_state);
    void cmatrixlusolvefast(/* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void cmatrixmixedsolvem(/* Complex */ ae_matrix* a, /* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void cmatrixmixedsolve(/* Complex */ ae_matrix* a, /* Complex */ ae_matrix* lua, /* Integer */ ae_vector* p, ae_int_t n, /* Complex */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_vector* x, ae_state *_state);
    void spdmatrixsolvem(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void _pexec_spdmatrixsolvem(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void spdmatrixsolvemfast(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_spdmatrixsolvemfast(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void spdmatrixsolve(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void _pexec_spdmatrixsolve(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void spdmatrixsolvefast(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void _pexec_spdmatrixsolvefast(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void spdmatrixcholeskysolvem(/* Real    */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void _pexec_spdmatrixcholeskysolvem(/* Real    */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_matrix* x, ae_state *_state);
    void spdmatrixcholeskysolvemfast(/* Real    */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_spdmatrixcholeskysolvemfast(/* Real    */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Real    */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void spdmatrixcholeskysolve(/* Real    */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void spdmatrixcholeskysolvefast(/* Real    */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void hpdmatrixsolvem(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void _pexec_hpdmatrixsolvem(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void hpdmatrixsolvemfast(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_hpdmatrixsolvemfast(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void hpdmatrixsolve(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_vector* x, ae_state *_state);
    void _pexec_hpdmatrixsolve(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_vector* x, ae_state *_state);
    void hpdmatrixsolvefast(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void _pexec_hpdmatrixsolvefast(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* b, ae_int_t* info, ae_state *_state); 
    void hpdmatrixcholeskysolvem(/* Complex */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void _pexec_hpdmatrixcholeskysolvem(/* Complex */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_matrix* x, ae_state *_state);
    void hpdmatrixcholeskysolvemfast(/* Complex */ ae_matrix* cha,
         ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void _pexec_hpdmatrixcholeskysolvemfast(/* Complex */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Complex */ ae_matrix* b, ae_int_t m, ae_int_t* info, ae_state *_state);
    void hpdmatrixcholeskysolve(/* Complex */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* b, ae_int_t* info, densesolverreport* rep, /* Complex */ ae_vector* x, ae_state *_state);
    void hpdmatrixcholeskysolvefast(/* Complex */ ae_matrix* cha, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* b, ae_int_t* info, ae_state *_state);
    void rmatrixsolvels(/* Real    */ ae_matrix* a, ae_int_t nrows, ae_int_t ncols, /* Real    */ ae_vector* b, double threshold, ae_int_t* info, densesolverlsreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void _pexec_rmatrixsolvels(/* Real    */ ae_matrix* a, ae_int_t nrows, ae_int_t ncols, /* Real    */ ae_vector* b, double threshold, ae_int_t* info, densesolverlsreport* rep, /* Real    */ ae_vector* x, ae_state *_state);
    void _densesolverreport_init(void* _p, ae_state *_state);
    void _densesolverreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _densesolverreport_clear(void* _p);
    void _densesolverreport_destroy(void* _p);
    void _densesolverlsreport_init(void* _p, ae_state *_state);
    void _densesolverlsreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _densesolverlsreport_clear(void* _p);
    void _densesolverlsreport_destroy(void* _p);
    void linlsqrcreate(ae_int_t m, ae_int_t n, linlsqrstate* state, ae_state *_state);
    void linlsqrsetb(linlsqrstate* state, /* Real    */ ae_vector* b, ae_state *_state);
    void linlsqrsetprecunit(linlsqrstate* state, ae_state *_state);
    void linlsqrsetprecdiag(linlsqrstate* state, ae_state *_state);
    void linlsqrsetlambdai(linlsqrstate* state, double lambdai, ae_state *_state);
    ae_bool linlsqriteration(linlsqrstate* state, ae_state *_state);
    void linlsqrsolvesparse(linlsqrstate* state, sparsematrix* a, /* Real    */ ae_vector* b, ae_state *_state);
    void linlsqrsetcond(linlsqrstate* state, double epsa, double epsb, ae_int_t maxits, ae_state *_state);
    void linlsqrresults(linlsqrstate* state, /* Real    */ ae_vector* x, linlsqrreport* rep, ae_state *_state);
    void linlsqrsetxrep(linlsqrstate* state, ae_bool needxrep, ae_state *_state);
    void linlsqrrestart(linlsqrstate* state, ae_state *_state);
    void _linlsqrstate_init(void* _p, ae_state *_state);
    void _linlsqrstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _linlsqrstate_clear(void* _p);
    void _linlsqrstate_destroy(void* _p);
    void _linlsqrreport_init(void* _p, ae_state *_state);
    void _linlsqrreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _linlsqrreport_clear(void* _p);
    void _linlsqrreport_destroy(void* _p);
    void lincgcreate(ae_int_t n, lincgstate* state, ae_state *_state);
    void lincgsetstartingpoint(lincgstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void lincgsetb(lincgstate* state, /* Real    */ ae_vector* b, ae_state *_state);
    void lincgsetprecunit(lincgstate* state, ae_state *_state);
    void lincgsetprecdiag(lincgstate* state, ae_state *_state);
    void lincgsetcond(lincgstate* state, double epsf, ae_int_t maxits, ae_state *_state);
    ae_bool lincgiteration(lincgstate* state, ae_state *_state);
    void lincgsolvesparse(lincgstate* state, sparsematrix* a, ae_bool isupper, /* Real    */ ae_vector* b, ae_state *_state);
    void lincgresults(lincgstate* state, /* Real    */ ae_vector* x, lincgreport* rep, ae_state *_state);
    void lincgsetrestartfreq(lincgstate* state, ae_int_t srf, ae_state *_state);
    void lincgsetrupdatefreq(lincgstate* state, ae_int_t freq, ae_state *_state);
    void lincgsetxrep(lincgstate* state, ae_bool needxrep, ae_state *_state);
    void lincgrestart(lincgstate* state, ae_state *_state);
    void _lincgstate_init(void* _p, ae_state *_state);
    void _lincgstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _lincgstate_clear(void* _p);
    void _lincgstate_destroy(void* _p);
    void _lincgreport_init(void* _p, ae_state *_state);
    void _lincgreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _lincgreport_clear(void* _p);
    void _lincgreport_destroy(void* _p);
    void nleqcreatelm(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, nleqstate* state, ae_state *_state);
    void nleqsetcond(nleqstate* state, double epsf, ae_int_t maxits, ae_state *_state);
    void nleqsetxrep(nleqstate* state, ae_bool needxrep, ae_state *_state);
    void nleqsetstpmax(nleqstate* state, double stpmax, ae_state *_state);
    ae_bool nleqiteration(nleqstate* state, ae_state *_state);
    void nleqresults(nleqstate* state, /* Real    */ ae_vector* x, nleqreport* rep, ae_state *_state);
    void nleqresultsbuf(nleqstate* state, /* Real    */ ae_vector* x, nleqreport* rep, ae_state *_state);
    void nleqrestartfrom(nleqstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void _nleqstate_init(void* _p, ae_state *_state);
    void _nleqstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _nleqstate_clear(void* _p);
    void _nleqstate_destroy(void* _p);
    void _nleqreport_init(void* _p, ae_state *_state);
    void _nleqreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _nleqreport_clear(void* _p);
    void _nleqreport_destroy(void* _p);
    void polynomialsolve(/* Real    */ ae_vector* a, ae_int_t n, /* Complex */ ae_vector* x, polynomialsolverreport* rep, ae_state *_state);
    void _polynomialsolverreport_init(void* _p, ae_state *_state);
    void _polynomialsolverreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _polynomialsolverreport_clear(void* _p);
    void _polynomialsolverreport_destroy(void* _p);
}