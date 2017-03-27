module alglib.linalg;
import alglib.ap;
import alglib.internal;
import alglib.misc;

struct sparsematrix
{
    ae_vector vals;
    ae_vector idx;
    ae_vector ridx;
    ae_vector didx;
    ae_vector uidx;
    ae_int_t matrixtype;
    ae_int_t m;
    ae_int_t n;
    ae_int_t nfree;
    ae_int_t ninitialized;
    ae_int_t tablesize;
}

struct sparsebuffers
{
    ae_vector d;
    ae_vector u;
    sparsematrix s;
}

struct matinvreport
{
    double r1;
    double rinf;
}

struct fblslincgstate
{
    double e1;
    double e2;
    ae_vector x;
    ae_vector ax;
    double xax;
    ae_int_t n;
    ae_vector rk;
    ae_vector rk1;
    ae_vector xk;
    ae_vector xk1;
    ae_vector pk;
    ae_vector pk1;
    ae_vector b;
    rcommstate rstate;
    ae_vector tmp2;
}

struct normestimatorstate
{
    ae_int_t n;
    ae_int_t m;
    ae_int_t nstart;
    ae_int_t nits;
    ae_int_t seedval;
    ae_vector x0;
    ae_vector x1;
    ae_vector t;
    ae_vector xbest;
    hqrndstate r;
    ae_vector x;
    ae_vector mv;
    ae_vector mtv;
    ae_bool needmv;
    ae_bool needmtv;
    double repnorm;
    rcommstate rstate;
}

/**
extern(C++,alglib):
{
    class _sparsematrix_owner
    {
        _sparsematrix_owner();
        _sparsematrix_owner(const _sparsematrix_owner &rhs);
        _sparsematrix_owner& operator=(const _sparsematrix_owner &rhs);
        virtual ~_sparsematrix_owner();
        alglib_impl::sparsematrix* c_ptr();
        alglib_impl::sparsematrix* c_ptr() const;
        protected:
            alglib_impl::sparsematrix *p_struct;
    }

    class sparsematrix :  _sparsematrix_owner
    {
        sparsematrix();
        sparsematrix(const sparsematrix &rhs);
        sparsematrix& operator=(const sparsematrix &rhs);
        virtual ~sparsematrix();

    }

    class _sparsebuffers_owner
    {
        _sparsebuffers_owner();
        _sparsebuffers_owner(const _sparsebuffers_owner &rhs);
        _sparsebuffers_owner& operator=(const _sparsebuffers_owner &rhs);
        virtual ~_sparsebuffers_owner();
        alglib_impl::sparsebuffers* c_ptr();
        alglib_impl::sparsebuffers* c_ptr() const;
        protected:
            alglib_impl::sparsebuffers *p_struct;
    }

    class sparsebuffers : public _sparsebuffers_owner
    {
        sparsebuffers();
        sparsebuffers(const sparsebuffers &rhs);
        sparsebuffers& operator=(const sparsebuffers &rhs);
        virtual ~sparsebuffers();
    }

    class _matinvreport_owner
    {
        _matinvreport_owner();
        _matinvreport_owner(const _matinvreport_owner &rhs);
        _matinvreport_owner& operator=(const _matinvreport_owner &rhs);
        virtual ~_matinvreport_owner();
        alglib_impl::matinvreport* c_ptr();
        alglib_impl::matinvreport* c_ptr() const;
        protected:
            alglib_impl::matinvreport *p_struct;
    }

    class matinvreport :  _matinvreport_owner
    {
        matinvreport();
        matinvreport(const matinvreport &rhs);
        matinvreport& operator=(const matinvreport &rhs);
        virtual ~matinvreport();
        double &r1;
        double &rinf;
    }


    class _normestimatorstate_owner
    {
        _normestimatorstate_owner();
        _normestimatorstate_owner(const _normestimatorstate_owner &rhs);
        _normestimatorstate_owner& operator=(const _normestimatorstate_owner &rhs);
        virtual ~_normestimatorstate_owner();
        alglib_impl::normestimatorstate* c_ptr();
        alglib_impl::normestimatorstate* c_ptr() const;
        protected:
            alglib_impl::normestimatorstate *p_struct;
    }

    class normestimatorstate : _normestimatorstate_owner
    {
        normestimatorstate();
        normestimatorstate(const normestimatorstate &rhs);
        normestimatorstate& operator=(const normestimatorstate &rhs);
        virtual ~normestimatorstate();

    }
    void cmatrixtranspose(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb);
    void rmatrixtranspose(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb);
    void rmatrixenforcesymmetricity(const real_2d_array &a, const ae_int_t n, const bool isupper);
    void cmatrixcopy(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_2d_array &b, const ae_int_t ib, const ae_int_t jb);
    void rmatrixcopy(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_2d_array &b, const ae_int_t ib, const ae_int_t jb);
    void cmatrixrank1(const ae_int_t m, const ae_int_t n, complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, complex_1d_array &u, const ae_int_t iu, complex_1d_array &v, const ae_int_t iv);
    void rmatrixrank1(const ae_int_t m, const ae_int_t n, real_2d_array &a, const ae_int_t ia, const ae_int_t ja, real_1d_array &u, const ae_int_t iu, real_1d_array &v, const ae_int_t iv);
    void cmatrixmv(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const complex_1d_array &x, const ae_int_t ix, complex_1d_array &y, const ae_int_t iy);
    void rmatrixmv(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t opa, const real_1d_array &x, const ae_int_t ix, real_1d_array &y, const ae_int_t iy);
    void cmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void smp_cmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void cmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void smp_cmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const complex_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const complex_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void rmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void smp_rmatrixrighttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void rmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void smp_rmatrixlefttrsm(const ae_int_t m, const ae_int_t n, const real_2d_array &a, const ae_int_t i1, const ae_int_t j1, const bool isupper, const bool isunit, const ae_int_t optype, const real_2d_array &x, const ae_int_t i2, const ae_int_t j2);
    void cmatrixherk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
    void smp_cmatrixherk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
    void rmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
    void smp_rmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
    void cmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const alglib::complex alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const complex_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const alglib::complex beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc);
    void smp_cmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const alglib::complex alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const complex_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const alglib::complex beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc);
    void rmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc);
    void smp_rmatrixgemm(const ae_int_t m, const ae_int_t n, const ae_int_t k, const double alpha, const real_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const real_2d_array &b, const ae_int_t ib, const ae_int_t jb, const ae_int_t optypeb, const double beta, const real_2d_array &c, const ae_int_t ic, const ae_int_t jc);
    void cmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
    void smp_cmatrixsyrk(const ae_int_t n, const ae_int_t k, const double alpha, const complex_2d_array &a, const ae_int_t ia, const ae_int_t ja, const ae_int_t optypea, const double beta, const complex_2d_array &c, const ae_int_t ic, const ae_int_t jc, const bool isupper);
    void rmatrixqr(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau);
    void smp_rmatrixqr(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau);
    void rmatrixlq(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau);
    void smp_rmatrixlq(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tau);
    void cmatrixqr(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau);
    void smp_cmatrixqr(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau);
    void cmatrixlq(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau);
    void smp_cmatrixlq(complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_1d_array &tau);
    void rmatrixqrunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qcolumns, real_2d_array &q);
    void smp_rmatrixqrunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qcolumns, real_2d_array &q);
    void rmatrixqrunpackr(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &r);
    void rmatrixlqunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qrows, real_2d_array &q);
    void smp_rmatrixlqunpackq(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const real_1d_array &tau, const ae_int_t qrows, real_2d_array &q);
    void rmatrixlqunpackl(const real_2d_array &a, const ae_int_t m, const ae_int_t n, real_2d_array &l);
    void cmatrixqrunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qcolumns, complex_2d_array &q);
    void smp_cmatrixqrunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qcolumns, complex_2d_array &q);
    void cmatrixqrunpackr(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &r);
    void cmatrixlqunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qrows, complex_2d_array &q);
    void smp_cmatrixlqunpackq(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, const complex_1d_array &tau, const ae_int_t qrows, complex_2d_array &q);
    void cmatrixlqunpackl(const complex_2d_array &a, const ae_int_t m, const ae_int_t n, complex_2d_array &l);
    void rmatrixbd(real_2d_array &a, const ae_int_t m, const ae_int_t n, real_1d_array &tauq, real_1d_array &taup);
    void rmatrixbdunpackq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, const ae_int_t qcolumns, real_2d_array &q);
    void rmatrixbdmultiplybyq(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &tauq, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose);
    void rmatrixbdunpackpt(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, const ae_int_t ptrows, real_2d_array &pt);
    void rmatrixbdmultiplybyp(const real_2d_array &qp, const ae_int_t m, const ae_int_t n, const real_1d_array &taup, real_2d_array &z, const ae_int_t zrows, const ae_int_t zcolumns, const bool fromtheright, const bool dotranspose);
    void rmatrixbdunpackdiagonals(const real_2d_array &b, const ae_int_t m, const ae_int_t n, bool &isupper, real_1d_array &d, real_1d_array &e);
    void rmatrixhessenberg(real_2d_array &a, const ae_int_t n, real_1d_array &tau);
    void rmatrixhessenbergunpackq(const real_2d_array &a, const ae_int_t n, const real_1d_array &tau, real_2d_array &q);
    void rmatrixhessenbergunpackh(const real_2d_array &a, const ae_int_t n, real_2d_array &h);
    void smatrixtd(real_2d_array &a, const ae_int_t n, const bool isupper, real_1d_array &tau, real_1d_array &d, real_1d_array &e);
    void smatrixtdunpackq(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &tau, real_2d_array &q);
    void hmatrixtd(complex_2d_array &a, const ae_int_t n, const bool isupper, complex_1d_array &tau, real_1d_array &d, real_1d_array &e);
    void hmatrixtdunpackq(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &tau, complex_2d_array &q);
    bool rmatrixbdsvd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const bool isupper, const bool isfractionalaccuracyrequired, real_2d_array &u, const ae_int_t nru, real_2d_array &c, const ae_int_t ncc, real_2d_array &vt, const ae_int_t ncvt);
    bool rmatrixsvd(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const ae_int_t uneeded, const ae_int_t vtneeded, const ae_int_t additionalmemory, real_1d_array &w, real_2d_array &u, real_2d_array &vt);
    bool smp_rmatrixsvd(const real_2d_array &a, const ae_int_t m, const ae_int_t n, const ae_int_t uneeded, const ae_int_t vtneeded, const ae_int_t additionalmemory, real_1d_array &w, real_2d_array &u, real_2d_array &vt);
    bool smatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, real_2d_array &z);
    bool smatrixevdr(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, real_2d_array &z);
    bool smatrixevdi(const real_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, real_2d_array &z);
    bool hmatrixevd(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, real_1d_array &d, complex_2d_array &z);
    bool hmatrixevdr(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const double b1, const double b2, ae_int_t &m, real_1d_array &w, complex_2d_array &z);
    bool hmatrixevdi(const complex_2d_array &a, const ae_int_t n, const ae_int_t zneeded, const bool isupper, const ae_int_t i1, const ae_int_t i2, real_1d_array &w, complex_2d_array &z);
    bool smatrixtdevd(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, real_2d_array &z);
    bool smatrixtdevdr(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const double a, const double b, ae_int_t &m, real_2d_array &z);
    bool smatrixtdevdi(real_1d_array &d, const real_1d_array &e, const ae_int_t n, const ae_int_t zneeded, const ae_int_t i1, const ae_int_t i2, real_2d_array &z);
    bool rmatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t vneeded, real_1d_array &wr, real_1d_array &wi, real_2d_array &vl, real_2d_array &vr);
    void rmatrixrndorthogonal(const ae_int_t n, real_2d_array &a);
    void rmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a);
    void cmatrixrndorthogonal(const ae_int_t n, complex_2d_array &a);
    void cmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a);
    void smatrixrndcond(const ae_int_t n, const double c, real_2d_array &a);
    void spdmatrixrndcond(const ae_int_t n, const double c, real_2d_array &a);
    void hmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a);
    void hpdmatrixrndcond(const ae_int_t n, const double c, complex_2d_array &a);
    void rmatrixrndorthogonalfromtheright(real_2d_array &a, const ae_int_t m, const ae_int_t n);
    void rmatrixrndorthogonalfromtheleft(real_2d_array &a, const ae_int_t m, const ae_int_t n);
    void cmatrixrndorthogonalfromtheleft(complex_2d_array &a, const ae_int_t m, const ae_int_t n);
    void smatrixrndmultiply(real_2d_array &a, const ae_int_t n);
    void hmatrixrndmultiply(complex_2d_array &a, const ae_int_t n);
    void sparsecreate(const ae_int_t m, const ae_int_t n, const ae_int_t k, sparsematrix &s);
    void sparsecreate(const ae_int_t m, const ae_int_t n, sparsematrix &s);
    void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const ae_int_t k, const sparsematrix &s);
    void sparsecreatebuf(const ae_int_t m, const ae_int_t n, const sparsematrix &s);
    void sparsecreatecrs(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, sparsematrix &s);
    void sparsecreatecrsbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &ner, const sparsematrix &s);
    void sparsecreatesks(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, sparsematrix &s);
    void sparsecreatesksbuf(const ae_int_t m, const ae_int_t n, const integer_1d_array &d, const integer_1d_array &u, const sparsematrix &s);
    void sparsecopy(const sparsematrix &s0, sparsematrix &s1);
    void sparsecopybuf(const sparsematrix &s0, const sparsematrix &s1);
    void sparseswap(const sparsematrix &s0, const sparsematrix &s1);
    void sparseadd(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v);
    void sparseset(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v);
    double sparsegetdiagonal(const sparsematrix &s, const ae_int_t i);
    void sparsemv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y);
    void sparsemtv(const sparsematrix &s, const real_1d_array &x, real_1d_array &y);
    void sparsemv2(const sparsematrix &s, const real_1d_array &x, real_1d_array &y0, real_1d_array &y1);
    void sparsesmv(const sparsematrix &s, const bool isupper, const real_1d_array &x, real_1d_array &y);
    double sparsevsmv(const sparsematrix &s, const bool isupper, const real_1d_array &x);
    void sparsemtm(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b);
    void sparsemm2(const sparsematrix &s, const real_2d_array &a, const ae_int_t k, real_2d_array &b0, real_2d_array &b1);
    void sparsesmm(const sparsematrix &s, const bool isupper, const real_2d_array &a, const ae_int_t k, real_2d_array &b);
    void sparsetrmv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x, real_1d_array &y);
    void sparsetrsv(const sparsematrix &s, const bool isupper, const bool isunit, const ae_int_t optype, const real_1d_array &x);
    void sparseresizematrix(const sparsematrix &s);
    bool sparseenumerate(const sparsematrix &s, ae_int_t &t0, ae_int_t &t1, ae_int_t &i, ae_int_t &j, double &v);
    bool sparserewriteexisting(const sparsematrix &s, const ae_int_t i, const ae_int_t j, const double v);
    void sparsegetrow(const sparsematrix &s, const ae_int_t i, real_1d_array &irow);
    void sparsegetcompressedrow(const sparsematrix &s, const ae_int_t i, integer_1d_array &colidx, real_1d_array &vals, ae_int_t &nzcnt);
    void sparsetransposesks(const sparsematrix &s);
    void sparsecopytobuf(const sparsematrix &s0, const ae_int_t fmt, const sparsematrix &s1);
    void sparseconverttohash(const sparsematrix &s);
    void sparsecopytohash(const sparsematrix &s0, sparsematrix &s1);
    void sparsecopytohashbuf(const sparsematrix &s0, const sparsematrix &s1);
    void sparseconverttocrs(const sparsematrix &s);
    void sparsecopytocrsbuf(const sparsematrix &s0, const sparsematrix &s1);
    void sparseconverttosks(const sparsematrix &s);
    void sparsecopytosks(const sparsematrix &s0, sparsematrix &s1);
    ae_int_t sparsegetmatrixtype(const sparsematrix &s);
    bool sparseishash(const sparsematrix &s);
    bool sparseiscrs(const sparsematrix &s);
    bool sparseissks(const sparsematrix &s);
    void sparsefree(sparsematrix &s);
    ae_int_t sparsegetnrows(const sparsematrix &s);
    ae_int_t sparsegetncols(const sparsematrix &s);
    ae_int_t sparsegetlowercount(const sparsematrix &s);
    void rmatrixlu(real_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots);
    void smp_rmatrixlu(real_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots);
    void cmatrixlu(complex_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots);
    void smp_cmatrixlu(complex_2d_array &a, const ae_int_t m, const ae_int_t n, integer_1d_array &pivots);
    bool hpdmatrixcholesky(complex_2d_array &a, const ae_int_t n, const bool isupper);
    bool smp_hpdmatrixcholesky(complex_2d_array &a, const ae_int_t n, const bool isupper);
    bool spdmatrixcholesky(real_2d_array &a, const ae_int_t n, const bool isupper);
    bool smp_spdmatrixcholesky(real_2d_array &a, const ae_int_t n, const bool isupper);
    void spdmatrixcholeskyupdateadd1(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u);
    void spdmatrixcholeskyupdatefix(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix);
    void spdmatrixcholeskyupdateadd1buf(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &u, real_1d_array &bufr);
    void spdmatrixcholeskyupdatefixbuf(const real_2d_array &a, const ae_int_t n, const bool isupper, const boolean_1d_array &fix, real_1d_array &bufr);
    bool sparsecholeskyskyline(const sparsematrix &a, const ae_int_t n, const bool isupper);
    double rmatrixrcond1(const real_2d_array &a, const ae_int_t n);
    double spdmatrixrcond(const real_2d_array &a, const ae_int_t n, const bool isupper);
    double rmatrixtrrcond1(const real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit);
    double hpdmatrixrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper);
    double cmatrixrcond1(const complex_2d_array &a, const ae_int_t n);
    double cmatrixrcondinf(const complex_2d_array &a, const ae_int_t n);
    double rmatrixlurcond1(const real_2d_array &lua, const ae_int_t n);
    double spdmatrixcholeskyrcond(const real_2d_array &a, const ae_int_t n, const bool isupper);
    double hpdmatrixcholeskyrcond(const complex_2d_array &a, const ae_int_t n, const bool isupper);
    double cmatrixlurcond1(const complex_2d_array &lua, const ae_int_t n);
    double cmatrixlurcondinf(const complex_2d_array &lua, const ae_int_t n);
    double cmatrixtrrcondinf(const complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit);
    void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void smp_rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep);
    void smp_rmatrixluinverse(real_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep);
    void rmatrixinverse(real_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void smp_rmatrixinverse(real_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void rmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
    void smp_rmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
    void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void smp_cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep);
    void smp_cmatrixluinverse(complex_2d_array &a, const integer_1d_array &pivots, ae_int_t &info, matinvreport &rep);
    void cmatrixinverse(complex_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void smp_cmatrixinverse(complex_2d_array &a, const ae_int_t n, ae_int_t &info, matinvreport &rep);
    void cmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
    void smp_cmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
    void spdmatrixcholeskyinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void smp_spdmatrixcholeskyinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void spdmatrixcholeskyinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
    void smp_spdmatrixcholeskyinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
    void spdmatrixinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void smp_spdmatrixinverse(real_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void spdmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
    void smp_spdmatrixinverse(real_2d_array &a, ae_int_t &info, matinvreport &rep);
    void hpdmatrixcholeskyinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void smp_hpdmatrixcholeskyinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void hpdmatrixcholeskyinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
    void smp_hpdmatrixcholeskyinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
    void hpdmatrixinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void smp_hpdmatrixinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, ae_int_t &info, matinvreport &rep);
    void hpdmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
    void smp_hpdmatrixinverse(complex_2d_array &a, ae_int_t &info, matinvreport &rep);
    void rmatrixtrinverse(real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep);
    void smp_rmatrixtrinverse(real_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep);
    void rmatrixtrinverse(real_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep);
    void smp_rmatrixtrinverse(real_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep);
    void cmatrixtrinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep);
    void smp_cmatrixtrinverse(complex_2d_array &a, const ae_int_t n, const bool isupper, const bool isunit, ae_int_t &info, matinvreport &rep);
    void cmatrixtrinverse(complex_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep);
    void smp_cmatrixtrinverse(complex_2d_array &a, const bool isupper, ae_int_t &info, matinvreport &rep);
    void normestimatorcreate(const ae_int_t m, const ae_int_t n, const ae_int_t nstart, const ae_int_t nits, normestimatorstate &state);
    void normestimatorsetseed(const normestimatorstate &state, const ae_int_t seedval);
    void normestimatorestimatesparse(const normestimatorstate &state, const sparsematrix &a);
    void normestimatorresults(const normestimatorstate &state, double &nrm);
    double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots, const ae_int_t n);
    double rmatrixludet(const real_2d_array &a, const integer_1d_array &pivots);
    double rmatrixdet(const real_2d_array &a, const ae_int_t n);
    double rmatrixdet(const real_2d_array &a);
    alglib::complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots, const ae_int_t n);
    alglib::complex cmatrixludet(const complex_2d_array &a, const integer_1d_array &pivots);
    alglib::complex cmatrixdet(const complex_2d_array &a, const ae_int_t n);
    alglib::complex cmatrixdet(const complex_2d_array &a);
    double spdmatrixcholeskydet(const real_2d_array &a, const ae_int_t n);
    double spdmatrixcholeskydet(const real_2d_array &a);
    double spdmatrixdet(const real_2d_array &a, const ae_int_t n, const bool isupper);
    double spdmatrixdet(const real_2d_array &a);
    bool smatrixgevd(const real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t zneeded, const ae_int_t problemtype, real_1d_array &d, real_2d_array &z);
    bool smatrixgevdreduce(real_2d_array &a, const ae_int_t n, const bool isuppera, const real_2d_array &b, const bool isupperb, const ae_int_t problemtype, real_2d_array &r, bool &isupperr);
    void rmatrixinvupdatesimple(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const ae_int_t updcolumn, const double updval);
    void rmatrixinvupdaterow(real_2d_array &inva, const ae_int_t n, const ae_int_t updrow, const real_1d_array &v);
    void rmatrixinvupdatecolumn(real_2d_array &inva, const ae_int_t n, const ae_int_t updcolumn, const real_1d_array &u);
    void rmatrixinvupdateuv(real_2d_array &inva, const ae_int_t n, const real_1d_array &u, const real_1d_array &v);
    bool rmatrixschur(real_2d_array &a, const ae_int_t n, real_2d_array &s);
}
*/

extern(C)
{
    void ablassplitlength(/* Real    */ ae_matrix* a, ae_int_t n, ae_int_t* n1, ae_int_t* n2, ae_state *_state);
    void ablascomplexsplitlength(/* Complex */ ae_matrix* a, ae_int_t n, ae_int_t* n1, ae_int_t* n2, ae_state *_state);
    ae_int_t ablasblocksize(/* Real    */ ae_matrix* a, ae_state *_state);
    ae_int_t ablascomplexblocksize(/* Complex */ ae_matrix* a, ae_state *_state);
    ae_int_t ablasmicroblocksize(ae_state *_state);
    void cmatrixtranspose(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, /* Complex */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_state *_state);
    void rmatrixtranspose(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, /* Real    */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_state *_state);
    void rmatrixenforcesymmetricity(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    void cmatrixcopy(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, /* Complex */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_state *_state);
    void rmatrixcopy(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, /* Real    */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_state *_state);
    void cmatrixrank1(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, /* Complex */ ae_vector* u, ae_int_t iu, /* Complex */ ae_vector* v, ae_int_t iv, ae_state *_state);
    void rmatrixrank1(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, /* Real    */ ae_vector* u, ae_int_t iu, /* Real    */ ae_vector* v, ae_int_t iv, ae_state *_state);
    void cmatrixmv(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t opa, /* Complex */ ae_vector* x, ae_int_t ix, /* Complex */ ae_vector* y, ae_int_t iy, ae_state *_state);
    void rmatrixmv(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t opa, /* Real    */ ae_vector* x, ae_int_t ix, /* Real    */ ae_vector* y, ae_int_t iy, ae_state *_state);
    void cmatrixrighttrsm(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Complex */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void _pexec_cmatrixrighttrsm(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Complex */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void cmatrixlefttrsm(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Complex */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void _pexec_cmatrixlefttrsm(ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Complex */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void rmatrixrighttrsm(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Real    */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void _pexec_rmatrixrighttrsm(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Real    */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void rmatrixlefttrsm(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Real    */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void _pexec_rmatrixlefttrsm(ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Real    */ ae_matrix* x, ae_int_t i2, ae_int_t j2, ae_state *_state);
    void cmatrixherk(ae_int_t n, ae_int_t k, double alpha, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, /* Complex */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_bool isupper, ae_state *_state);
    void _pexec_cmatrixherk(ae_int_t n, ae_int_t k, double alpha, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, /* Complex */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_bool isupper, ae_state *_state);
    void rmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, /* Real    */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_bool isupper, ae_state *_state);
    void _pexec_rmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, /* Real    */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_bool isupper, ae_state *_state);
    void cmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, /* Complex */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, ae_complex beta, /* Complex */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_state *_state);
    void _pexec_cmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, /* Complex */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, ae_complex beta, /* Complex */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_state *_state);
    void rmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, /* Real    */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, /* Real    */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_state *_state);
    void _pexec_rmatrixgemm(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, /* Real    */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, /* Real    */ ae_matrix* b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, /* Real    */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_state *_state);
    void cmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, /* Complex */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_bool isupper, ae_state *_state);
    void _pexec_cmatrixsyrk(ae_int_t n, ae_int_t k, double alpha, /* Complex */ ae_matrix* a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, /* Complex */ ae_matrix* c, ae_int_t ic, ae_int_t jc, ae_bool isupper, ae_state *_state);
    void rmatrixqr(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_state *_state);
    void _pexec_rmatrixqr(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_state *_state);
    void rmatrixlq(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_state *_state);
    void _pexec_rmatrixlq(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_state *_state);
    void cmatrixqr(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_state *_state);
    void _pexec_cmatrixqr(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_state *_state);
    void cmatrixlq(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_state *_state);
    void _pexec_cmatrixlq(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_state *_state);
    void rmatrixqrunpackq(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_int_t qcolumns, /* Real    */ ae_matrix* q, ae_state *_state);
    void _pexec_rmatrixqrunpackq(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_int_t qcolumns, /* Real    */ ae_matrix* q, ae_state *_state);
    void rmatrixqrunpackr(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* r, ae_state *_state);
    void rmatrixlqunpackq(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_int_t qrows, /* Real    */ ae_matrix* q, ae_state *_state);
    void _pexec_rmatrixlqunpackq(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tau, ae_int_t qrows, /* Real    */ ae_matrix* q, ae_state *_state);
    void rmatrixlqunpackl(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_matrix* l, ae_state *_state);
    void cmatrixqrunpackq(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_int_t qcolumns, /* Complex */ ae_matrix* q, ae_state *_state);
    void _pexec_cmatrixqrunpackq(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_int_t qcolumns, /* Complex */ ae_matrix* q, ae_state *_state);
    void cmatrixqrunpackr(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* r, ae_state *_state);
    void cmatrixlqunpackq(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_int_t qrows, /* Complex */ ae_matrix* q, ae_state *_state);
    void _pexec_cmatrixlqunpackq(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_vector* tau, ae_int_t qrows, /* Complex */ ae_matrix* q, ae_state *_state);
    void cmatrixlqunpackl(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Complex */ ae_matrix* l, ae_state *_state);
    void rmatrixqrbasecase(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* work, /* Real    */ ae_vector* t, /* Real    */ ae_vector* tau, ae_state *_state);
    void rmatrixlqbasecase(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* work, /* Real    */ ae_vector* t, /* Real    */ ae_vector* tau, ae_state *_state);
    void rmatrixbd(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tauq, /* Real    */ ae_vector* taup, ae_state *_state);
    void rmatrixbdunpackq(/* Real    */ ae_matrix* qp, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tauq, ae_int_t qcolumns, /* Real    */ ae_matrix* q, ae_state *_state);
    void rmatrixbdmultiplybyq(/* Real    */ ae_matrix* qp, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tauq, /* Real    */ ae_matrix* z, ae_int_t zrows, ae_int_t zcolumns, ae_bool fromtheright, ae_bool dotranspose, ae_state *_state);
    void rmatrixbdunpackpt(/* Real    */ ae_matrix* qp, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* taup, ae_int_t ptrows, /* Real    */ ae_matrix* pt, ae_state *_state);
    void rmatrixbdmultiplybyp(/* Real    */ ae_matrix* qp, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* taup, /* Real    */ ae_matrix* z, ae_int_t zrows, ae_int_t zcolumns, ae_bool fromtheright, ae_bool dotranspose, ae_state *_state);
    void rmatrixbdunpackdiagonals(/* Real    */ ae_matrix* b, ae_int_t m, ae_int_t n, ae_bool* isupper, /* Real    */ ae_vector* d, /* Real    */ ae_vector* e, ae_state *_state);
    void rmatrixhessenberg(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_vector* tau, ae_state *_state);
    void rmatrixhessenbergunpackq(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_vector* tau, /* Real    */ ae_matrix* q, ae_state *_state);
    void rmatrixhessenbergunpackh(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_matrix* h, ae_state *_state);
    void smatrixtd(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* tau, /* Real    */ ae_vector* d, /* Real    */ ae_vector* e, ae_state *_state);
    void smatrixtdunpackq(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* tau, /* Real    */ ae_matrix* q, ae_state *_state);
    void hmatrixtd(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* tau, /* Real    */ ae_vector* d, /* Real    */ ae_vector* e, ae_state *_state);
    void hmatrixtdunpackq(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Complex */ ae_vector* tau, /* Complex */ ae_matrix* q, ae_state *_state);
    ae_bool rmatrixbdsvd(/* Real    */ ae_vector* d, /* Real    */ ae_vector* e, ae_int_t n, ae_bool isupper, ae_bool isfractionalaccuracyrequired, /* Real    */ ae_matrix* u, ae_int_t nru, /* Real    */ ae_matrix* c, ae_int_t ncc, /* Real    */ ae_matrix* vt, ae_int_t ncvt, ae_state *_state);
    ae_bool bidiagonalsvddecomposition(/* Real    */ ae_vector* d, /* Real    */ ae_vector* e, ae_int_t n, ae_bool isupper, ae_bool isfractionalaccuracyrequired, /* Real    */ ae_matrix* u, ae_int_t nru, /* Real    */ ae_matrix* c, ae_int_t ncc, /* Real    */ ae_matrix* vt, ae_int_t ncvt, ae_state *_state);
    ae_bool rmatrixsvd(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, ae_int_t uneeded, ae_int_t vtneeded, ae_int_t additionalmemory, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* u, /* Real    */ ae_matrix* vt, ae_state *_state);
    ae_bool _pexec_rmatrixsvd(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, ae_int_t uneeded, ae_int_t vtneeded, ae_int_t additionalmemory, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* u, /* Real    */ ae_matrix* vt, ae_state *_state);
    ae_bool smatrixevd(/* Real    */ ae_matrix* a, ae_int_t n, ae_int_t zneeded, ae_bool isupper, /* Real    */ ae_vector* d, /* Real    */ ae_matrix* z, ae_state *_state);
    ae_bool smatrixevdr(/* Real    */ ae_matrix* a, ae_int_t n, ae_int_t zneeded, ae_bool isupper, double b1, double b2, ae_int_t* m, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* z, ae_state *_state);
    ae_bool smatrixevdi(/* Real    */ ae_matrix* a, ae_int_t n, ae_int_t zneeded, ae_bool isupper, ae_int_t i1, ae_int_t i2, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* z, ae_state *_state);
    ae_bool hmatrixevd(/* Complex */ ae_matrix* a, ae_int_t n, ae_int_t zneeded, ae_bool isupper, /* Real    */ ae_vector* d, /* Complex */ ae_matrix* z, ae_state *_state);
    ae_bool hmatrixevdr(/* Complex */ ae_matrix* a, ae_int_t n, ae_int_t zneeded, ae_bool isupper, double b1, double b2, ae_int_t* m, /* Real    */ ae_vector* w, /* Complex */ ae_matrix* z, ae_state *_state);
    ae_bool hmatrixevdi(/* Complex */ ae_matrix* a, ae_int_t n, ae_int_t zneeded, ae_bool isupper, ae_int_t i1, ae_int_t i2, /* Real    */ ae_vector* w, /* Complex */ ae_matrix* z, ae_state *_state);
    ae_bool smatrixtdevd(/* Real    */ ae_vector* d,
         /* Real    */ ae_vector* e,
         ae_int_t n,
         ae_int_t zneeded,
         /* Real    */ ae_matrix* z,
         ae_state *_state);
    ae_bool smatrixtdevdr(/* Real    */ ae_vector* d, /* Real    */ ae_vector* e, ae_int_t n, ae_int_t zneeded, double a, double b, ae_int_t* m, /* Real    */ ae_matrix* z, ae_state *_state);
    ae_bool smatrixtdevdi(/* Real    */ ae_vector* d, /* Real    */ ae_vector* e, ae_int_t n, ae_int_t zneeded, ae_int_t i1, ae_int_t i2, /* Real    */ ae_matrix* z, ae_state *_state);
    ae_bool rmatrixevd(/* Real    */ ae_matrix* a, ae_int_t n, ae_int_t vneeded, /* Real    */ ae_vector* wr, /* Real    */ ae_vector* wi, /* Real    */ ae_matrix* vl, /* Real    */ ae_matrix* vr, ae_state *_state);
    void rmatrixrndorthogonal(ae_int_t n, /* Real    */ ae_matrix* a, ae_state *_state);
    void rmatrixrndcond(ae_int_t n, double c, /* Real    */ ae_matrix* a, ae_state *_state);
    void cmatrixrndorthogonal(ae_int_t n, /* Complex */ ae_matrix* a, ae_state *_state);
    void cmatrixrndcond(ae_int_t n, double c, /* Complex */ ae_matrix* a, ae_state *_state);
    void smatrixrndcond(ae_int_t n, double c, /* Real    */ ae_matrix* a, ae_state *_state);
    void spdmatrixrndcond(ae_int_t n, double c, /* Real    */ ae_matrix* a, ae_state *_state);
    void hmatrixrndcond(ae_int_t n, double c, /* Complex */ ae_matrix* a, ae_state *_state);
    void hpdmatrixrndcond(ae_int_t n, double c, /* Complex */ ae_matrix* a, ae_state *_state);
    void rmatrixrndorthogonalfromtheright(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, ae_state *_state);
    void rmatrixrndorthogonalfromtheleft(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, ae_state *_state);
    void cmatrixrndorthogonalfromtheright(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, ae_state *_state);
    void cmatrixrndorthogonalfromtheleft(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, ae_state *_state);
    void smatrixrndmultiply(/* Real    */ ae_matrix* a, ae_int_t n, ae_state *_state);
    void hmatrixrndmultiply(/* Complex */ ae_matrix* a, ae_int_t n, ae_state *_state);
    void sparsecreate(ae_int_t m, ae_int_t n, ae_int_t k, sparsematrix* s, ae_state *_state);
    void sparsecreatebuf(ae_int_t m, ae_int_t n, ae_int_t k, sparsematrix* s, ae_state *_state);
    void sparsecreatecrs(ae_int_t m, ae_int_t n, /* Integer */ ae_vector* ner, sparsematrix* s, ae_state *_state);
    void sparsecreatecrsbuf(ae_int_t m, ae_int_t n, /* Integer */ ae_vector* ner, sparsematrix* s, ae_state *_state);
    void sparsecreatesks(ae_int_t m, ae_int_t n, /* Integer */ ae_vector* d, /* Integer */ ae_vector* u, sparsematrix* s, ae_state *_state);
    void sparsecreatesksbuf(ae_int_t m, ae_int_t n, /* Integer */ ae_vector* d, /* Integer */ ae_vector* u, sparsematrix* s, ae_state *_state);
    void sparsecopy(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparsecopybuf(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparseswap(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparseadd(sparsematrix* s, ae_int_t i, ae_int_t j, double v, ae_state *_state);
    void sparseset(sparsematrix* s, ae_int_t i, ae_int_t j, double v, ae_state *_state);
    double sparseget(sparsematrix* s, ae_int_t i, ae_int_t j, ae_state *_state);
    double sparsegetdiagonal(sparsematrix* s, ae_int_t i, ae_state *_state);
    void sparsemv(sparsematrix* s, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void sparsemtv(sparsematrix* s, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void sparsemv2(sparsematrix* s, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y0, /* Real    */ ae_vector* y1, ae_state *_state);
    void sparsesmv(sparsematrix* s, ae_bool isupper, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    double sparsevsmv(sparsematrix* s, ae_bool isupper, /* Real    */ ae_vector* x, ae_state *_state);
    void sparsemm(sparsematrix* s, /* Real    */ ae_matrix* a, ae_int_t k, /* Real    */ ae_matrix* b, ae_state *_state);
    void sparsemtm(sparsematrix* s, /* Real    */ ae_matrix* a, ae_int_t k, /* Real    */ ae_matrix* b, ae_state *_state);
    void sparsemm2(sparsematrix* s, /* Real    */ ae_matrix* a, ae_int_t k, /* Real    */ ae_matrix* b0, /* Real    */ ae_matrix* b1, ae_state *_state);
    void sparsesmm(sparsematrix* s, ae_bool isupper, /* Real    */ ae_matrix* a, ae_int_t k, /* Real    */ ae_matrix* b, ae_state *_state);
    void sparsetrmv(sparsematrix* s, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void sparsetrsv(sparsematrix* s, ae_bool isupper, ae_bool isunit, ae_int_t optype, /* Real    */ ae_vector* x, ae_state *_state);
    void sparseresizematrix(sparsematrix* s, ae_state *_state);
    double sparsegetaveragelengthofchain(sparsematrix* s, ae_state *_state);
    ae_bool sparseenumerate(sparsematrix* s, ae_int_t* t0, ae_int_t* t1, ae_int_t* i, ae_int_t* j, double* v, ae_state *_state);
    ae_bool sparserewriteexisting(sparsematrix* s, ae_int_t i, ae_int_t j, double v, ae_state *_state);
    void sparsegetrow(sparsematrix* s, ae_int_t i, /* Real    */ ae_vector* irow, ae_state *_state);
    void sparsegetcompressedrow(sparsematrix* s, ae_int_t i, /* Integer */ ae_vector* colidx, /* Real    */ ae_vector* vals, ae_int_t* nzcnt, ae_state *_state);
    void sparsetransposesks(sparsematrix* s, ae_state *_state);
    void sparseconvertto(sparsematrix* s0, ae_int_t fmt, ae_state *_state);
    void sparsecopytobuf(sparsematrix* s0, ae_int_t fmt, sparsematrix* s1, ae_state *_state);
    void sparseconverttohash(sparsematrix* s, ae_state *_state);
    void sparsecopytohash(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparsecopytohashbuf(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparseconverttocrs(sparsematrix* s, ae_state *_state);
    void sparsecopytocrs(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparsecopytocrsbuf(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparseconverttosks(sparsematrix* s, ae_state *_state);
    void sparsecopytosks(sparsematrix* s0, sparsematrix* s1, ae_state *_state);
    void sparsecopytosksbuf(sparsematrix* s0,
         sparsematrix* s1,
         ae_state *_state);
    ae_int_t sparsegetmatrixtype(sparsematrix* s, ae_state *_state);
    ae_bool sparseishash(sparsematrix* s, ae_state *_state);
    ae_bool sparseiscrs(sparsematrix* s, ae_state *_state);
    ae_bool sparseissks(sparsematrix* s, ae_state *_state);
    void sparsefree(sparsematrix* s, ae_state *_state);
    ae_int_t sparsegetnrows(sparsematrix* s, ae_state *_state);
    ae_int_t sparsegetncols(sparsematrix* s, ae_state *_state);
    ae_int_t sparsegetuppercount(sparsematrix* s, ae_state *_state);
    ae_int_t sparsegetlowercount(sparsematrix* s, ae_state *_state);
    void _sparsematrix_init(void* _p, ae_state *_state);
    void _sparsematrix_init_copy(void* _dst, void* _src, ae_state *_state);
    void _sparsematrix_clear(void* _p);
    void _sparsematrix_destroy(void* _p);
    void _sparsebuffers_init(void* _p, ae_state *_state);
    void _sparsebuffers_init_copy(void* _dst, void* _src, ae_state *_state);
    void _sparsebuffers_clear(void* _p);
    void _sparsebuffers_destroy(void* _p);
    void rmatrixlu(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    void _pexec_rmatrixlu(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    void cmatrixlu(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    void _pexec_cmatrixlu(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    ae_bool hpdmatrixcholesky(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    ae_bool _pexec_hpdmatrixcholesky(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    ae_bool spdmatrixcholesky(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    ae_bool _pexec_spdmatrixcholesky(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    void spdmatrixcholeskyupdateadd1(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* u, ae_state *_state);
    void spdmatrixcholeskyupdatefix(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Boolean */ ae_vector* fix, ae_state *_state);
    void spdmatrixcholeskyupdateadd1buf(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* u, /* Real    */ ae_vector* bufr, ae_state *_state);
    void spdmatrixcholeskyupdatefixbuf(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, /* Boolean */ ae_vector* fix, /* Real    */ ae_vector* bufr, ae_state *_state);
    ae_bool sparsecholeskyskyline(sparsematrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    ae_bool sparsecholeskyx(sparsematrix* a, ae_int_t n, ae_bool isupper, /* Integer */ ae_vector* p0, /* Integer */ ae_vector* p1, ae_int_t ordering, ae_int_t algo, ae_int_t fmt, sparsebuffers* buf, sparsematrix* c, ae_state *_state);
    void rmatrixlup(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    void cmatrixlup(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    void rmatrixplu(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    void cmatrixplu(/* Complex */ ae_matrix* a, ae_int_t m, ae_int_t n, /* Integer */ ae_vector* pivots, ae_state *_state);
    ae_bool spdmatrixcholeskyrec(/* Real    */ ae_matrix* a, ae_int_t offs, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* tmp, ae_state *_state);
    double rmatrixrcond1(/* Real    */ ae_matrix* a, ae_int_t n, ae_state *_state);
    double rmatrixrcondinf(/* Real    */ ae_matrix* a, ae_int_t n, ae_state *_state);
    double spdmatrixrcond(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    double rmatrixtrrcond1(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_state *_state);
    double rmatrixtrrcondinf(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_state *_state);
    double hpdmatrixrcond(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    double cmatrixrcond1(/* Complex */ ae_matrix* a, ae_int_t n, ae_state *_state);
    double cmatrixrcondinf(/* Complex */ ae_matrix* a, ae_int_t n, ae_state *_state);
    double rmatrixlurcond1(/* Real    */ ae_matrix* lua, ae_int_t n, ae_state *_state);
    double rmatrixlurcondinf(/* Real    */ ae_matrix* lua, ae_int_t n, ae_state *_state);
    double spdmatrixcholeskyrcond(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    double hpdmatrixcholeskyrcond(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    double cmatrixlurcond1(/* Complex */ ae_matrix* lua, ae_int_t n, ae_state *_state);
    double cmatrixlurcondinf(/* Complex */ ae_matrix* lua, ae_int_t n, ae_state *_state);
    double cmatrixtrrcond1(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_state *_state);
    double cmatrixtrrcondinf(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_state *_state);
    double rcondthreshold(ae_state *_state);
    void rmatrixluinverse(/* Real    */ ae_matrix* a, /* Integer */ ae_vector* pivots, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_rmatrixluinverse(/* Real    */ ae_matrix* a, /* Integer */ ae_vector* pivots, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void rmatrixinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_rmatrixinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void cmatrixluinverse(/* Complex */ ae_matrix* a, /* Integer */ ae_vector* pivots, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_cmatrixluinverse(/* Complex */ ae_matrix* a, /* Integer */ ae_vector* pivots, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void cmatrixinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_cmatrixinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void spdmatrixcholeskyinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_spdmatrixcholeskyinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void spdmatrixinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_spdmatrixinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void hpdmatrixcholeskyinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_hpdmatrixcholeskyinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void hpdmatrixinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_hpdmatrixinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void rmatrixtrinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_rmatrixtrinverse(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void cmatrixtrinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void _pexec_cmatrixtrinverse(/* Complex */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_bool isunit, ae_int_t* info, matinvreport* rep, ae_state *_state);
    void spdmatrixcholeskyinverserec(/* Real    */ ae_matrix* a, ae_int_t offs, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* tmp, ae_state *_state);
    void _matinvreport_init(void* _p, ae_state *_state);
    void _matinvreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _matinvreport_clear(void* _p);
    void _matinvreport_destroy(void* _p);
    void fblscholeskysolve(/* Real    */ ae_matrix* cha, double sqrtscalea, ae_int_t n, ae_bool isupper, /* Real    */ ae_vector* xb, /* Real    */ ae_vector* tmp, ae_state *_state);
    void fblssolvecgx(/* Real    */ ae_matrix* a, ae_int_t m, ae_int_t n, double alpha, /* Real    */ ae_vector* b, /* Real    */ ae_vector* x, /* Real    */ ae_vector* buf, ae_state *_state);
    void fblscgcreate(/* Real    */ ae_vector* x, /* Real    */ ae_vector* b, ae_int_t n, fblslincgstate* state, ae_state *_state);
    ae_bool fblscgiteration(fblslincgstate* state, ae_state *_state);
    void fblssolvels(/* Real    */ ae_matrix* a, /* Real    */ ae_vector* b, ae_int_t m, ae_int_t n, /* Real    */ ae_vector* tmp0, /* Real    */ ae_vector* tmp1, /* Real    */ ae_vector* tmp2, ae_state *_state);
    void _fblslincgstate_init(void* _p, ae_state *_state);
    void _fblslincgstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _fblslincgstate_clear(void* _p);
    void _fblslincgstate_destroy(void* _p);
    void normestimatorcreate(ae_int_t m, ae_int_t n, ae_int_t nstart, ae_int_t nits, normestimatorstate* state, ae_state *_state);
    void normestimatorsetseed(normestimatorstate* state, ae_int_t seedval, ae_state *_state);
    ae_bool normestimatoriteration(normestimatorstate* state, ae_state *_state);
    void normestimatorestimatesparse(normestimatorstate* state, sparsematrix* a, ae_state *_state);
    void normestimatorresults(normestimatorstate* state, double* nrm, ae_state *_state);
    void normestimatorrestart(normestimatorstate* state, ae_state *_state);
    void _normestimatorstate_init(void* _p, ae_state *_state);
    void _normestimatorstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _normestimatorstate_clear(void* _p);
    void _normestimatorstate_destroy(void* _p);
    double rmatrixludet(/* Real    */ ae_matrix* a, /* Integer */ ae_vector* pivots, ae_int_t n, ae_state *_state);
    double rmatrixdet(/* Real    */ ae_matrix* a, ae_int_t n, ae_state *_state);
    ae_complex cmatrixludet(/* Complex */ ae_matrix* a, /* Integer */ ae_vector* pivots, ae_int_t n, ae_state *_state);
    ae_complex cmatrixdet(/* Complex */ ae_matrix* a, ae_int_t n, ae_state *_state);
    double spdmatrixcholeskydet(/* Real    */ ae_matrix* a, ae_int_t n, ae_state *_state);
    double spdmatrixdet(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isupper, ae_state *_state);
    ae_bool smatrixgevd(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isuppera, /* Real    */ ae_matrix* b, ae_bool isupperb, ae_int_t zneeded, ae_int_t problemtype, /* Real    */ ae_vector* d, /* Real    */ ae_matrix* z, ae_state *_state);
    ae_bool smatrixgevdreduce(/* Real    */ ae_matrix* a, ae_int_t n, ae_bool isuppera, /* Real    */ ae_matrix* b, ae_bool isupperb, ae_int_t problemtype, /* Real    */ ae_matrix* r, ae_bool* isupperr, ae_state *_state);
    void rmatrixinvupdatesimple(/* Real    */ ae_matrix* inva, ae_int_t n, ae_int_t updrow, ae_int_t updcolumn, double updval, ae_state *_state);
    void rmatrixinvupdaterow(/* Real    */ ae_matrix* inva, ae_int_t n, ae_int_t updrow, /* Real    */ ae_vector* v, ae_state *_state); void rmatrixinvupdatecolumn(/* Real    */ ae_matrix* inva, ae_int_t n, ae_int_t updcolumn, /* Real    */ ae_vector* u, ae_state *_state);
    void rmatrixinvupdateuv(/* Real    */ ae_matrix* inva, ae_int_t n, /* Real    */ ae_vector* u, /* Real    */ ae_vector* v, ae_state *_state);
    ae_bool rmatrixschur(/* Real    */ ae_matrix* a, ae_int_t n, /* Real    */ ae_matrix* s, ae_state *_state);
}
