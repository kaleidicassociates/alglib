module alglib.interpolation;
import alglib.ap;
import alglib.internal;
import alglib.misc;
import alglib.linalg;
import alglib.solvers;
import alglib.optimization;
import alglib.specialfunctions;
import alglib.integration;

struct idwinterpolant
{
    ae_int_t n;
    ae_int_t nx;
    ae_int_t d;
    double r;
    ae_int_t nw;
    kdtree tree;
    ae_int_t modeltype;
    ae_matrix q;
    ae_vector xbuf;
    ae_vector tbuf;
    ae_vector rbuf;
    ae_matrix xybuf;
    ae_int_t debugsolverfailures;
    double debugworstrcond;
    double debugbestrcond;
}

struct barycentricinterpolant
{
    ae_int_t n;
    double sy;
    ae_vector x;
    ae_vector y;
    ae_vector w;
}

struct spline1dinterpolant
{
    ae_bool periodic;
    ae_int_t n;
    ae_int_t k;
    ae_int_t continuity;
    ae_vector x;
    ae_vector c;
}

struct polynomialfitreport
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
}

struct barycentricfitreport
{
    double taskrcond;
    ae_int_t dbest;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} 

struct spline1dfitreport
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
}

struct lsfitreport
{
    double taskrcond;
    ae_int_t iterationscount;
    ae_int_t varidx;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
    double wrmserror;
    ae_matrix covpar;
    ae_vector errpar;
    ae_vector errcurve;
    ae_vector noise;
    double r2;
}

struct lsfitstate
{
    ae_int_t optalgo;
    ae_int_t m;
    ae_int_t k;
    double epsf;
    double epsx;
    ae_int_t maxits;
    double stpmax;
    ae_bool xrep;
    ae_vector s;
    ae_vector bndl;
    ae_vector bndu;
    ae_matrix taskx;
    ae_vector tasky;
    ae_int_t npoints;
    ae_vector taskw;
    ae_int_t nweights;
    ae_int_t wkind;
    ae_int_t wits;
    double diffstep;
    double teststep;
    ae_bool xupdated;
    ae_bool needf;
    ae_bool needfg;
    ae_bool needfgh;
    ae_int_t pointindex;
    ae_vector x;
    ae_vector c;
    double f;
    ae_vector g;
    ae_matrix h;
    ae_vector wcur;
    ae_vector tmp;
    ae_vector tmpf;
    ae_matrix tmpjac;
    ae_matrix tmpjacw;
    double tmpnoise;
    matinvreport invrep;
    ae_int_t repiterationscount;
    ae_int_t repterminationtype;
    ae_int_t repvaridx;
    double reprmserror;
    double repavgerror;
    double repavgrelerror;
    double repmaxerror;
    double repwrmserror;
    lsfitreport rep;
    minlmstate optstate;
    minlmreport optrep;
    ae_int_t prevnpt;
    ae_int_t prevalgo;
    rcommstate rstate;
}

struct pspline2interpolant
{
    ae_int_t n;
    ae_bool periodic;
    ae_vector p;
    spline1dinterpolant x;
    spline1dinterpolant y;
}


struct pspline3interpolant
{
    ae_int_t n;
    ae_bool periodic;
    ae_vector p;
    spline1dinterpolant x;
    spline1dinterpolant y;
    spline1dinterpolant z;
} 

struct rbfmodel
{
    ae_int_t ny;
    ae_int_t nx;
    ae_int_t nc;
    ae_int_t nl;
    kdtree tree;
    ae_matrix xc;
    ae_matrix wr;
    double rmax;
    ae_matrix v;
    ae_int_t gridtype;
    ae_bool fixrad;
    double lambdav;
    double radvalue;
    double radzvalue;
    ae_int_t nlayers;
    ae_int_t aterm;
    ae_int_t algorithmtype;
    double epsort;
    double epserr;
    ae_int_t maxits;
    double h;
    ae_int_t n;
    ae_matrix x;
    ae_matrix y;
    ae_vector calcbufxcx;
    ae_matrix calcbufx;
    ae_vector calcbuftags;
}

struct rbfreport
{
    ae_int_t arows;
    ae_int_t acols;
    ae_int_t annz;
    ae_int_t iterationscount;
    ae_int_t nmv;
    ae_int_t terminationtype;
}

struct spline2dinterpolant
{
    ae_int_t k;
    ae_int_t stype;
    ae_int_t n;
    ae_int_t m;
    ae_int_t d;
    ae_vector x;
    ae_vector y;
    ae_vector f;
}

struct spline3dinterpolant
{
    ae_int_t k;
    ae_int_t stype;
    ae_int_t n;
    ae_int_t m;
    ae_int_t l;
    ae_int_t d;
    ae_vector x;
    ae_vector y;
    ae_vector z;
    ae_vector f;
}
/**
extern(C++,alglib)
{
    class _idwinterpolant_owner
    {
        _idwinterpolant_owner();
        _idwinterpolant_owner(const _idwinterpolant_owner &rhs);
        _idwinterpolant_owner& operator=(const _idwinterpolant_owner &rhs);
        virtual ~_idwinterpolant_owner();
        alglib_impl::idwinterpolant* c_ptr();
        alglib_impl::idwinterpolant* c_ptr() const;
    protected:
        alglib_impl::idwinterpolant *p_struct;
    }

    class idwinterpolant :  _idwinterpolant_owner
    {
        idwinterpolant();
        idwinterpolant(const idwinterpolant &rhs);
        idwinterpolant& operator=(const idwinterpolant &rhs);
        virtual ~idwinterpolant();
    }

    class _barycentricinterpolant_owner
    {
        _barycentricinterpolant_owner();
        _barycentricinterpolant_owner(const _barycentricinterpolant_owner &rhs);
        _barycentricinterpolant_owner& operator=(const _barycentricinterpolant_owner &rhs);
        virtual ~_barycentricinterpolant_owner();
        alglib_impl::barycentricinterpolant* c_ptr();
        alglib_impl::barycentricinterpolant* c_ptr() const;
        protected:
            alglib_impl::barycentricinterpolant *p_struct;
    }

    class barycentricinterpolant :  _barycentricinterpolant_owner
    {
        barycentricinterpolant();
        barycentricinterpolant(const barycentricinterpolant &rhs);
        barycentricinterpolant& operator=(const barycentricinterpolant &rhs);
        virtual ~barycentricinterpolant();
    }

    class _spline1dinterpolant_owner
    {
        _spline1dinterpolant_owner();
        _spline1dinterpolant_owner(const _spline1dinterpolant_owner &rhs);
        _spline1dinterpolant_owner& operator=(const _spline1dinterpolant_owner &rhs);
        virtual ~_spline1dinterpolant_owner();
        alglib_impl::spline1dinterpolant* c_ptr();
        alglib_impl::spline1dinterpolant* c_ptr() const;
        protected:
            alglib_impl::spline1dinterpolant *p_struct;
    }

    class spline1dinterpolant :  _spline1dinterpolant_owner
    {
        spline1dinterpolant();
        spline1dinterpolant(const spline1dinterpolant &rhs);
        spline1dinterpolant& operator=(const spline1dinterpolant &rhs);
        virtual ~spline1dinterpolant();
    }

    class _polynomialfitreport_owner
    {
        _polynomialfitreport_owner();
        _polynomialfitreport_owner(const _polynomialfitreport_owner &rhs);
        _polynomialfitreport_owner& operator=(const _polynomialfitreport_owner &rhs);
        virtual ~_polynomialfitreport_owner();
        alglib_impl::polynomialfitreport* c_ptr();
        alglib_impl::polynomialfitreport* c_ptr() const;
        protected:
            alglib_impl::polynomialfitreport *p_struct;
    }

    class polynomialfitreport :  _polynomialfitreport_owner
    {
        polynomialfitreport();
        polynomialfitreport(const polynomialfitreport &rhs);
        polynomialfitreport& operator=(const polynomialfitreport &rhs);
        virtual ~polynomialfitreport();
        double &taskrcond;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
        double &maxerror;
    }


    class _barycentricfitreport_owner
    {
        _barycentricfitreport_owner();
        _barycentricfitreport_owner(const _barycentricfitreport_owner &rhs);
        _barycentricfitreport_owner& operator=(const _barycentricfitreport_owner &rhs);
        virtual ~_barycentricfitreport_owner();
        alglib_impl::barycentricfitreport* c_ptr();
        alglib_impl::barycentricfitreport* c_ptr() const;
        protected:
            alglib_impl::barycentricfitreport *p_struct;
    }

    class barycentricfitreport :  _barycentricfitreport_owner
    {
        barycentricfitreport();
        barycentricfitreport(const barycentricfitreport &rhs);
        barycentricfitreport& operator=(const barycentricfitreport &rhs);
        virtual ~barycentricfitreport();
        double &taskrcond;
        ae_int_t &dbest;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
        double &maxerror;
    }

    class _spline1dfitreport_owner
    {
        _spline1dfitreport_owner();
        _spline1dfitreport_owner(const _spline1dfitreport_owner &rhs);
        _spline1dfitreport_owner& operator=(const _spline1dfitreport_owner &rhs);
        virtual ~_spline1dfitreport_owner();
        alglib_impl::spline1dfitreport* c_ptr();
        alglib_impl::spline1dfitreport* c_ptr() const;
        protected:
            alglib_impl::spline1dfitreport *p_struct;
    }

    class spline1dfitreport : _spline1dfitreport_owner
    {
        spline1dfitreport();
        spline1dfitreport(const spline1dfitreport &rhs);
        spline1dfitreport& operator=(const spline1dfitreport &rhs);
        virtual ~spline1dfitreport();
        double &taskrcond;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
        double &maxerror;
    }

    class _lsfitreport_owner
    {
        _lsfitreport_owner();
        _lsfitreport_owner(const _lsfitreport_owner &rhs);
        _lsfitreport_owner& operator=(const _lsfitreport_owner &rhs);
        virtual ~_lsfitreport_owner();
        alglib_impl::lsfitreport* c_ptr();
        alglib_impl::lsfitreport* c_ptr() const;
        protected:
            alglib_impl::lsfitreport *p_struct;
    }

    class lsfitreport :  _lsfitreport_owner
    {
        lsfitreport();
        lsfitreport(const lsfitreport &rhs);
        lsfitreport& operator=(const lsfitreport &rhs);
        virtual ~lsfitreport();
        double &taskrcond;
        ae_int_t &iterationscount;
        ae_int_t &varidx;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
        double &maxerror;
        double &wrmserror;
        real_2d_array covpar;
        real_1d_array errpar;
        real_1d_array errcurve;
        real_1d_array noise;
        double &r2;
    }


    class _lsfitstate_owner
    {
        _lsfitstate_owner();
        _lsfitstate_owner(const _lsfitstate_owner &rhs);
        _lsfitstate_owner& operator=(const _lsfitstate_owner &rhs);
        virtual ~_lsfitstate_owner();
        alglib_impl::lsfitstate* c_ptr();
        alglib_impl::lsfitstate* c_ptr() const;
        protected:
            alglib_impl::lsfitstate *p_struct;
    }

    class lsfitstate :  _lsfitstate_owner
    {
        lsfitstate();
        lsfitstate(const lsfitstate &rhs);
        lsfitstate& operator=(const lsfitstate &rhs);
        virtual ~lsfitstate();
        ae_bool &needf;
        ae_bool &needfg;
        ae_bool &needfgh;
        ae_bool &xupdated;
        real_1d_array c;
        double &f;
        real_1d_array g;
        real_2d_array h;
        real_1d_array x;

    }

    class _pspline2interpolant_owner
    {
        _pspline2interpolant_owner();
        _pspline2interpolant_owner(const _pspline2interpolant_owner &rhs);
        _pspline2interpolant_owner& operator=(const _pspline2interpolant_owner &rhs);
        virtual ~_pspline2interpolant_owner();
        alglib_impl::pspline2interpolant* c_ptr();
        alglib_impl::pspline2interpolant* c_ptr() const;
        protected:
            alglib_impl::pspline2interpolant *p_struct;
    }

    class pspline2interpolant :  _pspline2interpolant_owner
    {
        pspline2interpolant();
        pspline2interpolant(const pspline2interpolant &rhs);
        pspline2interpolant& operator=(const pspline2interpolant &rhs);
        virtual ~pspline2interpolant();
    }

    class _pspline3interpolant_owner
    {
        _pspline3interpolant_owner();
        _pspline3interpolant_owner(const _pspline3interpolant_owner &rhs);
        _pspline3interpolant_owner& operator=(const _pspline3interpolant_owner &rhs);
        virtual ~_pspline3interpolant_owner();
        alglib_impl::pspline3interpolant* c_ptr();
        alglib_impl::pspline3interpolant* c_ptr() const;
        protected:
            alglib_impl::pspline3interpolant *p_struct;
    }

    class pspline3interpolant :  _pspline3interpolant_owner
    {
        pspline3interpolant();
        pspline3interpolant(const pspline3interpolant &rhs);
        pspline3interpolant& operator=(const pspline3interpolant &rhs);
        virtual ~pspline3interpolant();
    }

    class _rbfmodel_owner
    {
        _rbfmodel_owner();
        _rbfmodel_owner(const _rbfmodel_owner &rhs);
        _rbfmodel_owner& operator=(const _rbfmodel_owner &rhs);
        virtual ~_rbfmodel_owner();
        alglib_impl::rbfmodel* c_ptr();
        alglib_impl::rbfmodel* c_ptr() const;
        protected:
            alglib_impl::rbfmodel *p_struct;
    }

    class rbfmodel :  _rbfmodel_owner
    {
        rbfmodel();
        rbfmodel(const rbfmodel &rhs);
        rbfmodel& operator=(const rbfmodel &rhs);
        virtual ~rbfmodel();
    }


    class _rbfreport_owner
    {
        _rbfreport_owner();
        _rbfreport_owner(const _rbfreport_owner &rhs);
        _rbfreport_owner& operator=(const _rbfreport_owner &rhs);
        virtual ~_rbfreport_owner();
        alglib_impl::rbfreport* c_ptr();
        alglib_impl::rbfreport* c_ptr() const;
        protected:
            alglib_impl::rbfreport *p_struct;
    }

    class rbfreport :  _rbfreport_owner
    {
        rbfreport();
        rbfreport(const rbfreport &rhs);
        rbfreport& operator=(const rbfreport &rhs);
        virtual ~rbfreport();
        ae_int_t &arows;
        ae_int_t &acols;
        ae_int_t &annz;
        ae_int_t &iterationscount;
        ae_int_t &nmv;
        ae_int_t &terminationtype;
    }

    class _spline2dinterpolant_owner
    {
        _spline2dinterpolant_owner();
        _spline2dinterpolant_owner(const _spline2dinterpolant_owner &rhs);
        _spline2dinterpolant_owner& operator=(const _spline2dinterpolant_owner &rhs);
        virtual ~_spline2dinterpolant_owner();
        alglib_impl::spline2dinterpolant* c_ptr();
        alglib_impl::spline2dinterpolant* c_ptr() const;
        protected:
            alglib_impl::spline2dinterpolant *p_struct;
    }

    class spline2dinterpolant :  _spline2dinterpolant_owner
    {
        spline2dinterpolant();
        spline2dinterpolant(const spline2dinterpolant &rhs);
        spline2dinterpolant& operator=(const spline2dinterpolant &rhs);
        virtual ~spline2dinterpolant();
    }

    class _spline3dinterpolant_owner
    {
        _spline3dinterpolant_owner();
        _spline3dinterpolant_owner(const _spline3dinterpolant_owner &rhs);
        _spline3dinterpolant_owner& operator=(const _spline3dinterpolant_owner &rhs);
        virtual ~_spline3dinterpolant_owner();
        alglib_impl::spline3dinterpolant* c_ptr();
        alglib_impl::spline3dinterpolant* c_ptr() const;
        protected:
            alglib_impl::spline3dinterpolant *p_struct;
    }

    class spline3dinterpolant :  _spline3dinterpolant_owner
    {
        spline3dinterpolant();
        spline3dinterpolant(const spline3dinterpolant &rhs);
        spline3dinterpolant& operator=(const spline3dinterpolant &rhs);
        virtual ~spline3dinterpolant();
    }

    double idwcalc(const idwinterpolant &z, const real_1d_array &x);
    void idwbuildmodifiedshepard(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t d, const ae_int_t nq, const ae_int_t nw, idwinterpolant &z);
    void idwbuildmodifiedshepardr(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const double r, idwinterpolant &z);
    void idwbuildnoisy(const real_2d_array &xy, const ae_int_t n, const ae_int_t nx, const ae_int_t d, const ae_int_t nq, const ae_int_t nw, idwinterpolant &z);
    double barycentriccalc(const barycentricinterpolant &b, const double t);
    void barycentricdiff1(const barycentricinterpolant &b, const double t, double &f, double &df);
    void barycentricdiff2(const barycentricinterpolant &b, const double t, double &f, double &df, double &d2f);
    void barycentriclintransx(const barycentricinterpolant &b, const double ca, const double cb);
    void barycentriclintransy(const barycentricinterpolant &b, const double ca, const double cb);
    void barycentricunpack(const barycentricinterpolant &b, ae_int_t &n, real_1d_array &x, real_1d_array &y, real_1d_array &w);
    void barycentricbuildxyw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, barycentricinterpolant &b);
    void barycentricbuildfloaterhormann(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t d, barycentricinterpolant &b);
    void polynomialbar2cheb(const barycentricinterpolant &p, const double a, const double b, real_1d_array &t);
    void polynomialcheb2bar(const real_1d_array &t, const ae_int_t n, const double a, const double b, barycentricinterpolant &p);
    void polynomialcheb2bar(const real_1d_array &t, const double a, const double b, barycentricinterpolant &p);
    void polynomialbar2pow(const barycentricinterpolant &p, const double c, const double s, real_1d_array &a);
    void polynomialbar2pow(const barycentricinterpolant &p, real_1d_array &a);
    void polynomialpow2bar(const real_1d_array &a, const ae_int_t n, const double c, const double s, barycentricinterpolant &p);
    void polynomialpow2bar(const real_1d_array &a, barycentricinterpolant &p);
    void polynomialbuild(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p);
    void polynomialbuild(const real_1d_array &x, const real_1d_array &y, barycentricinterpolant &p);
    void polynomialbuildeqdist(const double a, const double b, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p);
    void polynomialbuildeqdist(const double a, const double b, const real_1d_array &y, barycentricinterpolant &p);
    void polynomialbuildcheb1(const double a, const double b, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p);
    void polynomialbuildcheb1(const double a, const double b, const real_1d_array &y, barycentricinterpolant &p);
    void polynomialbuildcheb2(const double a, const double b, const real_1d_array &y, const ae_int_t n, barycentricinterpolant &p);
    void polynomialbuildcheb2(const double a, const double b, const real_1d_array &y, barycentricinterpolant &p);
    double polynomialcalceqdist(const double a, const double b, const real_1d_array &f, const ae_int_t n, const double t);
    double polynomialcalceqdist(const double a, const double b, const real_1d_array &f, const double t);
    double polynomialcalccheb1(const double a, const double b, const real_1d_array &f, const ae_int_t n, const double t);
    double polynomialcalccheb1(const double a, const double b, const real_1d_array &f, const double t);
    double polynomialcalccheb2(const double a, const double b, const real_1d_array &f, const ae_int_t n, const double t);
    double polynomialcalccheb2(const double a, const double b, const real_1d_array &f, const double t);
    void spline1dbuildlinear(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, spline1dinterpolant &c);
    void spline1dbuildlinear(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c);
    void spline1dbuildcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, spline1dinterpolant &c);
    void spline1dbuildcubic(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c);
    void spline1dgriddiffcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, real_1d_array &d);
    void spline1dgriddiffcubic(const real_1d_array &x, const real_1d_array &y, real_1d_array &d);
    void spline1dgriddiff2cubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, real_1d_array &d1, real_1d_array &d2);
    void spline1dgriddiff2cubic(const real_1d_array &x, const real_1d_array &y, real_1d_array &d1, real_1d_array &d2);
    void spline1dconvcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, const real_1d_array &x2, const ae_int_t n2, real_1d_array &y2);
    void spline1dconvcubic(const real_1d_array &x, const real_1d_array &y, const real_1d_array &x2, real_1d_array &y2);
    void spline1dconvdiffcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, const real_1d_array &x2, const ae_int_t n2, real_1d_array &y2, real_1d_array &d2);
    void spline1dconvdiffcubic(const real_1d_array &x, const real_1d_array &y, const real_1d_array &x2, real_1d_array &y2, real_1d_array &d2);
    void spline1dconvdiff2cubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundltype, const double boundl, const ae_int_t boundrtype, const double boundr, const real_1d_array &x2, const ae_int_t n2, real_1d_array &y2, real_1d_array &d2, real_1d_array &dd2);
    void spline1dconvdiff2cubic(const real_1d_array &x, const real_1d_array &y, const real_1d_array &x2, real_1d_array &y2, real_1d_array &d2, real_1d_array &dd2);
    void spline1dbuildcatmullrom(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t boundtype, const double tension, spline1dinterpolant &c);
    void spline1dbuildcatmullrom(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c);
    void spline1dbuildhermite(const real_1d_array &x, const real_1d_array &y, const real_1d_array &d, const ae_int_t n, spline1dinterpolant &c);
    void spline1dbuildhermite(const real_1d_array &x, const real_1d_array &y, const real_1d_array &d, spline1dinterpolant &c);
    void spline1dbuildakima(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, spline1dinterpolant &c);
    void spline1dbuildakima(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c);
    double spline1dcalc(const spline1dinterpolant &c, const double x);
    void spline1ddiff(const spline1dinterpolant &c, const double x, double &s, double &ds, double &d2s);
    void spline1dunpack(const spline1dinterpolant &c, ae_int_t &n, real_2d_array &tbl);
    void spline1dlintransx(const spline1dinterpolant &c, const double a, const double b);
    void spline1dlintransy(const spline1dinterpolant &c, const double a, const double b);
    double spline1dintegrate(const spline1dinterpolant &c, const double x);
    void spline1dbuildmonotone(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, spline1dinterpolant &c);
    void spline1dbuildmonotone(const real_1d_array &x, const real_1d_array &y, spline1dinterpolant &c);
    void lstfitpiecewiselinearrdpfixed(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, real_1d_array &x2, real_1d_array &y2, ae_int_t &nsections);
    void lstfitpiecewiselinearrdp(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double eps, real_1d_array &x2, real_1d_array &y2, ae_int_t &nsections);
    void polynomialfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    void smp_polynomialfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    void polynomialfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    void smp_polynomialfit(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    void polynomialfitwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    void smp_polynomialfitwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    void polynomialfitwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    void smp_polynomialfitwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, barycentricinterpolant &p, polynomialfitreport &rep);
    double logisticcalc4(const double x, const double a, const double b, const double c, const double d);
    double logisticcalc5(const double x, const double a, const double b, const double c, const double d, const double g);
    void logisticfit4(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, double &a, double &b, double &c, double &d, lsfitreport &rep);
    void logisticfit4ec(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double cnstrleft, const double cnstrright, double &a, double &b, double &c, double &d, lsfitreport &rep);
    void logisticfit5(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, double &a, double &b, double &c, double &d, double &g, lsfitreport &rep);
    void logisticfit5ec(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double cnstrleft, const double cnstrright, double &a, double &b, double &c, double &d, double &g, lsfitreport &rep);
    void logisticfit45x(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const double cnstrleft, const double cnstrright, const bool is4pl, const double lambdav, const double epsx, const ae_int_t rscnt, double &a, double &b, double &c, double &d, double &g, lsfitreport &rep);
    void barycentricfitfloaterhormannwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, barycentricinterpolant &b, barycentricfitreport &rep);
    void smp_barycentricfitfloaterhormannwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, barycentricinterpolant &b, barycentricfitreport &rep);
    void barycentricfitfloaterhormann(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, barycentricinterpolant &b, barycentricfitreport &rep);
    void smp_barycentricfitfloaterhormann(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, barycentricinterpolant &b, barycentricfitreport &rep);
    void spline1dfitpenalized(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitpenalized(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfitpenalized(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitpenalized(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfitpenalizedw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitpenalizedw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfitpenalizedw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitpenalizedw(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t m, const double rho, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfitcubicwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitcubicwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfitcubicwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitcubicwc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfithermitewc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfithermitewc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const ae_int_t n, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t k, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfithermitewc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfithermitewc(const real_1d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &xc, const real_1d_array &yc, const integer_1d_array &dc, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfitcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfitcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfitcubic(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfithermite(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfithermite(const real_1d_array &x, const real_1d_array &y, const ae_int_t n, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void spline1dfithermite(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void smp_spline1dfithermite(const real_1d_array &x, const real_1d_array &y, const ae_int_t m, ae_int_t &info, spline1dinterpolant &s, spline1dfitreport &rep);
    void lsfitlinearw(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const ae_int_t n, const ae_int_t m, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinearw(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const ae_int_t n, const ae_int_t m, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitlinearw(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinearw(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitlinearwc(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const real_2d_array &cmatrix, const ae_int_t n, const ae_int_t m, const ae_int_t k, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinearwc(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const real_2d_array &cmatrix, const ae_int_t n, const ae_int_t m, const ae_int_t k, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitlinearwc(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const real_2d_array &cmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinearwc(const real_1d_array &y, const real_1d_array &w, const real_2d_array &fmatrix, const real_2d_array &cmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitlinear(const real_1d_array &y, const real_2d_array &fmatrix, const ae_int_t n, const ae_int_t m, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinear(const real_1d_array &y, const real_2d_array &fmatrix, const ae_int_t n, const ae_int_t m, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitlinear(const real_1d_array &y, const real_2d_array &fmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinear(const real_1d_array &y, const real_2d_array &fmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitlinearc(const real_1d_array &y, const real_2d_array &fmatrix, const real_2d_array &cmatrix, const ae_int_t n, const ae_int_t m, const ae_int_t k, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinearc(const real_1d_array &y, const real_2d_array &fmatrix, const real_2d_array &cmatrix, const ae_int_t n, const ae_int_t m, const ae_int_t k, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitlinearc(const real_1d_array &y, const real_2d_array &fmatrix, const real_2d_array &cmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void smp_lsfitlinearc(const real_1d_array &y, const real_2d_array &fmatrix, const real_2d_array &cmatrix, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitcreatewf(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const double diffstep, lsfitstate &state);
    void lsfitcreatewf(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const double diffstep, lsfitstate &state);
    void lsfitcreatef(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const double diffstep, lsfitstate &state);
    void lsfitcreatef(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const double diffstep, lsfitstate &state);
    void lsfitcreatewfg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const bool cheapfg, lsfitstate &state);
    void lsfitcreatewfg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const bool cheapfg, lsfitstate &state);
    void lsfitcreatefg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, const bool cheapfg, lsfitstate &state);
    void lsfitcreatefg(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const bool cheapfg, lsfitstate &state);
    void lsfitcreatewfgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, lsfitstate &state);
    void lsfitcreatewfgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &w, const real_1d_array &c, lsfitstate &state);
    void lsfitcreatefgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, const ae_int_t n, const ae_int_t m, const ae_int_t k, lsfitstate &state);
    void lsfitcreatefgh(const real_2d_array &x, const real_1d_array &y, const real_1d_array &c, lsfitstate &state);
    void lsfitsetcond(const lsfitstate &state, const double epsf, const double epsx, const ae_int_t maxits);
    void lsfitsetstpmax(const lsfitstate &state, const double stpmax);
    void lsfitsetxrep(const lsfitstate &state, const bool needxrep);
    void lsfitsetscale(const lsfitstate &state, const real_1d_array &s);
    void lsfitsetbc(const lsfitstate &state, const real_1d_array &bndl, const real_1d_array &bndu);
    bool lsfititeration(const lsfitstate &state);
    void lsfitfit(lsfitstate &state,
        void* function(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr) func,
        void* function(const real_1d_array &c, double func, void *ptr) rep= null,
        void *ptr = null);
    void lsfitfit(lsfitstate &state,
        void* function(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr) func,
        void* function(const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) grad,
        void* function(const real_1d_array &c, double func, void *ptr) rep= null,
        void *ptr = null);
    void lsfitfit(lsfitstate &state,
        void* function (const real_1d_array &c, const real_1d_array &x, double &func, void *ptr) func,
        void* function(const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) grad,
        void* function(const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr) hess,
        void* function(const real_1d_array &c, double func, void *ptr) rep= null,
        void *ptr = null);
    void lsfitresults(const lsfitstate &state, ae_int_t &info, real_1d_array &c, lsfitreport &rep);
    void lsfitsetgradientcheck(const lsfitstate &state, const double teststep);
    void pspline2build(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline2interpolant &p);
    void pspline3build(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline3interpolant &p);
    void pspline2buildperiodic(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline2interpolant &p);
    void pspline3buildperiodic(const real_2d_array &xy, const ae_int_t n, const ae_int_t st, const ae_int_t pt, pspline3interpolant &p);
    void pspline2parametervalues(const pspline2interpolant &p, ae_int_t &n, real_1d_array &t);
    void pspline3parametervalues(const pspline3interpolant &p, ae_int_t &n, real_1d_array &t);
    void pspline2calc(const pspline2interpolant &p, const double t, double &x, double &y);
    void pspline3calc(const pspline3interpolant &p, const double t, double &x, double &y, double &z);
    void pspline2tangent(const pspline2interpolant &p, const double t, double &x, double &y);
    void pspline3tangent(const pspline3interpolant &p, const double t, double &x, double &y, double &z);
    void pspline2diff(const pspline2interpolant &p, const double t, double &x, double &dx, double &y, double &dy);
    void pspline3diff(const pspline3interpolant &p, const double t, double &x, double &dx, double &y, double &dy, double &z, double &dz);
    void pspline2diff2(const pspline2interpolant &p, const double t, double &x, double &dx, double &d2x, double &y, double &dy, double &d2y);
    void pspline3diff2(const pspline3interpolant &p, const double t, double &x, double &dx, double &d2x, double &y, double &dy, double &d2y, double &z, double &dz, double &d2z);
    double pspline2arclength(const pspline2interpolant &p, const double a, const double b);
    double pspline3arclength(const pspline3interpolant &p, const double a, const double b);
    void parametricrdpfixed(const real_2d_array &x, const ae_int_t n, const ae_int_t d, const ae_int_t stopm, const double stopeps, real_2d_array &x2, integer_1d_array &idx2, ae_int_t &nsections);
    void rbfserialize(rbfmodel &obj, std::string &s_out);
    void rbfunserialize(std::string &s_in, rbfmodel &obj);
    void rbfcreate(const ae_int_t nx, const ae_int_t ny, rbfmodel &s);
    void rbfsetpoints(const rbfmodel &s, const real_2d_array &xy, const ae_int_t n);
    void rbfsetpoints(const rbfmodel &s, const real_2d_array &xy);
    void rbfsetalgoqnn(const rbfmodel &s, const double q, const double z);
    void rbfsetalgoqnn(const rbfmodel &s);
    void rbfsetalgomultilayer(const rbfmodel &s, const double rbase, const ae_int_t nlayers, const double lambdav);
    void rbfsetalgomultilayer(const rbfmodel &s, const double rbase, const ae_int_t nlayers);
    void rbfsetlinterm(const rbfmodel &s);
    void rbfsetconstterm(const rbfmodel &s);
    void rbfsetzeroterm(const rbfmodel &s);
    void rbfbuildmodel(const rbfmodel &s, rbfreport &rep);
    double rbfcalc2(const rbfmodel &s, const double x0, const double x1);
    double rbfcalc3(const rbfmodel &s, const double x0, const double x1, const double x2);
    void rbfcalc(const rbfmodel &s, const real_1d_array &x, real_1d_array &y);
    void rbfcalcbuf(const rbfmodel &s, const real_1d_array &x, real_1d_array &y);
    void rbfgridcalc2(const rbfmodel &s, const real_1d_array &x0, const ae_int_t n0, const real_1d_array &x1, const ae_int_t n1, real_2d_array &y);
    void rbfunpack(const rbfmodel &s, ae_int_t &nx, ae_int_t &ny, real_2d_array &xwr, ae_int_t &nc, real_2d_array &v);
    double spline2dcalc(const spline2dinterpolant &c, const double x, const double y);
    void spline2ddiff(const spline2dinterpolant &c, const double x, const double y, double &f, double &fx, double &fy, double &fxy);
    void spline2dlintransxy(const spline2dinterpolant &c, const double ax, const double bx, const double ay, const double by);
    void spline2dlintransf(const spline2dinterpolant &c, const double a, const double b);
    void spline2dcopy(const spline2dinterpolant &c, spline2dinterpolant &cc);
    void spline2dresamplebicubic(const real_2d_array &a, const ae_int_t oldheight, const ae_int_t oldwidth, real_2d_array &b, const ae_int_t newheight, const ae_int_t newwidth);
    void spline2dresamplebilinear(const real_2d_array &a, const ae_int_t oldheight, const ae_int_t oldwidth, real_2d_array &b, const ae_int_t newheight, const ae_int_t newwidth);
    void spline2dbuildbilinearv(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, const real_1d_array &f, const ae_int_t d, spline2dinterpolant &c);
    void spline2dbuildbicubicv(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, const real_1d_array &f, const ae_int_t d, spline2dinterpolant &c);
    void spline2dcalcvbuf(const spline2dinterpolant &c, const double x, const double y, real_1d_array &f);
    void spline2dcalcv(const spline2dinterpolant &c, const double x, const double y, real_1d_array &f);
    void spline2dunpackv(const spline2dinterpolant &c, ae_int_t &m, ae_int_t &n, ae_int_t &d, real_2d_array &tbl);
    void spline2dbuildbilinear(const real_1d_array &x, const real_1d_array &y, const real_2d_array &f, const ae_int_t m, const ae_int_t n, spline2dinterpolant &c);
    void spline2dbuildbicubic(const real_1d_array &x, const real_1d_array &y, const real_2d_array &f, const ae_int_t m, const ae_int_t n, spline2dinterpolant &c);
    void spline2dunpack(const spline2dinterpolant &c, ae_int_t &m, ae_int_t &n, real_2d_array &tbl);
    double spline3dcalc(const spline3dinterpolant &c, const double x, const double y, const double z);
    void spline3dlintransxyz(const spline3dinterpolant &c, const double ax, const double bx, const double ay, const double by, const double az, const double bz);
    void spline3dlintransf(const spline3dinterpolant &c, const double a, const double b);
    void spline3dresampletrilinear(const real_1d_array &a, const ae_int_t oldzcount, const ae_int_t oldycount, const ae_int_t oldxcount, const ae_int_t newzcount, const ae_int_t newycount, const ae_int_t newxcount, real_1d_array &b);
    void spline3dbuildtrilinearv(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, const real_1d_array &z, const ae_int_t l, const real_1d_array &f, const ae_int_t d, spline3dinterpolant &c);
    void spline3dcalcvbuf(const spline3dinterpolant &c, const double x, const double y, const double z, real_1d_array &f);
    void spline3dcalcv(const spline3dinterpolant &c, const double x, const double y, const double z, real_1d_array &f);
    void spline3dunpackv(const spline3dinterpolant &c, ae_int_t &n, ae_int_t &m, ae_int_t &l, ae_int_t &d, ae_int_t &stype, real_2d_array &tbl);
}
*/

extern(C)
{
    double idwcalc(idwinterpolant* z, /* Real    */ ae_vector* x, ae_state *_state);
    void idwbuildmodifiedshepard(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t nx, ae_int_t d, ae_int_t nq, ae_int_t nw, idwinterpolant* z, ae_state *_state);
    void idwbuildmodifiedshepardr(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t nx, double r, idwinterpolant* z, ae_state *_state);
    void idwbuildnoisy(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t nx, ae_int_t d, ae_int_t nq, ae_int_t nw, idwinterpolant* z, ae_state *_state);
    void _idwinterpolant_init(void* _p, ae_state *_state);
    void _idwinterpolant_init_copy(void* _dst, void* _src, ae_state *_state);
    void _idwinterpolant_clear(void* _p);
    void _idwinterpolant_destroy(void* _p);
    double barycentriccalc(barycentricinterpolant* b, double t, ae_state *_state);
    void barycentricdiff1(barycentricinterpolant* b, double t, double* f, double* df, ae_state *_state);
    void barycentricdiff2(barycentricinterpolant* b, double t, double* f, double* df, double* d2f, ae_state *_state);
    void barycentriclintransx(barycentricinterpolant* b, double ca, double cb, ae_state *_state);
    void barycentriclintransy(barycentricinterpolant* b, double ca, double cb, ae_state *_state);
    void barycentricunpack(barycentricinterpolant* b, ae_int_t* n, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_state *_state);
    void barycentricbuildxyw(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, barycentricinterpolant* b, ae_state *_state);
    void barycentricbuildfloaterhormann(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t d, barycentricinterpolant* b, ae_state *_state);
    void barycentriccopy(barycentricinterpolant* b, barycentricinterpolant* b2, ae_state *_state);
    void _barycentricinterpolant_init(void* _p, ae_state *_state);
    void _barycentricinterpolant_init_copy(void* _dst, void* _src, ae_state *_state);
    void _barycentricinterpolant_clear(void* _p);
    void _barycentricinterpolant_destroy(void* _p);
    void polynomialbar2cheb(barycentricinterpolant* p, double a, double b, /* Real    */ ae_vector* t, ae_state *_state);
    void polynomialcheb2bar(/* Real    */ ae_vector* t, ae_int_t n, double a, double b, barycentricinterpolant* p, ae_state *_state);
    void polynomialbar2pow(barycentricinterpolant* p, double c, double s, /* Real    */ ae_vector* a, ae_state *_state);
    void polynomialpow2bar(/* Real    */ ae_vector* a, ae_int_t n, double c, double s, barycentricinterpolant* p, ae_state *_state);
    void polynomialbuild(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, barycentricinterpolant* p, ae_state *_state);
    void polynomialbuildeqdist(double a, double b, /* Real    */ ae_vector* y, ae_int_t n, barycentricinterpolant* p, ae_state *_state);
    void polynomialbuildcheb1(double a, double b, /* Real    */ ae_vector* y, ae_int_t n, barycentricinterpolant* p, ae_state *_state);
    void polynomialbuildcheb2(double a, double b, /* Real    */ ae_vector* y, ae_int_t n, barycentricinterpolant* p, ae_state *_state);
    double polynomialcalceqdist(double a, double b, /* Real    */ ae_vector* f, ae_int_t n, double t, ae_state *_state);
    double polynomialcalccheb1(double a, double b, /* Real    */ ae_vector* f, ae_int_t n, double t, ae_state *_state);
    double polynomialcalccheb2(double a, double b, /* Real    */ ae_vector* f, ae_int_t n, double t, ae_state *_state);
    void spline1dbuildlinear(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, spline1dinterpolant* c, ae_state *_state);
    void spline1dbuildcubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, spline1dinterpolant* c, ae_state *_state);
    void spline1dgriddiffcubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, /* Real    */ ae_vector* d, ae_state *_state);
    void spline1dgriddiff2cubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, /* Real    */ ae_vector* d1, /* Real    */ ae_vector* d2, ae_state *_state);
    void spline1dconvcubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, /* Real    */ ae_vector* x2, ae_int_t n2, /* Real    */ ae_vector* y2, ae_state *_state);
    void spline1dconvdiffcubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, /* Real    */ ae_vector* x2, ae_int_t n2, /* Real    */ ae_vector* y2, /* Real    */ ae_vector* d2, ae_state *_state);
    void spline1dconvdiff2cubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t boundltype, double boundl, ae_int_t boundrtype, double boundr, /* Real    */ ae_vector* x2, ae_int_t n2, /* Real    */ ae_vector* y2, /* Real    */ ae_vector* d2, /* Real    */ ae_vector* dd2, ae_state *_state);
    void spline1dbuildcatmullrom(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t boundtype, double tension, spline1dinterpolant* c, ae_state *_state);
    void spline1dbuildhermite(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* d, ae_int_t n, spline1dinterpolant* c, ae_state *_state);
    void spline1dbuildakima(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, spline1dinterpolant* c, ae_state *_state);
    double spline1dcalc(spline1dinterpolant* c, double x, ae_state *_state);
    void spline1ddiff(spline1dinterpolant* c, double x, double* s, double* ds, double* d2s, ae_state *_state);
    void spline1dcopy(spline1dinterpolant* c, spline1dinterpolant* cc, ae_state *_state);
    void spline1dunpack(spline1dinterpolant* c, ae_int_t* n, /* Real    */ ae_matrix* tbl, ae_state *_state);
    void spline1dlintransx(spline1dinterpolant* c, double a, double b, ae_state *_state);
    void spline1dlintransy(spline1dinterpolant* c, double a, double b, ae_state *_state);
    double spline1dintegrate(spline1dinterpolant* c, double x, ae_state *_state);
    void spline1dconvdiffinternal(/* Real    */ ae_vector* xold, /* Real    */ ae_vector* yold, /* Real    */ ae_vector* dold, ae_int_t n, /* Real    */ ae_vector* x2, ae_int_t n2, /* Real    */ ae_vector* y, ae_bool needy, /* Real    */ ae_vector* d1, ae_bool needd1, /* Real    */ ae_vector* d2, ae_bool needd2, ae_state *_state);
    void spline1drootsandextrema(spline1dinterpolant* c, /* Real    */ ae_vector* r, ae_int_t* nr, ae_bool* dr, /* Real    */ ae_vector* e, /* Integer */ ae_vector* et, ae_int_t* ne, ae_bool* de, ae_state *_state);
    void heapsortdpoints(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* d, ae_int_t n, ae_state *_state);
    void solvepolinom2(double p0, double m0, double p1, double m1, double* x0, double* x1, ae_int_t* nr, ae_state *_state);
    void solvecubicpolinom(double pa, double ma, double pb, double mb, double a, double b, double* x0, double* x1, double* x2, double* ex0, double* ex1, ae_int_t* nr, ae_int_t* ne, /* Real    */ ae_vector* tempdata, ae_state *_state);
    ae_int_t bisectmethod(double pa, double ma, double pb, double mb, double a, double b, double* x, ae_state *_state);
    void spline1dbuildmonotone(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, spline1dinterpolant* c, ae_state *_state);
    void _spline1dinterpolant_init(void* _p, ae_state *_state);
    void _spline1dinterpolant_init_copy(void* _dst, void* _src, ae_state *_state);
    void _spline1dinterpolant_clear(void* _p);
    void _spline1dinterpolant_destroy(void* _p);
    void lstfitpiecewiselinearrdpfixed(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x2, /* Real    */ ae_vector* y2, ae_int_t* nsections, ae_state *_state);
    void lstfitpiecewiselinearrdp(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, double eps, /* Real    */ ae_vector* x2, /* Real    */ ae_vector* y2, ae_int_t* nsections, ae_state *_state);
    void polynomialfit(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, barycentricinterpolant* p, polynomialfitreport* rep, ae_state *_state);
    void _pexec_polynomialfit(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, barycentricinterpolant* p, polynomialfitreport* rep, ae_state *_state);
    void polynomialfitwc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, barycentricinterpolant* p, polynomialfitreport* rep, ae_state *_state);
    void _pexec_polynomialfitwc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, barycentricinterpolant* p, polynomialfitreport* rep, ae_state *_state);
    double logisticcalc4(double x, double a, double b, double c, double d, ae_state *_state);
    double logisticcalc5(double x, double a, double b, double c, double d, double g, ae_state *_state);
    void logisticfit4(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, double* a, double* b, double* c, double* d, lsfitreport* rep, ae_state *_state);
    void logisticfit4ec(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, double cnstrleft, double cnstrright, double* a, double* b, double* c, double* d, lsfitreport* rep, ae_state *_state);
    void logisticfit5(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, double* a, double* b, double* c, double* d, double* g, lsfitreport* rep, ae_state *_state);
    void logisticfit5ec(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, double cnstrleft, double cnstrright, double* a, double* b, double* c, double* d, double* g, lsfitreport* rep, ae_state *_state);
    void logisticfit45x(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, double cnstrleft, double cnstrright, ae_bool is4pl, double lambdav, double epsx, ae_int_t rscnt, double* a, double* b, double* c, double* d, double* g, lsfitreport* rep, ae_state *_state);
    void barycentricfitfloaterhormannwc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, barycentricinterpolant* b, barycentricfitreport* rep, ae_state *_state);
    void _pexec_barycentricfitfloaterhormannwc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, barycentricinterpolant* b, barycentricfitreport* rep, ae_state *_state);
    void barycentricfitfloaterhormann(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, barycentricinterpolant* b, barycentricfitreport* rep, ae_state *_state);
    void _pexec_barycentricfitfloaterhormann(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, barycentricinterpolant* b, barycentricfitreport* rep, ae_state *_state);
    void spline1dfitpenalized(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, double rho, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void _pexec_spline1dfitpenalized(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, double rho, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void spline1dfitpenalizedw(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, ae_int_t m, double rho, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void _pexec_spline1dfitpenalizedw(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, ae_int_t m, double rho, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void spline1dfitcubicwc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void _pexec_spline1dfitcubicwc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void spline1dfithermitewc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void _pexec_spline1dfithermitewc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void spline1dfitcubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void _pexec_spline1dfitcubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void spline1dfithermite(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void _pexec_spline1dfithermite(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_int_t n, ae_int_t m, ae_int_t* info, spline1dinterpolant* s, spline1dfitreport* rep, ae_state *_state);
    void lsfitlinearw(/* Real    */ ae_vector* y, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* fmatrix, ae_int_t n, ae_int_t m, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void _pexec_lsfitlinearw(/* Real    */ ae_vector* y, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* fmatrix, ae_int_t n, ae_int_t m, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void lsfitlinearwc(/* Real    */ ae_vector* y, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* fmatrix, /* Real    */ ae_matrix* cmatrix, ae_int_t n, ae_int_t m, ae_int_t k, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void _pexec_lsfitlinearwc(/* Real    */ ae_vector* y, /* Real    */ ae_vector* w, /* Real    */ ae_matrix* fmatrix, /* Real    */ ae_matrix* cmatrix, ae_int_t n, ae_int_t m, ae_int_t k, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void lsfitlinear(/* Real    */ ae_vector* y, /* Real    */ ae_matrix* fmatrix, ae_int_t n, ae_int_t m, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void _pexec_lsfitlinear(/* Real    */ ae_vector* y, /* Real    */ ae_matrix* fmatrix, ae_int_t n, ae_int_t m, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void lsfitlinearc(/* Real    */ ae_vector* y, /* Real    */ ae_matrix* fmatrix, /* Real    */ ae_matrix* cmatrix, ae_int_t n, ae_int_t m, ae_int_t k, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void _pexec_lsfitlinearc(/* Real    */ ae_vector* y, /* Real    */ ae_matrix* fmatrix, /* Real    */ ae_matrix* cmatrix, ae_int_t n, ae_int_t m, ae_int_t k, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void lsfitcreatewf(/* Real    */ ae_matrix* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, /* Real    */ ae_vector* c, ae_int_t n, ae_int_t m, ae_int_t k, double diffstep, lsfitstate* state, ae_state *_state);
    void lsfitcreatef(/* Real    */ ae_matrix* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* c, ae_int_t n, ae_int_t m, ae_int_t k, double diffstep, lsfitstate* state, ae_state *_state);
    void lsfitcreatewfg(/* Real    */ ae_matrix* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, /* Real    */ ae_vector* c, ae_int_t n, ae_int_t m, ae_int_t k, ae_bool cheapfg, lsfitstate* state, ae_state *_state);
    void lsfitcreatefg(/* Real    */ ae_matrix* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* c, ae_int_t n, ae_int_t m, ae_int_t k, ae_bool cheapfg, lsfitstate* state, ae_state *_state);
    void lsfitcreatewfgh(/* Real    */ ae_matrix* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, /* Real    */ ae_vector* c, ae_int_t n, ae_int_t m, ae_int_t k, lsfitstate* state, ae_state *_state);
    void lsfitcreatefgh(/* Real    */ ae_matrix* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* c, ae_int_t n, ae_int_t m, ae_int_t k, lsfitstate* state, ae_state *_state);
    void lsfitsetcond(lsfitstate* state, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void lsfitsetstpmax(lsfitstate* state, double stpmax, ae_state *_state);
    void lsfitsetxrep(lsfitstate* state, ae_bool needxrep, ae_state *_state);
    void lsfitsetscale(lsfitstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void lsfitsetbc(lsfitstate* state, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    ae_bool lsfititeration(lsfitstate* state, ae_state *_state);
    void lsfitresults(lsfitstate* state, ae_int_t* info, /* Real    */ ae_vector* c, lsfitreport* rep, ae_state *_state);
    void lsfitsetgradientcheck(lsfitstate* state, double teststep, ae_state *_state);
    void lsfitscalexy(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_vector* w, ae_int_t n, /* Real    */ ae_vector* xc, /* Real    */ ae_vector* yc, /* Integer */ ae_vector* dc, ae_int_t k, double* xa, double* xb, double* sa, double* sb, /* Real    */ ae_vector* xoriginal, /* Real    */ ae_vector* yoriginal, ae_state *_state);
    void _polynomialfitreport_init(void* _p, ae_state *_state);
    void _polynomialfitreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _polynomialfitreport_clear(void* _p);
    void _polynomialfitreport_destroy(void* _p);
    void _barycentricfitreport_init(void* _p, ae_state *_state);
    void _barycentricfitreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _barycentricfitreport_clear(void* _p);
    void _barycentricfitreport_destroy(void* _p);
    void _spline1dfitreport_init(void* _p, ae_state *_state);
    void _spline1dfitreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _spline1dfitreport_clear(void* _p);
    void _spline1dfitreport_destroy(void* _p);
    void _lsfitreport_init(void* _p, ae_state *_state);
    void _lsfitreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _lsfitreport_clear(void* _p);
    void _lsfitreport_destroy(void* _p);
    void _lsfitstate_init(void* _p, ae_state *_state);
    void _lsfitstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _lsfitstate_clear(void* _p);
    void _lsfitstate_destroy(void* _p);
    void pspline2build(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline2interpolant* p, ae_state *_state);
    void pspline3build(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline3interpolant* p, ae_state *_state);
    void pspline2buildperiodic(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline2interpolant* p, ae_state *_state);
    void pspline3buildperiodic(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t st, ae_int_t pt, pspline3interpolant* p, ae_state *_state);
    void pspline2parametervalues(pspline2interpolant* p, ae_int_t* n, /* Real    */ ae_vector* t, ae_state *_state);
    void pspline3parametervalues(pspline3interpolant* p, ae_int_t* n, /* Real    */ ae_vector* t, ae_state *_state);
    void pspline2calc(pspline2interpolant* p, double t, double* x, double* y, ae_state *_state);
    void pspline3calc(pspline3interpolant* p, double t, double* x, double* y, double* z, ae_state *_state);
    void pspline2tangent(pspline2interpolant* p, double t, double* x, double* y, ae_state *_state);
    void pspline3tangent(pspline3interpolant* p, double t, double* x, double* y, double* z, ae_state *_state);
    void pspline2diff(pspline2interpolant* p, double t, double* x, double* dx, double* y, double* dy, ae_state *_state);
    void pspline3diff(pspline3interpolant* p, double t, double* x, double* dx, double* y, double* dy, double* z, double* dz, ae_state *_state);
    void pspline2diff2(pspline2interpolant* p, double t, double* x, double* dx, double* d2x, double* y, double* dy, double* d2y, ae_state *_state);
    void pspline3diff2(pspline3interpolant* p, double t, double* x, double* dx, double* d2x, double* y, double* dy, double* d2y, double* z, double* dz, double* d2z, ae_state *_state);
    double pspline2arclength(pspline2interpolant* p, double a, double b, ae_state *_state);
    double pspline3arclength(pspline3interpolant* p, double a, double b, ae_state *_state);
    void parametricrdpfixed(/* Real    */ ae_matrix* x, ae_int_t n, ae_int_t d, ae_int_t stopm, double stopeps, /* Real    */ ae_matrix* x2, /* Integer */ ae_vector* idx2, ae_int_t* nsections, ae_state *_state);
    void _pspline2interpolant_init(void* _p, ae_state *_state);
    void _pspline2interpolant_init_copy(void* _dst, void* _src, ae_state *_state);
    void _pspline2interpolant_clear(void* _p);
    void _pspline2interpolant_destroy(void* _p);
    void _pspline3interpolant_init(void* _p, ae_state *_state);
    void _pspline3interpolant_init_copy(void* _dst, void* _src, ae_state *_state);
    void _pspline3interpolant_clear(void* _p);
    void _pspline3interpolant_destroy(void* _p);
    void rbfcreate(ae_int_t nx, ae_int_t ny, rbfmodel* s, ae_state *_state);
    void rbfsetpoints(rbfmodel* s, /* Real    */ ae_matrix* xy, ae_int_t n, ae_state *_state);
    void rbfsetalgoqnn(rbfmodel* s, double q, double z, ae_state *_state);
    void rbfsetalgomultilayer(rbfmodel* s, double rbase, ae_int_t nlayers, double lambdav, ae_state *_state);
    void rbfsetlinterm(rbfmodel* s, ae_state *_state);
    void rbfsetconstterm(rbfmodel* s, ae_state *_state);
    void rbfsetzeroterm(rbfmodel* s, ae_state *_state);
    void rbfsetcond(rbfmodel* s, double epsort, double epserr, ae_int_t maxits, ae_state *_state);
    void rbfbuildmodel(rbfmodel* s, rbfreport* rep, ae_state *_state);
    double rbfcalc2(rbfmodel* s, double x0, double x1, ae_state *_state);
    double rbfcalc3(rbfmodel* s, double x0, double x1, double x2, ae_state *_state);
    void rbfcalc(rbfmodel* s, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void rbfcalcbuf(rbfmodel* s, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void rbfgridcalc2(rbfmodel* s, /* Real    */ ae_vector* x0, ae_int_t n0, /* Real    */ ae_vector* x1, ae_int_t n1, /* Real    */ ae_matrix* y, ae_state *_state);
    void rbfunpack(rbfmodel* s, ae_int_t* nx, ae_int_t* ny, /* Real    */ ae_matrix* xwr, ae_int_t* nc, /* Real    */ ae_matrix* v, ae_state *_state);
    void rbfalloc(ae_serializer* s, rbfmodel* model, ae_state *_state);
    void rbfserialize(ae_serializer* s, rbfmodel* model, ae_state *_state);
    void rbfunserialize(ae_serializer* s, rbfmodel* model, ae_state *_state);
    void _rbfmodel_init(void* _p, ae_state *_state);
    void _rbfmodel_init_copy(void* _dst, void* _src, ae_state *_state);
    void _rbfmodel_clear(void* _p);
    void _rbfmodel_destroy(void* _p);
    void _rbfreport_init(void* _p, ae_state *_state);
    void _rbfreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _rbfreport_clear(void* _p);
    void _rbfreport_destroy(void* _p);
    double spline2dcalc(spline2dinterpolant* c, double x, double y, ae_state *_state);
    void spline2ddiff(spline2dinterpolant* c, double x, double y, double* f, double* fx, double* fy, double* fxy, ae_state *_state);
    void spline2dlintransxy(spline2dinterpolant* c, double ax, double bx, double ay, double by, ae_state *_state);
    void spline2dlintransf(spline2dinterpolant* c, double a, double b, ae_state *_state);
    void spline2dcopy(spline2dinterpolant* c, spline2dinterpolant* cc, ae_state *_state);
    void spline2dresamplebicubic(/* Real    */ ae_matrix* a, ae_int_t oldheight, ae_int_t oldwidth, /* Real    */ ae_matrix* b, ae_int_t newheight, ae_int_t newwidth, ae_state *_state);
    void spline2dresamplebilinear(/* Real    */ ae_matrix* a, ae_int_t oldheight, ae_int_t oldwidth, /* Real    */ ae_matrix* b, ae_int_t newheight, ae_int_t newwidth, ae_state *_state);
    void spline2dbuildbilinearv(/* Real    */ ae_vector* x, ae_int_t n, /* Real    */ ae_vector* y, ae_int_t m, /* Real    */ ae_vector* f, ae_int_t d, spline2dinterpolant* c, ae_state *_state);
    void spline2dbuildbicubicv(/* Real    */ ae_vector* x, ae_int_t n, /* Real    */ ae_vector* y, ae_int_t m, /* Real    */ ae_vector* f, ae_int_t d, spline2dinterpolant* c, ae_state *_state);
    void spline2dcalcvbuf(spline2dinterpolant* c, double x, double y, /* Real    */ ae_vector* f, ae_state *_state);
    void spline2dcalcv(spline2dinterpolant* c, double x, double y, /* Real    */ ae_vector* f, ae_state *_state);
    void spline2dunpackv(spline2dinterpolant* c, ae_int_t* m, ae_int_t* n, ae_int_t* d, /* Real    */ ae_matrix* tbl, ae_state *_state);
    void spline2dbuildbilinear(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_matrix* f, ae_int_t m, ae_int_t n, spline2dinterpolant* c, ae_state *_state);
    void spline2dbuildbicubic(/* Real    */ ae_vector* x, /* Real    */ ae_vector* y, /* Real    */ ae_matrix* f, ae_int_t m, ae_int_t n, spline2dinterpolant* c, ae_state *_state);
    void spline2dunpack(spline2dinterpolant* c, ae_int_t* m, ae_int_t* n, /* Real    */ ae_matrix* tbl, ae_state *_state);
    void _spline2dinterpolant_init(void* _p, ae_state *_state);
    void _spline2dinterpolant_init_copy(void* _dst, void* _src, ae_state *_state);
    void _spline2dinterpolant_clear(void* _p);
    void _spline2dinterpolant_destroy(void* _p);
    double spline3dcalc(spline3dinterpolant* c, double x, double y, double z, ae_state *_state);
    void spline3dlintransxyz(spline3dinterpolant* c, double ax, double bx, double ay, double by, double az, double bz, ae_state *_state);
    void spline3dlintransf(spline3dinterpolant* c, double a, double b, ae_state *_state);
    void spline3dcopy(spline3dinterpolant* c, spline3dinterpolant* cc, ae_state *_state);
    void spline3dresampletrilinear(/* Real    */ ae_vector* a, ae_int_t oldzcount, ae_int_t oldycount, ae_int_t oldxcount, ae_int_t newzcount, ae_int_t newycount, ae_int_t newxcount, /* Real    */ ae_vector* b, ae_state *_state);
    void spline3dbuildtrilinearv(/* Real    */ ae_vector* x, ae_int_t n, /* Real    */ ae_vector* y, ae_int_t m, /* Real    */ ae_vector* z, ae_int_t l, /* Real    */ ae_vector* f, ae_int_t d, spline3dinterpolant* c, ae_state *_state);
    void spline3dcalcvbuf(spline3dinterpolant* c, double x, double y, double z, /* Real    */ ae_vector* f, ae_state *_state);
    void spline3dcalcv(spline3dinterpolant* c, double x, double y, double z, /* Real    */ ae_vector* f, ae_state *_state);
    void spline3dunpackv(spline3dinterpolant* c, ae_int_t* n, ae_int_t* m, ae_int_t* l, ae_int_t* d, ae_int_t* stype, /* Real    */ ae_matrix* tbl, ae_state *_state);
    void _spline3dinterpolant_init(void* _p, ae_state *_state);
    void _spline3dinterpolant_init_copy(void* _dst, void* _src, ae_state *_state);
    void _spline3dinterpolant_clear(void* _p);
    void _spline3dinterpolant_destroy(void* _p);
}

