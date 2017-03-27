module alglib.optimization;
import alglib.ap;
import alglib.internal;
import alglib.misc;
import alglib.linalg;
import alglib.solvers;

struct precbuflbfgs
{
    ae_vector norms;
    ae_vector alpha;
    ae_vector rho;
    ae_matrix yk;
    ae_vector idx;
    ae_vector bufa;
    ae_vector bufb;
} 

struct precbuflowrank
{
    ae_int_t n;
    ae_int_t k;
    ae_vector d;
    ae_matrix v;
    ae_vector bufc;
    ae_matrix bufz;
    ae_matrix bufw;
    ae_vector tmp;
}

struct convexquadraticmodel
{
    ae_int_t n;
    ae_int_t k;
    double alpha;
    double tau;
    double theta;
    ae_matrix a;
    ae_matrix q;
    ae_vector b;
    ae_vector r;
    ae_vector xc;
    ae_vector d;
    ae_vector activeset;
    ae_matrix tq2dense;
    ae_matrix tk2;
    ae_vector tq2diag;
    ae_vector tq1;
    ae_vector tk1;
    double tq0;
    double tk0;
    ae_vector txc;
    ae_vector tb;
    ae_int_t nfree;
    ae_int_t ecakind;
    ae_matrix ecadense;
    ae_matrix eq;
    ae_matrix eccm;
    ae_vector ecadiag;
    ae_vector eb;
    double ec;
    ae_vector tmp0;
    ae_vector tmp1;
    ae_vector tmpg;
    ae_matrix tmp2;
    ae_bool ismaintermchanged;
    ae_bool issecondarytermchanged;
    ae_bool islineartermchanged;
    ae_bool isactivesetchanged;
}

struct snnlssolver
{
    ae_int_t ns;
    ae_int_t nd;
    ae_int_t nr;
    ae_matrix densea;
    ae_vector b;
    ae_vector nnc;
    double debugflops;
    ae_int_t debugmaxinnerits;
    ae_vector xn;
    ae_vector xp;
    ae_matrix tmpz;
    ae_matrix tmpca;
    ae_matrix tmplq;
    ae_matrix trda;
    ae_vector trdd;
    ae_vector crb;
    ae_vector g;
    ae_vector d;
    ae_vector dx;
    ae_vector diagaa;
    ae_vector cb;
    ae_vector cx;
    ae_vector cborg;
    ae_vector tmpcholesky;
    ae_vector r;
    ae_vector regdiag;
    ae_vector tmp0;
    ae_vector tmp1;
    ae_vector tmp2;
    ae_vector rdtmprowmap;
}

struct sactiveset
{
    ae_int_t n;
    ae_int_t algostate;
    ae_vector xc;
    ae_bool hasxc;
    ae_vector s;
    ae_vector h;
    ae_vector activeset;
    ae_bool basisisready;
    ae_matrix sbasis;
    ae_matrix pbasis;
    ae_matrix ibasis;
    ae_int_t basissize;
    ae_bool constraintschanged;
    ae_vector hasbndl;
    ae_vector hasbndu;
    ae_vector bndl;
    ae_vector bndu;
    ae_matrix cleic;
    ae_int_t nec;
    ae_int_t nic;
    ae_vector mtx;
    ae_vector mtas;
    ae_vector cdtmp;
    ae_vector corrtmp;
    ae_vector unitdiagonal;
    snnlssolver solver;
    ae_vector scntmp;
    ae_vector tmp0;
    ae_vector tmpfeas;
    ae_matrix tmpm0;
    ae_vector rctmps;
    ae_vector rctmpg;
    ae_vector rctmprightpart;
    ae_matrix rctmpdense0;
    ae_matrix rctmpdense1;
    ae_vector rctmpisequality;
    ae_vector rctmpconstraintidx;
    ae_vector rctmplambdas;
    ae_matrix tmpbasis;
}

struct mincgstate
{
    ae_int_t n;
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
    double stpmax;
    double suggestedstep;
    ae_bool xrep;
    ae_bool drep;
    ae_int_t cgtype;
    ae_int_t prectype;
    ae_vector diagh;
    ae_vector diaghl2;
    ae_matrix vcorr;
    ae_int_t vcnt;
    ae_vector s;
    double diffstep;
    ae_int_t nfev;
    ae_int_t mcstage;
    ae_int_t k;
    ae_vector xk;
    ae_vector dk;
    ae_vector xn;
    ae_vector dn;
    ae_vector d;
    double fold;
    double stp;
    double curstpmax;
    ae_vector yk;
    double lastgoodstep;
    double lastscaledstep;
    ae_int_t mcinfo;
    ae_bool innerresetneeded;
    ae_bool terminationneeded;
    double trimthreshold;
    ae_int_t rstimer;
    ae_vector x;
    double f;
    ae_vector g;
    ae_bool needf;
    ae_bool needfg;
    ae_bool xupdated;
    ae_bool algpowerup;
    ae_bool lsstart;
    ae_bool lsend;
    ae_bool userterminationneeded;
    double teststep;
    rcommstate rstate;
    ae_int_t repiterationscount;
    ae_int_t repnfev;
    ae_int_t repvaridx;
    ae_int_t repterminationtype;
    ae_int_t debugrestartscount;
    linminstate lstate;
    double fbase;
    double fm2;
    double fm1;
    double fp1;
    double fp2;
    double betahs;
    double betady;
    ae_vector work0;
    ae_vector work1;
}

struct mincgreport
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t varidx;
    ae_int_t terminationtype;
}

struct minbleicstate
{
    ae_int_t nmain;
    ae_int_t nslack;
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
    ae_bool xrep;
    ae_bool drep;
    double stpmax;
    double diffstep;
    sactiveset sas;
    ae_vector s;
    ae_int_t prectype;
    ae_vector diagh;
    ae_vector x;
    double f;
    ae_vector g;
    ae_bool needf;
    ae_bool needfg;
    ae_bool xupdated;
    ae_bool lsstart;
    ae_bool steepestdescentstep;
    ae_bool boundedstep;
    ae_bool userterminationneeded;
    double teststep;
    rcommstate rstate;
    ae_vector ugc;
    ae_vector cgc;
    ae_vector xn;
    ae_vector ugn;
    ae_vector cgn;
    ae_vector xp;
    double fc;
    double fn;
    double fp;
    ae_vector d;
    ae_matrix cleic;
    ae_int_t nec;
    ae_int_t nic;
    double lastgoodstep;
    double lastscaledgoodstep;
    double maxscaledgrad;
    ae_vector hasbndl;
    ae_vector hasbndu;
    ae_vector bndl;
    ae_vector bndu;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repnfev;
    ae_int_t repvaridx;
    ae_int_t repterminationtype;
    double repdebugeqerr;
    double repdebugfs;
    double repdebugff;
    double repdebugdx;
    ae_int_t repdebugfeasqpits;
    ae_int_t repdebugfeasgpaits;
    ae_vector xstart;
    snnlssolver solver;
    double fbase;
    double fm2;
    double fm1;
    double fp1;
    double fp2;
    double xm1;
    double xp1;
    double gm1;
    double gp1;
    ae_int_t cidx;
    double cval;
    ae_vector tmpprec;
    ae_vector tmp0;
    ae_int_t nfev;
    ae_int_t mcstage;
    double stp;
    double curstpmax;
    double activationstep;
    ae_vector work;
    linminstate lstate;
    double trimthreshold;
    ae_int_t nonmonotoniccnt;
    ae_matrix bufyk;
    ae_matrix bufsk;
    ae_vector bufrho;
    ae_vector buftheta;
    ae_int_t bufsize;
}

struct minbleicreport
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t varidx;
    ae_int_t terminationtype;
    double debugeqerr;
    double debugfs;
    double debugff;
    double debugdx;
    ae_int_t debugfeasqpits;
    ae_int_t debugfeasgpaits;
    ae_int_t inneriterationscount;
    ae_int_t outeriterationscount;
}

struct minlbfgsstate
{
    ae_int_t n;
    ae_int_t m;
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
    ae_bool xrep;
    double stpmax;
    ae_vector s;
    double diffstep;
    ae_int_t nfev;
    ae_int_t mcstage;
    ae_int_t k;
    ae_int_t q;
    ae_int_t p;
    ae_vector rho;
    ae_matrix yk;
    ae_matrix sk;
    ae_vector xp;
    ae_vector theta;
    ae_vector d;
    double stp;
    ae_vector work;
    double fold;
    double trimthreshold;
    ae_int_t prectype;
    double gammak;
    ae_matrix denseh;
    ae_vector diagh;
    ae_vector precc;
    ae_vector precd;
    ae_matrix precw;
    ae_int_t preck;
    precbuflbfgs precbuf;
    precbuflowrank lowrankbuf;
    double fbase;
    double fm2;
    double fm1;
    double fp1;
    double fp2;
    ae_vector autobuf;
    ae_vector x;
    double f;
    ae_vector g;
    ae_bool needf;
    ae_bool needfg;
    ae_bool xupdated;
    ae_bool userterminationneeded;
    double teststep;
    rcommstate rstate;
    ae_int_t repiterationscount;
    ae_int_t repnfev;
    ae_int_t repvaridx;
    ae_int_t repterminationtype;
    linminstate lstate;
}

struct minlbfgsreport
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t varidx;
    ae_int_t terminationtype;
}

struct qqpsettings
{
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxouterits;
    ae_bool cgphase;
    ae_bool cnphase;
    ae_int_t cgminits;
    ae_int_t cgmaxits;
    ae_int_t cnmaxupdates;
    ae_int_t sparsesolver;
}

struct qqpbuffers
{
    ae_int_t n;
    ae_int_t nmain;
    ae_int_t nslack;
    ae_int_t nec;
    ae_int_t nic;
    ae_int_t akind;
    ae_matrix densea;
    sparsematrix sparsea;
    ae_bool sparseupper;
    double absamax;
    double absasum;
    double absasum2;
    ae_vector b;
    ae_vector bndl;
    ae_vector bndu;
    ae_vector havebndl;
    ae_vector havebndu;
    ae_matrix cleic;
    ae_vector xs;
    ae_vector gc;
    ae_vector xp;
    ae_vector dc;
    ae_vector dp;
    ae_vector cgc;
    ae_vector cgp;
    sactiveset sas;
    ae_vector activated;
    ae_int_t nfree;
    ae_int_t cnmodelage;
    ae_matrix densez;
    sparsematrix sparsecca;
    ae_vector yidx;
    ae_vector regdiag;
    ae_vector regx0;
    ae_vector tmpcn;
    ae_vector tmpcni;
    ae_vector tmpcnb;
    ae_vector tmp0;
    ae_vector stpbuf;
    sparsebuffers sbuf;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repncholesky;
    ae_int_t repncupdates;
}

struct qpbleicsettings
{
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
}

struct qpbleicbuffers
{
    minbleicstate solver;
    minbleicreport solverrep;
    ae_vector tmp0;
    ae_vector tmp1;
    ae_vector tmpi;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
}

struct qpcholeskysettings
{
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
}

struct qpcholeskybuffers
{
    sactiveset sas;
    ae_vector pg;
    ae_vector gc;
    ae_vector xs;
    ae_vector xn;
    ae_vector workbndl;
    ae_vector workbndu;
    ae_vector havebndl;
    ae_vector havebndu;
    ae_matrix workcleic;
    ae_vector rctmpg;
    ae_vector tmp0;
    ae_vector tmp1;
    ae_vector tmpb;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repncholesky;
}

struct minqpstate
{
    ae_int_t n;
    qqpsettings qqpsettingsuser;
    qqpsettings qqpsettingscurrent;
    qpbleicsettings qpbleicsettingsuser;
    qpbleicsettings qpbleicsettingscurrent;
    ae_int_t algokind;
    ae_int_t akind;
    convexquadraticmodel a;
    sparsematrix sparsea;
    ae_bool sparseaupper;
    double absamax;
    double absasum;
    double absasum2;
    ae_vector b;
    ae_vector bndl;
    ae_vector bndu;
    ae_vector s;
    ae_vector havebndl;
    ae_vector havebndu;
    ae_vector xorigin;
    ae_vector startx;
    ae_bool havex;
    ae_matrix cleic;
    ae_int_t nec;
    ae_int_t nic;
    ae_vector xs;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repncholesky;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    ae_vector tmp0;
    ae_bool qpbleicfirstcall;
    qpbleicbuffers qpbleicbuf;
    qqpbuffers qqpbuf;
    qpcholeskybuffers qpcholeskybuf;
    normestimatorstate estimator;
}

struct minqpreport
{
    ae_int_t inneriterationscount;
    ae_int_t outeriterationscount;
    ae_int_t nmv;
    ae_int_t ncholesky;
    ae_int_t terminationtype;
}

struct minlmstate
{
    ae_int_t n;
    ae_int_t m;
    double diffstep;
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
    ae_bool xrep;
    double stpmax;
    ae_int_t maxmodelage;
    ae_bool makeadditers;
    ae_vector x;
    double f;
    ae_vector fi;
    ae_matrix j;
    ae_matrix h;
    ae_vector g;
    ae_bool needf;
    ae_bool needfg;
    ae_bool needfgh;
    ae_bool needfij;
    ae_bool needfi;
    ae_bool xupdated;
    ae_bool userterminationneeded;
    ae_int_t algomode;
    ae_bool hasf;
    ae_bool hasfi;
    ae_bool hasg;
    ae_vector xbase;
    double fbase;
    ae_vector fibase;
    ae_vector gbase;
    ae_matrix quadraticmodel;
    ae_vector bndl;
    ae_vector bndu;
    ae_vector havebndl;
    ae_vector havebndu;
    ae_vector s;
    double lambdav;
    double nu;
    ae_int_t modelage;
    ae_vector xdir;
    ae_vector deltax;
    ae_vector deltaf;
    ae_bool deltaxready;
    ae_bool deltafready;
    double teststep;
    ae_int_t repiterationscount;
    ae_int_t repterminationtype;
    ae_int_t repfuncidx;
    ae_int_t repvaridx;
    ae_int_t repnfunc;
    ae_int_t repnjac;
    ae_int_t repngrad;
    ae_int_t repnhess;
    ae_int_t repncholesky;
    rcommstate rstate;
    ae_vector choleskybuf;
    ae_vector tmp0;
    double actualdecrease;
    double predicteddecrease;
    double xm1;
    double xp1;
    ae_vector fm1;
    ae_vector fp1;
    ae_vector fc1;
    ae_vector gm1;
    ae_vector gp1;
    ae_vector gc1;
    minlbfgsstate internalstate;
    minlbfgsreport internalrep;
    minqpstate qpstate;
    minqpreport qprep;
}

struct  minlmreport
{
    ae_int_t iterationscount;
    ae_int_t terminationtype;
    ae_int_t funcidx;
    ae_int_t varidx;
    ae_int_t nfunc;
    ae_int_t njac;
    ae_int_t ngrad;
    ae_int_t nhess;
    ae_int_t ncholesky;
}

struct minasastate
{
    ae_int_t n;
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
    ae_bool xrep;
    double stpmax;
    ae_int_t cgtype;
    ae_int_t k;
    ae_int_t nfev;
    ae_int_t mcstage;
    ae_vector bndl;
    ae_vector bndu;
    ae_int_t curalgo;
    ae_int_t acount;
    double mu;
    double finit;
    double dginit;
    ae_vector ak;
    ae_vector xk;
    ae_vector dk;
    ae_vector an;
    ae_vector xn;
    ae_vector dn;
    ae_vector d;
    double fold;
    double stp;
    ae_vector work;
    ae_vector yk;
    ae_vector gc;
    double laststep;
    ae_vector x;
    double f;
    ae_vector g;
    ae_bool needfg;
    ae_bool xupdated;
    rcommstate rstate;
    ae_int_t repiterationscount;
    ae_int_t repnfev;
    ae_int_t repterminationtype;
    ae_int_t debugrestartscount;
    linminstate lstate;
    double betahs;
    double betady;
} 

struct minasareport
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
    ae_int_t activeconstraints;
}

struct minnlcstate
{
    double stabilizingpoint;
    double initialinequalitymultiplier;
    ae_int_t solvertype;
    ae_int_t prectype;
    ae_int_t updatefreq;
    double rho;
    ae_int_t n;
    double epsg;
    double epsf;
    double epsx;
    ae_int_t maxits;
    ae_int_t aulitscnt;
    ae_bool xrep;
    double diffstep;
    double teststep;
    ae_vector s;
    ae_vector bndl;
    ae_vector bndu;
    ae_vector hasbndl;
    ae_vector hasbndu;
    ae_int_t nec;
    ae_int_t nic;
    ae_matrix cleic;
    ae_int_t ng;
    ae_int_t nh;
    ae_vector x;
    double f;
    ae_vector fi;
    ae_matrix j;
    ae_bool needfij;
    ae_bool needfi;
    ae_bool xupdated;
    rcommstate rstate;
    rcommstate rstateaul;
    ae_vector scaledbndl;
    ae_vector scaledbndu;
    ae_matrix scaledcleic;
    ae_vector xc;
    ae_vector xstart;
    ae_vector xbase;
    ae_vector fbase;
    ae_vector dfbase;
    ae_vector fm2;
    ae_vector fm1;
    ae_vector fp1;
    ae_vector fp2;
    ae_vector dfm1;
    ae_vector dfp1;
    ae_vector bufd;
    ae_vector bufc;
    ae_matrix bufw;
    ae_vector xk;
    ae_vector xk1;
    ae_vector gk;
    ae_vector gk1;
    double gammak;
    ae_bool xkpresent;
    minlbfgsstate auloptimizer;
    minlbfgsreport aulreport;
    ae_vector nubc;
    ae_vector nulc;
    ae_vector nunlc;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repnfev;
    ae_int_t repvaridx;
    ae_int_t repfuncidx;
    ae_int_t repterminationtype;
    ae_int_t repdbgphase0its;
}

struct minnlcreport
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t varidx;
    ae_int_t funcidx;
    ae_int_t terminationtype;
    ae_int_t dbgphase0its;
}

struct minnsqp
{
    double fc;
    double fn;
    ae_vector xc;
    ae_vector xn;
    ae_vector x0;
    ae_vector gc;
    ae_vector d;
    ae_matrix uh;
    ae_matrix ch;
    ae_matrix rk;
    ae_vector invutc;
    ae_vector tmp0;
    ae_vector tmpidx;
    ae_vector tmpd;
    ae_vector tmpc;
    ae_vector tmplambdas;
    ae_matrix tmpc2;
    ae_vector tmpb;
    snnlssolver nnls;
}

struct minnsstate
{
    ae_int_t solvertype;
    ae_int_t n;
    double epsx;
    ae_int_t maxits;
    ae_bool xrep;
    double diffstep;
    ae_vector s;
    ae_vector bndl;
    ae_vector bndu;
    ae_vector hasbndl;
    ae_vector hasbndu;
    ae_int_t nec;
    ae_int_t nic;
    ae_matrix cleic;
    ae_int_t ng;
    ae_int_t nh;
    ae_vector x;
    double f;
    ae_vector fi;
    ae_matrix j;
    ae_bool needfij;
    ae_bool needfi;
    ae_bool xupdated;
    rcommstate rstate;
    rcommstate rstateags;
    hqrndstate agsrs;
    double agsradius;
    ae_int_t agssamplesize;
    double agsraddecay;
    double agsalphadecay;
    double agsdecrease;
    double agsinitstp;
    double agsstattold;
    double agsshortstpabs;
    double agsshortstprel;
    double agsshortf;
    ae_int_t agsshortlimit;
    double agsrhononlinear;
    ae_int_t agsminupdate;
    ae_int_t agsmaxraddecays;
    ae_int_t agsmaxbacktrack;
    ae_int_t agsmaxbacktracknonfull;
    double agspenaltylevel;
    double agspenaltyincrease;
    ae_vector xstart;
    ae_vector xc;
    ae_vector xn;
    ae_vector grs;
    ae_vector d;
    ae_vector colmax;
    ae_vector diagh;
    ae_vector signmin;
    ae_vector signmax;
    ae_bool userterminationneeded;
    ae_vector scaledbndl;
    ae_vector scaledbndu;
    ae_matrix scaledcleic;
    ae_vector rholinear;
    ae_matrix samplex;
    ae_matrix samplegm;
    ae_matrix samplegmbc;
    ae_vector samplef;
    ae_vector samplef0;
    minnsqp nsqp;
    ae_vector tmp0;
    ae_vector tmp1;
    ae_matrix tmp2;
    ae_vector tmp3;
    ae_vector xbase;
    ae_vector fp;
    ae_vector fm;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repnfev;
    ae_int_t repvaridx;
    ae_int_t repfuncidx;
    ae_int_t repterminationtype;
    double replcerr;
    double repnlcerr;
    ae_int_t dbgncholesky;
}

struct minnsreport
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    double cerr;
    double lcerr;
    double nlcerr;
    ae_int_t terminationtype;
    ae_int_t varidx;
    ae_int_t funcidx;
}


/**
extern(C++,alglib)
{
    class _mincgstate_owner
    {
            _mincgstate_owner();
            _mincgstate_owner(const _mincgstate_owner &rhs);
            _mincgstate_owner& operator=(const _mincgstate_owner &rhs);
            virtual ~_mincgstate_owner();
            alglib_impl::mincgstate* c_ptr();
            alglib_impl::mincgstate* c_ptr() const;
            protected:
                alglib_impl::mincgstate *p_struct;
    }

    class mincgstate : _mincgstate_owner
    {
        mincgstate();
        mincgstate(const mincgstate &rhs);
        mincgstate& operator=(const mincgstate &rhs);
        virtual ~mincgstate();
        ae_bool &needf;
        ae_bool &needfg;
        ae_bool &xupdated;
        double &f;
        real_1d_array g;
        real_1d_array x;
    }

    class _mincgreport_owner
    {
        _mincgreport_owner();
        _mincgreport_owner(const _mincgreport_owner &rhs);
        _mincgreport_owner& operator=(const _mincgreport_owner &rhs);
        virtual ~_mincgreport_owner();
        alglib_impl::mincgreport* c_ptr();
        alglib_impl::mincgreport* c_ptr() const;
        protected:
            alglib_impl::mincgreport *p_struct;
    }

    class mincgreport :  _mincgreport_owner
    {
        mincgreport();
        mincgreport(const mincgreport &rhs);
        mincgreport& operator=(const mincgreport &rhs);
        virtual ~mincgreport();
        ae_int_t &iterationscount;
        ae_int_t &nfev;
        ae_int_t &varidx;
        ae_int_t &terminationtype;

    }

    class _minbleicstate_owner
    {
        _minbleicstate_owner();
        _minbleicstate_owner(const _minbleicstate_owner &rhs);
        _minbleicstate_owner& operator=(const _minbleicstate_owner &rhs);
        virtual ~_minbleicstate_owner();
        alglib_impl::minbleicstate* c_ptr();
        alglib_impl::minbleicstate* c_ptr() const;
        protected:
            alglib_impl::minbleicstate *p_struct;
    }
    class minbleicstate :  _minbleicstate_owner
    {
        minbleicstate();
        minbleicstate(const minbleicstate &rhs);
        minbleicstate& operator=(const minbleicstate &rhs);
        virtual ~minbleicstate();
        ae_bool &needf;
        ae_bool &needfg;
        ae_bool &xupdated;
        double &f;
        real_1d_array g;
        real_1d_array x;
    }

    class _minbleicreport_owner
    {
        _minbleicreport_owner();
        _minbleicreport_owner(const _minbleicreport_owner &rhs);
        _minbleicreport_owner& operator=(const _minbleicreport_owner &rhs);
        virtual ~_minbleicreport_owner();
        alglib_impl::minbleicreport* c_ptr();
        alglib_impl::minbleicreport* c_ptr() const;
        protected:
            alglib_impl::minbleicreport *p_struct;
    }

    class minbleicreport :  _minbleicreport_owner
    {
        minbleicreport();
        minbleicreport(const minbleicreport &rhs);
        minbleicreport& operator=(const minbleicreport &rhs);
        virtual ~minbleicreport();
        ae_int_t &iterationscount;
        ae_int_t &nfev;
        ae_int_t &varidx;
        ae_int_t &terminationtype;
        double &debugeqerr;
        double &debugfs;
        double &debugff;
        double &debugdx;
        ae_int_t &debugfeasqpits;
        ae_int_t &debugfeasgpaits;
        ae_int_t &inneriterationscount;
        ae_int_t &outeriterationscount;
    }

    class _minlbfgsstate_owner
    {
        _minlbfgsstate_owner();
        _minlbfgsstate_owner(const _minlbfgsstate_owner &rhs);
        _minlbfgsstate_owner& operator=(const _minlbfgsstate_owner &rhs);
        virtual ~_minlbfgsstate_owner();
        alglib_impl::minlbfgsstate* c_ptr();
        alglib_impl::minlbfgsstate* c_ptr() const;
        protected:
            alglib_impl::minlbfgsstate *p_struct;
    }

    class minlbfgsstate : _minlbfgsstate_owner
    {
        minlbfgsstate();
        minlbfgsstate(const minlbfgsstate &rhs);
        minlbfgsstate& operator=(const minlbfgsstate &rhs);
        virtual ~minlbfgsstate();
        ae_bool &needf;
        ae_bool &needfg;
        ae_bool &xupdated;
        double &f;
        real_1d_array g;
        real_1d_array x;
    }

    class _minlbfgsreport_owner
    {
        _minlbfgsreport_owner();
        _minlbfgsreport_owner(const _minlbfgsreport_owner &rhs);
        _minlbfgsreport_owner& operator=(const _minlbfgsreport_owner &rhs);
        virtual ~_minlbfgsreport_owner();
        alglib_impl::minlbfgsreport* c_ptr();
        alglib_impl::minlbfgsreport* c_ptr() const;
        protected:
            alglib_impl::minlbfgsreport *p_struct;
    }

    class minlbfgsreport : _minlbfgsreport_owner
    {
        minlbfgsreport();
        minlbfgsreport(const minlbfgsreport &rhs);
        minlbfgsreport& operator=(const minlbfgsreport &rhs);
        virtual ~minlbfgsreport();
        ae_int_t &iterationscount;
        ae_int_t &nfev;
        ae_int_t &varidx;
        ae_int_t &terminationtype;

    }
    class _minqpstate_owner
    {
        _minqpstate_owner();
        _minqpstate_owner(const _minqpstate_owner &rhs);
        _minqpstate_owner& operator=(const _minqpstate_owner &rhs);
        virtual ~_minqpstate_owner();
        alglib_impl::minqpstate* c_ptr();
        alglib_impl::minqpstate* c_ptr() const;
        protected:
            alglib_impl::minqpstate *p_struct;
    }

    class minqpstate : _minqpstate_owner
    {
        minqpstate();
        minqpstate(const minqpstate &rhs);
        minqpstate& operator=(const minqpstate &rhs);
        virtual ~minqpstate();

    }

    class _minqpreport_owner
    {
        _minqpreport_owner();
        _minqpreport_owner(const _minqpreport_owner &rhs);
        _minqpreport_owner& operator=(const _minqpreport_owner &rhs);
        virtual ~_minqpreport_owner();
        alglib_impl::minqpreport* c_ptr();
        alglib_impl::minqpreport* c_ptr() const;
        protected:
            alglib_impl::minqpreport *p_struct;
    }

    class minqpreport : _minqpreport_owner
    {
        minqpreport();
        minqpreport(const minqpreport &rhs);
        minqpreport& operator=(const minqpreport &rhs);
        virtual ~minqpreport();
        ae_int_t &inneriterationscount;
        ae_int_t &outeriterationscount;
        ae_int_t &nmv;
        ae_int_t &ncholesky;
        ae_int_t &terminationtype;
    }

    class _minlmstate_owner
    {
        _minlmstate_owner();
        _minlmstate_owner(const _minlmstate_owner &rhs);
        _minlmstate_owner& operator=(const _minlmstate_owner &rhs);
        virtual ~_minlmstate_owner();
        alglib_impl::minlmstate* c_ptr();
        alglib_impl::minlmstate* c_ptr() const;
        protected:
            alglib_impl::minlmstate *p_struct;
    }

    class minlmstate :  _minlmstate_owner
    {
        minlmstate();
        minlmstate(const minlmstate &rhs);
        minlmstate& operator=(const minlmstate &rhs);
        virtual ~minlmstate();
        ae_bool &needf;
        ae_bool &needfg;
        ae_bool &needfgh;
        ae_bool &needfi;
        ae_bool &needfij;
        ae_bool &xupdated;
        double &f;
        real_1d_array fi;
        real_1d_array g;
        real_2d_array h;
        real_2d_array j;
        real_1d_array x;
    }

    class _minlmreport_owner
    {
        _minlmreport_owner();
        _minlmreport_owner(const _minlmreport_owner &rhs);
        _minlmreport_owner& operator=(const _minlmreport_owner &rhs);
        virtual ~_minlmreport_owner();
        alglib_impl::minlmreport* c_ptr();
        alglib_impl::minlmreport* c_ptr() const;
        protected:
            alglib_impl::minlmreport *p_struct;
    }

    class minlmreport :  _minlmreport_owner
    {
        minlmreport();
        minlmreport(const minlmreport &rhs);
        minlmreport& operator=(const minlmreport &rhs);
        virtual ~minlmreport();
        ae_int_t &iterationscount;
        ae_int_t &terminationtype;
        ae_int_t &funcidx;
        ae_int_t &varidx;
        ae_int_t &nfunc;
        ae_int_t &njac;
        ae_int_t &ngrad;
        ae_int_t &nhess;
        ae_int_t &ncholesky;

    }

    class _minasastate_owner
    {
        _minasastate_owner();
        _minasastate_owner(const _minasastate_owner &rhs);
        _minasastate_owner& operator=(const _minasastate_owner &rhs);
        virtual ~_minasastate_owner();
        alglib_impl::minasastate* c_ptr();
        alglib_impl::minasastate* c_ptr() const;
        protected:
            alglib_impl::minasastate *p_struct;
    }

    class minasastate : _minasastate_owner
    {
        public:
            minasastate();
            minasastate(const minasastate &rhs);
            minasastate& operator=(const minasastate &rhs);
            virtual ~minasastate();
            ae_bool &needfg;
            ae_bool &xupdated;
            double &f;
            real_1d_array g;
            real_1d_array x;
    }

    class _minasareport_owner
    {
        public:
            _minasareport_owner();
            _minasareport_owner(const _minasareport_owner &rhs);
            _minasareport_owner& operator=(const _minasareport_owner &rhs);
            virtual ~_minasareport_owner();
            alglib_impl::minasareport* c_ptr();
            alglib_impl::minasareport* c_ptr() const;
        protected:
            alglib_impl::minasareport *p_struct;
    }

    class minasareport : public _minasareport_owner
    {
        public:
            minasareport();
            minasareport(const minasareport &rhs);
            minasareport& operator=(const minasareport &rhs);
            virtual ~minasareport();
            ae_int_t &iterationscount;
            ae_int_t &nfev;
            ae_int_t &terminationtype;
            ae_int_t &activeconstraints;
    }

    class _minnlcstate_owner
    {
        public:
            _minnlcstate_owner();
            _minnlcstate_owner(const _minnlcstate_owner &rhs);
            _minnlcstate_owner& operator=(const _minnlcstate_owner &rhs);
            virtual ~_minnlcstate_owner();
            alglib_impl::minnlcstate* c_ptr();
            alglib_impl::minnlcstate* c_ptr() const;
        protected:
            alglib_impl::minnlcstate *p_struct;
    }

    class minnlcstate : _minnlcstate_owner
    {
        public:
            minnlcstate();
            minnlcstate(const minnlcstate &rhs);
            minnlcstate& operator=(const minnlcstate &rhs);
            virtual ~minnlcstate();
            ae_bool &needfi;
            ae_bool &needfij;
            ae_bool &xupdated;
            double &f;
            real_1d_array fi;
            real_2d_array j;
            real_1d_array x;
    }
    class _minnlcreport_owner
    {
        public:
            _minnlcreport_owner();
            _minnlcreport_owner(const _minnlcreport_owner &rhs);
            _minnlcreport_owner& operator=(const _minnlcreport_owner &rhs);
            virtual ~_minnlcreport_owner();
            alglib_impl::minnlcreport* c_ptr();
            alglib_impl::minnlcreport* c_ptr() const;
        protected:
            alglib_impl::minnlcreport *p_struct;
        };
        class minnlcreport : public _minnlcreport_owner
        {
        public:
            minnlcreport();
            minnlcreport(const minnlcreport &rhs);
            minnlcreport& operator=(const minnlcreport &rhs);
            virtual ~minnlcreport();
            ae_int_t &iterationscount;
            ae_int_t &nfev;
            ae_int_t &varidx;
            ae_int_t &funcidx;
            ae_int_t &terminationtype;
            ae_int_t &dbgphase0its;
    }

    class _minnsstate_owner
    {
        public:
            _minnsstate_owner();
            _minnsstate_owner(const _minnsstate_owner &rhs);
            _minnsstate_owner& operator=(const _minnsstate_owner &rhs);
            virtual ~_minnsstate_owner();
            alglib_impl::minnsstate* c_ptr();
            alglib_impl::minnsstate* c_ptr() const;
        protected:
            alglib_impl::minnsstate *p_struct;
    }

    class minnsstate : public _minnsstate_owner
    {
    public:
        minnsstate();
        minnsstate(const minnsstate &rhs);
        minnsstate& operator=(const minnsstate &rhs);
        virtual ~minnsstate();
        ae_bool &needfi;
        ae_bool &needfij;
        ae_bool &xupdated;
        double &f;
        real_1d_array fi;
        real_2d_array j;
        real_1d_array x;

    }
    class _minnsreport_owner
    {
        public:
            _minnsreport_owner();
            _minnsreport_owner(const _minnsreport_owner &rhs);
            _minnsreport_owner& operator=(const _minnsreport_owner &rhs);
            //virtual ~_minnsreport_owner();
            alglib_impl::minnsreport* c_ptr();
            alglib_impl::minnsreport* c_ptr() const;
        //protected:
          //  alglib_impl::minnsreport *p_struct;
    }
    
    class minnsreport : _minnsreport_owner
    {
        //minnsreport();
        //minnsreport(const minnsreport &rhs);
        //minnsreport& operator=(const minnsreport &rhs);
        // virtual ~minnsreport();
        ae_int_t* iterationscount;
        ae_int_t* nfev;
        double* cerr;
        double* lcerr;
        double* nlcerr;
        ae_int_t* terminationtype;
        ae_int_t* varidx;
        ae_int_t* funcidx;
    }

    void mincgcreate(const ae_int_t n, const real_1d_array &x, mincgstate &state);
    void mincgcreate(const real_1d_array &x, mincgstate &state);
    void mincgcreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, mincgstate &state);
    void mincgcreatef(const real_1d_array &x, const double diffstep, mincgstate &state);
    void mincgsetcond(const mincgstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits);
    void mincgsetscale(const mincgstate &state, const real_1d_array &s);
    void mincgsetxrep(const mincgstate &state, const bool needxrep);
    void mincgsetcgtype(const mincgstate &state, const ae_int_t cgtype);
    void mincgsetstpmax(const mincgstate &state, const double stpmax);
    void mincgsuggeststep(const mincgstate &state, const double stp);
    void mincgsetprecdefault(const mincgstate &state);
    void mincgsetprecdiag(const mincgstate &state, const real_1d_array &d);
    void mincgsetprecscale(const mincgstate &state);
    bool mincgiteration(const mincgstate &state);
    void mincgoptimize(mincgstate &state,
        void* function(const real_1d_array &x, double &func, void *ptr) func,
        void* function (const real_1d_array &x, double func, void *ptr) rep = null,
        void *ptr = null);
    void mincgoptimize(mincgstate &state,
        void* function(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) grad,
        void* function(const real_1d_array &x, double func, void *ptr) rep= null,
        void *ptr = null);
    void mincgresults(const mincgstate &state, real_1d_array &x, mincgreport &rep);
    void mincgresultsbuf(const mincgstate &state, real_1d_array &x, mincgreport &rep);
    void mincgrestartfrom(const mincgstate &state, const real_1d_array &x);
    void mincgrequesttermination(const mincgstate &state);
    void mincgsetgradientcheck(const mincgstate &state, const double teststep);
    void minbleiccreate(const ae_int_t n, const real_1d_array &x, minbleicstate &state);
    void minbleiccreate(const real_1d_array &x, minbleicstate &state);
    void minbleiccreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minbleicstate &state);
    void minbleiccreatef(const real_1d_array &x, const double diffstep, minbleicstate &state);
    void minbleicsetbc(const minbleicstate &state, const real_1d_array &bndl, const real_1d_array &bndu);
    void minbleicsetlc(const minbleicstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
    void minbleicsetlc(const minbleicstate &state, const real_2d_array &c, const integer_1d_array &ct);
    void minbleicsetcond(const minbleicstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits);
    void minbleicsetscale(const minbleicstate &state, const real_1d_array &s);
    void minbleicsetprecdefault(const minbleicstate &state);
    void minbleicsetprecdiag(const minbleicstate &state, const real_1d_array &d);
    void minbleicsetprecscale(const minbleicstate &state);
    void minbleicsetxrep(const minbleicstate &state, const bool needxrep);
    void minbleicsetstpmax(const minbleicstate &state, const double stpmax);
    bool minbleiciteration(const minbleicstate &state);
    void minbleicoptimize(minbleicstate &state,
        void (*func)(const real_1d_array &x, double &func, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minbleicoptimize(minbleicstate &state,
        void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minbleicresults(const minbleicstate &state, real_1d_array &x, minbleicreport &rep);
    void minbleicresultsbuf(const minbleicstate &state, real_1d_array &x, minbleicreport &rep);
    void minbleicrestartfrom(const minbleicstate &state, const real_1d_array &x);
    void minbleicrequesttermination(const minbleicstate &state);
    void minbleicsetgradientcheck(const minbleicstate &state, const double teststep);
    void minlbfgscreate(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlbfgsstate &state);
    void minlbfgscreate(const ae_int_t m, const real_1d_array &x, minlbfgsstate &state);
    void minlbfgscreatef(const ae_int_t n, const ae_int_t m, const real_1d_array &x, const double diffstep, minlbfgsstate &state);
    void minlbfgscreatef(const ae_int_t m, const real_1d_array &x, const double diffstep, minlbfgsstate &state);
    void minlbfgssetcond(const minlbfgsstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits);
    void minlbfgssetxrep(const minlbfgsstate &state, const bool needxrep);
    void minlbfgssetstpmax(const minlbfgsstate &state, const double stpmax);
    void minlbfgssetscale(const minlbfgsstate &state, const real_1d_array &s);
    void minlbfgssetprecdefault(const minlbfgsstate &state);
    void minlbfgssetpreccholesky(const minlbfgsstate &state, const real_2d_array &p, const bool isupper);
    void minlbfgssetprecdiag(const minlbfgsstate &state, const real_1d_array &d);
    void minlbfgssetprecscale(const minlbfgsstate &state);
    bool minlbfgsiteration(const minlbfgsstate &state);
    void minlbfgsoptimize(minlbfgsstate &state,
        void (*func)(const real_1d_array &x, double &func, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minlbfgsoptimize(minlbfgsstate &state,
        void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minlbfgsresults(const minlbfgsstate &state, real_1d_array &x, minlbfgsreport &rep);
    void minlbfgsresultsbuf(const minlbfgsstate &state, real_1d_array &x, minlbfgsreport &rep);
    void minlbfgsrestartfrom(const minlbfgsstate &state, const real_1d_array &x);
    void minlbfgsrequesttermination(const minlbfgsstate &state);
    void minlbfgssetgradientcheck(const minlbfgsstate &state, const double teststep);
    void minqpcreate(const ae_int_t n, minqpstate &state);


    void minqpsetlinearterm(const minqpstate &state, const real_1d_array &b);
    void minqpsetquadraticterm(const minqpstate &state, const real_2d_array &a, const bool isupper);
    void minqpsetquadraticterm(const minqpstate &state, const real_2d_array &a);
    void minqpsetquadratictermsparse(const minqpstate &state, const sparsematrix &a, const bool isupper);
    void minqpsetstartingpoint(const minqpstate &state, const real_1d_array &x);
    void minqpsetorigin(const minqpstate &state, const real_1d_array &xorigin);
    void minqpsetscale(const minqpstate &state, const real_1d_array &s);
    void minqpsetalgocholesky(const minqpstate &state);
    void minqpsetalgobleic(const minqpstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits);
    void minqpsetalgoquickqp(const minqpstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxouterits, const bool usenewton);
    void minqpsetbc(const minqpstate &state, const real_1d_array &bndl, const real_1d_array &bndu);
    void minqpsetlc(const minqpstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
    void minqpsetlc(const minqpstate &state, const real_2d_array &c, const integer_1d_array &ct);
    void minqpoptimize(const minqpstate &state);
    void minqpresults(const minqpstate &state, real_1d_array &x, minqpreport &rep);
    void minqpresultsbuf(const minqpstate &state, real_1d_array &x, minqpreport &rep);
    void minlmcreatevj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmcreatevj(const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmcreatev(const ae_int_t n, const ae_int_t m, const real_1d_array &x, const double diffstep, minlmstate &state);
    void minlmcreatev(const ae_int_t m, const real_1d_array &x, const double diffstep, minlmstate &state);
    void minlmcreatefgh(const ae_int_t n, const real_1d_array &x, minlmstate &state);
    void minlmcreatefgh(const real_1d_array &x, minlmstate &state);
    void minlmsetcond(const minlmstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits);
    void minlmsetxrep(const minlmstate &state, const bool needxrep);
    void minlmsetstpmax(const minlmstate &state, const double stpmax);
    void minlmsetscale(const minlmstate &state, const real_1d_array &s);
    void minlmsetbc(const minlmstate &state, const real_1d_array &bndl, const real_1d_array &bndu);
    void minlmsetacctype(const minlmstate &state, const ae_int_t acctype);
    bool minlmiteration(const minlmstate &state);
    void minlmoptimize(minlmstate &state,
        void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minlmoptimize(minlmstate &state,
        void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr),
        void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minlmoptimize(minlmstate &state,
        void (*func)(const real_1d_array &x, double &func, void *ptr),
        void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
        void (*hess)(const real_1d_array &x, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minlmoptimize(minlmstate &state,
        void (*func)(const real_1d_array &x, double &func, void *ptr),
        void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minlmoptimize(minlmstate &state,
        void (*func)(const real_1d_array &x, double &func, void *ptr),
        void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
        void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minlmresults(const minlmstate &state, real_1d_array &x, minlmreport &rep);
    void minlmresultsbuf(const minlmstate &state, real_1d_array &x, minlmreport &rep);
    void minlmrestartfrom(const minlmstate &state, const real_1d_array &x);
    void minlmrequesttermination(const minlmstate &state);
    void minlmcreatevgj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmcreatevgj(const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmcreatefgj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmcreatefgj(const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmcreatefj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmcreatefj(const ae_int_t m, const real_1d_array &x, minlmstate &state);
    void minlmsetgradientcheck(const minlmstate &state, const double teststep);
    void minlbfgssetdefaultpreconditioner(const minlbfgsstate &state);
    void minlbfgssetcholeskypreconditioner(const minlbfgsstate &state, const real_2d_array &p, const bool isupper);
    void minbleicsetbarrierwidth(const minbleicstate &state, const double mu);
    void minbleicsetbarrierdecay(const minbleicstate &state, const double mudecay);
    void minasacreate(const ae_int_t n, const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu, minasastate &state);
    void minasacreate(const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu, minasastate &state);
    void minasasetcond(const minasastate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits);
    void minasasetxrep(const minasastate &state, const bool needxrep);
    void minasasetalgorithm(const minasastate &state, const ae_int_t algotype);
    void minasasetstpmax(const minasastate &state, const double stpmax);
    bool minasaiteration(const minasastate &state);
    void minasaoptimize(minasastate &state,
        void* function(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) grad,
        void* function(const real_1d_array &x, double func, void *ptr) rep= null,
        void *ptr = null);
    void minasaresults(const minasastate &state, real_1d_array &x, minasareport &rep);
    void minasaresultsbuf(const minasastate &state, real_1d_array &x, minasareport &rep);
    void minasarestartfrom(const minasastate &state, const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu);
    void minnlccreate(const ae_int_t n, const real_1d_array &x, minnlcstate &state);
    void minnlccreate(const real_1d_array &x, minnlcstate &state);
    void minnlccreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minnlcstate &state);
    void minnlccreatef(const real_1d_array &x, const double diffstep, minnlcstate &state);
    void minnlcsetbc(const minnlcstate &state, const real_1d_array &bndl, const real_1d_array &bndu);
    void minnlcsetlc(const minnlcstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
    void minnlcsetlc(const minnlcstate &state, const real_2d_array &c, const integer_1d_array &ct);
    void minnlcsetnlc(const minnlcstate &state, const ae_int_t nlec, const ae_int_t nlic);
    void minnlcsetcond(const minnlcstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits);
    void minnlcsetscale(const minnlcstate &state, const real_1d_array &s);
    void minnlcsetprecinexact(const minnlcstate &state);
    void minnlcsetprecnone(const minnlcstate &state);
    void minnlcsetalgoaul(const minnlcstate &state, const double rho, const ae_int_t itscnt);
    void minnlcsetxrep(const minnlcstate &state, const bool needxrep);
    bool minnlciteration(const minnlcstate &state);
    void minnlcoptimize(minnlcstate &state,
        void* function(const real_1d_array &x, real_1d_array &fi, void *ptr) fvec,
        void* function(const real_1d_array &x, double func, void *ptr) rep = null,
        void *ptr = null);
    void minnlcoptimize(minnlcstate &state,
        void* function(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr) jac,
        void* function(const real_1d_array &x, double func, void *ptr) rep= null,
        void *ptr = null);
    void minnlcresults(const minnlcstate &state, real_1d_array &x, minnlcreport &rep);
    void minnlcresultsbuf(const minnlcstate &state, real_1d_array &x, minnlcreport &rep);
    void minnlcrestartfrom(const minnlcstate &state, const real_1d_array &x);
    void minnlcsetgradientcheck(const minnlcstate &state, const double teststep);
    void minnscreate(const ae_int_t n, const real_1d_array &x, minnsstate &state);
    void minnscreate(const real_1d_array &x, minnsstate &state);
    void minnscreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minnsstate &state);
    void minnscreatef(const real_1d_array &x, const double diffstep, minnsstate &state);
    void minnssetbc(const minnsstate &state, const real_1d_array &bndl, const real_1d_array &bndu);
    void minnssetlc(const minnsstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
    void minnssetlc(const minnsstate &state, const real_2d_array &c, const integer_1d_array &ct);
    void minnssetnlc(const minnsstate &state, const ae_int_t nlec, const ae_int_t nlic);
    void minnssetcond(const minnsstate &state, const double epsx, const ae_int_t maxits);
    void minnssetscale(const minnsstate &state, const real_1d_array &s);
    void minnssetalgoags(const minnsstate &state, const double radius, const double penalty);
    void minnssetxrep(const minnsstate &state, const bool needxrep);
    void minnsrequesttermination(const minnsstate &state);
    bool minnsiteration(const minnsstate &state);
    void minnsoptimize(minnsstate &state,
        void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minnsoptimize(minnsstate &state,
        void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
        void  (*rep)(const real_1d_array &x, double func, void *ptr) = null,
        void *ptr = null);
    void minnsresults(const minnsstate &state, real_1d_array &x, minnsreport &rep);
    void minnsresultsbuf(const minnsstate &state, real_1d_array &x, minnsreport &rep);
    void minnsrestartfrom(const minnsstate &state, const real_1d_array &x);
}
*/

extern(C)
{
    void trimprepare(double f, double* threshold, ae_state *_state);
    void trimfunction(double* f, /* Real    */ ae_vector* g, ae_int_t n, double threshold, ae_state *_state);
    ae_bool enforceboundaryconstraints(/* Real    */ ae_vector* x, /* Real    */ ae_vector* bl, /* Boolean */ ae_vector* havebl, /* Real    */ ae_vector* bu, /* Boolean */ ae_vector* havebu, ae_int_t nmain, ae_int_t nslack, ae_state *_state);
    void projectgradientintobc(/* Real    */ ae_vector* x, /* Real    */ ae_vector* g, /* Real    */ ae_vector* bl, /* Boolean */ ae_vector* havebl, /* Real    */ ae_vector* bu, /* Boolean */ ae_vector* havebu, ae_int_t nmain, ae_int_t nslack, ae_state *_state);
    void calculatestepbound(/* Real    */ ae_vector* x, /* Real    */ ae_vector* d, double alpha, /* Real    */ ae_vector* bndl, /* Boolean */ ae_vector* havebndl, /* Real    */ ae_vector* bndu, /* Boolean */ ae_vector* havebndu, ae_int_t nmain, ae_int_t nslack, ae_int_t* variabletofreeze, double* valuetofreeze, double* maxsteplen, ae_state *_state);
    ae_int_t postprocessboundedstep(/* Real    */ ae_vector* x, /* Real    */ ae_vector* xprev, /* Real    */ ae_vector* bndl, /* Boolean */ ae_vector* havebndl, /* Real    */ ae_vector* bndu, /* Boolean */ ae_vector* havebndu, ae_int_t nmain, ae_int_t nslack, ae_int_t variabletofreeze, double valuetofreeze, double steptaken, double maxsteplen, ae_state *_state);
    void filterdirection(/* Real    */ ae_vector* d, /* Real    */ ae_vector* x, /* Real    */ ae_vector* bndl, /* Boolean */ ae_vector* havebndl, /* Real    */ ae_vector* bndu, /* Boolean */ ae_vector* havebndu, /* Real    */ ae_vector* s, ae_int_t nmain, ae_int_t nslack, double droptol, ae_state *_state);
    ae_int_t numberofchangedconstraints(/* Real    */ ae_vector* x, /* Real    */ ae_vector* xprev, /* Real    */ ae_vector* bndl, /* Boolean */ ae_vector* havebndl, /* Real    */ ae_vector* bndu, /* Boolean */ ae_vector* havebndu, ae_int_t nmain, ae_int_t nslack, ae_state *_state);
    ae_bool findfeasiblepoint(/* Real    */ ae_vector* x, /* Real    */ ae_vector* bndl, /* Boolean */ ae_vector* havebndl, /* Real    */ ae_vector* bndu, /* Boolean */ ae_vector* havebndu, ae_int_t nmain, ae_int_t nslack, /* Real    */ ae_matrix* ce, ae_int_t k, double epsi, ae_int_t* qpits, ae_int_t* gpaits, ae_state *_state);
    ae_bool derivativecheck(double f0, double df0, double f1, double df1, double f, double df, double width, ae_state *_state);
    void estimateparabolicmodel(double absasum, double absasum2, double mx, double mb, double md, double d1, double d2, ae_int_t* d1est, ae_int_t* d2est, ae_state *_state);
    void inexactlbfgspreconditioner(/* Real    */ ae_vector* s, ae_int_t n, /* Real    */ ae_vector* d, /* Real    */ ae_vector* c, /* Real    */ ae_matrix* w, ae_int_t k, precbuflbfgs* buf, ae_state *_state);
    void preparelowrankpreconditioner(/* Real    */ ae_vector* d, /* Real    */ ae_vector* c, /* Real    */ ae_matrix* w, ae_int_t n, ae_int_t k, precbuflowrank* buf, ae_state *_state);
    void applylowrankpreconditioner(/* Real    */ ae_vector* s, precbuflowrank* buf, ae_state *_state);
    void _precbuflbfgs_init(void* _p, ae_state *_state);
    void _precbuflbfgs_init_copy(void* _dst, void* _src, ae_state *_state);
    void _precbuflbfgs_clear(void* _p);
    void _precbuflbfgs_destroy(void* _p);
    void _precbuflowrank_init(void* _p, ae_state *_state);
    void _precbuflowrank_init_copy(void* _dst, void* _src, ae_state *_state);
    void _precbuflowrank_clear(void* _p);
    void _precbuflowrank_destroy(void* _p);
    void cqminit(ae_int_t n, convexquadraticmodel* s, ae_state *_state);
    void cqmseta(convexquadraticmodel* s, /* Real    */ ae_matrix* a, ae_bool isupper, double alpha, ae_state *_state);
    void cqmgeta(convexquadraticmodel* s, /* Real    */ ae_matrix* a, ae_state *_state);
    void cqmrewritedensediagonal(convexquadraticmodel* s, /* Real    */ ae_vector* z, ae_state *_state);
    void cqmsetd(convexquadraticmodel* s, /* Real    */ ae_vector* d, double tau, ae_state *_state);
    void cqmdropa(convexquadraticmodel* s, ae_state *_state);
    void cqmsetb(convexquadraticmodel* s, /* Real    */ ae_vector* b, ae_state *_state);
    void cqmsetq(convexquadraticmodel* s, /* Real    */ ae_matrix* q, /* Real    */ ae_vector* r, ae_int_t k, double theta, ae_state *_state);
    void cqmsetactiveset(convexquadraticmodel* s, /* Real    */ ae_vector* x, /* Boolean */ ae_vector* activeset, ae_state *_state);
    double cqmeval(convexquadraticmodel* s, /* Real    */ ae_vector* x, ae_state *_state);
    void cqmevalx(convexquadraticmodel* s, /* Real    */ ae_vector* x, double* r, double* noise, ae_state *_state);
    void cqmgradunconstrained(convexquadraticmodel* s, /* Real    */ ae_vector* x, /* Real    */ ae_vector* g, ae_state *_state);
    double cqmxtadx2(convexquadraticmodel* s, /* Real    */ ae_vector* x, ae_state *_state);
    void cqmadx(convexquadraticmodel* s, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    ae_bool cqmconstrainedoptimum(convexquadraticmodel* s, /* Real    */ ae_vector* x, ae_state *_state);
    void cqmscalevector(convexquadraticmodel* s, /* Real    */ ae_vector* x, ae_state *_state);
    double cqmdebugconstrainedevalt(convexquadraticmodel* s, /* Real    */ ae_vector* x, ae_state *_state);
    double cqmdebugconstrainedevale(convexquadraticmodel* s, /* Real    */ ae_vector* x, ae_state *_state);
    void _convexquadraticmodel_init(void* _p, ae_state *_state);
    void _convexquadraticmodel_init_copy(void* _dst, void* _src, ae_state *_state);
    void _convexquadraticmodel_clear(void* _p);
    void _convexquadraticmodel_destroy(void* _p);
    void snnlsinit(ae_int_t nsmax, ae_int_t ndmax, ae_int_t nrmax, snnlssolver* s, ae_state *_state);
    void snnlssetproblem(snnlssolver* s, /* Real    */ ae_matrix* a, /* Real    */ ae_vector* b, ae_int_t ns, ae_int_t nd, ae_int_t nr, ae_state *_state);
    void snnlsdropnnc(snnlssolver* s, ae_int_t idx, ae_state *_state);
    void snnlssolve(snnlssolver* s, /* Real    */ ae_vector* x, ae_state *_state);
    void _snnlssolver_init(void* _p, ae_state *_state);
    void _snnlssolver_init_copy(void* _dst, void* _src, ae_state *_state);
    void _snnlssolver_clear(void* _p);
    void _snnlssolver_destroy(void* _p);
    void sasinit(ae_int_t n, sactiveset* s, ae_state *_state);
    void sassetscale(sactiveset* state, /* Real    */ ae_vector* s, ae_state *_state);
    void sassetprecdiag(sactiveset* state, /* Real    */ ae_vector* d, ae_state *_state);
    void sassetbc(sactiveset* state, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    void sassetlc(sactiveset* state, /* Real    */ ae_matrix* c, /* Integer */ ae_vector* ct, ae_int_t k, ae_state *_state);
    void sassetlcx(sactiveset* state, /* Real    */ ae_matrix* cleic, ae_int_t nec, ae_int_t nic, ae_state *_state);
    ae_bool sasstartoptimization(sactiveset* state, /* Real    */ ae_vector* x, ae_state *_state);
    void sasexploredirection(sactiveset* state, /* Real    */ ae_vector* d, double* stpmax, ae_int_t* cidx, double* vval, ae_state *_state);
    ae_int_t sasmoveto(sactiveset* state, /* Real    */ ae_vector* xn, ae_bool needact, ae_int_t cidx, double cval, ae_state *_state);
    void sasimmediateactivation(sactiveset* state, ae_int_t cidx, double cval, ae_state *_state);
    void sasconstraineddescent(sactiveset* state, /* Real    */ ae_vector* g, /* Real    */ ae_vector* d, ae_state *_state);
    void sasconstraineddescentprec(sactiveset* state, /* Real    */ ae_vector* g, /* Real    */ ae_vector* d, ae_state *_state);
    void sasconstraineddirection(sactiveset* state, /* Real    */ ae_vector* d, ae_state *_state);
    void sasconstraineddirectionprec(sactiveset* state, /* Real    */ ae_vector* d, ae_state *_state);
    void sascorrection(sactiveset* state, /* Real    */ ae_vector* x, double* penalty, ae_state *_state);
    double sasactivelcpenalty1(sactiveset* state, /* Real    */ ae_vector* x, ae_state *_state);
    double sasscaledconstrainednorm(sactiveset* state, /* Real    */ ae_vector* d, ae_state *_state);
    void sasstopoptimization(sactiveset* state, ae_state *_state);
    void sasreactivateconstraints(sactiveset* state, /* Real    */ ae_vector* gc, ae_state *_state);
    void sasreactivateconstraintsprec(sactiveset* state, /* Real    */ ae_vector* gc, ae_state *_state);
    void sasrebuildbasis(sactiveset* state, ae_state *_state);
    void _sactiveset_init(void* _p, ae_state *_state);
    void _sactiveset_init_copy(void* _dst, void* _src, ae_state *_state);
    void _sactiveset_clear(void* _p);
    void _sactiveset_destroy(void* _p);
    void mincgcreate(ae_int_t n, /* Real    */ ae_vector* x, mincgstate* state, ae_state *_state);
    void mincgcreatef(ae_int_t n, /* Real    */ ae_vector* x, double diffstep, mincgstate* state, ae_state *_state);
    void mincgsetcond(mincgstate* state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void mincgsetscale(mincgstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void mincgsetxrep(mincgstate* state, ae_bool needxrep, ae_state *_state);
    void mincgsetdrep(mincgstate* state, ae_bool needdrep, ae_state *_state);
    void mincgsetcgtype(mincgstate* state, ae_int_t cgtype, ae_state *_state);
    void mincgsetstpmax(mincgstate* state, double stpmax, ae_state *_state);
    void mincgsuggeststep(mincgstate* state, double stp, ae_state *_state);
    double mincglastgoodstep(mincgstate* state, ae_state *_state);
    void mincgsetprecdefault(mincgstate* state, ae_state *_state);
    void mincgsetprecdiag(mincgstate* state, /* Real    */ ae_vector* d, ae_state *_state);
    void mincgsetprecscale(mincgstate* state, ae_state *_state);
    ae_bool mincgiteration(mincgstate* state, ae_state *_state);
    void mincgresults(mincgstate* state, /* Real    */ ae_vector* x, mincgreport* rep, ae_state *_state);
    void mincgresultsbuf(mincgstate* state, /* Real    */ ae_vector* x, mincgreport* rep, ae_state *_state);
    void mincgrestartfrom(mincgstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void mincgrequesttermination(mincgstate* state, ae_state *_state);
    void mincgsetprecdiagfast(mincgstate* state, /* Real    */ ae_vector* d, ae_state *_state);
    void mincgsetpreclowrankfast(mincgstate* state, /* Real    */ ae_vector* d1, /* Real    */ ae_vector* c, /* Real    */ ae_matrix* v, ae_int_t vcnt, ae_state *_state);
    void mincgsetprecvarpart(mincgstate* state, /* Real    */ ae_vector* d2, ae_state *_state);
    void mincgsetgradientcheck(mincgstate* state, double teststep, ae_state *_state);
    void _mincgstate_init(void* _p, ae_state *_state);
    void _mincgstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mincgstate_clear(void* _p);
    void _mincgstate_destroy(void* _p);
    void _mincgreport_init(void* _p, ae_state *_state);
    void _mincgreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mincgreport_clear(void* _p);
    void _mincgreport_destroy(void* _p);
    void minbleiccreate(ae_int_t n, /* Real    */ ae_vector* x, minbleicstate* state, ae_state *_state);
    void minbleiccreatef(ae_int_t n, /* Real    */ ae_vector* x, double diffstep, minbleicstate* state, ae_state *_state);
    void minbleicsetbc(minbleicstate* state, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    void minbleicsetlc(minbleicstate* state, /* Real    */ ae_matrix* c, /* Integer */ ae_vector* ct, ae_int_t k, ae_state *_state);
    void minbleicsetcond(minbleicstate* state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void minbleicsetscale(minbleicstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void minbleicsetprecdefault(minbleicstate* state, ae_state *_state);
    void minbleicsetprecdiag(minbleicstate* state, /* Real    */ ae_vector* d, ae_state *_state);
    void minbleicsetprecscale(minbleicstate* state, ae_state *_state);
    void minbleicsetxrep(minbleicstate* state, ae_bool needxrep, ae_state *_state);
    void minbleicsetdrep(minbleicstate* state, ae_bool needdrep, ae_state *_state);
    void minbleicsetstpmax(minbleicstate* state, double stpmax, ae_state *_state);
    ae_bool minbleiciteration(minbleicstate* state, ae_state *_state);
    void minbleicresults(minbleicstate* state, /* Real    */ ae_vector* x, minbleicreport* rep, ae_state *_state);
    void minbleicresultsbuf(minbleicstate* state, /* Real    */ ae_vector* x, minbleicreport* rep, ae_state *_state);
    void minbleicrestartfrom(minbleicstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void minbleicrequesttermination(minbleicstate* state, ae_state *_state);
    void minbleicemergencytermination(minbleicstate* state, ae_state *_state);
    void minbleicsetgradientcheck(minbleicstate* state, double teststep, ae_state *_state);
    void _minbleicstate_init(void* _p, ae_state *_state);
    void _minbleicstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minbleicstate_clear(void* _p);
    void _minbleicstate_destroy(void* _p);
    void _minbleicreport_init(void* _p, ae_state *_state);
    void _minbleicreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minbleicreport_clear(void* _p);
    void _minbleicreport_destroy(void* _p);
    void minlbfgscreate(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, minlbfgsstate* state, ae_state *_state);
    void minlbfgscreatef(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, double diffstep, minlbfgsstate* state, ae_state *_state);
    void minlbfgssetcond(minlbfgsstate* state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void minlbfgssetxrep(minlbfgsstate* state, ae_bool needxrep, ae_state *_state);
    void minlbfgssetstpmax(minlbfgsstate* state, double stpmax, ae_state *_state);
    void minlbfgssetscale(minlbfgsstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void minlbfgscreatex(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, ae_int_t flags, double diffstep, minlbfgsstate* state, ae_state *_state);
    void minlbfgssetprecdefault(minlbfgsstate* state, ae_state *_state);
    void minlbfgssetpreccholesky(minlbfgsstate* state, /* Real    */ ae_matrix* p, ae_bool isupper, ae_state *_state);
    void minlbfgssetprecdiag(minlbfgsstate* state, /* Real    */ ae_vector* d, ae_state *_state);
    void minlbfgssetprecscale(minlbfgsstate* state, ae_state *_state);
    void minlbfgssetprecrankklbfgsfast(minlbfgsstate* state, /* Real    */ ae_vector* d, /* Real    */ ae_vector* c, /* Real    */ ae_matrix* w, ae_int_t cnt, ae_state *_state);
    void minlbfgssetpreclowrankexact(minlbfgsstate* state, /* Real    */ ae_vector* d, /* Real    */ ae_vector* c, /* Real    */ ae_matrix* w, ae_int_t cnt, ae_state *_state);
    ae_bool minlbfgsiteration(minlbfgsstate* state, ae_state *_state);
    void minlbfgsresults(minlbfgsstate* state, /* Real    */ ae_vector* x, minlbfgsreport* rep, ae_state *_state);
    void minlbfgsresultsbuf(minlbfgsstate* state, /* Real    */ ae_vector* x, minlbfgsreport* rep, ae_state *_state);
    void minlbfgsrestartfrom(minlbfgsstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void minlbfgsrequesttermination(minlbfgsstate* state, ae_state *_state);
    void minlbfgssetgradientcheck(minlbfgsstate* state, double teststep, ae_state *_state);
    void _minlbfgsstate_init(void* _p, ae_state *_state);
    void _minlbfgsstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minlbfgsstate_clear(void* _p);
    void _minlbfgsstate_destroy(void* _p);
    void _minlbfgsreport_init(void* _p, ae_state *_state);
    void _minlbfgsreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minlbfgsreport_clear(void* _p);
    void _minlbfgsreport_destroy(void* _p);
    void qqploaddefaults(ae_int_t nmain, qqpsettings* s, ae_state *_state);
    void qqpcopysettings(qqpsettings* src, qqpsettings* dst, ae_state *_state);
    void qqpoptimize(convexquadraticmodel* ac, sparsematrix* sparseac, ae_int_t akind, ae_bool sparseupper, /* Real    */ ae_vector* bc, /* Real    */ ae_vector* bndlc, /* Real    */ ae_vector* bnduc, /* Real    */ ae_vector* sc, /* Real    */ ae_vector* xoriginc, ae_int_t nc, /* Real    */ ae_matrix* cleicc, ae_int_t nec, ae_int_t nic, qqpsettings* settings, qqpbuffers* sstate, /* Real    */ ae_vector* xs, ae_int_t* terminationtype, ae_state *_state);
    void _qqpsettings_init(void* _p, ae_state *_state);
    void _qqpsettings_init_copy(void* _dst, void* _src, ae_state *_state);
    void _qqpsettings_clear(void* _p);
    void _qqpsettings_destroy(void* _p);
    void _qqpbuffers_init(void* _p, ae_state *_state);
    void _qqpbuffers_init_copy(void* _dst, void* _src, ae_state *_state);
    void _qqpbuffers_clear(void* _p);
    void _qqpbuffers_destroy(void* _p);
    void qpbleicloaddefaults(ae_int_t nmain, qpbleicsettings* s, ae_state *_state);
    void qpbleiccopysettings(qpbleicsettings* src, qpbleicsettings* dst, ae_state *_state);
    void qpbleicoptimize(convexquadraticmodel* a, sparsematrix* sparsea, ae_int_t akind, ae_bool sparseaupper, double absasum, double absasum2, /* Real    */ ae_vector* b, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, /* Real    */ ae_vector* s, /* Real    */ ae_vector* xorigin, ae_int_t n, /* Real    */ ae_matrix* cleic, ae_int_t nec, ae_int_t nic, qpbleicsettings* settings, qpbleicbuffers* sstate, ae_bool* firstcall, /* Real    */ ae_vector* xs, ae_int_t* terminationtype, ae_state *_state);
    void _qpbleicsettings_init(void* _p, ae_state *_state);
    void _qpbleicsettings_init_copy(void* _dst, void* _src, ae_state *_state);
    void _qpbleicsettings_clear(void* _p);
    void _qpbleicsettings_destroy(void* _p);
    void _qpbleicbuffers_init(void* _p, ae_state *_state);
    void _qpbleicbuffers_init_copy(void* _dst, void* _src, ae_state *_state);
    void _qpbleicbuffers_clear(void* _p);
    void _qpbleicbuffers_destroy(void* _p);
    void qpcholeskyloaddefaults(ae_int_t nmain, qpcholeskysettings* s, ae_state *_state);
    void qpcholeskycopysettings(qpcholeskysettings* src, qpcholeskysettings* dst, ae_state *_state);
    void qpcholeskyoptimize(convexquadraticmodel* a, double anorm, /* Real    */ ae_vector* b, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, /* Real    */ ae_vector* s, /* Real    */ ae_vector* xorigin, ae_int_t n, /* Real    */ ae_matrix* cleic, ae_int_t nec, ae_int_t nic, qpcholeskybuffers* sstate, /* Real    */ ae_vector* xsc, ae_int_t* terminationtype, ae_state *_state);
    void _qpcholeskysettings_init(void* _p, ae_state *_state);
    void _qpcholeskysettings_init_copy(void* _dst, void* _src, ae_state *_state);
    void _qpcholeskysettings_clear(void* _p);
    void _qpcholeskysettings_destroy(void* _p);
    void _qpcholeskybuffers_init(void* _p, ae_state *_state);
    void _qpcholeskybuffers_init_copy(void* _dst, void* _src, ae_state *_state);
    void _qpcholeskybuffers_clear(void* _p);
    void _qpcholeskybuffers_destroy(void* _p);
    void minqpcreate(ae_int_t n, minqpstate* state, ae_state *_state);
    void minqpsetlinearterm(minqpstate* state, /* Real    */ ae_vector* b, ae_state *_state);
    void minqpsetquadraticterm(minqpstate* state, /* Real    */ ae_matrix* a, ae_bool isupper, ae_state *_state);
    void minqpsetquadratictermsparse(minqpstate* state, sparsematrix* a, ae_bool isupper, ae_state *_state);
    void minqpsetstartingpoint(minqpstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void minqpsetorigin(minqpstate* state, /* Real    */ ae_vector* xorigin, ae_state *_state);
    void minqpsetscale(minqpstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void minqpsetalgocholesky(minqpstate* state, ae_state *_state);
    void minqpsetalgobleic(minqpstate* state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void minqpsetalgoquickqp(minqpstate* state, double epsg, double epsf, double epsx, ae_int_t maxouterits, ae_bool usenewton, ae_state *_state);
    void minqpsetbc(minqpstate* state, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    void minqpsetlc(minqpstate* state, /* Real    */ ae_matrix* c, /* Integer */ ae_vector* ct, ae_int_t k, ae_state *_state);
    void minqpoptimize(minqpstate* state, ae_state *_state);
    void minqpresults(minqpstate* state, /* Real    */ ae_vector* x, minqpreport* rep, ae_state *_state);
    void minqpresultsbuf(minqpstate* state, /* Real    */ ae_vector* x, minqpreport* rep, ae_state *_state);
    void minqpsetlineartermfast(minqpstate* state, /* Real    */ ae_vector* b, ae_state *_state);
    void minqpsetquadratictermfast(minqpstate* state, /* Real    */ ae_matrix* a, ae_bool isupper, double s, ae_state *_state);
    void minqprewritediagonal(minqpstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void minqpsetstartingpointfast(minqpstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void minqpsetoriginfast(minqpstate* state, /* Real    */ ae_vector* xorigin, ae_state *_state);
    void _minqpstate_init(void* _p, ae_state *_state);
    void _minqpstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minqpstate_clear(void* _p);
    void _minqpstate_destroy(void* _p);
    void _minqpreport_init(void* _p, ae_state *_state);
    void _minqpreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minqpreport_clear(void* _p);
    void _minqpreport_destroy(void* _p);
    void minlmcreatevj(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, minlmstate* state, ae_state *_state);
    void minlmcreatev(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, double diffstep, minlmstate* state, ae_state *_state);
    void minlmcreatefgh(ae_int_t n, /* Real    */ ae_vector* x, minlmstate* state, ae_state *_state);
    void minlmsetcond(minlmstate* state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void minlmsetxrep(minlmstate* state, ae_bool needxrep, ae_state *_state);
    void minlmsetstpmax(minlmstate* state, double stpmax, ae_state *_state);
    void minlmsetscale(minlmstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void minlmsetbc(minlmstate* state, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    void minlmsetacctype(minlmstate* state, ae_int_t acctype, ae_state *_state);
    ae_bool minlmiteration(minlmstate* state, ae_state *_state);
    void minlmresults(minlmstate* state, /* Real    */ ae_vector* x, minlmreport* rep, ae_state *_state);
    void minlmresultsbuf(minlmstate* state, /* Real    */ ae_vector* x, minlmreport* rep, ae_state *_state);
    void minlmrestartfrom(minlmstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void minlmrequesttermination(minlmstate* state, ae_state *_state);
    void minlmcreatevgj(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, minlmstate* state, ae_state *_state);
    void minlmcreatefgj(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, minlmstate* state, ae_state *_state);
    void minlmcreatefj(ae_int_t n, ae_int_t m, /* Real    */ ae_vector* x, minlmstate* state, ae_state *_state);
    void minlmsetgradientcheck(minlmstate* state, double teststep, ae_state *_state);
    void _minlmstate_init(void* _p, ae_state *_state);
    void _minlmstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minlmstate_clear(void* _p);
    void _minlmstate_destroy(void* _p);
    void _minlmreport_init(void* _p, ae_state *_state);
    void _minlmreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minlmreport_clear(void* _p);
    void _minlmreport_destroy(void* _p);
    void minlbfgssetdefaultpreconditioner(minlbfgsstate* state, ae_state *_state);
    void minlbfgssetcholeskypreconditioner(minlbfgsstate* state, /* Real    */ ae_matrix* p, ae_bool isupper, ae_state *_state);
    void minbleicsetbarrierwidth(minbleicstate* state, double mu, ae_state *_state);
    void minbleicsetbarrierdecay(minbleicstate* state, double mudecay, ae_state *_state);
    void minasacreate(ae_int_t n, /* Real    */ ae_vector* x, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, minasastate* state, ae_state *_state);
    void minasasetcond(minasastate* state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void minasasetxrep(minasastate* state, ae_bool needxrep, ae_state *_state);
    void minasasetalgorithm(minasastate* state, ae_int_t algotype, ae_state *_state);
    void minasasetstpmax(minasastate* state, double stpmax, ae_state *_state);
    ae_bool minasaiteration(minasastate* state, ae_state *_state);
    void minasaresults(minasastate* state, /* Real    */ ae_vector* x, minasareport* rep, ae_state *_state);
    void minasaresultsbuf(minasastate* state, /* Real    */ ae_vector* x, minasareport* rep, ae_state *_state);
    void minasarestartfrom(minasastate* state, /* Real    */ ae_vector* x, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    void _minasastate_init(void* _p, ae_state *_state);
    void _minasastate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minasastate_clear(void* _p);
    void _minasastate_destroy(void* _p);
    void _minasareport_init(void* _p, ae_state *_state);
    void _minasareport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minasareport_clear(void* _p);
    void _minasareport_destroy(void* _p);
    void minnlccreate(ae_int_t n, /* Real    */ ae_vector* x, minnlcstate* state, ae_state *_state);
    void minnlccreatef(ae_int_t n, /* Real    */ ae_vector* x, double diffstep, minnlcstate* state, ae_state *_state);
    void minnlcsetbc(minnlcstate* state, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    void minnlcsetlc(minnlcstate* state, /* Real    */ ae_matrix* c, /* Integer */ ae_vector* ct, ae_int_t k, ae_state *_state);
    void minnlcsetnlc(minnlcstate* state, ae_int_t nlec, ae_int_t nlic, ae_state *_state);
    void minnlcsetcond(minnlcstate* state, double epsg, double epsf, double epsx, ae_int_t maxits, ae_state *_state);
    void minnlcsetscale(minnlcstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void minnlcsetprecinexact(minnlcstate* state, ae_state *_state);
    void minnlcsetprecexactlowrank(minnlcstate* state, ae_int_t updatefreq, ae_state *_state);
    void minnlcsetprecnone(minnlcstate* state, ae_state *_state);
    void minnlcsetalgoaul(minnlcstate* state, double rho, ae_int_t itscnt, ae_state *_state);
    void minnlcsetxrep(minnlcstate* state, ae_bool needxrep, ae_state *_state);
    ae_bool minnlciteration(minnlcstate* state, ae_state *_state);
    void minnlcresults(minnlcstate* state, /* Real    */ ae_vector* x, minnlcreport* rep, ae_state *_state);
    void minnlcresultsbuf(minnlcstate* state, /* Real    */ ae_vector* x, minnlcreport* rep, ae_state *_state);
    void minnlcrestartfrom(minnlcstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void minnlcsetgradientcheck(minnlcstate* state, double teststep, ae_state *_state);
    void minnlcequalitypenaltyfunction(double alpha, double* f, double* df, double* d2f, ae_state *_state);
    void minnlcinequalitypenaltyfunction(double alpha, double stabilizingpoint, double* f, double* df, double* d2f, ae_state *_state);
    void minnlcinequalityshiftfunction(double alpha, double* f, double* df, double* d2f, ae_state *_state);
    void _minnlcstate_init(void* _p, ae_state *_state);
    void _minnlcstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minnlcstate_clear(void* _p);
    void _minnlcstate_destroy(void* _p);
    void _minnlcreport_init(void* _p, ae_state *_state);
    void _minnlcreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minnlcreport_clear(void* _p);
    void _minnlcreport_destroy(void* _p);
    void minnscreate(ae_int_t n, /* Real    */ ae_vector* x, minnsstate* state, ae_state *_state);
    void minnscreatef(ae_int_t n, /* Real    */ ae_vector* x, double diffstep, minnsstate* state, ae_state *_state);
    void minnssetbc(minnsstate* state, /* Real    */ ae_vector* bndl, /* Real    */ ae_vector* bndu, ae_state *_state);
    void minnssetlc(minnsstate* state, /* Real    */ ae_matrix* c, /* Integer */ ae_vector* ct, ae_int_t k, ae_state *_state);
    void minnssetnlc(minnsstate* state, ae_int_t nlec, ae_int_t nlic, ae_state *_state);
    void minnssetcond(minnsstate* state, double epsx, ae_int_t maxits, ae_state *_state);
    void minnssetscale(minnsstate* state, /* Real    */ ae_vector* s, ae_state *_state);
    void minnssetalgoags(minnsstate* state, double radius, double penalty, ae_state *_state);
    void minnssetxrep(minnsstate* state, ae_bool needxrep, ae_state *_state);
    void minnsrequesttermination(minnsstate* state, ae_state *_state);
    ae_bool minnsiteration(minnsstate* state, ae_state *_state);
    void minnsresults(minnsstate* state, /* Real    */ ae_vector* x, minnsreport* rep, ae_state *_state);
    void minnsresultsbuf(minnsstate* state, /* Real    */ ae_vector* x, minnsreport* rep, ae_state *_state);
    void minnsrestartfrom(minnsstate* state, /* Real    */ ae_vector* x, ae_state *_state);
    void _minnsqp_init(void* _p, ae_state *_state);
    void _minnsqp_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minnsqp_clear(void* _p);
    void _minnsqp_destroy(void* _p);
    void _minnsstate_init(void* _p, ae_state *_state);
    void _minnsstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minnsstate_clear(void* _p);
    void _minnsstate_destroy(void* _p);
    void _minnsreport_init(void* _p, ae_state *_state);
    void _minnsreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _minnsreport_clear(void* _p);
    void _minnsreport_destroy(void* _p);
}

