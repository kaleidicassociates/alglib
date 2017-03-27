module alglib.dataanalysis;
import alglib.ap;
import alglib.internal;
import alglib.linalg;
import alglib.statistics;
import alglib.misc;
import alglib.specialfunctions;
import alglib.solvers;
import alglib.optimization;

struct cvreport
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
}

struct kmeansbuffers
{
    ae_matrix ct;
    ae_matrix ctbest;
    ae_vector xycbest;
    ae_vector xycprev;
    ae_vector d2;
    ae_vector csizes;
    apbuffers initbuf;
    ae_shared_pool updatepool;
}

struct clusterizerstate
{
    ae_int_t npoints;
    ae_int_t nfeatures;
    ae_int_t disttype;
    ae_matrix xy;
    ae_matrix d;
    ae_int_t ahcalgo;
    ae_int_t kmeansrestarts;
    ae_int_t kmeansmaxits;
    ae_int_t kmeansinitalgo;
    ae_bool kmeansdbgnoits;
    ae_matrix tmpd;
    apbuffers distbuf;
    kmeansbuffers kmeanstmp;
}

struct ahcreport
{
    ae_int_t terminationtype;
    ae_int_t npoints;
    ae_vector p;
    ae_matrix z;
    ae_matrix pz;
    ae_matrix pm;
    ae_vector mergedist;
}

struct kmeansreport
{
    ae_int_t npoints;
    ae_int_t nfeatures;
    ae_int_t terminationtype;
    ae_int_t iterationscount;
    double energy;
    ae_int_t k;
    ae_matrix c;
    ae_vector cidx;
}

struct decisionforest
{
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t ntrees;
    ae_int_t bufsize;
    ae_vector trees;
}

struct dfreport
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double oobrelclserror;
    double oobavgce;
    double oobrmserror;
    double oobavgerror;
    double oobavgrelerror;
} 

struct dfinternalbuffers
{
    ae_vector treebuf;
    ae_vector idxbuf;
    ae_vector tmpbufr;
    ae_vector tmpbufr2;
    ae_vector tmpbufi;
    ae_vector classibuf;
    ae_vector sortrbuf;
    ae_vector sortrbuf2;
    ae_vector sortibuf;
    ae_vector varpool;
    ae_vector evsbin;
    ae_vector evssplits;
}

struct linearmodel
{
    ae_vector w;
} 

struct lrreport
{
    ae_matrix c;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double cvrmserror;
    double cvavgerror;
    double cvavgrelerror;
    ae_int_t ncvdefects;
    ae_vector cvdefects;
}

struct modelerrors
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
}

struct smlpgrad
{
    double f;
    ae_vector g;
}

struct multilayerperceptron
{
    ae_int_t hlnetworktype;
    ae_int_t hlnormtype;
    ae_vector hllayersizes;
    ae_vector hlconnections;
    ae_vector hlneurons;
    ae_vector structinfo;
    ae_vector weights;
    ae_vector columnmeans;
    ae_vector columnsigmas;
    ae_vector neurons;
    ae_vector dfdnet;
    ae_vector derror;
    ae_vector x;
    ae_vector y;
    ae_matrix xy;
    ae_vector xyrow;
    ae_vector nwbuf;
    ae_vector integerbuf;
    modelerrors err;
    ae_vector rndbuf;
    ae_shared_pool buf;
    ae_shared_pool gradbuf;
    ae_matrix dummydxy;
    sparsematrix dummysxy;
    ae_vector dummyidx;
    ae_shared_pool dummypool;
}

struct logitmodel
{
    ae_vector w;
}

struct logitmcstate
{
    ae_bool brackt;
    ae_bool stage1;
    ae_int_t infoc;
    double dg;
    double dgm;
    double dginit;
    double dgtest;
    double dgx;
    double dgxm;
    double dgy;
    double dgym;
    double finit;
    double ftest1;
    double fm;
    double fx;
    double fxm;
    double fy;
    double fym;
    double stx;
    double sty;
    double stmin;
    double stmax;
    double width;
    double width1;
    double xtrapf;
}

struct mnlreport
{
    ae_int_t ngrad;
    ae_int_t nhess;
}

struct mcpdstate
{
    ae_int_t n;
    ae_vector states;
    ae_int_t npairs;
    ae_matrix data;
    ae_matrix ec;
    ae_matrix bndl;
    ae_matrix bndu;
    ae_matrix c;
    ae_vector ct;
    ae_int_t ccnt;
    ae_vector pw;
    ae_matrix priorp;
    double regterm;
    minbleicstate bs;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repnfev;
    ae_int_t repterminationtype;
    minbleicreport br;
    ae_vector tmpp;
    ae_vector effectivew;
    ae_vector effectivebndl;
    ae_vector effectivebndu;
    ae_matrix effectivec;
    ae_vector effectivect;
    ae_vector h;
    ae_matrix p;
}

struct mcpdreport
{
    ae_int_t inneriterationscount;
    ae_int_t outeriterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
}

struct mlpensemble
{
    ae_int_t ensemblesize;
    ae_vector weights;
    ae_vector columnmeans;
    ae_vector columnsigmas;
    multilayerperceptron network;
    ae_vector y;
}

struct mlpreport
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
    ae_int_t ngrad;
    ae_int_t nhess;
    ae_int_t ncholesky;
}

struct mlpcvreport
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
}

struct smlptrnsession
{
    ae_vector bestparameters;
    double bestrmserror;
    ae_bool randomizenetwork;
    multilayerperceptron network;
    minlbfgsstate optimizer;
    minlbfgsreport optimizerrep;
    ae_vector wbuf0;
    ae_vector wbuf1;
    ae_vector allminibatches;
    ae_vector currentminibatch;
    rcommstate rstate;
    ae_int_t algoused;
    ae_int_t minibatchsize;
    hqrndstate generator;
}

struct mlpetrnsession
{
    ae_vector trnsubset;
    ae_vector valsubset;
    ae_shared_pool mlpsessions;
    mlpreport mlprep;
    multilayerperceptron network;
}

struct mlptrainer
{
    ae_int_t nin;
    ae_int_t nout;
    ae_bool rcpar;
    ae_int_t lbfgsfactor;
    double decay;
    double wstep;
    ae_int_t maxits;
    ae_int_t datatype;
    ae_int_t npoints;
    ae_matrix densexy;
    sparsematrix sparsexy;
    smlptrnsession session;
    ae_int_t ngradbatch;
    ae_vector subset;
    ae_int_t subsetsize;
    ae_vector valsubset;
    ae_int_t valsubsetsize;
    ae_int_t algokind;
    ae_int_t minibatchsize;
}

struct mlpparallelizationcv
{
    multilayerperceptron network;
    mlpreport rep;
    ae_vector subset;
    ae_int_t subsetsize;
    ae_vector xyrow;
    ae_vector y;
    ae_int_t ngrad;
    ae_shared_pool trnpool;
}

/**
extern(C++,alglib)
{

    class _clusterizerstate_owner
    {
        _clusterizerstate_owner();
        _clusterizerstate_owner(const _clusterizerstate_owner &rhs);
        _clusterizerstate_owner& operator=(const _clusterizerstate_owner &rhs);
        virtual ~_clusterizerstate_owner();
        alglib_impl::clusterizerstate* c_ptr();
        alglib_impl::clusterizerstate* c_ptr() const;
        protected:
            alglib_impl::clusterizerstate *p_struct;
    }

    class clusterizerstate : _clusterizerstate_owner
    {
        clusterizerstate();
        clusterizerstate(const clusterizerstate &rhs);
        clusterizerstate& operator=(const clusterizerstate &rhs);
        virtual ~clusterizerstate();
    }
    class _ahcreport_owner
    {
        _ahcreport_owner();
        _ahcreport_owner(const _ahcreport_owner &rhs);
        _ahcreport_owner& operator=(const _ahcreport_owner &rhs);
        virtual ~_ahcreport_owner();
        alglib_impl::ahcreport* c_ptr();
        alglib_impl::ahcreport* c_ptr() const;
        protected:
            alglib_impl::ahcreport *p_struct;
    }

    class ahcreport : public _ahcreport_owner
    {
        ahcreport();
        ahcreport(const ahcreport &rhs);
        ahcreport& operator=(const ahcreport &rhs);
        virtual ~ahcreport();
        ae_int_t &terminationtype;
        ae_int_t &npoints;
        integer_1d_array p;
        integer_2d_array z;
        integer_2d_array pz;
        integer_2d_array pm;
        real_1d_array mergedist;
    }


    class _kmeansreport_owner
    {
        _kmeansreport_owner();
        _kmeansreport_owner(const _kmeansreport_owner &rhs);
        _kmeansreport_owner& operator=(const _kmeansreport_owner &rhs);
        virtual ~_kmeansreport_owner();
        alglib_impl::kmeansreport* c_ptr();
        alglib_impl::kmeansreport* c_ptr() const;
        protected:
            alglib_impl::kmeansreport *p_struct;
    }

    class kmeansreport : _kmeansreport_owner
    {
        kmeansreport();
        kmeansreport(const kmeansreport &rhs);
        kmeansreport& operator=(const kmeansreport &rhs);
        virtual ~kmeansreport();
        ae_int_t &npoints;
        ae_int_t &nfeatures;
        ae_int_t &terminationtype;
        ae_int_t &iterationscount;
        double &energy;
        ae_int_t &k;
        real_2d_array c;
        integer_1d_array cidx;
    }

    class _decisionforest_owner
    {
        _decisionforest_owner();
        _decisionforest_owner(const _decisionforest_owner &rhs);
        _decisionforest_owner& operator=(const _decisionforest_owner &rhs);
        virtual ~_decisionforest_owner();
        alglib_impl::decisionforest* c_ptr();
        alglib_impl::decisionforest* c_ptr() const;
        protected:
            alglib_impl::decisionforest *p_struct;
    }

    class decisionforest : _decisionforest_owner
    {
        decisionforest();
        decisionforest(const decisionforest &rhs);
        decisionforest& operator=(const decisionforest &rhs);
        virtual ~decisionforest();
    }

    class _dfreport_owner
    {
        _dfreport_owner();
        _dfreport_owner(const _dfreport_owner &rhs);
        _dfreport_owner& operator=(const _dfreport_owner &rhs);
        virtual ~_dfreport_owner();
        alglib_impl::dfreport* c_ptr();
        alglib_impl::dfreport* c_ptr() const;
        protected:
            alglib_impl::dfreport *p_struct;
    }

    class dfreport :  _dfreport_owner
    {
        dfreport();
        dfreport(const dfreport &rhs);
        dfreport& operator=(const dfreport &rhs);
        virtual ~dfreport();
        double &relclserror;
        double &avgce;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
        double &oobrelclserror;
        double &oobavgce;
        double &oobrmserror;
        double &oobavgerror;
        double &oobavgrelerror;
    }

    class _linearmodel_owner
    {
        _linearmodel_owner();
        _linearmodel_owner(const _linearmodel_owner &rhs);
        _linearmodel_owner& operator=(const _linearmodel_owner &rhs);
        virtual ~_linearmodel_owner();
        alglib_impl::linearmodel* c_ptr();
        alglib_impl::linearmodel* c_ptr() const;
        protected:
            alglib_impl::linearmodel *p_struct;
    }

    class linearmodel :  _linearmodel_owner
    {
        linearmodel();
        linearmodel(const linearmodel &rhs);
        linearmodel& operator=(const linearmodel &rhs);
        virtual ~linearmodel();
    }

    class _lrreport_owner
    {
        _lrreport_owner();
        _lrreport_owner(const _lrreport_owner &rhs);
        _lrreport_owner& operator=(const _lrreport_owner &rhs);
        virtual ~_lrreport_owner();
        alglib_impl::lrreport* c_ptr();
        alglib_impl::lrreport* c_ptr() const;
        protected:
            alglib_impl::lrreport *p_struct;
    }

    class lrreport :  _lrreport_owner
    {
        lrreport();
        lrreport(const lrreport &rhs);
        lrreport& operator=(const lrreport &rhs);
        virtual ~lrreport();
        real_2d_array c;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
        double &cvrmserror;
        double &cvavgerror;
        double &cvavgrelerror;
        ae_int_t &ncvdefects;
        integer_1d_array cvdefects;
    }

    class _modelerrors_owner
    {
        _modelerrors_owner();
        _modelerrors_owner(const _modelerrors_owner &rhs);
        _modelerrors_owner& operator=(const _modelerrors_owner &rhs);
        virtual ~_modelerrors_owner();
        alglib_impl::modelerrors* c_ptr();
        alglib_impl::modelerrors* c_ptr() const;
        protected:
            alglib_impl::modelerrors *p_struct;
    }

    class modelerrors :  _modelerrors_owner
    {
        modelerrors();
        modelerrors(const modelerrors &rhs);
        modelerrors& operator=(const modelerrors &rhs);
        virtual ~modelerrors();
        double &relclserror;
        double &avgce;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
    }

    class _multilayerperceptron_owner
    {
        _multilayerperceptron_owner();
        _multilayerperceptron_owner(const _multilayerperceptron_owner &rhs);
        _multilayerperceptron_owner& operator=(const _multilayerperceptron_owner &rhs);
        virtual ~_multilayerperceptron_owner();
        alglib_impl::multilayerperceptron* c_ptr();
        alglib_impl::multilayerperceptron* c_ptr() const;
        protected:
            alglib_impl::multilayerperceptron *p_struct;
    }

    class multilayerperceptron :  _multilayerperceptron_owner
    {
        multilayerperceptron();
        multilayerperceptron(const multilayerperceptron &rhs);
        multilayerperceptron& operator=(const multilayerperceptron &rhs);
        virtual ~multilayerperceptron();
    }
    
    class _logitmodel_owner
    {
        _logitmodel_owner();
        _logitmodel_owner(const _logitmodel_owner &rhs);
        _logitmodel_owner& operator=(const _logitmodel_owner &rhs);
        virtual ~_logitmodel_owner();
        alglib_impl::logitmodel* c_ptr();
        alglib_impl::logitmodel* c_ptr() const;
        protected:
            alglib_impl::logitmodel *p_struct;
    };
    class logitmodel : public _logitmodel_owner
    {
        logitmodel();
        logitmodel(const logitmodel &rhs);
        logitmodel& operator=(const logitmodel &rhs);
        virtual ~logitmodel();
    }

    class _mnlreport_owner
    {
        _mnlreport_owner();
        _mnlreport_owner(const _mnlreport_owner &rhs);
        _mnlreport_owner& operator=(const _mnlreport_owner &rhs);
        virtual ~_mnlreport_owner();
        alglib_impl::mnlreport* c_ptr();
        alglib_impl::mnlreport* c_ptr() const;
        protected:
            alglib_impl::mnlreport *p_struct;
    }

    class mnlreport :  _mnlreport_owner
    {
        mnlreport();
        mnlreport(const mnlreport &rhs);
        mnlreport& operator=(const mnlreport &rhs);
        virtual ~mnlreport();
        ae_int_t &ngrad;
        ae_int_t &nhess;
    }

    class _mcpdstate_owner
    {
        _mcpdstate_owner();
        _mcpdstate_owner(const _mcpdstate_owner &rhs);
        _mcpdstate_owner& operator=(const _mcpdstate_owner &rhs);
        virtual ~_mcpdstate_owner();
        alglib_impl::mcpdstate* c_ptr();
        alglib_impl::mcpdstate* c_ptr() const;
        protected:
            alglib_impl::mcpdstate *p_struct;
    }

    class mcpdstate :  _mcpdstate_owner
    {
        mcpdstate();
        mcpdstate(const mcpdstate &rhs);
        mcpdstate& operator=(const mcpdstate &rhs);
        virtual ~mcpdstate();
    }

    class _mcpdreport_owner
    {
        _mcpdreport_owner();
        _mcpdreport_owner(const _mcpdreport_owner &rhs);
        _mcpdreport_owner& operator=(const _mcpdreport_owner &rhs);
        virtual ~_mcpdreport_owner();
        alglib_impl::mcpdreport* c_ptr();
        alglib_impl::mcpdreport* c_ptr() const;
        protected:
            alglib_impl::mcpdreport *p_struct;
    }

    class mcpdreport : _mcpdreport_owner
    {
        mcpdreport();
        mcpdreport(const mcpdreport &rhs);
        mcpdreport& operator=(const mcpdreport &rhs);
        virtual ~mcpdreport();
        ae_int_t &inneriterationscount;
        ae_int_t &outeriterationscount;
        ae_int_t &nfev;
        ae_int_t &terminationtype;
    }

    class _mlpensemble_owner
    {
        _mlpensemble_owner();
        _mlpensemble_owner(const _mlpensemble_owner &rhs);
        _mlpensemble_owner& operator=(const _mlpensemble_owner &rhs);
        virtual ~_mlpensemble_owner();
        alglib_impl::mlpensemble* c_ptr();
        alglib_impl::mlpensemble* c_ptr() const;
        protected:
            alglib_impl::mlpensemble *p_struct;
    }

    class mlpensemble :  _mlpensemble_owner
    {
        mlpensemble();
        mlpensemble(const mlpensemble &rhs);
        mlpensemble& operator=(const mlpensemble &rhs);
        virtual ~mlpensemble();
    }

    class _mlpreport_owner
    {
        _mlpreport_owner();
        _mlpreport_owner(const _mlpreport_owner &rhs);
        _mlpreport_owner& operator=(const _mlpreport_owner &rhs);
        virtual ~_mlpreport_owner();
        alglib_impl::mlpreport* c_ptr();
        alglib_impl::mlpreport* c_ptr() const;
        protected:
            alglib_impl::mlpreport *p_struct;
    }

    class mlpreport : _mlpreport_owner
    {
        mlpreport();
        mlpreport(const mlpreport &rhs);
        mlpreport& operator=(const mlpreport &rhs);
        virtual ~mlpreport();
        double &relclserror;
        double &avgce;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
        ae_int_t &ngrad;
        ae_int_t &nhess;
        ae_int_t &ncholesky;
    }


    class _mlpcvreport_owner
    {
        _mlpcvreport_owner();
        _mlpcvreport_owner(const _mlpcvreport_owner &rhs);
        _mlpcvreport_owner& operator=(const _mlpcvreport_owner &rhs);
        virtual ~_mlpcvreport_owner();
        alglib_impl::mlpcvreport* c_ptr();
        alglib_impl::mlpcvreport* c_ptr() const;
        protected:
            alglib_impl::mlpcvreport *p_struct;
    }

    class mlpcvreport :  _mlpcvreport_owner
    {
        mlpcvreport();
        mlpcvreport(const mlpcvreport &rhs);
        mlpcvreport& operator=(const mlpcvreport &rhs);
        virtual ~mlpcvreport();
        double &relclserror;
        double &avgce;
        double &rmserror;
        double &avgerror;
        double &avgrelerror;
    }

    class _mlptrainer_owner
    {
        _mlptrainer_owner();
        _mlptrainer_owner(const _mlptrainer_owner &rhs);
        _mlptrainer_owner& operator=(const _mlptrainer_owner &rhs);
        virtual ~_mlptrainer_owner();
        alglib_impl::mlptrainer* c_ptr();
        alglib_impl::mlptrainer* c_ptr() const;
        protected:
            alglib_impl::mlptrainer *p_struct;
    }

    class mlptrainer :  _mlptrainer_owner
    {
        mlptrainer();
        mlptrainer(const mlptrainer &rhs);
        mlptrainer& operator=(const mlptrainer &rhs);
        virtual ~mlptrainer();
    }
    
    void dsoptimalsplit2(const real_1d_array &a, const integer_1d_array &c, const ae_int_t n, ae_int_t &info, double &threshold, double &pal, double &pbl, double &par, double &pbr, double &cve);
    void dsoptimalsplit2fast(real_1d_array &a, integer_1d_array &c, integer_1d_array &tiesbuf, integer_1d_array &cntbuf, real_1d_array &bufr, integer_1d_array &bufi, const ae_int_t n, const ae_int_t nc, const double alpha, ae_int_t &info, double &threshold, double &rms, double &cvrms);
    void clusterizercreate(clusterizerstate &s);
    void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype);
    void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t disttype);
    void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const ae_int_t npoints, const bool isupper);
    void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const bool isupper);
    void clusterizersetahcalgo(const clusterizerstate &s, const ae_int_t algo);
    void clusterizersetkmeanslimits(const clusterizerstate &s, const ae_int_t restarts, const ae_int_t maxits);
    void clusterizersetkmeansinit(const clusterizerstate &s, const ae_int_t initalgo);
    void clusterizerrunahc(const clusterizerstate &s, ahcreport &rep);
    void smp_clusterizerrunahc(const clusterizerstate &s, ahcreport &rep);
    void clusterizerrunkmeans(const clusterizerstate &s, const ae_int_t k, kmeansreport &rep);
    void smp_clusterizerrunkmeans(const clusterizerstate &s, const ae_int_t k, kmeansreport &rep);
    void clusterizergetdistances(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, real_2d_array &d);
    void smp_clusterizergetdistances(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, real_2d_array &d);
    void clusterizergetkclusters(const ahcreport &rep, const ae_int_t k, integer_1d_array &cidx, integer_1d_array &cz);
    void clusterizerseparatedbydist(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);
    void clusterizerseparatedbycorr(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);
    void kmeansgenerate(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t k, const ae_int_t restarts, ae_int_t &info, real_2d_array &c, integer_1d_array &xyc);
    void dfserialize(decisionforest &obj, std::string &s_out);
    void dfunserialize(std::string &s_in, decisionforest &obj);
    void dfbuildrandomdecisionforest(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);
    void dfbuildrandomdecisionforestx1(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const ae_int_t nrndvars, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);
    void dfprocess(const decisionforest &df, const real_1d_array &x, real_1d_array &y);
    void dfprocessi(const decisionforest &df, const real_1d_array &x, real_1d_array &y);
    double dfrelclserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
    double dfavgce(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
    double dfrmserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
    double dfavgerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
    double dfavgrelerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);
    void lrbuild(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
    void lrbuilds(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
    void lrbuildzs(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
    void lrbuildz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);
    void lrunpack(const linearmodel &lm, real_1d_array &v, ae_int_t &nvars);
    void lrpack(const real_1d_array &v, const ae_int_t nvars, linearmodel &lm);
    double lrprocess(const linearmodel &lm, const real_1d_array &x);
    double lrrmserror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    double lravgerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    double lravgrelerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    void filtersma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
    void filtersma(real_1d_array &x, const ae_int_t k);
    void filterema(real_1d_array &x, const ae_int_t n, const double alpha);
    void filterema(real_1d_array &x, const double alpha);
    void filterlrma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
    void filterlrma(real_1d_array &x, const ae_int_t k);
    void fisherlda(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_1d_array &w);
    void fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w);
    void smp_fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w);
    void mlpserialize(multilayerperceptron &obj, std::string &s_out);
    void mlpunserialize(std::string &s_in, multilayerperceptron &obj);
    void mlpcreate0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);
    void mlpcreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);
    void mlpcreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);
    void mlpcreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);
    void mlpcreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);
    void mlpcreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);
    void mlpcreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);
    void mlpcreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);
    void mlpcreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);
    void mlpcreatec0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);
    void mlpcreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);
    void mlpcreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);
    void mlpcopy(const multilayerperceptron &network1, multilayerperceptron &network2);
    void mlpcopytunableparameters(const multilayerperceptron &network1, const multilayerperceptron &network2);
    void mlprandomize(const multilayerperceptron &network);
    void mlprandomizefull(const multilayerperceptron &network);
    void mlpinitpreprocessor(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);
    void mlpproperties(const multilayerperceptron &network, ae_int_t &nin, ae_int_t &nout, ae_int_t &wcount);
    ae_int_t mlpgetinputscount(const multilayerperceptron &network);
    ae_int_t mlpgetoutputscount(const multilayerperceptron &network);
    ae_int_t mlpgetweightscount(const multilayerperceptron &network);
    bool mlpissoftmax(const multilayerperceptron &network);
    ae_int_t mlpgetlayerscount(const multilayerperceptron &network);
    ae_int_t mlpgetlayersize(const multilayerperceptron &network, const ae_int_t k);
    void mlpgetinputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);
    void mlpgetoutputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);
    void mlpgetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, ae_int_t &fkind, double &threshold);
    double mlpgetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1);
    void mlpsetinputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);
    void mlpsetoutputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);
    void mlpsetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, const ae_int_t fkind, const double threshold);
    void mlpsetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const double w);
    void mlpactivationfunction(const double net, const ae_int_t k, double &f, double &df, double &d2f);
    void mlpprocess(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);
    void mlpprocessi(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);
    double mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double smp_mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double mlperrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double smp_mlperrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double mlperrorn(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);
    ae_int_t mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    ae_int_t smp_mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double smp_mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double mlprelclserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double smp_mlprelclserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double smp_mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double mlpavgcesparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double smp_mlpavgcesparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double smp_mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double mlprmserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double smp_mlprmserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double smp_mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double mlpavgerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double smp_mlpavgerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double smp_mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
    double mlpavgrelerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    double smp_mlpavgrelerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
    void mlpgrad(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);
    void mlpgradn(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);
    void mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
    void smp_mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
    void mlpgradbatchsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
    void smp_mlpgradbatchsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
    void mlpgradbatchsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
    void smp_mlpgradbatchsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
    void mlpgradbatchsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
    void smp_mlpgradbatchsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
    void mlpgradnbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
    void mlphessiannbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);
    void mlphessianbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);
    void mlpallerrorssubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
    void smp_mlpallerrorssubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
    void mlpallerrorssparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
    void smp_mlpallerrorssparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
    double mlperrorsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
    double smp_mlperrorsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
    double mlperrorsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
    double smp_mlperrorsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
    void mnltrainh(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, logitmodel &lm, mnlreport &rep);
    void mnlprocess(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);
    void mnlprocessi(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);
    void mnlunpack(const logitmodel &lm, real_2d_array &a, ae_int_t &nvars, ae_int_t &nclasses);
    void mnlpack(const real_2d_array &a, const ae_int_t nvars, const ae_int_t nclasses, logitmodel &lm);
    double mnlavgce(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    double mnlrelclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    double mnlrmserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    double mnlavgerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    double mnlavgrelerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t ssize);
    ae_int_t mnlclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);
    void mcpdcreate(const ae_int_t n, mcpdstate &s);
    void mcpdcreateentry(const ae_int_t n, const ae_int_t entrystate, mcpdstate &s);
    void mcpdcreateexit(const ae_int_t n, const ae_int_t exitstate, mcpdstate &s);
    void mcpdcreateentryexit(const ae_int_t n, const ae_int_t entrystate, const ae_int_t exitstate, mcpdstate &s);
    void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const ae_int_t k);
    void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy);
    void mcpdsetec(const mcpdstate &s, const real_2d_array &ec);
    void mcpdaddec(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double c);
    void mcpdsetbc(const mcpdstate &s, const real_2d_array &bndl, const real_2d_array &bndu);
    void mcpdaddbc(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double bndl, const double bndu);
    void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
    void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct);
    void mcpdsettikhonovregularizer(const mcpdstate &s, const double v);
    void mcpdsetprior(const mcpdstate &s, const real_2d_array &pp);
    void mcpdsetpredictionweights(const mcpdstate &s, const real_1d_array &pw);
    void mcpdsolve(const mcpdstate &s);
    void mcpdresults(const mcpdstate &s, real_2d_array &p, mcpdreport &rep);
    void mlpeserialize(mlpensemble &obj, std::string &s_out);
    void mlpeunserialize(std::string &s_in, mlpensemble &obj);
    void mlpecreate0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreatec0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlpecreatefromnetwork(const multilayerperceptron &network, const ae_int_t ensemblesize, mlpensemble &ensemble);
    void mlperandomize(const mlpensemble &ensemble);
    void mlpeproperties(const mlpensemble &ensemble, ae_int_t &nin, ae_int_t &nout);
    bool mlpeissoftmax(const mlpensemble &ensemble);
    void mlpeprocess(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);
    void mlpeprocessi(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);
    double mlperelclserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
    double mlpeavgce(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
    double mlpermserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
    double mlpeavgerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
    double mlpeavgrelerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);
    void mlptrainlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);
    void mlptrainlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep);
    void mlptraines(const multilayerperceptron &network, const real_2d_array &trnxy, const ae_int_t trnsize, const real_2d_array &valxy, const ae_int_t valsize, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);
    void mlpkfoldcvlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);
    void mlpkfoldcvlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);
    void mlpkfoldcv(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, const ae_int_t foldscount, mlpreport &rep);
    void smp_mlpkfoldcv(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, const ae_int_t foldscount, mlpreport &rep);
    void mlpcreatetrainercls(const ae_int_t nin, const ae_int_t nclasses, mlptrainer &s);
    void mlpsetdataset(const mlptrainer &s, const real_2d_array &xy, const ae_int_t npoints);
    void mlpsetsparsedataset(const mlptrainer &s, const sparsematrix &xy, const ae_int_t npoints);
    void mlpsetcond(const mlptrainer &s, const double wstep, const ae_int_t maxits);
    void mlpsetalgobatch(const mlptrainer &s);
    void mlptrainnetwork(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, mlpreport &rep);
    void smp_mlptrainnetwork(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, mlpreport &rep);
    void mlpstarttraining(const mlptrainer &s, const multilayerperceptron &network, const bool randomstart);
    bool mlpcontinuetraining(const mlptrainer &s, const multilayerperceptron &network);
    bool smp_mlpcontinuetraining(const mlptrainer &s, const multilayerperceptron &network);
    void mlpebagginglm(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);
    void mlpetraines(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);
    void mlptrainensemblees(const mlptrainer &s, const mlpensemble &ensemble, const ae_int_t nrestarts, mlpreport &rep);
    void smp_mlptrainensemblees(const mlptrainer &s, const mlpensemble &ensemble, const ae_int_t nrestarts, mlpreport &rep);
    void pcabuildbasis(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, real_1d_array &s2, real_2d_array &v);
}
*/

extern(C)
{
    void dserrallocate(ae_int_t nclasses, /* Real    */ ae_vector* buf, ae_state *_state);
    void dserraccumulate(/* Real    */ ae_vector* buf, /* Real    */ ae_vector* y, /* Real    */ ae_vector* desiredy, ae_state *_state);
    void dserrfinish(/* Real    */ ae_vector* buf, ae_state *_state);
    void dsnormalize(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t* info, /* Real    */ ae_vector* means, /* Real    */ ae_vector* sigmas, ae_state *_state);
    void dsnormalizec(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t* info, /* Real    */ ae_vector* means, /* Real    */ ae_vector* sigmas, ae_state *_state);
    double dsgetmeanmindistance(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_state *_state);
    void dstie(/* Real    */ ae_vector* a, ae_int_t n, /* Integer */ ae_vector* ties, ae_int_t* tiecount, /* Integer */ ae_vector* p1, /* Integer */ ae_vector* p2, ae_state *_state);
    void dstiefasti(/* Real    */ ae_vector* a, /* Integer */ ae_vector* b, ae_int_t n, /* Integer */ ae_vector* ties, ae_int_t* tiecount, /* Real    */ ae_vector* bufr, /* Integer */ ae_vector* bufi, ae_state *_state);
    void dsoptimalsplit2(/* Real    */ ae_vector* a, /* Integer */ ae_vector* c, ae_int_t n, ae_int_t* info, double* threshold, double* pal, double* pbl, double* par, double* pbr, double* cve, ae_state *_state);
    void dsoptimalsplit2fast(/* Real    */ ae_vector* a, /* Integer */ ae_vector* c, /* Integer */ ae_vector* tiesbuf, /* Integer */ ae_vector* cntbuf, /* Real    */ ae_vector* bufr, /* Integer */ ae_vector* bufi, ae_int_t n, ae_int_t nc, double alpha, ae_int_t* info, double* threshold, double* rms, double* cvrms, ae_state *_state);
    void dssplitk(/* Real    */ ae_vector* a, /* Integer */ ae_vector* c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t* info, /* Real    */ ae_vector* thresholds, ae_int_t* ni, double* cve, ae_state *_state);
    void dsoptimalsplitk(/* Real    */ ae_vector* a, /* Integer */ ae_vector* c, ae_int_t n, ae_int_t nc, ae_int_t kmax, ae_int_t* info, /* Real    */ ae_vector* thresholds, ae_int_t* ni, double* cve, ae_state *_state);
    void _cvreport_init(void* _p, ae_state *_state);
    void _cvreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _cvreport_clear(void* _p);
    void _cvreport_destroy(void* _p);
    void clusterizercreate(clusterizerstate* s, ae_state *_state);
    void clusterizersetpoints(clusterizerstate* s, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, ae_state *_state);
    void clusterizersetdistances(clusterizerstate* s, /* Real    */ ae_matrix* d, ae_int_t npoints, ae_bool isupper, ae_state *_state);
    void clusterizersetahcalgo(clusterizerstate* s, ae_int_t algo, ae_state *_state);
    void clusterizersetkmeanslimits(clusterizerstate* s, ae_int_t restarts, ae_int_t maxits, ae_state *_state);
    void clusterizersetkmeansinit(clusterizerstate* s, ae_int_t initalgo, ae_state *_state);
    void clusterizerrunahc(clusterizerstate* s, ahcreport* rep, ae_state *_state);
    void _pexec_clusterizerrunahc(clusterizerstate* s, ahcreport* rep, ae_state *_state);
    void clusterizerrunkmeans(clusterizerstate* s, ae_int_t k, kmeansreport* rep, ae_state *_state);
    void _pexec_clusterizerrunkmeans(clusterizerstate* s, ae_int_t k, kmeansreport* rep, ae_state *_state);
    void clusterizergetdistances(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, /* Real    */ ae_matrix* d, ae_state *_state);
    void _pexec_clusterizergetdistances(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, /* Real    */ ae_matrix* d, ae_state *_state);
    void clusterizergetdistancesbuf(apbuffers* buf, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nfeatures, ae_int_t disttype, /* Real    */ ae_matrix* d, ae_state *_state);
    void clusterizergetkclusters(ahcreport* rep, ae_int_t k, /* Integer */ ae_vector* cidx, /* Integer */ ae_vector* cz, ae_state *_state);
    void clusterizerseparatedbydist(ahcreport* rep, double r, ae_int_t* k, /* Integer */ ae_vector* cidx, /* Integer */ ae_vector* cz, ae_state *_state);
    void clusterizerseparatedbycorr(ahcreport* rep, double r, ae_int_t* k, /* Integer */ ae_vector* cidx, /* Integer */ ae_vector* cz, ae_state *_state);
    void kmeansinitbuf(kmeansbuffers* buf, ae_state *_state);
    void kmeansgenerateinternal(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t initalgo, ae_int_t maxits, ae_int_t restarts, ae_bool kmeansdbgnoits, ae_int_t* info, ae_int_t* iterationscount, /* Real    */ ae_matrix* ccol, ae_bool needccol, /* Real    */ ae_matrix* crow, ae_bool needcrow, /* Integer */ ae_vector* xyc, double* energy, kmeansbuffers* buf, ae_state *_state);
    void kmeansupdatedistances(/* Real    */ ae_matrix* xy, ae_int_t idx0, ae_int_t idx1, ae_int_t nvars, /* Real    */ ae_matrix* ct, ae_int_t cidx0, ae_int_t cidx1, /* Integer */ ae_vector* xyc, /* Real    */ ae_vector* xydist2, ae_shared_pool* bufferpool, ae_state *_state);
    void _kmeansbuffers_init(void* _p, ae_state *_state);
    void _kmeansbuffers_init_copy(void* _dst, void* _src, ae_state *_state);
    void _kmeansbuffers_clear(void* _p);
    void _kmeansbuffers_destroy(void* _p);
    void _clusterizerstate_init(void* _p, ae_state *_state);
    void _clusterizerstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _clusterizerstate_clear(void* _p);
    void _clusterizerstate_destroy(void* _p);
    void _ahcreport_init(void* _p, ae_state *_state);
    void _ahcreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _ahcreport_clear(void* _p);
    void _ahcreport_destroy(void* _p);
    void _kmeansreport_init(void* _p, ae_state *_state);
    void _kmeansreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _kmeansreport_clear(void* _p);
    void _kmeansreport_destroy(void* _p);
    void kmeansgenerate(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t k, ae_int_t restarts, ae_int_t* info, /* Real    */ ae_matrix* c, /* Integer */ ae_vector* xyc, ae_state *_state);
    void dfbuildrandomdecisionforest(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, double r, ae_int_t* info, decisionforest* df, dfreport* rep, ae_state *_state);
    void dfbuildrandomdecisionforestx1(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t nrndvars, double r, ae_int_t* info, decisionforest* df, dfreport* rep, ae_state *_state);
    void dfbuildinternal(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t ntrees, ae_int_t samplesize, ae_int_t nfeatures, ae_int_t flags, ae_int_t* info, decisionforest* df, dfreport* rep, ae_state *_state);
    void dfprocess(decisionforest* df, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void dfprocessi(decisionforest* df, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    double dfrelclserror(decisionforest* df, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double dfavgce(decisionforest* df, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double dfrmserror(decisionforest* df, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double dfavgerror(decisionforest* df, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double dfavgrelerror(decisionforest* df, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    void dfcopy(decisionforest* df1, decisionforest* df2, ae_state *_state);
    void dfalloc(ae_serializer* s, decisionforest* forest, ae_state *_state);
    void dfserialize(ae_serializer* s, decisionforest* forest, ae_state *_state);
    void dfunserialize(ae_serializer* s, decisionforest* forest, ae_state *_state);
    void _decisionforest_init(void* _p, ae_state *_state);
    void _decisionforest_init_copy(void* _dst, void* _src, ae_state *_state);
    void _decisionforest_clear(void* _p);
    void _decisionforest_destroy(void* _p);
    void _dfreport_init(void* _p, ae_state *_state);
    void _dfreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _dfreport_clear(void* _p);
    void _dfreport_destroy(void* _p);
    void _dfinternalbuffers_init(void* _p, ae_state *_state);
    void _dfinternalbuffers_init_copy(void* _dst, void* _src, ae_state *_state);
    void _dfinternalbuffers_clear(void* _p);
    void _dfinternalbuffers_destroy(void* _p);
    void lrbuild(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t* info, linearmodel* lm, lrreport* ar, ae_state *_state);
    void lrbuilds(/* Real    */ ae_matrix* xy, /* Real    */ ae_vector* s, ae_int_t npoints, ae_int_t nvars, ae_int_t* info, linearmodel* lm, lrreport* ar, ae_state *_state);
    void lrbuildzs(/* Real    */ ae_matrix* xy, /* Real    */ ae_vector* s, ae_int_t npoints, ae_int_t nvars, ae_int_t* info, linearmodel* lm, lrreport* ar, ae_state *_state);
    void lrbuildz(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t* info, linearmodel* lm, lrreport* ar, ae_state *_state);
    void lrunpack(linearmodel* lm, /* Real    */ ae_vector* v, ae_int_t* nvars, ae_state *_state);
    void lrpack(/* Real    */ ae_vector* v, ae_int_t nvars, linearmodel* lm, ae_state *_state);
    double lrprocess(linearmodel* lm, /* Real    */ ae_vector* x, ae_state *_state);
    double lrrmserror(linearmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double lravgerror(linearmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double lravgrelerror(linearmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    void lrcopy(linearmodel* lm1, linearmodel* lm2, ae_state *_state);
    void lrlines(/* Real    */ ae_matrix* xy, /* Real    */ ae_vector* s, ae_int_t n, ae_int_t* info, double* a, double* b, double* vara, double* varb, double* covab, double* corrab, double* p, ae_state *_state);
    void lrline(/* Real    */ ae_matrix* xy, ae_int_t n, ae_int_t* info, double* a, double* b, ae_state *_state);
    void _linearmodel_init(void* _p, ae_state *_state);
    void _linearmodel_init_copy(void* _dst, void* _src, ae_state *_state);
    void _linearmodel_clear(void* _p);
    void _linearmodel_destroy(void* _p);
    void _lrreport_init(void* _p, ae_state *_state);
    void _lrreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _lrreport_clear(void* _p);
    void _lrreport_destroy(void* _p);
    void filtersma(/* Real    */ ae_vector* x, ae_int_t n, ae_int_t k, ae_state *_state);
    void filterema(/* Real    */ ae_vector* x, ae_int_t n, double alpha, ae_state *_state);
    void filterlrma(/* Real    */ ae_vector* x, ae_int_t n, ae_int_t k, ae_state *_state);
    void fisherlda(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t* info, /* Real    */ ae_vector* w, ae_state *_state);
    void fisherldan(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t* info, /* Real    */ ae_matrix* w, ae_state *_state);
    void _pexec_fisherldan(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t* info, /* Real    */ ae_matrix* w, ae_state *_state);
    ae_int_t mlpgradsplitcost(ae_state *_state);
    ae_int_t mlpgradsplitsize(ae_state *_state);
    void mlpcreate0(ae_int_t nin, ae_int_t nout, multilayerperceptron* network, ae_state *_state);
    void mlpcreate1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, multilayerperceptron* network, ae_state *_state);
    void mlpcreate2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, multilayerperceptron* network, ae_state *_state);
    void mlpcreateb0(ae_int_t nin, ae_int_t nout, double b, double d, multilayerperceptron* network, ae_state *_state);
    void mlpcreateb1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double b, double d, multilayerperceptron* network, ae_state *_state);
    void mlpcreateb2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double b, double d, multilayerperceptron* network, ae_state *_state);
    void mlpcreater0(ae_int_t nin, ae_int_t nout, double a, double b, multilayerperceptron* network, ae_state *_state);
    void mlpcreater1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double a, double b, multilayerperceptron* network, ae_state *_state);
    void mlpcreater2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double a, double b, multilayerperceptron* network, ae_state *_state);
    void mlpcreatec0(ae_int_t nin, ae_int_t nout, multilayerperceptron* network, ae_state *_state);
    void mlpcreatec1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, multilayerperceptron* network, ae_state *_state);
    void mlpcreatec2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, multilayerperceptron* network, ae_state *_state);
    void mlpcopy(multilayerperceptron* network1, multilayerperceptron* network2, ae_state *_state);
    void mlpcopyshared(multilayerperceptron* network1, multilayerperceptron* network2, ae_state *_state);
    ae_bool mlpsamearchitecture(multilayerperceptron* network1, multilayerperceptron* network2, ae_state *_state);
    void mlpcopytunableparameters(multilayerperceptron* network1, multilayerperceptron* network2, ae_state *_state);
    void mlpexporttunableparameters(multilayerperceptron* network, /* Real    */ ae_vector* p, ae_int_t* pcount, ae_state *_state);
    void mlpimporttunableparameters(multilayerperceptron* network, /* Real    */ ae_vector* p, ae_state *_state);
    void mlpserializeold(multilayerperceptron* network, /* Real    */ ae_vector* ra, ae_int_t* rlen, ae_state *_state);
    void mlpunserializeold(/* Real    */ ae_vector* ra, multilayerperceptron* network, ae_state *_state);
    void mlprandomize(multilayerperceptron* network, ae_state *_state);
    void mlprandomizefull(multilayerperceptron* network, ae_state *_state);
    void mlpinitpreprocessor(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t ssize, ae_state *_state);
    void mlpinitpreprocessorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t ssize, ae_state *_state);
    void mlpinitpreprocessorsubset(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t setsize, /* Integer */ ae_vector* idx, ae_int_t subsetsize, ae_state *_state);
    void mlpinitpreprocessorsparsesubset(multilayerperceptron* network, sparsematrix* xy, ae_int_t setsize, /* Integer */ ae_vector* idx, ae_int_t subsetsize, ae_state *_state);
    void mlpproperties(multilayerperceptron* network, ae_int_t* nin, ae_int_t* nout, ae_int_t* wcount, ae_state *_state);
    ae_int_t mlpntotal(multilayerperceptron* network, ae_state *_state);
    ae_int_t mlpgetinputscount(multilayerperceptron* network, ae_state *_state);
    ae_int_t mlpgetoutputscount(multilayerperceptron* network, ae_state *_state);
    ae_int_t mlpgetweightscount(multilayerperceptron* network, ae_state *_state);
    ae_bool mlpissoftmax(multilayerperceptron* network, ae_state *_state);
    ae_int_t mlpgetlayerscount(multilayerperceptron* network, ae_state *_state);
    ae_int_t mlpgetlayersize(multilayerperceptron* network, ae_int_t k, ae_state *_state);
    void mlpgetinputscaling(multilayerperceptron* network, ae_int_t i, double* mean, double* sigma, ae_state *_state);
    void mlpgetoutputscaling(multilayerperceptron* network, ae_int_t i, double* mean, double* sigma, ae_state *_state);
    void mlpgetneuroninfo(multilayerperceptron* network, ae_int_t k, ae_int_t i, ae_int_t* fkind, double* threshold, ae_state *_state);
    double mlpgetweight(multilayerperceptron* network, ae_int_t k0, ae_int_t i0, ae_int_t k1, ae_int_t i1, ae_state *_state);
    void mlpsetinputscaling(multilayerperceptron* network, ae_int_t i, double mean, double sigma, ae_state *_state);
    void mlpsetoutputscaling(multilayerperceptron* network, ae_int_t i, double mean, double sigma, ae_state *_state);
    void mlpsetneuroninfo(multilayerperceptron* network, ae_int_t k, ae_int_t i, ae_int_t fkind, double threshold, ae_state *_state);
    void mlpsetweight(multilayerperceptron* network, ae_int_t k0, ae_int_t i0, ae_int_t k1, ae_int_t i1, double w, ae_state *_state);
    void mlpactivationfunction(double net, ae_int_t k, double* f, double* df, double* d2f, ae_state *_state);
    void mlpprocess(multilayerperceptron* network, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void mlpprocessi(multilayerperceptron* network, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    double mlperror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlperror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlperrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlperrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double mlperrorn(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t ssize, ae_state *_state);
    ae_int_t mlpclserror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    ae_int_t _pexec_mlpclserror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlprelclserror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlprelclserror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlprelclserrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlprelclserrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpavgce(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlpavgce(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpavgcesparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlpavgcesparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double mlprmserror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlprmserror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlprmserrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlprmserrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpavgerror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlpavgerror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpavgerrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlpavgerrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpavgrelerror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlpavgrelerror(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpavgrelerrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    double _pexec_mlpavgrelerrorsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    void mlpgrad(multilayerperceptron* network, /* Real    */ ae_vector* x, /* Real    */ ae_vector* desiredy, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void mlpgradn(multilayerperceptron* network, /* Real    */ ae_vector* x, /* Real    */ ae_vector* desiredy, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void mlpgradbatch(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t ssize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void _pexec_mlpgradbatch(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t ssize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void mlpgradbatchsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t ssize, double* e, /* Real    */ ae_vector* grad, ae_state *_state); void _pexec_mlpgradbatchsparse(multilayerperceptron* network, sparsematrix* xy, ae_int_t ssize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void mlpgradbatchsubset(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t setsize, /* Integer */ ae_vector* idx, ae_int_t subsetsize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void _pexec_mlpgradbatchsubset(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t setsize, /* Integer */ ae_vector* idx, ae_int_t subsetsize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void mlpgradbatchsparsesubset(multilayerperceptron* network, sparsematrix* xy, ae_int_t setsize, /* Integer */ ae_vector* idx, ae_int_t subsetsize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void _pexec_mlpgradbatchsparsesubset(multilayerperceptron* network, sparsematrix* xy, ae_int_t setsize, /* Integer */ ae_vector* idx, ae_int_t subsetsize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void mlpgradbatchx(multilayerperceptron* network, /* Real    */ ae_matrix* densexy, sparsematrix* sparsexy, ae_int_t datasetsize, ae_int_t datasettype, /* Integer */ ae_vector* idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool* buf, ae_shared_pool* gradbuf, ae_state *_state);
    void mlpgradnbatch(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t ssize, double* e, /* Real    */ ae_vector* grad, ae_state *_state);
    void mlphessiannbatch(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t ssize, double* e, /* Real    */ ae_vector* grad, /* Real    */ ae_matrix* h, ae_state *_state);
    void mlphessianbatch(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t ssize, double* e, /* Real    */ ae_vector* grad, /* Real    */ ae_matrix* h, ae_state *_state);
    void mlpinternalprocessvector(/* Integer */ ae_vector* structinfo, /* Real    */ ae_vector* weights, /* Real    */ ae_vector* columnmeans, /* Real    */ ae_vector* columnsigmas, /* Real    */ ae_vector* neurons, /* Real    */ ae_vector* dfdnet, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void mlpalloc(ae_serializer* s, multilayerperceptron* network, ae_state *_state);
    void mlpserialize(ae_serializer* s, multilayerperceptron* network, ae_state *_state);
    void mlpunserialize(ae_serializer* s, multilayerperceptron* network, ae_state *_state);
    void mlpallerrorssubset(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, modelerrors* rep, ae_state *_state);
    void _pexec_mlpallerrorssubset(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, modelerrors* rep, ae_state *_state);
    void mlpallerrorssparsesubset(multilayerperceptron* network, sparsematrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, modelerrors* rep, ae_state *_state);
    void _pexec_mlpallerrorssparsesubset(multilayerperceptron* network, sparsematrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, modelerrors* rep, ae_state *_state);
    double mlperrorsubset(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, ae_state *_state);
    double _pexec_mlperrorsubset(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, ae_state *_state);
    double mlperrorsparsesubset(multilayerperceptron* network, sparsematrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, ae_state *_state);
    double _pexec_mlperrorsparsesubset(multilayerperceptron* network, sparsematrix* xy, ae_int_t setsize, /* Integer */ ae_vector* subset, ae_int_t subsetsize, ae_state *_state);
    void mlpallerrorsx(multilayerperceptron* network, /* Real    */ ae_matrix* densexy, sparsematrix* sparsexy, ae_int_t datasetsize, ae_int_t datasettype, /* Integer */ ae_vector* idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool* buf, modelerrors* rep, ae_state *_state);
    void _modelerrors_init(void* _p, ae_state *_state);
    void _modelerrors_init_copy(void* _dst, void* _src, ae_state *_state);
    void _modelerrors_clear(void* _p);
    void _modelerrors_destroy(void* _p);
    void _smlpgrad_init(void* _p, ae_state *_state);
    void _smlpgrad_init_copy(void* _dst, void* _src, ae_state *_state);
    void _smlpgrad_clear(void* _p);
    void _smlpgrad_destroy(void* _p);
    void _multilayerperceptron_init(void* _p, ae_state *_state);
    void _multilayerperceptron_init_copy(void* _dst, void* _src, ae_state *_state);
    void _multilayerperceptron_clear(void* _p);
    void _multilayerperceptron_destroy(void* _p);
    void mnltrainh(/* Real    */ ae_matrix* xy, ae_int_t npoints, ae_int_t nvars, ae_int_t nclasses, ae_int_t* info, logitmodel* lm, mnlreport* rep, ae_state *_state);
    void mnlprocess(logitmodel* lm, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void mnlprocessi(logitmodel* lm, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void mnlunpack(logitmodel* lm, /* Real    */ ae_matrix* a, ae_int_t* nvars, ae_int_t* nclasses, ae_state *_state);
    void mnlpack(/* Real    */ ae_matrix* a, ae_int_t nvars, ae_int_t nclasses, logitmodel* lm, ae_state *_state);
    void mnlcopy(logitmodel* lm1, logitmodel* lm2, ae_state *_state);
    double mnlavgce(logitmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mnlrelclserror(logitmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mnlrmserror(logitmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mnlavgerror(logitmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mnlavgrelerror(logitmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t ssize, ae_state *_state);
    ae_int_t mnlclserror(logitmodel* lm, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    void _logitmodel_init(void* _p, ae_state *_state);
    void _logitmodel_init_copy(void* _dst, void* _src, ae_state *_state);
    void _logitmodel_clear(void* _p);
    void _logitmodel_destroy(void* _p);
    void _logitmcstate_init(void* _p, ae_state *_state);
    void _logitmcstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _logitmcstate_clear(void* _p);
    void _logitmcstate_destroy(void* _p);
    void _mnlreport_init(void* _p, ae_state *_state);
    void _mnlreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mnlreport_clear(void* _p);
    void _mnlreport_destroy(void* _p);
    void mcpdcreate(ae_int_t n, mcpdstate* s, ae_state *_state);
    void mcpdcreateentry(ae_int_t n, ae_int_t entrystate, mcpdstate* s, ae_state *_state);
    void mcpdcreateexit(ae_int_t n, ae_int_t exitstate, mcpdstate* s, ae_state *_state);
    void mcpdcreateentryexit(ae_int_t n, ae_int_t entrystate, ae_int_t exitstate, mcpdstate* s, ae_state *_state);
    void mcpdaddtrack(mcpdstate* s, /* Real    */ ae_matrix* xy, ae_int_t k, ae_state *_state);
    void mcpdsetec(mcpdstate* s, /* Real    */ ae_matrix* ec, ae_state *_state);
    void mcpdaddec(mcpdstate* s, ae_int_t i, ae_int_t j, double c, ae_state *_state);
    void mcpdsetbc(mcpdstate* s, /* Real    */ ae_matrix* bndl, /* Real    */ ae_matrix* bndu, ae_state *_state);
    void mcpdaddbc(mcpdstate* s, ae_int_t i, ae_int_t j, double bndl, double bndu, ae_state *_state);
    void mcpdsetlc(mcpdstate* s, /* Real    */ ae_matrix* c, /* Integer */ ae_vector* ct, ae_int_t k, ae_state *_state);
    void mcpdsettikhonovregularizer(mcpdstate* s, double v, ae_state *_state);
    void mcpdsetprior(mcpdstate* s, /* Real    */ ae_matrix* pp, ae_state *_state);
    void mcpdsetpredictionweights(mcpdstate* s, /* Real    */ ae_vector* pw, ae_state *_state);
    void mcpdsolve(mcpdstate* s, ae_state *_state);
    void mcpdresults(mcpdstate* s, /* Real    */ ae_matrix* p, mcpdreport* rep, ae_state *_state);
    void _mcpdstate_init(void* _p, ae_state *_state);
    void _mcpdstate_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mcpdstate_clear(void* _p);
    void _mcpdstate_destroy(void* _p);
    void _mcpdreport_init(void* _p, ae_state *_state);
    void _mcpdreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mcpdreport_clear(void* _p);
    void _mcpdreport_destroy(void* _p);
    void mlpecreate0(ae_int_t nin, ae_int_t nout, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreate1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreate2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreateb0(ae_int_t nin, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreateb1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreateb2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double b, double d, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreater0(ae_int_t nin, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreater1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreater2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, double a, double b, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreatec0(ae_int_t nin, ae_int_t nout, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreatec1(ae_int_t nin, ae_int_t nhid, ae_int_t nout, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreatec2(ae_int_t nin, ae_int_t nhid1, ae_int_t nhid2, ae_int_t nout, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecreatefromnetwork(multilayerperceptron* network, ae_int_t ensemblesize, mlpensemble* ensemble, ae_state *_state);
    void mlpecopy(mlpensemble* ensemble1, mlpensemble* ensemble2, ae_state *_state);
    void mlperandomize(mlpensemble* ensemble, ae_state *_state);
    void mlpeproperties(mlpensemble* ensemble, ae_int_t* nin, ae_int_t* nout, ae_state *_state);
    ae_bool mlpeissoftmax(mlpensemble* ensemble, ae_state *_state);
    void mlpeprocess(mlpensemble* ensemble, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void mlpeprocessi(mlpensemble* ensemble, /* Real    */ ae_vector* x, /* Real    */ ae_vector* y, ae_state *_state);
    void mlpeallerrorsx(mlpensemble* ensemble, /* Real    */ ae_matrix* densexy, sparsematrix* sparsexy, ae_int_t datasetsize, ae_int_t datasettype, /* Integer */ ae_vector* idx, ae_int_t subset0, ae_int_t subset1, ae_int_t subsettype, ae_shared_pool* buf, modelerrors* rep, ae_state *_state);
    void mlpeallerrorssparse(mlpensemble* ensemble, sparsematrix* xy, ae_int_t npoints, double* relcls, double* avgce, double* rms, double* avg, double* avgrel, ae_state *_state);
    double mlperelclserror(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpeavgce(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpermserror(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpeavgerror(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    double mlpeavgrelerror(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    void mlpealloc(ae_serializer* s, mlpensemble* ensemble, ae_state *_state);
    void mlpeserialize(ae_serializer* s, mlpensemble* ensemble, ae_state *_state);
    void mlpeunserialize(ae_serializer* s, mlpensemble* ensemble, ae_state *_state);
    void _mlpensemble_init(void* _p, ae_state *_state);
    void _mlpensemble_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mlpensemble_clear(void* _p);
    void _mlpensemble_destroy(void* _p);
    void mlptrainlm(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t* info, mlpreport* rep, ae_state *_state);
    void mlptrainlbfgs(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t* info, mlpreport* rep, ae_state *_state);
    void mlptraines(multilayerperceptron* network, /* Real    */ ae_matrix* trnxy, ae_int_t trnsize, /* Real    */ ae_matrix* valxy, ae_int_t valsize, double decay, ae_int_t restarts, ae_int_t* info, mlpreport* rep, ae_state *_state);
    void mlpkfoldcvlbfgs(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t foldscount, ae_int_t* info, mlpreport* rep, mlpcvreport* cvrep, ae_state *_state);
    void mlpkfoldcvlm(multilayerperceptron* network, /* Real    */ ae_matrix* xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t foldscount, ae_int_t* info, mlpreport* rep, mlpcvreport* cvrep, ae_state *_state);
    void mlpkfoldcv(mlptrainer* s, multilayerperceptron* network, ae_int_t nrestarts, ae_int_t foldscount, mlpreport* rep, ae_state *_state);
    void _pexec_mlpkfoldcv(mlptrainer* s, multilayerperceptron* network, ae_int_t nrestarts, ae_int_t foldscount, mlpreport* rep, ae_state *_state);
    void mlpcreatetrainer(ae_int_t nin, ae_int_t nout, mlptrainer* s, ae_state *_state);
    void mlpcreatetrainercls(ae_int_t nin, ae_int_t nclasses, mlptrainer* s, ae_state *_state);
    void mlpsetdataset(mlptrainer* s, /* Real    */ ae_matrix* xy, ae_int_t npoints, ae_state *_state);
    void mlpsetsparsedataset(mlptrainer* s, sparsematrix* xy, ae_int_t npoints, ae_state *_state);
    void mlpsetdecay(mlptrainer* s, double decay, ae_state *_state);
    void mlpsetcond(mlptrainer* s, double wstep, ae_int_t maxits, ae_state *_state);
    void mlpsetalgobatch(mlptrainer* s, ae_state *_state);
    void mlptrainnetwork(mlptrainer* s, multilayerperceptron* network, ae_int_t nrestarts, mlpreport* rep, ae_state *_state);
    void _pexec_mlptrainnetwork(mlptrainer* s, multilayerperceptron* network, ae_int_t nrestarts, mlpreport* rep, ae_state *_state);
    void mlpstarttraining(mlptrainer* s, multilayerperceptron* network, ae_bool randomstart, ae_state *_state);
    ae_bool mlpcontinuetraining(mlptrainer* s, multilayerperceptron* network, ae_state *_state);
    ae_bool _pexec_mlpcontinuetraining(mlptrainer* s, multilayerperceptron* network, ae_state *_state);
    void mlpebagginglm(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t* info, mlpreport* rep, mlpcvreport* ooberrors, ae_state *_state);
    void mlpebagginglbfgs(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, double decay, ae_int_t restarts, double wstep, ae_int_t maxits, ae_int_t* info, mlpreport* rep, mlpcvreport* ooberrors, ae_state *_state);
    void mlpetraines(mlpensemble* ensemble, /* Real    */ ae_matrix* xy, ae_int_t npoints, double decay, ae_int_t restarts, ae_int_t* info, mlpreport* rep, ae_state *_state);
    void mlptrainensemblees(mlptrainer* s, mlpensemble* ensemble, ae_int_t nrestarts, mlpreport* rep, ae_state *_state);
    void _pexec_mlptrainensemblees(mlptrainer* s, mlpensemble* ensemble, ae_int_t nrestarts, mlpreport* rep, ae_state *_state);
    void _mlpreport_init(void* _p, ae_state *_state);
    void _mlpreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mlpreport_clear(void* _p);
    void _mlpreport_destroy(void* _p);
    void _mlpcvreport_init(void* _p, ae_state *_state);
    void _mlpcvreport_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mlpcvreport_clear(void* _p);
    void _mlpcvreport_destroy(void* _p);
    void _smlptrnsession_init(void* _p, ae_state *_state);
    void _smlptrnsession_init_copy(void* _dst, void* _src, ae_state *_state);
    void _smlptrnsession_clear(void* _p);
    void _smlptrnsession_destroy(void* _p);
    void _mlpetrnsession_init(void* _p, ae_state *_state);
    void _mlpetrnsession_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mlpetrnsession_clear(void* _p);
    void _mlpetrnsession_destroy(void* _p);
    void _mlptrainer_init(void* _p, ae_state *_state);
    void _mlptrainer_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mlptrainer_clear(void* _p);
    void _mlptrainer_destroy(void* _p);
    void _mlpparallelizationcv_init(void* _p, ae_state *_state);
    void _mlpparallelizationcv_init_copy(void* _dst, void* _src, ae_state *_state);
    void _mlpparallelizationcv_clear(void* _p);
    void _mlpparallelizationcv_destroy(void* _p);
    void pcabuildbasis(/* Real    */ ae_matrix* x, ae_int_t npoints, ae_int_t nvars, ae_int_t* info, /* Real    */ ae_vector* s2, /* Real    */ ae_matrix* v, ae_state *_state);
}

