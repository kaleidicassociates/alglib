module alglib.ap;

alias ae_int_t = int;
alias ae_int32_t = int;
alias ae_int64_t = long;
alias ae_datatype = ae_int_t;

/+
    AE_INT32_T ae_int32_t;
    int32_t ae_int32_t;
    __int32 ae_int32_t;
    int ae_int32_t;

    AE_INT64_T ae_int64_t;
    int64_t ae_int64_t;
    __int64 ae_int64_t;
    ptrdiff_t ae_int_t;
+/

enum AE_OS = AE_POSIX;
enum AE_COMPILER = AE_GNUC; // AE_MSVC
enum AE_UNKNOWN = 0;
enum AE_MSVC = 1;
enum AE_GNUC = 2;
enum AE_SUNC = 3;
enum AE_INTEL = 1;
enum AE_SPARC = 2;
enum AE_WINDOWS = 1;
enum AE_POSIX = 2;
enum AE_LOCK_ALIGNMENT = 16;


version(AE_USE_CPP_BOOL)
{
    alias ae_bool = bool;
    enum ae_true = true;
    enum ae_false = false;
}
else
{
    alias ae_bool = char;
    enum ae_true = 1;
    enum ae_false = 0;
}

struct ae_complex
{
    double x, y;
}

enum ae_error_type
{
    ERR_OK = 0,
    ERR_OUT_OF_MEMORY = 1,
    ERR_XARRAY_TOO_LARGE = 2,
    ERR_ASSERTION_FAILED = 3
}


enum
{
    OWN_CALLER=1,
    OWN_AE=2
}

enum
{
    ACT_UNCHANGED=1,
    ACT_SAME_LOCATION=2,
    ACT_NEW_LOCATION=3
}

enum
{
    DT_BOOL=1,
    DT_INT=2,
    DT_REAL=3,
    DT_COMPLEX=4
}

enum
{
    CPU_SSE2=1
}

struct x_string
{
    align(8) ae_int64_t     owner;
    align(8)  ae_int64_t     last_action;
    align(8)  char *ptr;
}

struct x_vector
{
    align(8)  ae_int64_t     cnt;
    align(8)  ae_int64_t     datatype;
    align(8)  ae_int64_t     owner;
    align(8)  ae_int64_t     last_action;
    align(8)  void *ptr;
}

struct x_matrix
{
    align(8)  ae_int64_t     rows;
    align(8)  ae_int64_t     cols;
    align(8)  ae_int64_t     stride;
    align(8)  ae_int64_t     datatype;
    align(8)  ae_int64_t     owner;
    align(8)  ae_int64_t     last_action;
    align(8)  void *ptr;
}

struct ae_dyn_block
{
    shared ae_dyn_block* p_next;
    /* void *deallocator; */
    void* function(void*) deallocator;
    shared(void*) ptr;
}

struct ae_frame
{
    ae_dyn_block db_marker;
}

struct ae_state
{
    ae_int_t endianness;
    double v_nan;
    double v_posinf;
    double v_neginf;
    shared(ae_dyn_block)*  p_top_block;
    ae_dyn_block last_block;
    version(AE_USE_CPP_ERROR_HANDLING)
    {
    }
    else
    {
        // shared(jmp_buf)*  break_jump;
    }
    shared(ae_error_type) last_error;
    shared(const(char))* error_msg;
    void *worker_thread;
    void *parent_task;
    extern(C) void *function(void*) thread_exception_handler;
}

struct ae_serializer
{
    ae_int_t mode;
    ae_int_t entries_needed;
    ae_int_t entries_saved;
    ae_int_t bytes_asked;
    ae_int_t bytes_written;

//#ifdef AE_USE_CPP_SERIALIZATION
//    std::string     *out_cppstr;
//#endif
    char *out_str;
    const(char) *in_str;
} 

void* function(void*) ae_deallocator;

struct ae_vector
{
    ae_int_t cnt;
    ae_datatype datatype;
    ae_bool is_attached;
    ae_dyn_block data;
    union Ptr
    {
        void *p_ptr;
        ae_bool *p_bool;
        ae_int_t *p_int;
        double *p_double;
        ae_complex *p_complex;
    }
    Ptr ptr;
}

struct ae_matrix
{
    ae_int_t rows;
    ae_int_t cols;
    ae_int_t stride;
    ae_datatype datatype;
    ae_bool is_attached;
   
    ae_dyn_block data;
    union Ptr
    {
        void *p_ptr;
        void **pp_void;
        ae_bool **pp_bool;
        ae_int_t **pp_int;
        double **pp_double;
        ae_complex **pp_complex;
    }
    Ptr ptr;
}

struct ae_smart_ptr
{
    void **subscriber;
    void *ptr;
    ae_bool is_owner;
    ae_bool is_dynamic;
    void* function(void*) destroy;
    ae_dyn_block frame_entry;
}

struct ae_lock
{
    void *ptr;
}

struct ae_shared_pool_entry
{
    shared(void*) obj;
    shared(void*)  next_entry;
}

struct ae_shared_pool
{
    ae_lock pool_lock;
    shared(void*)  seed_object;
    shared(ae_shared_pool_entry)* recycled_objects;
    shared(ae_shared_pool_entry)*  recycled_entries;
    shared(ae_shared_pool_entry)* enumeration_counter;
    ae_int_t                size_of_object;
    void* function(void* dst, ae_state* state) init;
    void* function(void* dst, void* src, ae_state* state) init_copy;
    void* function(void* ptr) destroy;
    ae_dyn_block frame_entry;
}

ae_int_t ae_misalignment(const void *ptr, size_t alignment);
void* ae_align(void *ptr, size_t alignment);
void* aligned_malloc(size_t size, size_t alignment);
void  aligned_free(void *block);

void* ae_malloc(size_t size, ae_state *state);
void  ae_free(void *p);
ae_int_t ae_sizeof(ae_datatype datatype);
void ae_touch_ptr(void *p);

void ae_state_init(ae_state *state);
void ae_state_clear(ae_state *state);
version(AE_USE_CPP_ERROR_HANDLING)
{}
else
{
    // void ae_state_set_break_jump(ae_state *state, jmp_buf *buf);
}
void ae_break(ae_state *state, ae_error_type error_type, const(char)* msg);

void ae_frame_make(ae_state *state, ae_frame *tmp);
void ae_frame_leave(ae_state *state);

void ae_db_attach(ae_dyn_block *block, ae_state *state);
ae_bool ae_db_malloc(ae_dyn_block *block, ae_int_t size, ae_state *state, ae_bool make_automatic);
ae_bool ae_db_realloc(ae_dyn_block *block, ae_int_t size, ae_state *state);
void ae_db_free(ae_dyn_block *block);
void ae_db_swap(ae_dyn_block *block1, ae_dyn_block *block2);

void ae_vector_init(ae_vector *dst, ae_int_t size, ae_datatype datatype, ae_state *state);
void ae_vector_init_copy(ae_vector *dst, ae_vector *src, ae_state *state);
void ae_vector_init_from_x(ae_vector *dst, x_vector *src, ae_state *state);
void ae_vector_attach_to_x(ae_vector *dst, x_vector *src, ae_state *state);
ae_bool ae_vector_set_length(ae_vector *dst, ae_int_t newsize, ae_state *state);
void ae_vector_clear(ae_vector *dst);
void ae_vector_destroy(ae_vector *dst);
void ae_swap_vectors(ae_vector *vec1, ae_vector *vec2);
void ae_matrix_init(ae_matrix *dst, ae_int_t rows, ae_int_t cols, ae_datatype datatype, ae_state *state);
void ae_matrix_init_copy(ae_matrix *dst, ae_matrix *src, ae_state *state);
void ae_matrix_init_from_x(ae_matrix *dst, x_matrix *src, ae_state *state);
void ae_matrix_attach_to_x(ae_matrix *dst, x_matrix *src, ae_state *state);
ae_bool ae_matrix_set_length(ae_matrix *dst, ae_int_t rows, ae_int_t cols, ae_state *state);
void ae_matrix_clear(ae_matrix *dst);
void ae_matrix_destroy(ae_matrix *dst);
void ae_swap_matrices(ae_matrix *mat1, ae_matrix *mat2);
void ae_smart_ptr_init(ae_smart_ptr *dst, void **subscriber, ae_state *state);
void ae_smart_ptr_clear(void *_dst); /* accepts ae_smart_ptr* */
void ae_smart_ptr_destroy(void *_dst);
void ae_smart_ptr_assign(ae_smart_ptr *dst, void *new_ptr, ae_bool is_owner, ae_bool is_dynamic, void* function(void*) destroy);
void ae_smart_ptr_release(ae_smart_ptr *dst);
void ae_yield();
void ae_init_lock(ae_lock *lock);
void ae_acquire_lock(ae_lock *lock);
void ae_release_lock(ae_lock *lock);
void ae_free_lock(ae_lock *lock);
void ae_shared_pool_init(void *_dst, ae_state *state);
void ae_shared_pool_init_copy(void *_dst, void *_src, ae_state *state);
void ae_shared_pool_clear(void *dst);
void ae_shared_pool_destroy(void *dst);
ae_bool ae_shared_pool_is_initialized(void *_dst);
void ae_shared_pool_set_seed(ae_shared_pool  *dst, void *seed_object, ae_int_t size_of_object,
    void* function(void* dst, ae_state* state) init, void* function(void* dst, void* src, ae_state* state) init_copy,
    void* function(void* ptr) destroy, ae_state *state);
void ae_shared_pool_retrieve(ae_shared_pool  *pool, ae_smart_ptr    *pptr, ae_state        *state);
void ae_shared_pool_recycle(ae_shared_pool  *pool, ae_smart_ptr    *pptr, ae_state        *state);
void ae_shared_pool_clear_recycled(ae_shared_pool  *pool, ae_state        *state);
void ae_shared_pool_first_recycled(ae_shared_pool  *pool, ae_smart_ptr    *pptr, ae_state        *state);
void ae_shared_pool_next_recycled(ae_shared_pool  *pool, ae_smart_ptr    *pptr, ae_state        *state);
void ae_shared_pool_reset(ae_shared_pool  *pool, ae_state        *state);
void ae_x_set_vector(x_vector *dst, ae_vector *src, ae_state *state);
void ae_x_set_matrix(x_matrix *dst, ae_matrix *src, ae_state *state);
void ae_x_attach_to_vector(x_vector *dst, ae_vector *src);
void ae_x_attach_to_matrix(x_matrix *dst, ae_matrix *src);
void x_vector_clear(x_vector *dst);
ae_bool x_is_symmetric(x_matrix *a);
ae_bool x_is_hermitian(x_matrix *a);
ae_bool x_force_symmetric(x_matrix *a);
ae_bool x_force_hermitian(x_matrix *a);
ae_bool ae_is_symmetric(ae_matrix *a);
ae_bool ae_is_hermitian(ae_matrix *a);
ae_bool ae_force_symmetric(ae_matrix *a);
ae_bool ae_force_hermitian(ae_matrix *a);
void ae_serializer_init(ae_serializer *serializer);
void ae_serializer_clear(ae_serializer *serializer);
void ae_serializer_alloc_start(ae_serializer *serializer);
void ae_serializer_alloc_entry(ae_serializer *serializer);
ae_int_t ae_serializer_get_alloc_size(ae_serializer *serializer);

extern(C++,alglib)
{
    version(AE_USE_CPP_SERIALIZATION)
    {
        void ae_serializer_sstart_str(ae_serializer *serializer, std.string *buf);
        void ae_serializer_ustart_str(ae_serializer *serializer, const std.string *buf);
    }
    void ae_serializer_sstart_str(ae_serializer *serializer, char *buf);
    void ae_serializer_ustart_str(ae_serializer *serializer, const(char)* buf);
    void ae_serializer_serialize_bool(ae_serializer *serializer, ae_bool v, ae_state *state);
    void ae_serializer_serialize_int(ae_serializer *serializer, ae_int_t v, ae_state *state);
    void ae_serializer_serialize_double(ae_serializer *serializer, double v, ae_state *state);
    void ae_serializer_unserialize_bool(ae_serializer *serializer, ae_bool *v, ae_state *state);
    void ae_serializer_unserialize_int(ae_serializer *serializer, ae_int_t *v, ae_state *state);
    void ae_serializer_unserialize_double(ae_serializer *serializer, double *v, ae_state *state);
    void ae_serializer_stop(ae_serializer *serializer);
    void ae_assert(ae_bool cond, const(char)* msg, ae_state *state);
    ae_int_t ae_cpuid();
    ae_bool ae_fp_eq(double v1, double v2);
    ae_bool ae_fp_neq(double v1, double v2);
    ae_bool ae_fp_less(double v1, double v2);
    ae_bool ae_fp_less_eq(double v1, double v2);
    ae_bool ae_fp_greater(double v1, double v2);
    ae_bool ae_fp_greater_eq(double v1, double v2);
    ae_bool ae_isfinite_stateless(double x, ae_int_t endianness);
    ae_bool ae_isnan_stateless(double x,    ae_int_t endianness);
    ae_bool ae_isinf_stateless(double x,    ae_int_t endianness);
    ae_bool ae_isposinf_stateless(double x, ae_int_t endianness);
    ae_bool ae_isneginf_stateless(double x, ae_int_t endianness);
    ae_int_t ae_get_endianness();
    ae_bool ae_isfinite(double x,ae_state *state);
    ae_bool ae_isnan(double x,   ae_state *state);
    ae_bool ae_isinf(double x,   ae_state *state);
    ae_bool ae_isposinf(double x,ae_state *state);
    ae_bool ae_isneginf(double x,ae_state *state);
    double   ae_fabs(double x,   ae_state *state);
    ae_int_t ae_iabs(ae_int_t x, ae_state *state);
    double   ae_sqr(double x,    ae_state *state);
    double   ae_sqrt(double x,   ae_state *state);
    ae_int_t ae_sign(double x,   ae_state *state);
    ae_int_t ae_round(double x,  ae_state *state);
    ae_int_t ae_trunc(double x,  ae_state *state);
    ae_int_t ae_ifloor(double x, ae_state *state);
    ae_int_t ae_iceil(double x,  ae_state *state);
    ae_int_t ae_maxint(ae_int_t m1, ae_int_t m2, ae_state *state);
    ae_int_t ae_minint(ae_int_t m1, ae_int_t m2, ae_state *state);
    double   ae_maxreal(double m1, double m2, ae_state *state);
    double   ae_minreal(double m1, double m2, ae_state *state);
    double   ae_randomreal(ae_state *state);
    ae_int_t ae_randominteger(ae_int_t maxv, ae_state *state);
    double   ae_sin(double x, ae_state *state);
    double   ae_cos(double x, ae_state *state);
    double   ae_tan(double x, ae_state *state);
    double   ae_sinh(double x, ae_state *state);
    double   ae_cosh(double x, ae_state *state);
    double   ae_tanh(double x, ae_state *state);
    double   ae_asin(double x, ae_state *state);
    double   ae_acos(double x, ae_state *state);
    double   ae_atan(double x, ae_state *state);
    double   ae_atan2(double y, double x, ae_state *state);
    double   ae_log(double x, ae_state *state);
    double   ae_pow(double x, double y, ae_state *state);
    double   ae_exp(double x, ae_state *state);
    ae_complex ae_complex_from_i(ae_int_t v);
    ae_complex ae_complex_from_d(double v);
    ae_complex ae_c_neg(ae_complex lhs);
    ae_bool ae_c_eq(ae_complex lhs,       ae_complex rhs);
    ae_bool ae_c_neq(ae_complex lhs,      ae_complex rhs);
    ae_complex ae_c_add(ae_complex lhs,   ae_complex rhs);
    ae_complex ae_c_mul(ae_complex lhs,   ae_complex rhs);
    ae_complex ae_c_sub(ae_complex lhs,   ae_complex rhs);
    ae_complex ae_c_div(ae_complex lhs,   ae_complex rhs);
    ae_bool ae_c_eq_d(ae_complex lhs,     double rhs);
    ae_bool ae_c_neq_d(ae_complex lhs,    double rhs);
    ae_complex ae_c_add_d(ae_complex lhs, double rhs);
    ae_complex ae_c_mul_d(ae_complex lhs, double rhs);
    ae_complex ae_c_sub_d(ae_complex lhs, double rhs);
    ae_complex ae_c_d_sub(double lhs,     ae_complex rhs);
    ae_complex ae_c_div_d(ae_complex lhs, double rhs);
    ae_complex ae_c_d_div(double lhs,   ae_complex rhs);
    ae_complex ae_c_conj(ae_complex lhs, ae_state *state);
    ae_complex ae_c_sqr(ae_complex lhs, ae_state *state);
    double     ae_c_abs(ae_complex z, ae_state *state);
    ae_complex ae_v_cdotproduct(const ae_complex *v0, ae_int_t stride0, const(char)* conj0, const ae_complex *v1, ae_int_t stride1, const(char)* conj1, ae_int_t n);
    void ae_v_cmove(ae_complex *vdst,    ae_int_t stride_dst, const ae_complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void ae_v_cmoveneg(ae_complex *vdst, ae_int_t stride_dst, const ae_complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void ae_v_cmoved(ae_complex *vdst,   ae_int_t stride_dst, const ae_complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, double alpha);
    void ae_v_cmovec(ae_complex *vdst,   ae_int_t stride_dst, const ae_complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, ae_complex alpha);
    void ae_v_cadd(ae_complex *vdst,     ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void ae_v_caddd(ae_complex *vdst,    ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, double alpha);
    void ae_v_caddc(ae_complex *vdst,    ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, ae_complex alpha);
    void ae_v_csub(ae_complex *vdst,     ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void ae_v_csubd(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, double alpha);
    void ae_v_csubc(ae_complex *vdst, ae_int_t stride_dst, const ae_complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, ae_complex alpha);
    void ae_v_cmuld(ae_complex *vdst, ae_int_t stride_dst, ae_int_t n, double alpha);
    void ae_v_cmulc(ae_complex *vdst, ae_int_t stride_dst, ae_int_t n, ae_complex alpha);
    double ae_v_dotproduct(const double *v0, ae_int_t stride0, const double *v1, ae_int_t stride1, ae_int_t n);
    void ae_v_move(double *vdst,    ae_int_t stride_dst, const double* vsrc,  ae_int_t stride_src, ae_int_t n);
    void ae_v_moveneg(double *vdst, ae_int_t stride_dst, const double* vsrc,  ae_int_t stride_src, ae_int_t n);
    void ae_v_moved(double *vdst,   ae_int_t stride_dst, const double* vsrc,  ae_int_t stride_src, ae_int_t n, double alpha);
    void ae_v_add(double *vdst,     ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n);
    void ae_v_addd(double *vdst,    ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n, double alpha);
    void ae_v_sub(double *vdst,     ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n);
    void ae_v_subd(double *vdst,    ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n, double alpha);
    void ae_v_muld(double *vdst,  ae_int_t stride_dst, ae_int_t n, double alpha);
    ae_int_t ae_v_len(ae_int_t a, ae_int_t b);

    /*
    extern const double ae_machineepsilon;
    extern const double ae_maxrealnumber;
    extern const double ae_minrealnumber;
    extern const double ae_pi;
    */
    enum ae_machineepsilon=5E-16;
    enum ae_maxrealnumber=1E300;
    enum ae_minrealnumber=1E-300;
    enum ae_pi=3.1415926535897932384626433832795;

    struct rcommstate
    {
        int stage;
        ae_vector ia;
        ae_vector ba;
        ae_vector ra;
        ae_vector ca;
    }
    void _rcommstate_init(rcommstate* p, ae_state *_state);
    void _rcommstate_init_copy(rcommstate* dst, rcommstate* src, ae_state *_state);
    void _rcommstate_clear(rcommstate* p);
    void _rcommstate_destroy(rcommstate* p);
    extern(C) __gshared ae_int64_t _alloc_counter;
    extern(C) __gshared ae_bool    _use_alloc_counter;
    version(AE_DEBUG4WINDOWS)
    {
        void flushconsole(T)(T s) { stdout.flush; }
        void tickcount(T)(T s) { _tickcount(); }
        int _tickcount();
    }
    else version(AE_DEBUG4POSIX)
    {
        void flushconsole(T)(T s) { stdout.flush; }
        void tickcount(T)(T s) { _tickcount(); }
        int _tickcount();
    }


    //alglib_impl.ae_int_t ae_int_t;
    class complex;
    ae_int_t vlen(ae_int_t n1, ae_int_t n2);

    class ap_error
    {
        // std.string msg;    
        //ap_error();
        //ap_error(const(char)* s);
        static void make_assertion(bool bClause);
        static void make_assertion(bool bClause, const(char)* p_msg);
    }

}
/**
    class complex
    {
        complex();
        complex(const double &_x);
        complex(const double &_x, const double &_y);
        complex(const complex &z);

        complex& operator= (const double& v);
        complex& operator+=(const double& v);
        complex& operator-=(const double& v);
        complex& operator*=(const double& v);
        complex& operator/=(const double& v);

        complex& operator= (const complex& z);
        complex& operator+=(const complex& z);
        complex& operator-=(const complex& z);
        complex& operator*=(const complex& z);
        complex& operator/=(const complex& z);

        alglib_impl::ae_complex*       c_ptr();
        const alglib_impl::ae_complex* c_ptr() const;
        
        std::string tostring(int dps) const;

        double x, y;
        const alglib::complex operator/(const alglib::complex& lhs, const alglib::complex& rhs);
        const bool operator==(const alglib::complex& lhs, const alglib::complex& rhs);
        const bool operator!=(const alglib::complex& lhs, const alglib::complex& rhs);
        const alglib::complex operator+(const alglib::complex& lhs);
        const alglib::complex operator-(const alglib::complex& lhs);
        const alglib::complex operator+(const alglib::complex& lhs, const alglib::complex& rhs);
        const alglib::complex operator+(const alglib::complex& lhs, const double& rhs);
        const alglib::complex operator+(const double& lhs, const alglib::complex& rhs);
        const alglib::complex operator-(const alglib::complex& lhs, const alglib::complex& rhs);
        const alglib::complex operator-(const alglib::complex& lhs, const double& rhs);
        const alglib::complex operator-(const double& lhs, const alglib::complex& rhs);
        const alglib::complex operator*(const alglib::complex& lhs, const alglib::complex& rhs);
        const alglib::complex operator*(const alglib::complex& lhs, const double& rhs);
        const alglib::complex operator*(const double& lhs, const alglib::complex& rhs);
        const alglib::complex operator/(const alglib::complex& lhs, const alglib::complex& rhs);
        const alglib::complex operator/(const double& lhs, const alglib::complex& rhs);
        const alglib::complex operator/(const alglib::complex& lhs, const double& rhs);
        alglib::complex conj(const alglib::complex &z);
        alglib::complex csqr(const alglib::complex &z);
        alglib::complex vdotproduct(const alglib::complex *v0, ae_int_t stride0, const(char)* conj0, const alglib::complex *v1, ae_int_t stride1, const(char)* conj1, ae_int_t n);
        alglib::complex vdotproduct(const alglib::complex *v1, const alglib::complex *v2, ae_int_t N);
        double abscomplex(const alglib::complex &z);
    }

    void setnworkers(alglib::ae_int_t nworkers);

    double vdotproduct(const double *v0, ae_int_t stride0, const double *v1, ae_int_t stride1, ae_int_t n);
    double vdotproduct(const double *v1, const double *v2, ae_int_t N);

    
    void vmove(double *vdst,  ae_int_t stride_dst, const double* vsrc,  ae_int_t stride_src, ae_int_t n);
    void vmove(double *vdst, const double* vsrc, ae_int_t N);

    void vmove(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void vmove(alglib::complex *vdst, const alglib::complex* vsrc, ae_int_t N);

    void vmoveneg(double *vdst,  ae_int_t stride_dst, const double* vsrc,  ae_int_t stride_src, ae_int_t n);
    void vmoveneg(double *vdst, const double *vsrc, ae_int_t N);

    void vmoveneg(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void vmoveneg(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N);

    void vmove(double *vdst,  ae_int_t stride_dst, const double* vsrc,  ae_int_t stride_src, ae_int_t n, double alpha);
    void vmove(double *vdst, const double *vsrc, ae_int_t N, double alpha);

    void vmove(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, double alpha);
    void vmove(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N, double alpha);

    void vmove(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex* vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, alglib::complex alpha);
    void vmove(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N, alglib::complex alpha);

    void vadd(double *vdst,  ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n);
    void vadd(double *vdst, const double *vsrc, ae_int_t N);

    void vadd(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void vadd(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N);

    void vadd(double *vdst,  ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n, double alpha);
    void vadd(double *vdst, const double *vsrc, ae_int_t N, double alpha);

    void vadd(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, double alpha);
    void vadd(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N, double alpha);

    void vadd(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, alglib::complex alpha);
    void vadd(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N, alglib::complex alpha);

    void vsub(double *vdst,  ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n);
    void vsub(double *vdst, const double *vsrc, ae_int_t N);

    void vsub(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n);
    void vsub(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N);

    void vsub(double *vdst,  ae_int_t stride_dst, const double *vsrc,  ae_int_t stride_src, ae_int_t n, double alpha);
    void vsub(double *vdst, const double *vsrc, ae_int_t N, double alpha);

    void vsub(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, double alpha);
    void vsub(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N, double alpha);

    void vsub(alglib::complex *vdst, ae_int_t stride_dst, const alglib::complex *vsrc, ae_int_t stride_src, const(char)* conj_src, ae_int_t n, alglib::complex alpha);
    void vsub(alglib::complex *vdst, const alglib::complex *vsrc, ae_int_t N, alglib::complex alpha);

    void vmul(double *vdst,  ae_int_t stride_dst, ae_int_t n, double alpha);
    void vmul(double *vdst, ae_int_t N, double alpha);

    void vmul(alglib::complex *vdst, ae_int_t stride_dst, ae_int_t n, double alpha);
    void vmul(alglib::complex *vdst, ae_int_t N, double alpha);

    void vmul(alglib::complex *vdst, ae_int_t stride_dst, ae_int_t n, alglib::complex alpha);
    void vmul(alglib::complex *vdst, ae_int_t N, alglib::complex alpha);

    class ae_vector_wrapper
    {
        ae_vector_wrapper();
        virtual ~ae_vector_wrapper();

        void setlength(ae_int_t iLen);
        ae_int_t length() const;

        void attach_to(alglib_impl::ae_vector *ptr);
        void allocate_own(ae_int_t size, alglib_impl::ae_datatype datatype);
        const alglib_impl::ae_vector* c_ptr() const;
        alglib_impl::ae_vector* c_ptr();
        private:
            ae_vector_wrapper(const ae_vector_wrapper &rhs);
            const ae_vector_wrapper& operator=(const ae_vector_wrapper &rhs);
        protected:
            void create(const ae_vector_wrapper &rhs);
            void create(const(char)* s, alglib_impl::ae_datatype datatype);
            void assign(const ae_vector_wrapper &rhs);      
            alglib_impl::ae_vector *p_vec;
            alglib_impl::ae_vector vec;
    }

    class boolean_1d_array : ae_vector_wrapper
    {
        boolean_1d_array();
        boolean_1d_array(const(char)* s);
        boolean_1d_array(const boolean_1d_array &rhs);
        boolean_1d_array(alglib_impl::ae_vector *p);
        const boolean_1d_array& operator=(const boolean_1d_array &rhs);
        virtual ~boolean_1d_array() ;

        const ae_bool& operator()(ae_int_t i) const;
        ae_bool& operator()(ae_int_t i);

        const ae_bool& operator[](ae_int_t i) const;
        ae_bool& operator[](ae_int_t i);

        void setcontent(ae_int_t iLen, const bool *pContent );
        ae_bool* getcontent();
        const ae_bool* getcontent() const;

        std::string tostring() const;
    }

    class integer_1d_array :  ae_vector_wrapper
    {
        integer_1d_array();
        integer_1d_array(const(char)* s);
        integer_1d_array(const integer_1d_array &rhs);
        integer_1d_array(alglib_impl::ae_vector *p);
        const integer_1d_array& operator=(const integer_1d_array &rhs);
        virtual ~integer_1d_array();

        const ae_int_t& operator()(ae_int_t i) const;
        ae_int_t& operator()(ae_int_t i);

        const ae_int_t& operator[](ae_int_t i) const;
        ae_int_t& operator[](ae_int_t i);

        void setcontent(ae_int_t iLen, const ae_int_t *pContent );

        ae_int_t* getcontent();
        const ae_int_t* getcontent() const;

        std::string tostring() const;
    }

    class real_1d_array : ae_vector_wrapper
    {
        real_1d_array();
        real_1d_array(const(char)* s);
        real_1d_array(const real_1d_array &rhs);
        real_1d_array(alglib_impl::ae_vector *p);
        const real_1d_array& operator=(const real_1d_array &rhs);
        virtual ~real_1d_array();

        const double& operator()(ae_int_t i) const;
        double& operator()(ae_int_t i);

        const double& operator[](ae_int_t i) const;
        double& operator[](ae_int_t i);

        void setcontent(ae_int_t iLen, const double *pContent );
        double* getcontent();
        const double* getcontent() const;

        std::string tostring(int dps) const;
    }

    class complex_1d_array : ae_vector_wrapper
    {
        complex_1d_array();
        complex_1d_array(const(char)* s);
        complex_1d_array(const complex_1d_array &rhs);
        complex_1d_array(alglib_impl::ae_vector *p);
        const complex_1d_array& operator=(const complex_1d_array &rhs);
        virtual ~complex_1d_array();

        const alglib::complex& operator()(ae_int_t i) const;
        alglib::complex& operator()(ae_int_t i);

        const alglib::complex& operator[](ae_int_t i) const;
        alglib::complex& operator[](ae_int_t i);

        void setcontent(ae_int_t iLen, const alglib::complex *pContent );
        alglib::complex* getcontent();
        const alglib::complex* getcontent() const;

        std::string tostring(int dps) const;
    }

    class ae_matrix_wrapper
    {
        ae_matrix_wrapper();
        virtual ~ae_matrix_wrapper();
        const ae_matrix_wrapper& operator=(const ae_matrix_wrapper &rhs);

        void setlength(ae_int_t rows, ae_int_t cols);
        ae_int_t rows() const;
        ae_int_t cols() const;
        bool isempty() const;
    	ae_int_t getstride() const;

        void attach_to(alglib_impl::ae_matrix *ptr);
        void allocate_own(ae_int_t rows, ae_int_t cols, alglib_impl::ae_datatype datatype);
        const alglib_impl::ae_matrix* c_ptr() const;
        alglib_impl::ae_matrix* c_ptr();
        private:
            ae_matrix_wrapper(const ae_matrix_wrapper &rhs);
        protected:
            void create(const ae_matrix_wrapper &rhs);
            void create(const(char)* s, alglib_impl::ae_datatype datatype);
            void assign(const ae_matrix_wrapper &rhs);       
            alglib_impl::ae_matrix *p_mat;
            alglib_impl::ae_matrix mat;
    }

    class boolean_2d_array :  ae_matrix_wrapper
    {
        boolean_2d_array();
        boolean_2d_array(const boolean_2d_array &rhs);
        boolean_2d_array(alglib_impl::ae_matrix *p);
        boolean_2d_array(const(char)* s);
        virtual ~boolean_2d_array();

        const ae_bool& operator()(ae_int_t i, ae_int_t j) const;
        ae_bool& operator()(ae_int_t i, ae_int_t j);

        const ae_bool* operator[](ae_int_t i) const;
        ae_bool* operator[](ae_int_t i);
        
        void setcontent(ae_int_t irows, ae_int_t icols, const bool *pContent );
        
        std::string tostring() const ;
    }

    class integer_2d_array :  ae_matrix_wrapper
    {
        integer_2d_array();
        integer_2d_array(const integer_2d_array &rhs);
        integer_2d_array(alglib_impl::ae_matrix *p);
        integer_2d_array(const(char)* s);
        virtual ~integer_2d_array();

        const ae_int_t& operator()(ae_int_t i, ae_int_t j) const;
        ae_int_t& operator()(ae_int_t i, ae_int_t j);

        const ae_int_t* operator[](ae_int_t i) const;
        ae_int_t* operator[](ae_int_t i);

        void setcontent(ae_int_t irows, ae_int_t icols, const ae_int_t *pContent );
        
        std::string tostring() const;
    }

    class real_2d_array : ae_matrix_wrapper
    {
        real_2d_array();
        real_2d_array(const real_2d_array &rhs);
        real_2d_array(alglib_impl::ae_matrix *p);
        real_2d_array(const(char)* s);
        virtual ~real_2d_array();

        const double& operator()(ae_int_t i, ae_int_t j) const;
        double& operator()(ae_int_t i, ae_int_t j);

        const double* operator[](ae_int_t i) const;
        double* operator[](ae_int_t i);

        void setcontent(ae_int_t irows, ae_int_t icols, const double *pContent );

        std::string tostring(int dps) const;
    }

    class complex_2d_array : ae_matrix_wrapper
    {
        complex_2d_array();
        complex_2d_array(const complex_2d_array &rhs);
        complex_2d_array(alglib_impl::ae_matrix *p);
        complex_2d_array(const(char)* s);
        virtual ~complex_2d_array();

        const alglib::complex& operator()(ae_int_t i, ae_int_t j) const;
        alglib::complex& operator()(ae_int_t i, ae_int_t j);

        const alglib::complex* operator[](ae_int_t i) const;
        alglib::complex* operator[](ae_int_t i);

        void setcontent(ae_int_t irows, ae_int_t icols, const alglib::complex *pContent );

        std::string tostring(int dps) const;
    }

    void read_csv(const(char)* filename, char separator, int flags, alglib::real_2d_array &out);
}
*/

extern(C)
{
    extern __gshared double machineepsilon;
    extern __gshared double maxrealnumber;
    extern __gshared double minrealnumber;
    extern __gshared double fp_nan;
    extern __gshared double fp_posinf;
    extern __gshared double fp_neginf;
    extern __gshared ae_int_t endianness;
//    extern __gshared int CSV_DEFAULT = 0x0;
    //extern __gshared int CSV_SKIP_HEADERS = 0x1;

    int sign(double x);
    double randomreal();
    ae_int_t randominteger(ae_int_t maxv);
    int round(double x);
    int trunc(double x);
    int ifloor(double x);
    int iceil(double x);
    double pi();
    double sqr(double x);
    int maxint(int m1, int m2);
    int minint(int m1, int m2);
    double maxreal(double m1, double m2);
    double minreal(double m1, double m2);

    bool fp_eq(double v1, double v2);
    bool fp_neq(double v1, double v2);
    bool fp_less(double v1, double v2);
    bool fp_less_eq(double v1, double v2);
    bool fp_greater(double v1, double v2);
    bool fp_greater_eq(double v1, double v2);

    bool fp_isnan(double x);
    bool fp_isposinf(double x);
    bool fp_isneginf(double x);
    bool fp_isinf(double x);
    bool fp_isfinite(double x);

    version=ALGLIB_INTERCEPTS_ABLAS;
    void _ialglib_vzero(ae_int_t n, double *p, ae_int_t stride);
    void _ialglib_vzero_complex(ae_int_t n, ae_complex *p, ae_int_t stride);
    void _ialglib_vcopy(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb);
    void _ialglib_vcopy_complex(ae_int_t n, const ae_complex *a, ae_int_t stridea, double *b, ae_int_t strideb, const(char)* conj);
    void _ialglib_vcopy_dcomplex(ae_int_t n, const double *a, ae_int_t stridea, double *b, ae_int_t strideb, const(char)* conj);
    void _ialglib_mcopyblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_int_t stride, double *b);
    void _ialglib_mcopyunblock(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, double *b, ae_int_t stride);
    void _ialglib_mcopyblock_complex(ae_int_t m, ae_int_t n, const ae_complex *a, ae_int_t op, ae_int_t stride, double *b);
    void _ialglib_mcopyunblock_complex(ae_int_t m, ae_int_t n, const double *a, ae_int_t op, ae_complex* b, ae_int_t stride);
    ae_bool _ialglib_i_rmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc);
    ae_bool _ialglib_i_cmatrixgemmf(ae_int_t m, ae_int_t n, ae_int_t k, ae_complex alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, ae_matrix *b, ae_int_t ib, ae_int_t jb, ae_int_t optypeb, ae_complex beta, ae_matrix *c, ae_int_t ic, ae_int_t jc);
    ae_bool _ialglib_i_cmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
    ae_bool _ialglib_i_rmatrixrighttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
    ae_bool _ialglib_i_cmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
    ae_bool _ialglib_i_rmatrixlefttrsmf(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t i1, ae_int_t j1, ae_bool isupper, ae_bool isunit, ae_int_t optype, ae_matrix *x, ae_int_t i2, ae_int_t j2);
    ae_bool _ialglib_i_cmatrixherkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, ae_bool isupper);
    ae_bool _ialglib_i_rmatrixsyrkf(ae_int_t n, ae_int_t k, double alpha, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_int_t optypea, double beta, ae_matrix *c, ae_int_t ic, ae_int_t jc, ae_bool isupper);
    ae_bool _ialglib_i_cmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
    ae_bool _ialglib_i_rmatrixrank1f(ae_int_t m, ae_int_t n, ae_matrix *a, ae_int_t ia, ae_int_t ja, ae_vector *u, ae_int_t uoffs, ae_vector *v, ae_int_t voffs);
}