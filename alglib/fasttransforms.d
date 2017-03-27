module alglib.fasttransforms;
import alglib.ap;
import alglib.internal;

/**
extern(C++,alglib)
{
  void fftc1d(complex_1d_array &a, const ae_int_t n);
  void fftc1d(complex_1d_array &a);
  void fftc1dinv(complex_1d_array &a, const ae_int_t n);
  void fftc1dinv(complex_1d_array &a);
  void fftr1d(const real_1d_array &a, const ae_int_t n, complex_1d_array &f);
  void fftr1d(const real_1d_array &a, complex_1d_array &f);
  void fftr1dinv(const complex_1d_array &f, const ae_int_t n, real_1d_array &a);
  void fftr1dinv(const complex_1d_array &f, real_1d_array &a);
  void convc1d(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
  void convc1dinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
  void convc1dcircular(const complex_1d_array &s, const ae_int_t m, const complex_1d_array &r, const ae_int_t n, complex_1d_array &c);
  void convc1dcircularinv(const complex_1d_array &a, const ae_int_t m, const complex_1d_array &b, const ae_int_t n, complex_1d_array &r);
  void convr1d(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
  void convr1dinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
  void convr1dcircular(const real_1d_array &s, const ae_int_t m, const real_1d_array &r, const ae_int_t n, real_1d_array &c);
  void convr1dcircularinv(const real_1d_array &a, const ae_int_t m, const real_1d_array &b, const ae_int_t n, real_1d_array &r);
  void corrc1d(const complex_1d_array &signal, const ae_int_t n, const complex_1d_array &pattern, const ae_int_t m, complex_1d_array &r);
  void corrc1dcircular(const complex_1d_array &signal, const ae_int_t m, const complex_1d_array &pattern, const ae_int_t n, complex_1d_array &c);
  void corrr1d(const real_1d_array &signal, const ae_int_t n, const real_1d_array &pattern, const ae_int_t m, real_1d_array &r);
  void corrr1dcircular(const real_1d_array &signal, const ae_int_t m, const real_1d_array &pattern, const ae_int_t n, real_1d_array &c);
  void fhtr1d(real_1d_array &a, const ae_int_t n);
  void fhtr1dinv(real_1d_array &a, const ae_int_t n);
}
*/

extern(C)
{
  void fftc1d(/* Complex */ ae_vector* a, ae_int_t n, ae_state *_state);
  void fftc1dinv(/* Complex */ ae_vector* a, ae_int_t n, ae_state *_state);
  void fftr1d(/* Real    */ ae_vector* a, ae_int_t n, /* Complex */ ae_vector* f, ae_state *_state);
  void fftr1dinv(/* Complex */ ae_vector* f, ae_int_t n, /* Real    */ ae_vector* a, ae_state *_state);
  void fftr1dinternaleven(/* Real    */ ae_vector* a, ae_int_t n, /* Real    */ ae_vector* buf, fasttransformplan* plan, ae_state *_state);
  void fftr1dinvinternaleven(/* Real    */ ae_vector* a, ae_int_t n, /* Real    */ ae_vector* buf, fasttransformplan* plan, ae_state *_state);
  void convc1d(/* Complex */ ae_vector* a, ae_int_t m, /* Complex */ ae_vector* b, ae_int_t n, /* Complex */ ae_vector* r, ae_state *_state);
  void convc1dinv(/* Complex */ ae_vector* a, ae_int_t m, /* Complex */ ae_vector* b, ae_int_t n, /* Complex */ ae_vector* r, ae_state *_state);
  void convc1dcircular(/* Complex */ ae_vector* s, ae_int_t m, /* Complex */ ae_vector* r, ae_int_t n, /* Complex */ ae_vector* c, ae_state *_state);
  void convc1dcircularinv(/* Complex */ ae_vector* a, ae_int_t m, /* Complex */ ae_vector* b, ae_int_t n, /* Complex */ ae_vector* r, ae_state *_state);
  void convr1d(/* Real    */ ae_vector* a, ae_int_t m, /* Real    */ ae_vector* b, ae_int_t n, /* Real    */ ae_vector* r, ae_state *_state);
  void convr1dinv(/* Real    */ ae_vector* a, ae_int_t m, /* Real    */ ae_vector* b, ae_int_t n, /* Real    */ ae_vector* r, ae_state *_state);
  void convr1dcircular(/* Real    */ ae_vector* s, ae_int_t m, /* Real    */ ae_vector* r, ae_int_t n, /* Real    */ ae_vector* c, ae_state *_state);
  void convr1dcircularinv(/* Real    */ ae_vector* a, ae_int_t m, /* Real    */ ae_vector* b, ae_int_t n, /* Real    */ ae_vector* r, ae_state *_state);
  void convc1dx(/* Complex */ ae_vector* a, ae_int_t m, /* Complex */ ae_vector* b, ae_int_t n, ae_bool circular, ae_int_t alg, ae_int_t q, /* Complex */ ae_vector* r, ae_state *_state);
  void convr1dx(/* Real    */ ae_vector* a, ae_int_t m, /* Real    */ ae_vector* b, ae_int_t n, ae_bool circular, ae_int_t alg, ae_int_t q, /* Real    */ ae_vector* r, ae_state *_state);
  void corrc1d(/* Complex */ ae_vector* signal, ae_int_t n, /* Complex */ ae_vector* pattern, ae_int_t m, /* Complex */ ae_vector* r, ae_state *_state);
  void corrc1dcircular(/* Complex */ ae_vector* signal, ae_int_t m, /* Complex */ ae_vector* pattern, ae_int_t n, /* Complex */ ae_vector* c, ae_state *_state);
  void corrr1d(/* Real    */ ae_vector* signal, ae_int_t n, /* Real    */ ae_vector* pattern, ae_int_t m, /* Real    */ ae_vector* r, ae_state *_state);
  void corrr1dcircular(/* Real    */ ae_vector* signal, ae_int_t m, /* Real    */ ae_vector* pattern, ae_int_t n, /* Real    */ ae_vector* c, ae_state *_state);
  void fhtr1d(/* Real    */ ae_vector* a, ae_int_t n, ae_state *_state);
  void fhtr1dinv(/* Real    */ ae_vector* a, ae_int_t n, ae_state *_state);
}