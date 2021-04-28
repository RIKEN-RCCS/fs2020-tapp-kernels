#ifndef TOOLS_H
#define TOOLS_H

#define MALLOCAREA	20000000

#define malloc		sim_malloc
#define posix_memalign	sim_memalign
#define memcpy		sim_memcpy
#define free		sim_free
#define printf		sim_printf
#define exit  		sim_exit

#if defined(__cplusplus)
extern "C" {
#endif

  int sim_memalign( void **out_p, int    w, size_t size);

  void *sim_malloc( size_t size );

  void *sim_memcpy(void *s1, const void *s2, size_t n);

  void sim_free( void *p );

  int sim_printf(const char *fmt, ...);

  int sim_exit(int in);

#ifdef __cplusplus
}
#endif

#endif
