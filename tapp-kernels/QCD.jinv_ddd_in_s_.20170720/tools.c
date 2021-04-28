#include <stddef.h>
#include "tools.h"

extern unsigned char *malloc_p;
extern size_t         malloc_len;

  int sim_memalign( void **out_p, int    w, size_t size)
  {
    unsigned char *p = malloc_p;
    size_t pad;
    pad =(size_t)p % w;
    p += (w - pad);
    if (malloc_len+size+w-pad > MALLOCAREA ) return 1;
    malloc_len += ( size + (w - pad ) );
    malloc_p += (size + (w - pad ) );
    *out_p = p;
    return 0;
  }

  void *sim_malloc( size_t size )
  {
    unsigned char *p = malloc_p;
    if (malloc_len+size > MALLOCAREA ) return NULL;
    malloc_len += size;
    malloc_p += size;
    return p;
  }

  void *sim_memcpy(void *s1, const void *s2, size_t n)
  {
    char        *p1 = (char *)s1;
    const char  *p2 = (const char *)s2;
  
    while (n-- > 0) {
        *p1 = *p2;
        p1++;
        p2++;
    }
    return (s1);
  }

  void sim_free( void *p )
  {
    return;
  }

  int sim_printf(const char *fmt, ...)
  {
    return 0;
  }

  int sim_exit(int in)
  {
    return 0;
  }
