#pragma once

#include<stddef.h>

extern "C" {
  int printf(const char *format, ...);
  void exit(int status);
//void *malloc(size_t size);
  void *my_malloc(size_t size);
  void free(void *ptr);
  int posix_memalign(void **memptr, size_t alignment, size_t size);
//int my_memalign(void **memptr, size_t alignment, size_t size);
  void *memcpy(void *dest, const void *src, size_t n);
}

template<class type>
class complex;

template<class type>
complex<type> operator*(const complex<type>& a, const complex<type>& b);

template<class type>
complex<type> operator/(const complex<type>& a, const complex<type>& b);

template<class type>
class complex {
  type r;
  type i;
public:
  complex(type r, type i)
    : r(r), i(i) {}
  complex(type r)
    : r(r), i(0) {}
  complex(){}
  type real() {
    return r;
  }
  type imag() {
    return i;
  }
  friend complex<type> operator* <>(const complex<type>& a, const complex<type>& b);
  friend complex<type> operator/ <>(const complex<type>& a, const complex<type>& b);
};

template<class type>
complex<type> operator*(const complex<type>& a, const complex<type>& b) {
  return complex<type>(a.r*b.r - a.i*b.i,
                       a.r*b.i + a.i*b.r);
}

template<class type>
complex<type> operator/(const complex<type>& a, const complex<type>& b) {
  type tmp = b.r*b.r + b.i*b.i;
  return complex<type>((a.r*b.r - a.i*b.i)/tmp,
                       (a.r*b.i + a.i*b.r)/tmp);
}
