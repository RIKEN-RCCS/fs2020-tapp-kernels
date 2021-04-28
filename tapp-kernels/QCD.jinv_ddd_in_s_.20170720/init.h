#include"report.h"

template<class ForwardIt>
void fill_rand(ForwardIt first, ForwardIt last)
{
  for (; first != last; ++first) {
    *first = get_rand(-0.5, 0.5);
  }
}
