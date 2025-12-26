#ifndef READER_HPP
#define READER_HPP

#include "model.hpp"
#include <random>

extern std::mt19937_64 rng;
extern bool allowIntra;
const int INF = 0x3f3f3f3f;

int uniform(int l, int r);


PCInstance read_pcinstance();

extern bool allowIntra;

#endif
