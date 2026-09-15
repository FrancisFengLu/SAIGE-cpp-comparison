// Prints the VR candidate marker indices saige-null draws when no VR bypass
// file is present (SAIGE_step1_fast.cpp: mt19937(20200814), 1000 draws of
// uniform_int_distribution<int>(0, M-1), then unique). Built by make_cases.py
// with the same compiler as saige-null so the libstdc++ distribution matches.
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
int main(int argc, char** argv) {
  if (argc < 2) { std::cerr << "usage: vrdraw M [seed]\n"; return 2; }
  const int M = std::atoi(argv[1]);
  const unsigned seed = argc > 2 ? (unsigned)std::atoi(argv[2]) : 20200814u;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(0, M - 1);
  std::set<int> s;
  for (int i = 0; i < 1000; ++i) s.insert(dist(gen));
  for (int x : s) std::cout << x << "\n";
}
