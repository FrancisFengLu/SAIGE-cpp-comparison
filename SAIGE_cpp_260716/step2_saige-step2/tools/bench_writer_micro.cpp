// Writer microbenchmark: where does "output write" actually go?
// Loads a real step-2 quantitative output file, then re-emits the same rows
// with several strategies, timing format-into-memory and the write() separately.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <fcntl.h>
#include <unistd.h>

static double now(){ using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count(); }

struct Rows {
  std::vector<std::string> chr, pos, mid, ref, alt, pval;
  std::vector<double> ac, af, miss, beta, se, tstat, var;
  std::vector<uint32_t> N;
  size_t n() const { return chr.size(); }
};

static void load(const char* path, Rows& R, size_t maxrows){
  std::ifstream in(path);
  std::string line; std::getline(in, line); // header
  while (R.chr.size() < maxrows && std::getline(in, line)){
    char* s = &line[0]; char* f[14]; int k=0;
    f[k++]=s; for(char* p=s; *p; ++p) if(*p=='\t'){ *p=0; if(k<14) f[k++]=p+1; }
    if(k!=14) continue;
    R.chr.emplace_back(f[0]); R.pos.emplace_back(f[1]); R.mid.emplace_back(f[2]);
    R.ref.emplace_back(f[3]); R.alt.emplace_back(f[4]);
    R.ac.push_back(strtod(f[5],0)); R.af.push_back(strtod(f[6],0));
    R.miss.push_back(strtod(f[7],0)); R.beta.push_back(strtod(f[8],0));
    R.se.push_back(strtod(f[9],0)); R.tstat.push_back(strtod(f[10],0));
    R.var.push_back(strtod(f[11],0)); R.pval.emplace_back(f[12]);
    R.N.push_back((uint32_t)strtoul(f[13],0,10));
  }
}

// ---- strategy 1: exactly what writeOutfile_single does today ----
static void ostream_rows(std::ostream& o, const Rows& R){
  for(size_t k=0;k<R.n();k++){
    o << R.chr[k]; o << "\t"; o << R.pos[k]; o << "\t"; o << R.mid[k]; o << "\t";
    o << R.ref[k]; o << "\t"; o << R.alt[k]; o << "\t";
    o << R.ac[k]; o << "\t"; o << R.af[k]; o << "\t"; o << R.miss[k]; o << "\t";
    o << R.beta[k]; o << "\t"; o << R.se[k]; o << "\t"; o << R.tstat[k]; o << "\t";
    o << R.var[k]; o << "\t"; o << R.pval[k]; o << "\t"; o << R.N[k]; o << "\n";
  }
}

// ---- strategy 2: snprintf into one big buffer ----
static inline char* putd(char* p, double v){ return p + snprintf(p, 32, "%.6g", v); }
static void snprintf_rows(std::string& buf, const Rows& R){
  buf.clear();
  char tmp[512];
  for(size_t k=0;k<R.n();k++){
    char* p = tmp;
    memcpy(p, R.chr[k].data(), R.chr[k].size()); p+=R.chr[k].size(); *p++='\t';
    memcpy(p, R.pos[k].data(), R.pos[k].size()); p+=R.pos[k].size(); *p++='\t';
    memcpy(p, R.mid[k].data(), R.mid[k].size()); p+=R.mid[k].size(); *p++='\t';
    memcpy(p, R.ref[k].data(), R.ref[k].size()); p+=R.ref[k].size(); *p++='\t';
    memcpy(p, R.alt[k].data(), R.alt[k].size()); p+=R.alt[k].size(); *p++='\t';
    p=putd(p,R.ac[k]); *p++='\t'; p=putd(p,R.af[k]); *p++='\t'; p=putd(p,R.miss[k]); *p++='\t';
    p=putd(p,R.beta[k]); *p++='\t'; p=putd(p,R.se[k]); *p++='\t';
    p=putd(p,R.tstat[k]); *p++='\t'; p=putd(p,R.var[k]); *p++='\t';
    memcpy(p, R.pval[k].data(), R.pval[k].size()); p+=R.pval[k].size(); *p++='\t';
    p += snprintf(p, 16, "%u", R.N[k]); *p++='\n';
    buf.append(tmp, p-tmp);
  }
}

// ---- strategy 3: memcpy only (no number formatting at all): the IO+copy floor
static void memcpy_rows(std::string& buf, const std::vector<std::string>& pre){
  buf.clear();
  for(size_t k=0;k<pre.size();k++) buf.append(pre[k]);
}

static double write_all(const char* path, const char* data, size_t n){
  double t0=now();
  int fd = open(path, O_WRONLY|O_CREAT|O_TRUNC, 0644);
  size_t off=0; while(off<n){ ssize_t w=write(fd,data+off,n-off); if(w<=0) break; off+=w; }
  if(strcmp(path,"/dev/null")) fsync(fd);
  close(fd);
  return now()-t0;
}

int main(int argc,char**argv){
  const char* src = argv[1];
  size_t maxrows = argc>2 ? strtoull(argv[2],0,10) : 1000000;
  const char* outdir = argc>3 ? argv[3] : ".";
  Rows R; double t0=now(); load(src,R,maxrows);
  fprintf(stderr,"loaded %zu rows in %.2f s\n", R.n(), now()-t0);

  std::string p1 = std::string(outdir)+"/w1.txt";
  std::string p2 = std::string(outdir)+"/w2.txt";
  std::string p3 = std::string(outdir)+"/w3.txt";

  // 1a ostream -> file (today's writer)
  { double t=now(); std::ofstream o(p1); ostream_rows(o,R); o.flush(); o.close();
    printf("1a ostream -> file        %8.3f s\n", now()-t); }
  // 1b ostream -> /dev/null (formatting only, no disk)
  { double t=now(); std::ofstream o("/dev/null"); ostream_rows(o,R); o.flush();
    printf("1b ostream -> /dev/null   %8.3f s   (= pure CPU of today's writer)\n", now()-t); }
  // 2 snprintf into buffer, then one write
  { std::string buf; buf.reserve(200ull*1024*1024);
    double t=now(); snprintf_rows(buf,R); double tf=now()-t;
    double tw = write_all(p2.c_str(), buf.data(), buf.size());
    double tn = write_all("/dev/null", buf.data(), buf.size());
    printf("2  snprintf format        %8.3f s   write(file) %.3f s   write(/dev/null) %.3f s   bytes %zu\n",
           tf, tw, tn, buf.size()); }
  // 3 memcpy of pre-built rows: the copy+IO floor with zero formatting
  { std::vector<std::string> pre; pre.reserve(R.n());
    { std::ifstream in(p2); std::string l; while(std::getline(in,l)) pre.push_back(l+"\n"); }
    std::string buf; buf.reserve(200ull*1024*1024);
    double t=now(); memcpy_rows(buf,pre); double tf=now()-t;
    double tw = write_all(p3.c_str(), buf.data(), buf.size());
    printf("3  memcpy only (no fmt)   %8.3f s   write(file) %.3f s   bytes %zu\n", tf, tw, buf.size()); }
  // 4 per-field cost of one %.6g
  { double acc=0; char b[32]; double t=now();
    for(size_t k=0;k<R.n();k++){ acc += snprintf(b,32,"%.6g",R.beta[k]); }
    printf("4  snprintf(%%.6g) x %zu  %8.3f s  (%.0f ns/field) [%.0f]\n", R.n(), now()-t,
           1e9*(now()-t)/R.n(), acc); }
  // 5 same through ostream
  { std::ofstream o("/dev/null"); double t=now();
    for(size_t k=0;k<R.n();k++) o << R.beta[k];
    o.flush(); printf("5  ostream<<double x %zu %8.3f s  (%.0f ns/field)\n", R.n(), now()-t,
           1e9*(now()-t)/R.n()); }
  // byte identity of 1a vs 2
  { char cmd[1024]; snprintf(cmd,1024,"cmp %s %s && echo '   [ok] snprintf output is byte-identical to ostream'", p1.c_str(), p2.c_str());
    if(system(cmd)) printf("   [FAIL] snprintf output DIFFERS from ostream\n"); }
  return 0;
}
