// sgs_fidelity -- what sgsPrecision: fp32 would cost, counted exactly, from an
// fp64 .sgs run.
//
//   g++ -std=c++17 -O2 -fopenmp tools/sgs_fidelity.cpp -o tools/sgs_fidelity
//
// The fp32 file is by construction the fp64 file with every floating column put
// through (float), so the fp64 file alone determines both texts: for each value
// this prints the field the converter would emit from the double and the field
// it would emit from (double)(float)double, and counts where they differ. Only
// rows the text writer actually emits (pval != "NA") are counted.
//
//   sgs_fidelity -m markers.sgs y1.txt.sgs y2.txt.sgs ...
//
// Cross-checked field for field against a real fp32 run at P=8 by
// check_sgs_roundtrip.sh; the two agree exactly, which is what lets the P=128
// numbers in out_fast.hpp be read off the fp64 file instead of 24 GB of text.

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../sgs_format.hpp"

using namespace SAIGE::sgs;

static bool slurp(const std::string& path, std::vector<char>& buf) {
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) { fprintf(stderr, "cannot open %s: %s\n", path.c_str(), strerror(errno)); return false; }
    struct stat st; if (fstat(fd, &st) != 0) { ::close(fd); return false; }
    buf.resize((size_t)st.st_size);
    size_t off = 0;
    while (off < buf.size()) {
        ssize_t r = ::read(fd, buf.data() + off, buf.size() - off);
        if (r <= 0) { ::close(fd); fprintf(stderr, "short read %s\n", path.c_str()); return false; }
        off += (size_t)r;
    }
    ::close(fd); return true;
}

template <typename T>
static void get_col_pod(Reader& r, std::vector<T>& v, size_t n) {
    v.assign(n, T());
    uint8_t e = r.u8();
    if (e == E_CONST) { const uint8_t* q = r.take(sizeof(T)); if (!q) return;
        T x; memcpy(&x, q, sizeof(T)); for (size_t i=0;i<n;i++) v[i]=x; }
    else if (e == E_RAW) { const uint8_t* q = r.take(n*sizeof(T)); if (q) memcpy(v.data(), q, n*sizeof(T)); }
    else r.bad = true;
}
static void get_col_f64(Reader& r, std::vector<double>& v, size_t n) {
    v.assign(n, 0.0);
    uint8_t e = r.u8();
    if (e == E_CONST) { const uint8_t* q=r.take(8); if(!q) return; double x; memcpy(&x,q,8); for(size_t i=0;i<n;i++) v[i]=x; }
    else if (e == E_RAW) { const uint8_t* q=r.take(n*8); if(q) memcpy(v.data(),q,n*8); }
    else if (e == E_CONST32) { const uint8_t* q=r.take(4); if(!q) return; float x; memcpy(&x,q,4); for(size_t i=0;i<n;i++) v[i]=x; }
    else if (e == E_RAW32) { const uint8_t* q=r.take(n*4); if(!q) return; const float* f=(const float*)(const void*)q; for(size_t i=0;i<n;i++) v[i]=f[i]; }
    else r.bad = true;
}
static void get_col_str(Reader& r, std::vector<std::string>& v, size_t n) {
    v.assign(n, std::string());
    uint8_t e = r.u8();
    if (e == E_CONST) { std::string x=r.sstr(); for(size_t i=0;i<n;i++) v[i]=x; }
    else if (e == E_RAW) { for(size_t i=0;i<n;i++) v[i]=r.sstr(); }
    else r.bad = true;
}
// doubles + which rows were exceptions (kept verbatim, so fp32 cannot touch them)
static void get_col_pval(Reader& r, std::vector<double>& d, std::vector<uint8_t>& exc, size_t n) {
    d.assign(n, 0.0); exc.assign(n, 0);
    uint8_t e = r.u8();
    if (e != E_PVAL && e != E_PVAL32) { r.bad = true; return; }
    const size_t w = (e == E_PVAL) ? 8 : 4;
    const uint8_t* q = r.take(n*w); if (!q) return;
    if (w == 8) memcpy(d.data(), q, n*8);
    else { const float* f=(const float*)(const void*)q; for(size_t i=0;i<n;i++) d[i]=f[i]; }
    uint32_t nexc = r.u32();
    for (uint32_t k=0;k<nexc;k++) { uint32_t idx=r.u32(); (void)r.sstr(); if (idx<n) exc[idx]=1; }
}

struct MarkerBlock {
    uint32_t nRows = 0;
    std::vector<std::string> chr,pos,mid,ref,alt;
    std::vector<double> ac,af,info;
};

// --------------------------------------------------------------- stats ----
struct ColStat {
    const char* name = "";
    unsigned long long n = 0, diff = 0;
    double maxAbs = 0, maxRel = 0;          // between the two PRINTED values
    double at64 = 0, at32 = 0;              // the pair that produced maxRel
};
struct Stats {
    ColStat g[8];                           // the %.6g columns
    unsigned long long pN = 0, pDiff = 0, pExc = 0;
    double pMaxDlog10 = 0, pAt64 = 0, pAt32 = 0;
    unsigned long long cross5e8 = 0;
    double pMinSeen = 1.0;
};

static inline void tally(ColStat& C, double v) {
    char a[48], b[48];
    snprintf(a, sizeof a, "%.6g", v);
    snprintf(b, sizeof b, "%.6g", (double)(float)v);
    C.n++;
    if (strcmp(a, b) != 0) {
        C.diff++;
        const double x = strtod(a, nullptr), y = strtod(b, nullptr);
        const double ad = fabs(x - y);
        if (ad > C.maxAbs) C.maxAbs = ad;
        const double rd = (x != 0.0) ? ad / fabs(x) : (ad > 0 ? 1e300 : 0.0);
        if (rd > C.maxRel) { C.maxRel = rd; C.at64 = x; C.at32 = y; }
    }
}

int main(int argc, char** argv) {
    std::string markerPath; std::vector<std::string> inputs; int jobs = 8;
    for (int i=1;i<argc;i++) { std::string a=argv[i];
        if (a=="-m"&&i+1<argc) markerPath=argv[++i];
        else if (a=="-j"&&i+1<argc) jobs=atoi(argv[++i]);
        else inputs.push_back(a); }
    if (inputs.empty()||markerPath.empty()) { fprintf(stderr,"usage: sgsfid -m markers.sgs trait.sgs...\n"); return 2; }

    std::vector<char> mbuf;
    if (!slurp(markerPath, mbuf)) return 1;
    std::vector<MarkerBlock> mblocks;
    {
        Reader r(mbuf.data(), mbuf.size());
        const uint8_t* mg = r.take(8);
        if (!mg||memcmp(mg,MAGIC_MARKER,8)!=0) { fprintf(stderr,"not a marker file\n"); return 1; }
        r.u32(); r.u32();
        while (!r.bad && (size_t)(r.end-r.p)>4) {
            uint32_t peek; memcpy(&peek,r.p,4); if (peek==END_MAGIC) break;
            MarkerBlock B;
            uint32_t m=r.u32(); if (m!=BLK_MAGIC) { fprintf(stderr,"corrupt marker block\n"); return 1; }
            B.nRows=r.u32(); const size_t n=B.nRows;
            get_col_str(r,B.chr,n); get_col_str(r,B.pos,n); get_col_str(r,B.mid,n);
            get_col_str(r,B.ref,n); get_col_str(r,B.alt,n);
            get_col_f64(r,B.ac,n);  get_col_f64(r,B.af,n);  get_col_f64(r,B.info,n);
            if (r.bad) { fprintf(stderr,"corrupt marker block\n"); return 1; }
            mblocks.push_back(B);
        }
    }

    const char* gname[8] = {"AC_Allele2","AF_Allele2","MissingRate/info","BETA","SE","Tstat","var","(other %.6g)"};
    Stats total; for (int i=0;i<8;i++) total.g[i].name = gname[i];

    int rc = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(jobs)
#endif
    for (long ii=0; ii<(long)inputs.size(); ii++) {
        Stats S; for (int i=0;i<8;i++) S.g[i].name = gname[i];
        std::vector<char> tbuf;
        if (!slurp(inputs[(size_t)ii], tbuf)) { rc=1; continue; }
        Reader r(tbuf.data(), tbuf.size());
        const uint8_t* mg = r.take(8);
        if (!mg||memcmp(mg,MAGIC_TRAIT,8)!=0) { fprintf(stderr,"%s not a trait file\n",inputs[ii].c_str()); rc=1; continue; }
        r.u32(); r.u32(); r.str(); r.str(); r.str(); r.u8(); r.u8(); r.str(); r.str();
        uint32_t nc=r.u32(); std::vector<uint8_t> cols(nc);
        for (uint32_t i=0;i<nc;i++) cols[i]=r.u8();

        std::vector<double> Beta,seB,Tst,vT,Bc,sec,Tc,vc,AFca,AFct,Nch,Nche,Ncth,Nctt,ac,af,info;
        std::vector<double> pv, pvNA, pvC, pvNAC;
        std::vector<uint8_t> pvE, pvNAE, pvCE, pvNACE;
        std::vector<char> isSPA; std::vector<uint32_t> Nca,Nct,Nq;

        size_t blk=0;
        while (!r.bad && (size_t)(r.end-r.p)>4) {
            uint32_t peek; memcpy(&peek,r.p,4); if (peek==END_MAGIC) break;
            if (blk>=mblocks.size()) { fprintf(stderr,"block overrun\n"); rc=1; break; }
            const MarkerBlock& MB = mblocks[blk++];
            if (r.u32()!=BLK_MAGIC) { rc=1; break; }
            const size_t n = r.u32();
            const uint32_t flags = r.u32();
            std::vector<uint8_t> present;
            if (flags & F_HAS_PRESENT) { present.assign(n,1); const uint8_t* q=r.take(n); if(q) memcpy(present.data(),q,n); }
            ac=MB.ac; af=MB.af; info=MB.info;
            if (flags & F_OVERRIDE_AC)   get_col_f64(r,ac,n);
            if (flags & F_OVERRIDE_AF)   get_col_f64(r,af,n);
            if (flags & F_OVERRIDE_MISS) get_col_f64(r,info,n);
            for (uint32_t i=0;i<nc;i++) switch (cols[i]) {
                case C_BETA:   get_col_f64(r,Beta,n); break;
                case C_SE:     get_col_f64(r,seB,n);  break;
                case C_TSTAT:  get_col_f64(r,Tst,n);  break;
                case C_VAR:    get_col_f64(r,vT,n);   break;
                case C_PVAL:   get_col_pval(r,pv,pvE,n); break;
                case C_PVALNA: get_col_pval(r,pvNA,pvNAE,n); break;
                case C_ISSPA:  get_col_pod(r,isSPA,n); break;
                case C_BETA_C: get_col_f64(r,Bc,n); break;
                case C_SE_C:   get_col_f64(r,sec,n); break;
                case C_TSTAT_C:get_col_f64(r,Tc,n); break;
                case C_VAR_C:  get_col_f64(r,vc,n); break;
                case C_PVAL_C: get_col_pval(r,pvC,pvCE,n); break;
                case C_PVALNA_C: get_col_pval(r,pvNAC,pvNACE,n); break;
                case C_AFCASE: get_col_f64(r,AFca,n); break;
                case C_AFCTRL: get_col_f64(r,AFct,n); break;
                case C_NCASE:  get_col_pod(r,Nca,n); break;
                case C_NCTRL:  get_col_pod(r,Nct,n); break;
                case C_NCASEHOM: get_col_f64(r,Nch,n); break;
                case C_NCASEHET: get_col_f64(r,Nche,n); break;
                case C_NCTRLHOM: get_col_f64(r,Ncth,n); break;
                case C_NCTRLHET: get_col_f64(r,Nctt,n); break;
                case C_N:      get_col_pod(r,Nq,n); break;
                default: rc=1;
            }
            if (r.bad) { rc=1; break; }

            for (size_t k=0;k<n;k++) {
                if (!present.empty() && !present[k]) continue;   // pval "NA": row not emitted
                tally(S.g[0], ac[k]);  tally(S.g[1], af[k]);  tally(S.g[2], info[k]);
                if (!Beta.empty()) tally(S.g[3], Beta[k]);
                if (!seB.empty())  tally(S.g[4], seB[k]);
                if (!Tst.empty())  tally(S.g[5], Tst[k]);
                if (!vT.empty())   tally(S.g[6], vT[k]);
                if (!AFca.empty()) tally(S.g[7], AFca[k]);
                if (!AFct.empty()) tally(S.g[7], AFct[k]);
                if (!pv.empty()) {
                    S.pN++;
                    if (pvE[k]) { S.pExc++; continue; }
                    char a[64], b[64];
                    snprintf(a,sizeof a,"%.6E", pv[k]);
                    snprintf(b,sizeof b,"%.6E", (double)(float)pv[k]);
                    const double x = strtod(a,nullptr), y = strtod(b,nullptr);
                    if (x < S.pMinSeen) S.pMinSeen = x;
                    if ((x < 5e-8) != (y < 5e-8)) S.cross5e8++;
                    if (strcmp(a,b)!=0) {
                        S.pDiff++;
                        if (x>0 && y>0) {
                            const double d = fabs(-log10(x) + log10(y));
                            if (d > S.pMaxDlog10) { S.pMaxDlog10=d; S.pAt64=x; S.pAt32=y; }
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            for (int i=0;i<8;i++) {
                total.g[i].n += S.g[i].n; total.g[i].diff += S.g[i].diff;
                if (S.g[i].maxAbs > total.g[i].maxAbs) total.g[i].maxAbs = S.g[i].maxAbs;
                if (S.g[i].maxRel > total.g[i].maxRel) { total.g[i].maxRel=S.g[i].maxRel; total.g[i].at64=S.g[i].at64; total.g[i].at32=S.g[i].at32; }
            }
            total.pN += S.pN; total.pDiff += S.pDiff; total.pExc += S.pExc;
            total.cross5e8 += S.cross5e8;
            if (S.pMaxDlog10 > total.pMaxDlog10) { total.pMaxDlog10=S.pMaxDlog10; total.pAt64=S.pAt64; total.pAt32=S.pAt32; }
            if (S.pMinSeen < total.pMinSeen) total.pMinSeen = S.pMinSeen;
        }
    }

    printf("traits: %zu\n", inputs.size());
    printf("%-20s %14s %14s %9s   %-12s %-12s %s\n",
           "column","fields","text differs","frac","max |dx|","max |dx/x|","worst pair (f64 -> f32)");
    unsigned long long gN=0,gD=0;
    for (int i=0;i<8;i++) {
        if (total.g[i].n==0) continue;
        gN += total.g[i].n; gD += total.g[i].diff;
        printf("%-20s %14llu %14llu %8.4f%%   %-12.3g %-12.3g %.7g -> %.7g\n",
               total.g[i].name, total.g[i].n, total.g[i].diff,
               100.0*(double)total.g[i].diff/(double)total.g[i].n,
               total.g[i].maxAbs, total.g[i].maxRel, total.g[i].at64, total.g[i].at32);
    }
    printf("%-20s %14llu %14llu %8.4f%%\n", "ALL %.6g fields", gN, gD, 100.0*(double)gD/(double)gN);
    printf("%-20s %14llu %14llu %8.4f%%   max |d(-log10 p)| = %.3g  at p %.7E -> %.7E\n",
           "p.value (%.6E)", total.pN, total.pDiff,
           100.0*(double)total.pDiff/(double)(total.pN?total.pN:1),
           total.pMaxDlog10, total.pAt64, total.pAt32);
    printf("p-values kept verbatim (exception list, fp32 cannot touch them): %llu\n", total.pExc);
    printf("p-values that cross the 5e-8 threshold: %llu   (smallest p seen %.6E)\n",
           total.cross5e8, total.pMinSeen);
    return rc;
}
