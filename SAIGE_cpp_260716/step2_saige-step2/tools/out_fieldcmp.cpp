// out_fieldcmp -- field-by-field comparison of two step-2 text outputs that
// have the same rows.
//
//   g++ -std=c++17 -O2 tools/out_fieldcmp.cpp -o tools/out_fieldcmp This is the fp32 round-trip check: `cmp` answers yes/no, and
// under sgsPrecision: fp32 the answer is always no, so what is wanted is how
// many fields and how far off.
//
//   out_fieldcmp ref.txt test.txt
//
// Prints, per column, the number of differing fields and the largest absolute
// and relative gap; for a column named p.value* also the largest gap in
// -log10(p) and the number of values that cross 5e-8.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

struct Col {
    std::string name;
    unsigned long long n=0, diff=0;
    double maxAbs=0, maxRel=0, a64=0, a32=0;
    double maxDlog10=0, p64=0, p32=0;
    unsigned long long cross=0;
    bool isP=false;
};

static bool getline_buf(FILE* f, std::vector<char>& b) {
    b.clear(); int c;
    while ((c=getc_unlocked(f))!=EOF) { if (c=='\n') { b.push_back(0); return true; } b.push_back((char)c); }
    if (b.empty()) return false;
    b.push_back(0); return true;
}
static void split(std::vector<char>& b, std::vector<char*>& out) {
    out.clear(); char* p=b.data(); out.push_back(p);
    for (char* q=p; *q; q++) if (*q=='\t') { *q=0; out.push_back(q+1); }
}

int main(int argc, char** argv) {
    if (argc != 3) { fprintf(stderr,"usage: fieldcmp ref.txt test.txt\n"); return 2; }
    FILE* fa=fopen(argv[1],"rb"); FILE* fb=fopen(argv[2],"rb");
    if (!fa||!fb) { fprintf(stderr,"cannot open inputs\n"); return 1; }
    setvbuf(fa,nullptr,_IOFBF,1<<22); setvbuf(fb,nullptr,_IOFBF,1<<22);
    std::vector<char> la,lb; std::vector<char*> va,vb;
    if (!getline_buf(fa,la)||!getline_buf(fb,lb)) { fprintf(stderr,"empty\n"); return 1; }
    if (strcmp(la.data(),lb.data())!=0) { fprintf(stderr,"HEADER DIFFERS\n"); return 1; }
    split(la,va);
    std::vector<Col> C(va.size());
    for (size_t i=0;i<va.size();i++) { C[i].name=va[i]; C[i].isP = (strncmp(va[i],"p.value",7)==0); }

    unsigned long long rows=0, rowsDiff=0;
    while (true) {
        bool oa=getline_buf(fa,la), ob=getline_buf(fb,lb);
        if (!oa&&!ob) break;
        if (oa!=ob) { fprintf(stderr,"ROW COUNT DIFFERS\n"); return 1; }
        rows++;
        if (strcmp(la.data(),lb.data())==0) continue;
        rowsDiff++;
        std::vector<char> ca(la), cb(lb);
        split(ca,va); split(cb,vb);
        if (va.size()!=vb.size()) { fprintf(stderr,"FIELD COUNT DIFFERS at row %llu\n",rows); return 1; }
        for (size_t i=0;i<va.size() && i<C.size();i++) {
            if (strcmp(va[i],vb[i])==0) continue;
            C[i].diff++;
            const double x=strtod(va[i],nullptr), y=strtod(vb[i],nullptr);
            const double ad=fabs(x-y);
            if (ad>C[i].maxAbs) C[i].maxAbs=ad;
            const double rd = (x!=0.0)? ad/fabs(x) : (ad>0?1e300:0.0);
            if (rd>C[i].maxRel) { C[i].maxRel=rd; C[i].a64=x; C[i].a32=y; }
            if (C[i].isP && x>0 && y>0) {
                const double d=fabs(-log10(x)+log10(y));
                if (d>C[i].maxDlog10) { C[i].maxDlog10=d; C[i].p64=x; C[i].p32=y; }
                if ((x<5e-8)!=(y<5e-8)) C[i].cross++;
            }
        }
    }
    for (size_t i=0;i<C.size();i++) C[i].n=rows;
    printf("rows %llu, rows that differ somewhere %llu (%.4f%%)\n", rows, rowsDiff,
           100.0*(double)rowsDiff/(double)(rows?rows:1));
    printf("%-14s %12s %9s  %-11s %-11s %s\n","column","fields diff","frac","max |dx|","max |dx/x|","worst (ref -> test)");
    unsigned long long tot=0;
    for (size_t i=0;i<C.size();i++) {
        if (C[i].diff==0) { printf("%-14s %12llu %8.4f%%\n", C[i].name.c_str(), 0ULL, 0.0); continue; }
        tot += C[i].diff;
        printf("%-14s %12llu %8.4f%%  %-11.3g %-11.3g %.7g -> %.7g",
               C[i].name.c_str(), C[i].diff, 100.0*(double)C[i].diff/(double)rows,
               C[i].maxAbs, C[i].maxRel, C[i].a64, C[i].a32);
        if (C[i].isP) printf("   max |d(-log10 p)| %.3g, cross 5e-8: %llu", C[i].maxDlog10, C[i].cross);
        printf("\n");
    }
    printf("total differing fields: %llu\n", tot);
    return 0;
}
