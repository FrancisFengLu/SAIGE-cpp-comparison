// firth_ref.cpp -- CPU reference for the Firth-on-GPU prototype.
//
// For every (marker, trait) pair in pairs.tsv it rebuilds the Firth inputs the
// way saige_test.cpp does (GVec from the .bed through the per-marker
// code->dosage table, gtilde = G - XXVX_inv * (XV * G) over the carriers only,
// exactly getadjGFast), then runs three fits on the same gtilde / y / offset:
//   legacy : fast_logistf_fit_simple copied verbatim from saige_test.cpp
//            (maxit 50, maxstep 15), plus iteration / singular diagnostics
//            appended after the loop; evaluation order untouched.
//   m15    : collaborator fit_two_parameter_moment with max_step 15
//            (the moment form of the legacy iteration, no line search since
//            line_search_trigger is left at 1.0 -- see firth_moment.hpp).
//   m1     : the same with max_step 1 (what the collaborator package ships).
// Build: see Makefile in this directory. Single-threaded by design.
#include <armadillo>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "firth_moment.hpp"

using SAIGE::firth::FitOptions; using SAIGE::firth::FitResult; using SAIGE::firth::FitFailure;

// ---- verbatim from saige_test.cpp (SAIGEClass::fast_logistf_fit_simple), made a
// free function; the only additions are the three diagnostic out-params at the end.
static void fast_logistf_fit_simple(arma::mat & x, arma::vec & y, arma::vec & offset, bool firth,
        arma::vec init, int maxit, int maxstep, int maxhs, double lconv, double gconv, double xconv,
        double & beta_G, double & sebeta_G, bool & isfirthconverge,
        int & d_iter, bool & d_singular, double & d_alpha, bool & d_strict){
  isfirthconverge = false;
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  thread_local arma::vec ws_pi_0, ws_pi, ws_wpi, ws_W2, ws_h, ws_ypih, ws_xcol;
  thread_local arma::mat ws_XW2, ws_Q, ws_XX_XW2;
  ws_pi_0 = -x * beta - offset;
  ws_pi_0 = arma::exp(ws_pi_0) + 1;
  ws_pi = 1/ws_pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  d_singular = false; d_strict = false;
  while(iter <= maxit){
        beta_old = beta;
        ws_wpi = ws_pi % (1 - ws_pi);
        ws_W2 = arma::sqrt(ws_wpi);
        ws_XW2.set_size(n, k);
        for(int j = 0; j < k; j++){
                ws_XW2.col(j) = x.col(j) % ws_W2;
        }
        arma::mat R;
        arma::qr_econ(ws_Q, R, ws_XW2);
        ws_h = ws_Q % ws_Q * oneVec;
        arma::vec U_star(2, arma::fill::zeros);
        if(firth){
                ws_ypih = (y - ws_pi) + (ws_h % (0.5 - ws_pi));
        }else{
                ws_ypih = (y - ws_pi);
        }
        U_star = x.t() * ws_ypih;
        ws_XX_XW2.set_size(n, k);
        for(int j = 0; j < k; j++){
                ws_xcol = x.col(j);
                ws_XX_XW2.col(j) = ws_xcol % ws_W2;
        }
        arma::mat XX_Fisher = ws_XX_XW2.t() * (ws_XX_XW2);
        bool isinv = arma::inv_sympd (XX_covs, XX_Fisher);
	if(!isinv){
                d_singular = true;
                break;
        }
        arma::vec delta = XX_covs * U_star;
        double mx = arma::max(arma::abs(delta))/maxstep;
        if(mx > 1){
                delta = delta/mx;
        }
        evals = evals + 1;
        iter = iter + 1;
        beta = beta + delta;
        ws_pi_0 = -x * beta - offset;
        ws_pi_0 = arma::exp(ws_pi_0) + 1;
        ws_pi = 1/ws_pi_0;
        if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		isfirthconverge = true;
                d_strict = (arma::max(arma::abs(delta)) <= xconv) && (abs(U_star).is_zero(gconv));
                break;
        }
  }
        d_iter = iter; d_alpha = beta(0);
        arma::mat var;
        if(XX_covs.has_nan()){
                var = XX_covs;
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = (XX_covs.n_elem == 4) ? sqrt(XX_covs(1,1)) : arma::datum::nan; // guard: inv_sympd failure resets XX_covs to 0x0
        }
}
// ---- end verbatim

static std::vector<double> load_npy_f64(const std::string& path, size_t& rows, size_t& cols) {
    std::ifstream f(path, std::ios::binary); if (!f) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(1); }
    char magic[6]; f.read(magic, 6); unsigned char ver[2]; f.read((char*)ver, 2);
    uint32_t hlen = 0; if (ver[0] == 1) { uint16_t h; f.read((char*)&h, 2); hlen = h; } else { f.read((char*)&hlen, 4); }
    std::string hdr(hlen, ' '); f.read(&hdr[0], hlen);
    if (hdr.find("'<f8'") == std::string::npos || hdr.find("'fortran_order': False") == std::string::npos) { fprintf(stderr, "npy: need C-order float64\n"); exit(1); }
    size_t a = hdr.find("'shape': (") + 10; size_t b = hdr.find(")", a); std::string sh = hdr.substr(a, b - a);
    rows = cols = 1; std::stringstream ss(sh); std::string tok; std::vector<size_t> dims;
    while (std::getline(ss, tok, ',')) { if (tok.find_first_of("0123456789") != std::string::npos) dims.push_back(std::stoul(tok)); }
    rows = dims[0]; cols = dims.size() > 1 ? dims[1] : 1;
    std::vector<double> v(rows * cols); f.read((char*)v.data(), v.size() * 8); return v;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: firth_ref <bed> <N> <pairdir> [out.tsv] [maxpairs]\n"); return 1; }
    std::string bedp = argv[1]; int N = atoi(argv[2]); std::string dir = argv[3];
    std::string outp = argc > 4 ? argv[4] : dir + "/cpu_ref.tsv"; long maxpairs = argc > 5 ? atol(argv[5]) : -1;
    const size_t B = (N + 3) / 4;
    std::ifstream bed(bedp, std::ios::binary); if (!bed) { fprintf(stderr, "no bed\n"); return 1; }
    size_t lr, lc; std::vector<double> lut = load_npy_f64(dir + "/lut.npy", lr, lc);
    // traits
    struct Trait { std::string name, mdir; arma::vec y, off; arma::mat XV, XXVX_inv; };
    std::vector<Trait> traits;
    { std::ifstream tf(dir + "/traits.tsv"); std::string l;
      while (std::getline(tf, l)) { std::stringstream ss(l); Trait t; int i; ss >> i >> t.name >> t.mdir;
        t.y.load(t.mdir + "/y.arma", arma::arma_binary); t.off.load(t.mdir + "/offset.arma", arma::arma_binary);
        t.XV.load(t.mdir + "/XV.arma", arma::arma_binary); t.XXVX_inv.load(t.mdir + "/XXVX_inv.arma", arma::arma_binary);
        if ((int)t.y.n_elem != N) { fprintf(stderr, "trait %s: y has %llu elems, N=%d\n", t.name.c_str(), (unsigned long long)t.y.n_elem, N); return 1; }
        traits.push_back(std::move(t)); } }
    FILE* out = fopen(outp.c_str(), "w");
    fprintf(out, "pair\ttrait\tmarker\tnnz\tbeta_leg\tse_leg\tflag_leg\titer_leg\tsing_leg\tstrict_leg\talpha_leg\t"
                 "beta_m15\tse_m15\tflag_m15\tconv_m15\titer_m15\tfail_m15\t"
                 "beta_m1\tse_m1\tflag_m1\tconv_m1\titer_m1\tfail_m1\tus_leg\tus_m15\tus_m1\n");
    std::ifstream pf(dir + "/pairs.tsv"); std::string line; std::getline(pf, line);
    arma::vec G(N), gt; std::vector<unsigned char> buf(B); long npairs = 0;
    double tot_leg = 0, tot_m15 = 0, tot_m1 = 0;
    while (std::getline(pf, line)) {
        if (maxpairs >= 0 && npairs >= maxpairs) break;
        std::stringstream ss(line); long pid; int ti; long m; ss >> pid >> ti >> m;
        bed.seekg(3 + (std::streamoff)m * B); bed.read((char*)buf.data(), B);
        arma::uvec iIndex; std::vector<arma::uword> nz; nz.reserve(N);
        for (int i = 0; i < N; i++) { int code = (buf[i >> 2] >> (2 * (i & 3))) & 3; G[i] = lut[m * 4 + code]; if (G[i] != 0) nz.push_back(i); }
        iIndex = arma::uvec(nz);
        Trait& T = traits[ti];
        // getadjGFast verbatim
        arma::vec m_XVG(T.XV.n_rows, arma::fill::zeros);
        for (unsigned int i = 0; i < iIndex.n_elem; i++) m_XVG += T.XV.col(iIndex(i)) * G(iIndex(i));
        gt = G - T.XXVX_inv * m_XVG;
        // legacy
        arma::mat x(N, 2); x.col(0).ones(); x.col(1) = gt; arma::vec init(2, arma::fill::zeros);
        double bl, sl; bool fl; int il; bool sg, st; double al;
        auto t0 = std::chrono::steady_clock::now();
        fast_logistf_fit_simple(x, T.y, T.off, true, init, 50, 15, 15, 1e-5, 1e-5, 1e-5, bl, sl, fl, il, sg, al, st);
        auto t1 = std::chrono::steady_clock::now();
        FitOptions o15; o15.max_iterations = 50; o15.max_step = 15.0; o15.gradient_tolerance = 1e-5; o15.parameter_tolerance = 1e-5;
        FitResult r15 = SAIGE::firth::fit_two_parameter_moment(gt, T.y, T.off, true, init, o15);
        auto t2 = std::chrono::steady_clock::now();
        FitOptions o1 = o15; o1.max_step = 1.0;
        FitResult r1 = SAIGE::firth::fit_two_parameter_moment(gt, T.y, T.off, true, init, o1);
        auto t3 = std::chrono::steady_clock::now();
        double us1 = std::chrono::duration<double, std::micro>(t1 - t0).count(), us2 = std::chrono::duration<double, std::micro>(t2 - t1).count(), us3 = std::chrono::duration<double, std::micro>(t3 - t2).count();
        tot_leg += us1; tot_m15 += us2; tot_m1 += us3;
        fprintf(out, "%ld\t%d\t%ld\t%llu\t%.17g\t%.17g\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\n",
            pid, ti, m, (unsigned long long)iIndex.n_elem, bl, sl, (int)fl, il, (int)sg, (int)st, al,
            r15.beta, r15.se_beta, (int)(r15.converged || r15.hit_iteration_limit), (int)r15.converged, r15.iterations, (int)r15.failure,
            r1.beta, r1.se_beta, (int)(r1.converged || r1.hit_iteration_limit), (int)r1.converged, r1.iterations, (int)r1.failure, us1, us2, us3);
        npairs++;
    }
    fclose(out);
    fprintf(stderr, "%ld pairs; mean per fit: legacy %.1f us, moment15 %.1f us, moment1 %.1f us\n", npairs, tot_leg / npairs, tot_m15 / npairs, tot_m1 / npairs);
    return 0;
}
