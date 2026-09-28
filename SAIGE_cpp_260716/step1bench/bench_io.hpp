// bench_io.hpp -- readers for the two inputs this harness takes.
//
// Deliberately minimal and self-contained: no yaml, no R, no Rcpp. The point
// of the harness is that the only interesting variable is the solver, so the
// input path is kept as short and as obviously correct as possible.
#ifndef STEP1BENCH_BENCH_IO_HPP
#define STEP1BENCH_BENCH_IO_HPP

#include <armadillo>
#include <string>
#include <vector>

namespace step1bench {

// Sparse GRM read from <prefix>.mtx (MatrixMarket coordinate real symmetric,
// 1-based) plus <prefix>.ids (one sample id per line, in matrix order).
//
// `loc` carries BOTH triangles, matching SAIGE's locationMat: gen_sp_Sigma
// walks every stored entry and only touches the diagonal, and BlockSigma's
// partition assumes both triangles are present ("kept exactly as they appear
// in locationMat/valueVec", block_sigma.hpp). An off-diagonal entry (i,j) in
// the file therefore becomes two columns of `loc`.
struct SparseGRM {
    int n = 0;
    arma::umat loc;      // 2 x nnz
    arma::vec  val;      // nnz
    std::vector<std::string> ids;
};

SparseGRM read_sparse_grm(const std::string& prefix);

// Phenotype table: header line, tab or whitespace separated, one id column.
struct PhenoTable {
    std::vector<std::string> header;
    std::vector<std::string> iid;
    arma::mat values;                       // rows x (header.size()-1)
    int colIndex(const std::string& name) const;   // -1 if absent
};

PhenoTable read_pheno(const std::string& path, const std::string& iidCol);

// Align the phenotype rows to GRM order and assemble y and X (intercept first,
// then the named covariates). Throws with a readable message if a GRM sample
// has no phenotype row, or if any used value is non-finite -- subsetting would
// change the GRM's indexing, so it is refused rather than done silently.
void assemble_design(const SparseGRM& grm, const PhenoTable& ph,
                     const std::string& yCol,
                     const std::vector<std::string>& covarCols,
                     arma::fvec& y, arma::fmat& X);

std::vector<std::string> split_commas(const std::string& s);

}  // namespace step1bench

#endif  // STEP1BENCH_BENCH_IO_HPP
