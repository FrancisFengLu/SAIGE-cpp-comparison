#include "bench_io.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace step1bench {

std::vector<std::string> split_commas(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') { if (!cur.empty()) out.push_back(cur); cur.clear(); }
        else if (!std::isspace((unsigned char)c)) cur.push_back(c);
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

SparseGRM read_sparse_grm(const std::string& prefix) {
    SparseGRM g;

    // --- ids ------------------------------------------------------------
    {
        std::ifstream f(prefix + ".ids");
        if (!f) throw std::runtime_error("cannot open " + prefix + ".ids");
        std::string line;
        while (std::getline(f, line)) {
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
            if (!line.empty()) g.ids.push_back(line);
        }
    }

    // --- matrix ---------------------------------------------------------
    std::ifstream f(prefix + ".mtx");
    if (!f) throw std::runtime_error("cannot open " + prefix + ".mtx");
    std::string line;
    bool symmetric = false, sawBanner = false;
    long long nrow = 0, ncol = 0, nnz = 0;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '%') {
            if (!sawBanner) {
                sawBanner = true;
                if (line.find("symmetric") != std::string::npos) symmetric = true;
            }
            continue;
        }
        std::istringstream is(line);
        is >> nrow >> ncol >> nnz;
        break;
    }
    if (nrow <= 0 || nrow != ncol)
        throw std::runtime_error(prefix + ".mtx: bad or non-square dimensions");
    g.n = (int)nrow;
    if (!g.ids.empty() && (long long)g.ids.size() != nrow)
        throw std::runtime_error(prefix + ": .ids has " + std::to_string(g.ids.size()) +
                                 " lines but .mtx is " + std::to_string(nrow) + " x " +
                                 std::to_string(nrow));

    std::vector<arma::uword> ri, ci;
    std::vector<double> vv;
    ri.reserve((size_t)nnz * 2); ci.reserve((size_t)nnz * 2); vv.reserve((size_t)nnz * 2);
    long long read = 0;
    long long i, j; double x;
    while (read < nnz && (f >> i >> j >> x)) {
        read++;
        if (i < 1 || j < 1 || i > nrow || j > nrow)
            throw std::runtime_error(prefix + ".mtx: index out of range at entry " +
                                     std::to_string(read));
        ri.push_back((arma::uword)(i - 1)); ci.push_back((arma::uword)(j - 1)); vv.push_back(x);
        if (symmetric && i != j) {
            ri.push_back((arma::uword)(j - 1)); ci.push_back((arma::uword)(i - 1)); vv.push_back(x);
        }
    }
    if (read != nnz)
        throw std::runtime_error(prefix + ".mtx: header says " + std::to_string(nnz) +
                                 " entries, found " + std::to_string(read));

    g.loc.set_size(2, ri.size());
    for (size_t k = 0; k < ri.size(); k++) { g.loc(0, k) = ri[k]; g.loc(1, k) = ci[k]; }
    g.val = arma::vec(vv.data(), vv.size());
    return g;
}

int PhenoTable::colIndex(const std::string& name) const {
    for (size_t k = 1; k < header.size(); k++)
        if (header[k] == name) return (int)k - 1;
    return -1;
}

PhenoTable read_pheno(const std::string& path, const std::string& iidCol) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);

    auto split_ws = [](const std::string& s) {
        std::vector<std::string> out;
        std::string cur;
        for (char c : s) {
            if (c == '\t' || c == ' ' || c == ',') { out.push_back(cur); cur.clear(); }
            else if (c != '\r') cur.push_back(c);
        }
        out.push_back(cur);
        return out;
    };

    std::string line;
    if (!std::getline(f, line)) throw std::runtime_error(path + ": empty file");
    std::vector<std::string> head = split_ws(line);

    int idc = -1;
    for (size_t k = 0; k < head.size(); k++) if (head[k] == iidCol) { idc = (int)k; break; }
    if (idc < 0) throw std::runtime_error(path + ": no column named '" + iidCol + "'");

    PhenoTable t;
    t.header.push_back(iidCol);
    for (size_t k = 0; k < head.size(); k++) if ((int)k != idc) t.header.push_back(head[k]);
    const size_t ncol = t.header.size() - 1;

    std::vector<double> flat;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> tok = split_ws(line);
        if (tok.size() != head.size()) continue;      // ragged line: skip
        t.iid.push_back(tok[(size_t)idc]);
        for (size_t k = 0; k < tok.size(); k++) {
            if ((int)k == idc) continue;
            const std::string& s = tok[k];
            double v;
            if (s.empty() || s == "NA" || s == "NaN" || s == ".")
                v = std::numeric_limits<double>::quiet_NaN();
            else {
                char* end = nullptr;
                v = std::strtod(s.c_str(), &end);
                if (end == s.c_str()) v = std::numeric_limits<double>::quiet_NaN();
            }
            flat.push_back(v);
        }
    }
    const size_t nrow = t.iid.size();
    t.values.set_size(nrow, ncol);
    for (size_t r = 0; r < nrow; r++)
        for (size_t c = 0; c < ncol; c++) t.values(r, c) = flat[r * ncol + c];
    return t;
}

void assemble_design(const SparseGRM& grm, const PhenoTable& ph,
                     const std::string& yCol,
                     const std::vector<std::string>& covarCols,
                     arma::fvec& y, arma::fmat& X) {
    const int yc = ph.colIndex(yCol);
    if (yc < 0) throw std::runtime_error("phenotype column '" + yCol + "' not found");
    std::vector<int> cc;
    for (const auto& c : covarCols) {
        const int k = ph.colIndex(c);
        if (k < 0) throw std::runtime_error("covariate column '" + c + "' not found");
        cc.push_back(k);
    }

    std::unordered_map<std::string, size_t> pos;
    pos.reserve(ph.iid.size() * 2);
    for (size_t r = 0; r < ph.iid.size(); r++) pos.emplace(ph.iid[r], r);

    const int n = grm.n;
    y.set_size((arma::uword)n);
    X.set_size((arma::uword)n, (arma::uword)(1 + cc.size()));

    for (int i = 0; i < n; i++) {
        const std::string& id = grm.ids.empty() ? std::string() : grm.ids[(size_t)i];
        auto it = pos.find(id);
        if (it == pos.end())
            throw std::runtime_error("GRM sample '" + id + "' has no phenotype row. "
                "Subsetting is refused: it would renumber the GRM and the three solvers "
                "would no longer be comparing the same matrix.");
        const size_t r = it->second;
        const double yv = ph.values(r, (arma::uword)yc);
        if (!std::isfinite(yv))
            throw std::runtime_error("sample '" + id + "' has a non-finite value in '" +
                                     yCol + "'; missing phenotypes are not handled.");
        y((arma::uword)i) = (float)yv;
        X((arma::uword)i, 0) = 1.0f;
        for (size_t k = 0; k < cc.size(); k++) {
            const double cv = ph.values(r, (arma::uword)cc[k]);
            if (!std::isfinite(cv))
                throw std::runtime_error("sample '" + id + "' has a non-finite covariate in '" +
                                         covarCols[k] + "'.");
            X((arma::uword)i, (arma::uword)(k + 1)) = (float)cv;
        }
    }
}

}  // namespace step1bench
