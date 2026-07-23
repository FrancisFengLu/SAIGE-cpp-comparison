Step-1 null fitting can run two ways:
  - via the C++ saige-null binary (configs in step1_saige-null/, see project docs), OR
  - via R/Docker SAIGE (wzhou88/saige:1.3.0 step1_fitNULLGLMM.R) when you need a
    full/dense GRM at UKB scale (the C++ dense-GRM path is bugged at large N -> K~=I;
    use SPARSE GRM in C++, or fit the dense-GRM null in R and convert the .rda to the
    C++ .arma null format).
For SAIGE-GENE+ region tests you want a SPARSE-fit / cate-VR null.
