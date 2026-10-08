#!/bin/bash
# Multi-thread matrix used to verify the no-R commit (aac50ac1): GPU dense cases at
# nthreads 4/8/8 vs main nthreads=1; CPU dense N=50k at 8 (twice), 4, and main at 8;
# in-place sparse GRM build vs main. Paths are this machine's.
S=/opt/saige/worktrees/step1-rng/SAIGE_cpp_260716/step1_saige-null/tests/no_r/run_vs_ref.sh
C=/home/francisfenglu4/miniforge3/envs/saige-build
BASE=/opt/saige/worktrees/step1-rng-base/SAIGE_cpp_260716/step1_saige-null/saige-null
E="env -i HOME=$HOME PATH=/usr/bin:/bin"
EB="env -i HOME=$HOME PATH=/usr/bin:/bin R_HOME=$C/lib/R LD_LIBRARY_PATH=$C/lib KEEP_R_HOME=1 BIN=$BASE"
G="m_b_full_gpu_mt m_q_full_gpu_mt"
$E NTHR=4 TAG=t4 bash $S $G
$E NTHR=8 TAG=t8a bash $S $G
$E NTHR=8 TAG=t8b bash $S $G
$EB NTHR=8 TAG=main_t8 bash $S $G
$E NTHR=8 TAG=t8b CMPTAG=t8a bash $S m_q_full_cpu
$E NTHR=4 TAG=t4 CMPTAG=t8a bash $S m_q_full_cpu
$EB NTHR=8 TAG=main_t8 CMPTAG=t8a bash $S m_q_full_cpu
$EB NTHR=1 TAG=main_t1 bash $S s_b_sgrm_build
$E NTHR=1 TAG=t1 CMPTAG=main_t1 bash $S s_b_sgrm_build
$E RCPP_PARALLEL_NUM_THREADS=8 NTHR=1 TAG=env8 CMPTAG=main_t1 bash $S s_b_sgrm_build
echo ALLDONE
