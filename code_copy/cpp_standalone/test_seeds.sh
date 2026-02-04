#!/bin/bash

# Script to test seeds 1-20 for R vs C++ SAIGE comparison
# Results will be saved to results.csv

RESULTS_FILE="/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_13_work/iteration_test/cpp_code/seed_test_results.csv"

# Initialize results file
echo "Seed,R_tau2,Cpp_tau2" > "$RESULTS_FILE"

for SEED in {1..20}; do
    echo "======================================"
    echo "Testing Seed $SEED"
    echo "======================================"

    # Step 1: Modify R code - line 5321 (set_seed)
    echo "Step 1: Modifying R code set_seed to $SEED..."
    sed -i '' "5321s/set_seed([0-9]*)/set_seed($SEED)/" /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/src/SAIGE_fitGLMM_fast.cpp

    # Step 2: Modify R code - line 5398 (output path)
    echo "Step 2: Modifying R code output path to seed$SEED.csv..."
    sed -i '' "5398s/random_vectors_seed[0-9]*.csv/random_vectors_seed$SEED.csv/" /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/src/SAIGE_fitGLMM_fast.cpp
    sed -i '' "5400s/random_vectors_seed[0-9]*.csv/random_vectors_seed$SEED.csv/" /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/src/SAIGE_fitGLMM_fast.cpp

    # Step 3: Recompile R
    echo "Step 3: Recompiling R SAIGE..."
    cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE
    rm -rf /Users/francis/Library/R/arm64/4.4/library/00LOCK-SAIGE
    ~/.pixi/bin/pixi run R CMD INSTALL --preclean . > /tmp/r_install_seed$SEED.log 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: R compilation failed for seed $SEED"
        tail -10 /tmp/r_install_seed$SEED.log
        exit 1
    fi
    echo "R compilation successful"
    tail -3 /tmp/r_install_seed$SEED.log

    # Step 4: Run R test
    echo "Step 4: Running R test..."
    cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE
    R_OUTPUT=$(~/.pixi/bin/pixi run Rscript /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/jan_14_comparison/Sparse_GRM_test_allmarkers/test_r_sparse_nocutoff.R 2>&1)
    echo "$R_OUTPUT" > /tmp/r_test_seed$SEED.log

    # Extract R tau[2] from "Tau:" line (third field after [1])
    R_TAU2=$(echo "$R_OUTPUT" | grep "^Tau:" -A 1 | tail -1 | awk '{print $3}')
    if [ -z "$R_TAU2" ]; then
        echo "ERROR: Could not extract R tau[2] for seed $SEED"
        echo "$R_OUTPUT" | grep -E "Tau:|CONVERGED" | tail -5
        exit 1
    fi
    echo "R tau[2] = $R_TAU2"

    # Step 5: Modify C++ code - line 5427 (csv path)
    echo "Step 5: Modifying C++ code to load seed$SEED.csv..."
    sed -i '' "5427s/random_vectors_seed[0-9]*.csv/random_vectors_seed$SEED.csv/" /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_13_work/iteration_test/cpp_code/SAIGE_step1_fast.cpp

    # Step 6: Recompile C++
    echo "Step 6: Recompiling C++..."
    cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_13_work/iteration_test/cpp_code
    make clean > /dev/null 2>&1
    make > /tmp/cpp_compile_seed$SEED.log 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: C++ compilation failed for seed $SEED"
        tail -10 /tmp/cpp_compile_seed$SEED.log
        exit 1
    fi
    echo "C++ compilation successful"
    tail -2 /tmp/cpp_compile_seed$SEED.log

    # Step 7: Run C++ test
    echo "Step 7: Running C++ test..."
    cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_13_work/iteration_test/cpp_code
    rm -f /tmp/r_Au_vec_*.txt
    CPP_OUTPUT=$(./saige-null -c /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/jan_14_comparison/Sparse_GRM_test_allmarkers/config_sparse_128k_nocov.yaml 2>&1)
    echo "$CPP_OUTPUT" > /tmp/cpp_test_seed$SEED.log

    # Extract C++ tau[2] from last tau_after_update before CONVERGED
    CPP_TAU2=$(echo "$CPP_OUTPUT" | grep "tau_after_update" | tail -1 | grep -oE '[0-9.e+-]+\]' | sed 's/\]//')
    if [ -z "$CPP_TAU2" ]; then
        echo "ERROR: Could not extract C++ tau[2] for seed $SEED"
        echo "$CPP_OUTPUT" | grep -E "tau_after_update|CONVERGED" | tail -5
        exit 1
    fi
    echo "C++ tau[2] = $CPP_TAU2"

    # Save results
    echo "$SEED,$R_TAU2,$CPP_TAU2" >> "$RESULTS_FILE"

    echo "Seed $SEED completed successfully!"
    echo ""
done

echo "======================================"
echo "All seeds tested! Results saved to:"
echo "$RESULTS_FILE"
echo "======================================"
