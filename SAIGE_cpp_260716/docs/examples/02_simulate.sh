#!/bin/bash
# 02_simulate.sh -- simulated test data for the examples (no real data needed).
#   $WORK/data/geno.{bed,bim,fam}     5,000 samples x 5,000 markers, 1% missing calls,
#                                     markers split over chromosomes 1 and 2
#   $WORK/data/geno.{pgen,pvar,psam}  the same genotypes as a hard-call PGEN
#   $WORK/data/geno.bgen + .sample    the same genotypes as BGEN 1.2 (8-bit)
#   $WORK/data/geno.vcf.gz            the same genotypes as a VCF with a DS field
#   $WORK/data/dosage.{pgen,pvar,psam} 1,000 markers with fractional dosages (same sample IDs)
#   $WORK/data/pheno.txt              IID, 4 binary traits, 2 quantitative traits, 2 covariates
set -euo pipefail
source "$(dirname "$0")/env.sh"
D=$WORK/data
mkdir -p "$D"
cd "$D"

# Genotypes. plink2 --dummy puts every marker on chromosome 1; the awk step moves
# the second half to chromosome 2 so that LOCO has two chromosomes to work with.
$PLINK2 --dummy 5000 5000 0.01 acgt --seed 1 --make-bed --out geno0 > /dev/null
awk 'BEGIN{OFS="\t"} {if (NR > 2500) $1 = 2; $2 = "snp" NR; $4 = 1000 * NR; print}' geno0.bim > geno.bim
mv geno0.bed geno.bed; mv geno0.fam geno.fam; rm -f geno0.*

# The same genotypes in the other formats step 2 reads.
$PLINK2 --bfile geno --make-pgen --out geno > /dev/null
$PLINK2 --bfile geno --export bgen-1.2 bits=8 --out geno > /dev/null
$PLINK2 --bfile geno --export vcf bgz vcf-dosage=DS-force --out geno > /dev/null
# A PGEN with fractional dosages (stays on the CPU in step 2).
$PLINK2 --dummy 5000 1000 0.01 acgt dosage-freq=0.3 --seed 2 --make-pgen --out dosage > /dev/null

# Phenotypes: one row per sample, tab separated, "NA" = missing.
# b1..b4 binary (0 = control, 1 = case), q1 q2 quantitative, x1 x2 covariates.
# b4 is missing for about 10% of the samples, so it has its own sample set.
awk 'BEGIN{srand(7); OFS="\t"; print "IID","b1","b2","b3","b4","q1","q2","x1","x2"}
     function gauss(){ return sqrt(-2*log(1-rand()))*cos(6.283185307*rand()) }
     { x1 = gauss(); x2 = (rand() < 0.5) ? 1 : 0
       b1 = (rand() < 0.10 + 0.03*x1) ? 1 : 0
       b2 = (rand() < 0.30) ? 1 : 0
       b3 = (rand() < 0.05) ? 1 : 0
       b4 = (rand() < 0.10) ? "NA" : ((rand() < 0.20) ? 1 : 0)
       q1 = 0.5*x1 + gauss(); q2 = 0.3*x2 + gauss()
       printf "%s\t%d\t%d\t%d\t%s\t%.4f\t%.4f\t%.4f\t%d\n", $2, b1, b2, b3, b4, q1, q2, x1, x2 }' geno.fam > pheno.txt

ls -l "$D"
head -3 pheno.txt
