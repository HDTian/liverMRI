#!/bin/bash
# =============================================================================
# Author: Haodong Tian
# Description: Submits PWAS/dosage-based association analysis jobs to the
#              compute cluster.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================
#$ -cwd
#$ -j y
#$ -o /home/user/Project_liver/PWAS/logs
#$ -l h_vmem=32G
#$ -l h_rt=24:00:00
#$ -t 1-22

set -euo pipefail

chr=${SGE_TASK_ID}

PWAS_DIR=/home/user/Project_liver/PWAS
mkdir -p ${PWAS_DIR}/logs
mkdir -p ${PWAS_DIR}/geno_raw
mkdir -p ${PWAS_DIR}/bgi

PLINK2=/path/to/server/tools/plink2
BGENIX=${PWAS_DIR}/bgen.tgz/build/apps/bgenix

SAMPLE_FILE=/path/to/server/ukb_imp_chr15_v3_s487395.sample  # UKB sample file (all individuals)
KEEP_FILE=${PWAS_DIR}/target_samples.txt                      # restrict to individuals of interest

BGEN_FILE=/path/to/server/ukb/imputed_v3/ukb_imp_chr${chr}_v3.bgen  # UKB v3 imputed genotypes split by chromosome
VAR_FILE=${PWAS_DIR}/target_variants_chr${chr}.txt           # target variant IDs split by chromosome

TMP_BGEN=${PWAS_DIR}/geno_raw/tmp_chr${chr}.bgen  # temporary BGEN extracted by bgenix (used instead of full BGEN to avoid OOM)
BGI_FILE=${PWAS_DIR}/bgi/ukb_imp_chr${chr}_v3.bgen.bgi  # BGEN index file for fast positional extraction

OUT_PREFIX=${PWAS_DIR}/geno_raw/pwas_chr${chr}

if [ ! -s "${VAR_FILE}" ]; then
  echo "chr${chr}: skip (no variant file)"
  exit 0
fi

## Build BGEN index (.bgi) if not already present
## The .bgi file is a positional index; it is independent of the target variants
if [ ! -f "${BGI_FILE}" ]; then
    ${BGENIX} -g ${BGEN_FILE} -i ${BGI_FILE} -index
fi

## Extract target variants using bgenix positional ranges, producing a temporary BGEN
## This avoids PLINK2 scanning the full chromosome BGEN, which would cause out-of-memory errors
RANGE_ARGS=""
while read rg; do
    RANGE_ARGS="${RANGE_ARGS} -incl-range ${rg}"
done < <(awk -F: '{chr=$1; if (chr ~ /^[0-9]+$/ && chr < 10) chr=sprintf("%02d", chr); print chr":"$2"-"$2}' ${VAR_FILE} | sort -u)
# Note: Broad server UKB chromosomes 1-9 are zero-padded (01, 02, ..., 09)
eval ${BGENIX} -g ${BGEN_FILE} -i ${BGI_FILE} ${RANGE_ARGS} ">" ${TMP_BGEN}

## Convert temporary BGEN to dosage format (additive coding), retaining only target individuals
${PLINK2} \
  --memory 12000 \
  --bgen ${TMP_BGEN} ref-first \
  --sample ${SAMPLE_FILE} \
  --keep ${KEEP_FILE} \
  --export A \
  --out ${OUT_PREFIX}

rm ${TMP_BGEN}  # remove temporary BGEN after dosage extraction
echo "chr${chr}: done"
