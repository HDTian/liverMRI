#!/bin/bash
# =============================================================================
# Author: Haodong Tian
# Description: Computes LD matrices from UKB individual-level genotype data
#              using PLINK2 for use in colocalization analysis.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================
#$ -cwd
#$ -j y
#$ -o /home/user/Project_liver/PWAS_Coloc/logs
#$ -l h_vmem=32G
#$ -l h_rt=4:00:00
#$ -t 1-28

set -x

COLoc_DIR=/home/user/Project_liver/PWAS_Coloc
SHARED_DIR=${COLoc_DIR}/shared_variant_lists
SAMPLE_FILE=/path/to/server/ukb_imp_chr15_v3_s487395.sample  # UKB sample file (all individuals)
KEEP_FILE=${COLoc_DIR}/UKB_MRI_unrelated_EUR.txt              # restrict to unrelated EUR MRI individuals
PLINK2=/home/user/software/PLINK2/plink2
BGENIX=/home/user/Project_liver/PWAS/bgen.tgz/build/apps/bgenix

mkdir -p ${COLoc_DIR}/logs ${COLoc_DIR}/LDmatrix ${COLoc_DIR}/bgi ${COLoc_DIR}/tmp_variant_lists

mapfile -t FILES < <(ls ${SHARED_DIR}/*_shared_variants.txt)
VAR_INFO_FILE=${FILES[$((SGE_TASK_ID-1))]}
prot=$(basename ${VAR_INFO_FILE} _shared_variants.txt)
OUT_PREFIX=${COLoc_DIR}/LDmatrix/${prot}

if [ -f "${OUT_PREFIX}.pvar" ] && [ -f "${OUT_PREFIX}.vcor1" -o -f "${OUT_PREFIX}.unphased.vcor1" ]; then echo "Skipping ${prot}"; exit 0; fi
echo "Processing ${prot}"

chr=$(awk 'NR==2{print $1}' ${VAR_INFO_FILE})
if [ -z "${chr}" ]; then echo "No chr found: ${prot}"; exit 1; fi

# Extract Broad_rsID for exact variant filtering with PLINK2 --extract
VAR_FILE=${COLoc_DIR}/tmp_variant_lists/${prot}_Broad_rsID.txt
awk 'NR>1 && $8!="NA" && $8!="" {print $8}' ${VAR_INFO_FILE} | sort -u > ${VAR_FILE}
if [ ! -s "${VAR_FILE}" ]; then echo "${prot}: no Broad_rsID"; exit 0; fi

# Extract CHR:BP positions for bgenix coarse region extraction
BROAD_ID_FILE=${COLoc_DIR}/tmp_variant_lists/${prot}_Broad_ID.txt
awk 'NR>1 && $1!="NA" && $2!="NA" {print $1":"$2}' ${VAR_INFO_FILE} | sort -u > ${BROAD_ID_FILE}
if [ ! -s "${BROAD_ID_FILE}" ]; then echo "${prot}: no Broad_ID"; exit 0; fi

BGEN_FILE=/path/to/server/ukb/imputed_v3/ukb_imp_chr${chr}_v3.bgen
BGI_FILE=${COLoc_DIR}/bgi/ukb_imp_chr${chr}_v3.bgen.bgi
TMP_BGEN=${COLoc_DIR}/LDmatrix/${prot}.tmp.bgen

if [ ! -f "${BGI_FILE}" ]; then ${BGENIX} -g ${BGEN_FILE} -i ${BGI_FILE} -index; fi

# Extract a temporary BGEN spanning the cis-region using bgenix
# Note: Broad server UKB chromosomes 1-9 are zero-padded (01, 02, ..., 09)
chr_pad=$(printf "%02d" ${chr})
bp_min=$(awk -F':' '{print $2}' ${BROAD_ID_FILE} | sort -n | head -1)
bp_max=$(awk -F':' '{print $2}' ${BROAD_ID_FILE} | sort -n | tail -1)
${BGENIX} -g ${BGEN_FILE} -i ${BGI_FILE} -incl-range ${chr_pad}:${bp_min}-${bp_max} > ${TMP_BGEN}

# Convert temporary BGEN to PLINK2 pfile format (retains .pvar with REF/ALT information)
# ref-first ensures dosage is computed relative to the REF allele, which is required for
# correct allele alignment when harmonizing betas in the colocalization analysis
${PLINK2} \
  --memory 24000 \
  --bgen ${TMP_BGEN} ref-first \
  --sample ${SAMPLE_FILE} \
  --keep ${KEEP_FILE} \
  --extract ${VAR_FILE} \
  --rm-dup force-first \
  --make-pgen \
  --out ${OUT_PREFIX}

# Compute LD correlation matrix from the pfile; ref-based ensures consistent allele direction
${PLINK2} \
  --memory 24000 \
  --pfile ${OUT_PREFIX} \
  --r-unphased square ref-based \
  --out ${OUT_PREFIX}

rm -f ${TMP_BGEN} ${OUT_PREFIX}.pgen ${OUT_PREFIX}.psam  # retain .pvar and LD matrix files
echo "${prot}: done"
