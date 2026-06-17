#!/bin/bash -l
# =============================================================================
# Author: Haodong Tian
# Description: REGENIE step 2 — performs per-chromosome association testing
#              for all 200 liver MRI radiomics features against imputed UKB
#              genotypes, using step 1 predictions as offsets. Runs as a
#              22-job array (one per autosome).
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================
#$ -wd /path/to/server/project/GWAS_regenie
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -pe smp 2 -R y -binding linear:4
#$ -l h_vmem=4G
#$ -o logs/Step2_new.$JOB_ID.out
#$ -e logs/Step2_new.$JOB_ID.err
#$ -t 1-22

source /path/to/server/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2


wdir=/path/to/server/project/GWAS_regenie
mkdir -p ${wdir}/Step2_new

phenofile=${wdir}/my_liver_UKB_more_covar_new.txt
covarfile=${wdir}/my_liver_UKB_more_covar_new.txt
covarCol=sex,age,age2,sex_age,sex_age2,PC{1:10}
catcovar=genotyping_array


chr=${SGE_TASK_ID}

/path/to/server/tools/regenie_v3.2.8.gz_x86_64_Centos7_mkl \
  --step 2 \
  --bgen /path/to/ukb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
  --sample /path/to/server/resources/UKBiobank/v3data/ukb_imp_chr15_v3.sample \
  --extract /path/to/server/project/GWAS_regenie/1000G/1000G_QC.snplist \
  --phenoFile ${phenofile} \
  --phenoColList shape_Elongation,shape_Flatness,shape_LeastAxisLength,shape_MajorAxisLength,shape_Maximum2DDiameterColumn,shape_Maximum2DDiameterRow,shape_Maximum2DDiameterSlice,shape_Maximum3DDiameter,shape_MeshVolume,shape_MinorAxisLength,shape_Sphericity,shape_SurfaceArea,shape_SurfaceVolumeRatio,shape_VoxelVolume,firstorder_10Percentile,firstorder_90Percentile,firstorder_Energy,firstorder_Entropy,firstorder_InterquartileRange,firstorder_Kurtosis,firstorder_Maximum,firstorder_MeanAbsoluteDeviation,firstorder_Mean,firstorder_Median,firstorder_Minimum,firstorder_Range,firstorder_RobustMeanAbsoluteDeviation,firstorder_RootMeanSquared,firstorder_Skewness,firstorder_TotalEnergy,firstorder_Uniformity,firstorder_Variance,glcm_Autocorrelation,glcm_ClusterProminence,glcm_ClusterShade,glcm_ClusterTendency,glcm_Contrast,glcm_Correlation,glcm_DifferenceAverage,glcm_DifferenceEntropy,glcm_DifferenceVariance,glcm_Id,glcm_Idm,glcm_Idmn,glcm_Idn,glcm_Imc1,glcm_Imc2,glcm_InverseVariance,glcm_JointAverage,glcm_JointEnergy,glcm_JointEntropy,glcm_MCC,glcm_MaximumProbability,glcm_SumAverage,glcm_SumEntropy,glcm_SumSquares,gldm_DependenceEntropy,gldm_DependenceNonUniformity,gldm_DependenceNonUniformityNormalized,gldm_DependenceVariance,gldm_GrayLevelNonUniformity,gldm_GrayLevelVariance,gldm_HighGrayLevelEmphasis,gldm_LargeDependenceEmphasis,gldm_LargeDependenceHighGrayLevelEmphasis,gldm_LargeDependenceLowGrayLevelEmphasis,gldm_LowGrayLevelEmphasis,gldm_SmallDependenceEmphasis,gldm_SmallDependenceHighGrayLevelEmphasis,gldm_SmallDependenceLowGrayLevelEmphasis,glrlm_GrayLevelNonUniformity,glrlm_GrayLevelNonUniformityNormalized,glrlm_GrayLevelVariance,glrlm_HighGrayLevelRunEmphasis,glrlm_LongRunEmphasis,glrlm_LongRunHighGrayLevelEmphasis,glrlm_LongRunLowGrayLevelEmphasis,glrlm_LowGrayLevelRunEmphasis,glrlm_RunEntropy,glrlm_RunLengthNonUniformity,glrlm_RunLengthNonUniformityNormalized,glrlm_RunPercentage,glrlm_RunVariance,glrlm_ShortRunEmphasis,glrlm_ShortRunHighGrayLevelEmphasis,glrlm_ShortRunLowGrayLevelEmphasis,glszm_GrayLevelNonUniformity,glszm_GrayLevelNonUniformityNormalized,glszm_GrayLevelVariance,glszm_HighGrayLevelZoneEmphasis,glszm_LargeAreaEmphasis,glszm_LargeAreaHighGrayLevelEmphasis,glszm_LargeAreaLowGrayLevelEmphasis,glszm_LowGrayLevelZoneEmphasis,glszm_SizeZoneNonUniformity,glszm_SizeZoneNonUniformityNormalized,glszm_SmallAreaEmphasis,glszm_SmallAreaHighGrayLevelEmphasis,glszm_SmallAreaLowGrayLevelEmphasis,glszm_ZoneEntropy,glszm_ZonePercentage,glszm_ZoneVariance,ngtdm_Busyness,ngtdm_Coarseness,ngtdm_Complexity,ngtdm_Contrast,ngtdm_Strength,firstorder_10Percentile_inp,firstorder_90Percentile_inp,firstorder_Energy_inp,firstorder_Entropy_inp,firstorder_InterquartileRange_inp,firstorder_Kurtosis_inp,firstorder_Maximum_inp,firstorder_MeanAbsoluteDeviation_inp,firstorder_Mean_inp,firstorder_Median_inp,firstorder_Minimum_inp,firstorder_Range_inp,firstorder_RobustMeanAbsoluteDeviation_inp,firstorder_RootMeanSquared_inp,firstorder_Skewness_inp,firstorder_TotalEnergy_inp,firstorder_Uniformity_inp,firstorder_Variance_inp,glcm_Autocorrelation_inp,glcm_ClusterProminence_inp,glcm_ClusterShade_inp,glcm_ClusterTendency_inp,glcm_Contrast_inp,glcm_Correlation_inp,glcm_DifferenceAverage_inp,glcm_DifferenceEntropy_inp,glcm_DifferenceVariance_inp,glcm_Id_inp,glcm_Idm_inp,glcm_Idmn_inp,glcm_Idc1_inp,glcm_Imc2_inp,glcm_InverseVariance_inp,glcm_JointAverage_inp,glcm_JointEnergy_inp,glcm_JointEntropy_inp,glcm_MCC_inp,glcm_MaximumProbability_inp,glcm_SumAverage_inp,glcm_SumEntropy_inp,glcm_SumSquares_inp,gldm_DependenceEntropy_inp,gldm_DependenceNonUniformity_inp,gldm_DependenceNonUniformityNormalized_inp,gldm_DependenceVariance_inp,gldm_GrayLevelNonUniformity_inp,gldm_GrayLevelVariance_inp,gldm_HighGrayLevelEmphasis_inp,gldm_LargeDependenceEmphasis_inp,gldm_LargeDependenceHighGrayLevelEmphasis_inp,gldm_LargeDependenceLowGrayLevelEmphasis_inp,gldm_LowGrayLevelEmphasis_inp,gldm_SmallDependenceEmphasis_inp,gldm_SmallDependenceHighGrayLevelEmphasis_inp,gldm_SmallDependenceLowGrayLevelEmphasis_inp,glrlm_GrayLevelNonUniformity_inp,glrlm_GrayLevelNonUniformityNormalized_inp,glrlm_GrayLevelVariance_inp,glrlm_HighGrayLevelRunEmphasis_inp,glrlm_LongRunEmphasis_inp,glrlm_LongRunHighGrayLevelEmphasis_inp,glrlm_LongRunLowGrayLevelEmphasis_inp,glrlm_LowGrayLevelRunEmphasis_inp,glrlm_RunEntropy_inp,glrlm_RunLengthNonUniformity_inp,glrlm_RunLengthNonUniformityNormalized_inp,glrlm_RunPercentage_inp,glrlm_RunVariance_inp,glrlm_ShortRunEmphasis_inp,glrlm_ShortRunHighGrayLevelEmphasis_inp,glrlm_ShortRunLowGrayLevelEmphasis_inp,glszm_GrayLevelNonUniformity_inp,glszm_GrayLevelNonUniformityNormalized_inp,glszm_GrayLevelVariance_inp,glszm_HighGrayLevelZoneEmphasis_inp,glszm_LargeAreaEmphasis_inp,glszm_LargeAreaHighGrayLevelEmphasis_inp,glszm_LargeAreaLowGrayLevelEmphasis_inp,glszm_LowGrayLevelZoneEmphasis_inp,glszm_SizeZoneNonUniformity_inp,glszm_SizeZoneNonUniformityNormalized_inp,glszm_SmallAreaEmphasis_inp,glszm_SmallAreaHighGrayLevelEmphasis_inp,glszm_SmallAreaLowGrayLevelEmphasis_inp,glszm_ZoneEntropy_inp,glszm_ZonePercentage_inp,glszm_ZoneVariance_inp,ngtdm_Busyness_inp,ngtdm_Coarseness_inp,ngtdm_Complexity_inp,ngtdm_Contrast_inp,ngtdm_Strength_inp \
  --apply-rint \
  --covarFile ${covarfile} \
  --covarColList ${covarCol} \
  --catCovarList ${catcovar} \
  --pred ${wdir}/Step1_new/liver_pred.list \
  --bsize 200 \
  --out ${wdir}/Step2_new/liver_chr${chr}
