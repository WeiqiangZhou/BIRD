####Convert CEL file format
The format of CEL files may need to be converted first using apt-cel-convert from Affymetrix Power Tools (APT) which can be obtained from http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx#1_2
```
apt-cel-convert --format xda --cel-files convert_file_list.txt -i
```
####Run GeneBASE
The Affymetrix tissue panel data is recommended to be combined with the user exon array data, which can be obtained from http://www.affymetrix.com/support/technical/sample_data/exon_array_data.affx

The required library files can be obtained from http://www.affymetrix.com/support/technical/byproduct.affx?product=huexon-st

To run GeneBASE:
```
ProbeEffects -par PAR_file.txt
```
####Example of PAR_file.txt
```
[log]
logfile	Log_file.txt

[exon_annotation]
probeset_annotation	HuEx-1_0-st-v2.na33.1.hg19.probeset.csv
pgf_file	HuEx-1_0-st-v2.r2.pgf
clf_file	HuEx-1_0-st-v2.r2.clf
num_cells	6553600
cell_dim	2560

[exon_data]
folder	/CEL_file_directory/
exon_cel_files	User_exon_array.CEL,huex_wta_breast_A.CEL,huex_wta_breast_B.CEL,huex_wta_breast_C.CEL,huex_wta_cerebellum_A.CEL,huex_wta_cerebellum_B.CEL,huex_wta_cerebellum_C.CEL,huex_wta_heart_A.CEL,huex_wta_heart_B.CEL,huex_wta_heart_C.CEL,huex_wta_kidney_A.CEL,huex_wta_kidney_B.CEL,huex_wta_kidney_C.CEL,huex_wta_liver_A.CEL,huex_wta_liver_B.CEL,huex_wta_liver_C.CEL,huex_wta_muscle_A.CEL,huex_wta_muscle_B.CEL,huex_wta_muscle_C.CEL,huex_wta_pancreas_A.CEL,huex_wta_pancreas_B.CEL,huex_wta_pancreas_C.CEL,huex_wta_prostate_A.CEL,huex_wta_prostate_B.CEL,huex_wta_prostate_C.CEL,huex_wta_spleen_A.CEL,huex_wta_spleen_B.CEL,huex_wta_spleen_C.CEL,huex_wta_testes_A.CEL,huex_wta_testes_B.CEL,huex_wta_testes_C.CEL,huex_wta_thyroid_A.CEL,huex_wta_thyroid_B.CEL,huex_wta_thyroid_C.CEL

[output]
output_model_fit	output_file.txt

[model]
array_type	exon
method	mat
summarize_expression	true
summary_method	selection
mat_training_probe_type	background
normalization_method	quantile
```
