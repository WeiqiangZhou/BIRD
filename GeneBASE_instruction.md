####Convert CEL file format
The format of CEL files may need to be converted first using apt-cel-convert from Affymetrix Power Tools (APT) which can be obtained from http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx#1_2
```
apt-cel-convert --format xda --cel-files convert_file_list.txt -i
```
####Run GeneBASE
The Affymetrix tissue panel data is recommended to be combined with the user exon array data, which can be obtained from http://www.affymetrix.com/support/technical/sample_data/exon_array_data.affx
To run GeneBASE:
```
ProbeEffects -par PAR_file.txt
```
####Example of PAR_file.txt
