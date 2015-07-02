####Convert CEL file format
The format of CEL files may need to be converted first using apt-cel-convert from Affymetrix Power Tools (APT) which can be obtained from http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx#1_2
```
apt-cel-convert --format xda --cel-files convert_file_list.txt -i
```
####Run GeneBASE
