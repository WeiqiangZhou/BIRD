#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <stdlib.h>

struct Exondata
{
	std::vector<std::string> sample_name;
	std::vector<std::string> TC_name;
	std::vector<std::vector<double> > data;
};

void ReleaseExondata(Exondata);
int CheckTCid(char **, std::vector<std::string>, int);
int ReadinExon(char [255], Exondata *);
void ShellSort(double *, int *, int);
void QuantileNorm(double *, double *, int);
void StandardizeRow(double **, double *, double *, int, int);
void StandardizeRow_r(double **, double *, double *, int, int);
void ClusterMean(double **, double **, int *, int, int, int);
int ReadinModel(char [255], double *, double *, double *, double **, double *, double *, int **, char **, int *, char **, int, int, int);
int ReadPar(char[255], int &, int &, int &, int &, int &);
void Regression(double **, double **, double **, int **, int, int, int);
int WriteWIG(double **, char **, std::vector<std::string>, char *, int, int, int, int);
void help_info();
