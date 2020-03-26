#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <unordered_map>

struct Exondata
{
	std::vector<std::string> sample_name;
	std::vector<std::string> TC_name;
	std::vector<std::vector<double> > data;
};

void ReleaseExondata(Exondata);
int CheckTCid(char **, std::vector<std::string>, int);
int MatchExon_sep(char **, Exondata *, int, int *, std::string);
int MatchExon(char **, Exondata *, int, int *);
int ReadinExon(char [255], Exondata *);
int WriteExpr(double **, std::vector<std::string>, char **, char *, int, int);
void ShellSort(double *, int *, int);
void GetRank(double *, double *,int);
void QuantileNorm(double *, double *, int);
void StandardizeRow(double **, double *, double *, int, int);
void StandardizeRow_r(double **, double *, double *, int, int);
void ClusterMean(double **, double **, int *, int, int, int);
int ReadinModel(char[255], double *, double *, double *, double **, double *, double *, int **, char **, int *, char **, int, int, int, double **, int **, double **, double **, double **, int **, int **, int **, int, int, int);
int ReadPar(char[255], int &, int &, int &, int &, int &, int &, int &, int &);
void Regression(double **, double **, double **, int **, int, int, int);
void ModelAverage(double **, double **, double **, double **, double **, int **, int, int);
int WriteWIG(double **, char **, std::vector<std::string>, char *, int, int, int, int, double, int);
void help_info();
