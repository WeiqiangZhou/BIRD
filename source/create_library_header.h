#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int ReadParam(char *, int *, char **);
int ReadVector_d(char *, double *);
int ReadVector_i(char *, int *);
int ReadMatrix_d(char *, double **);
int ReadMatrix_i(char *, int **);
int ReadName(char *, char **);
