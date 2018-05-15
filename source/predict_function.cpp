#include "predict_header.h"

void ReleaseExondata(Exondata target)
{
	std::vector<std::string>().swap(target.sample_name);
	std::vector<std::string>().swap(target.TC_name);
	std::vector<std::vector<double> >().swap(target.data);
}

//check the input data format
int CheckTCid(char **lib_TC, std::vector<std::string> in_TC, int Length)  
{
	for (int i = 0; i < Length; i++)
	{
		if (strcmp(lib_TC[i], in_TC[i].c_str()) != 0)
		{
			return 1;
		}
	}

	return 0;
}

//read in gene expression data
int ReadinExon(char filename[255], Exondata *indata)  
{
	std::ifstream infile(filename);
	std::string line_temp;
	std::string line_token;

	if (infile.is_open())
	{
		std::cout << "Reading input file " << filename << std::endl;
		std::getline(infile, line_temp);
                if(line_temp[line_temp.length()-1]=='\r') {line_temp[line_temp.length()-1] = '\0';}     //handle windows file
		std::stringstream each_line(line_temp);
		std::getline(each_line, line_token, '\t');          //ignore the first string 'TranscriptCluster'
		while (std::getline(each_line, line_token, '\t'))
		{
			if (line_token.empty()){ break;}
			indata->sample_name.push_back(line_token);
		}

		while (!infile.eof())
		{
			std::getline(infile, line_temp);
                        if(line_temp[line_temp.length()-1]=='\r') {line_temp[line_temp.length()-1] = '\0';}
			if (line_temp.empty()){ break; }
			std::stringstream each_line(line_temp);
			std::getline(each_line, line_token, '\t');
			indata->TC_name.push_back(line_token);
			std::vector<double> data_temp;
			while (std::getline(each_line, line_token, '\t'))
			{
				if (line_token.empty()){ break;}
				data_temp.push_back(atof(line_token.c_str()));
			}
			indata->data.push_back(data_temp);
		}

		infile.close();
	}
	else
	{
		std::cout << "Error! File " << filename << " not found!" << std::endl;
		return 1;
	}

	return 0;

}

//shell sorting algorithm
void ShellSort(double *num, int *index, int numLength)    
{
	int i, j, increment, temp_idx;
	double temp;

	for (increment = numLength / 2; increment > 0; increment /= 2)
	{
		for (i = increment; i<numLength; i++)
		{
			temp = num[i];
			temp_idx = index[i];
			for (j = i; j >= increment; j -= increment)
			{

				if (temp > num[j - increment])
				{
					num[j] = num[j - increment];
					index[j] = index[j - increment];
				}
				else
				{
					break;
				}
			}
			num[j] = temp;
			index[j] = temp_idx;
		}
	}
	return;
}

//quantile normalizaiton with a vector of known quantile
void QuantileNorm(double *indata, double *quantile, int dataLength) 
{
	int *index = new int[dataLength];
	double *indata_copy = new double[dataLength];
	int j;
	for (int i = 0; i < dataLength; i++)
	{
		index[i] = i;
		indata_copy[i] = indata[i];
	}

	ShellSort(indata_copy, index, dataLength);
	indata[index[0]] = quantile[0];
	j = 0;
	for (int i = 1; i < dataLength; i++)
	{
		if (indata_copy[i] == indata_copy[i - 1])
		{
			indata[index[i]] = quantile[j];
		}
		else
		{
			j++;
			indata[index[i]] = quantile[j];
		}
	}

	delete[] index;
	delete[] indata_copy;
	return;
}

//stadardization
void StandardizeRow(double **data_in, double *Mean, double *SD, int Length, int sample_size)  
{
	for (int i = 0; i < Length; i++)
	{
		if (SD[i] != 0)
		{
			for (int j = 0; j < sample_size; j++)
			{
				data_in[j][i] = (data_in[j][i] - Mean[i]) / SD[i];
			}
		}
		else
		{
			for (int j = 0; j < sample_size; j++)
			{
				data_in[j][i] = 0;
			}
		}

	}
	return;
}

//Reverse standardization
void StandardizeRow_r(double **data_in, double *Mean, double *SD, int Length, int sample_size)  
{
	for (int i = 0; i < Length; i++)
	{
		for (int j = 0; j < sample_size; j++)
		{
			data_in[j][i] = data_in[j][i] * SD[i] + Mean[i];
		}

	}
	return;
}

//calculate average value within each gene cluster
void ClusterMean(double **data_matrix, double **data_mean, int *cluster_idx, int p_length, int c_length, int sample_size)
{
	double *data_sum = new double[sample_size];
	int data_count;

	for (int i = 0; i < c_length; i++)
	{
		for (int m = 0; m < sample_size; m++)
		{
			data_sum[m] = 0;
		}
		data_count = 0;
		for (int j = 0; j < p_length; j++)
		{
			if (cluster_idx[j] == i + 1)                      
			{
				for (int k = 0; k < sample_size; k++)
				{
					data_sum[k] += data_matrix[k][j];
				}
				data_count++;
			}
		}
		for (int m = 0; m < sample_size; m++)
		{
			data_mean[m][i] = data_sum[m] / (double)data_count;
		}

	}

	delete[] data_sum;
	return;
}

//read in the pre-built prediction model file
int ReadinModel(char filename[255], double *quantile_in, double *exon_mean, double *exon_sd, double **coef, double *DNase_mean, double *DNase_sd, int **pre_idx, char **TC_id, int *cluster_idx, char **select_loci, int p_length, int var_length, int loci_length, double **dis_matrix, int **DH_cluster, double **DH_coef1, double **DH_coef2, double **DH_coef3, int **DH_pre_idx1, int **DH_pre_idx2, int **DH_pre_idx3, int DH_num1, int DH_num2, int DH_num3)
{
	int dump[8];
	FILE *pFile;
	pFile = fopen(filename, "rb");
	if (pFile != NULL)
	{
		std::cout << "Reading library file " << filename << std::endl;
		fread(dump, sizeof(int), 8, pFile);
		fread(quantile_in, sizeof(double), p_length, pFile);
		fread(exon_mean, sizeof(double), p_length, pFile);
		fread(exon_sd, sizeof(double), p_length, pFile);
		for (int k = 0; k < var_length; k++)
		{
			fread(coef[k], sizeof(double), loci_length, pFile);
		}
		fread(DNase_mean, sizeof(double), loci_length, pFile);
		fread(DNase_sd, sizeof(double), loci_length, pFile);
		for (int k = 0; k < var_length; k++)
		{
			fread(pre_idx[k], sizeof(int), loci_length, pFile);
		}
		fread(cluster_idx, sizeof(int), p_length, pFile);
		for (int k = 0; k < loci_length; k++)
		{
			fread(select_loci[k], sizeof(char), 30, pFile);
		}
		for (int k = 0; k < p_length; k++)
		{
			fread(TC_id[k], sizeof(char), 30, pFile);
		}
		for (int k = 0; k < 3; k++)
		{
			fread(dis_matrix[k], sizeof(double), loci_length, pFile);
		}
		for (int k = 0; k < 3; k++)
		{
			fread(DH_cluster[k], sizeof(int), loci_length, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_coef1[k], sizeof(double), DH_num1, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_coef2[k], sizeof(double), DH_num2, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_coef3[k], sizeof(double), DH_num3, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_pre_idx1[k], sizeof(int), DH_num1, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_pre_idx2[k], sizeof(int), DH_num2, pFile);
		}
		for (int k = 0; k < var_length; k++)
		{
			fread(DH_pre_idx3[k], sizeof(int), DH_num3, pFile);
		}
	}
	else
	{
		std::cout << "Error! File " << filename << " not found!" << std::endl;
		return 1;
	}

	fclose(pFile);
	return 0;
}

//read the parameters used in the prediction model
int ReadPar(char filename[255], int &loci_size, int &predictor_size, int &cluster_size, int &bin_size, int &var_size, int &DH_num1, int &DH_num2, int &DH_num3)
{
	FILE *pFile;
	int param[8];
	pFile = fopen(filename, "rb");
	if (pFile != NULL)
	{
		std::cout << "Reading model parameter..." << std::endl;
		fread(param, sizeof(int), 8, pFile);
		loci_size = param[0];
		predictor_size = param[1];
		cluster_size = param[2];
		bin_size = param[3];
		var_size = param[4];
		DH_num1 = param[5];
		DH_num2 = param[6];
		DH_num3 = param[7];
		fclose(pFile);

	}
	else
	{
		std::cout << "Error! File " << filename << " not found!" << std::endl;
		return 1;
	}
	return 0;
}

//regression according to known coefficients and predictor indexes
void Regression(double **predictor, double **output, double **coef, int **predictor_idx, int var_length, int loci_length, int sample_size)
{
	for (int i = 0; i < loci_length; i++)
	{
		for (int j = 0; j < sample_size; j++)
		{
			output[j][i] = 0;
			for (int k = 0; k < var_length; k++)
			{
				output[j][i] += coef[k][i] * predictor[j][predictor_idx[k][i] - 1];               //index from 1 in library file.
			}

		}
	}
	return;
}

//model average from different level of prediction
void ModelAverage(double **output, double **DH_pre1, double **DH_pre2, double **DH_pre3, double **dis_matrix, int **DH_cluster, int loci_length, int sample_size)
{
	double weight;
	for (int i = 0; i < loci_length; i++)
	{
		for (int j = 0; j < sample_size; j++)
		{
			weight = dis_matrix[0][i] + dis_matrix[1][i] + dis_matrix[2][i] + 1;
			output[j][i] = (output[j][i] + DH_pre1[j][DH_cluster[0][i] - 1] * dis_matrix[0][i] + DH_pre2[j][DH_cluster[1][i] - 1] * dis_matrix[1][i] + DH_pre3[j][DH_cluster[2][i] - 1] * dis_matrix[2][i]) / weight ;               //index from 1 in library file.
		}
	}
	return;
}

//write ouput file in the format of data matrix or wig file
int WriteWIG(double **data_out, char **select_idx, std::vector<std::string> outname, char *outfile, int bin_size, int loci_length, int sample_size, int flag, double up_bound)
{
	char *part;
	char temp_name[255];
	char chrname[10];
	int startsite, count;
	double value;
	FILE * pFile;

	//standard output
	if (flag != 1)
	{
		pFile = fopen(outfile, "w");
		if (pFile != NULL)
		{
			fprintf(pFile, "Chromosome\tStart\tEnd\t");
			for (int i = 0; i < sample_size-1; i++)
			{
				fprintf(pFile, "%s\t", outname[i].c_str());
			}
			fprintf(pFile, "%s\n", outname[sample_size-1].c_str());
			for (int i = 0; i < loci_length; i++)
			{
				fprintf(pFile, "%s\t", select_idx[i]);
				for (int j = 0; j < sample_size-1; j++)
				{
					value = data_out[j][i];
					if (value < 0)
					{
						value = 0;
					}
					else if (value > up_bound)
					{
						value = up_bound;
					}
					fprintf(pFile, "%lf\t", value);   //report original output value, log2(x+1) transformed.
				}

				value = data_out[sample_size-1][i];
				if (value < 0)
				{
					value = 0;
				}
				else if (value > up_bound)
				{
					value = up_bound;
				}
				fprintf(pFile, "%lf\n", value);

			}
			fclose(pFile);
		}
		else
		{
			std::cout << "Error! Cannot write file " << outfile << std::endl;
			return 1;
		}


	}
	else    //WIG output
	{
		std::vector<std::vector<int> > location;
		std::vector<int>::iterator it;
		for (int i = 0; i < 23; i++)
		{
			location.push_back(std::vector<int>());
		}

		for (int i = 0; i < loci_length; i++)
		{
			part = strtok(select_idx[i], "\t\0");
			if (strcmp(part, "chrX") == 0)
			{
				part = strtok(NULL, "\t\0");
				startsite = atoi(part);
				location[22].push_back(startsite);
			}
			else
			{
				for (int j = 0; j < 22; j++)
				{
					sprintf(chrname, "chr%d", j + 1);
					if (strcmp(part, chrname) == 0)
					{
						part = strtok(NULL, "\t\0");
						startsite = atoi(part);
						location[j].push_back(startsite);
						break;
					}
				}
			}
		}

		for (int i = 0; i < sample_size; i++)
		{
			count = 0;
			strcpy(temp_name, outfile);
			strcat(temp_name, ".");
			strcat(temp_name, outname[i].c_str());
			strcat(temp_name, ".wig");
			std::cout << "Writing file " << temp_name << std::endl;
			pFile = fopen(temp_name, "w");
			if (pFile != NULL)
			{
				fprintf(pFile, "track\ttype=wiggle_0\tname=%s\tvisibility=full\tautoScale=off\tmaxHeightPixels=100:50:10\tviewLimits=0.0:100.0\tyLineOnOff=off\n", outname[i].c_str());
				//print chr1-chr22
				for (int j = 0; j < 22; j++)
				{
					fprintf(pFile, "variableStep\tchrom=chr%d\tspan=%d\n", j + 1, bin_size);
					for (it = location[j].begin(); it != location[j].end(); it++)
					{
						value = pow(2, data_out[i][count]) - 1;
						if (value < 0)
						{
							value = 0;
						}
						else if (value > pow(2,up_bound))
						{
							value = pow(2,up_bound);
						}
						count++;
						fprintf(pFile, "%d\t%lf\n", (*it), value);        //report value is limited from 0 to 10000.
					}
				}
				//print chrX
				fprintf(pFile, "variableStep\tchrom=chrX\tspan=%d\n", bin_size);
				for (it = location[22].begin(); it != location[22].end(); it++)
				{
					value = pow(2, data_out[i][count]) - 1;
					if (value < 0)
					{
						value = 0;
					}
					else if (value > pow(2,up_bound))
					{
						value = pow(2,up_bound);
					}
					count++;
					fprintf(pFile, "%d\t%lf\n", (*it), value);
				}

				fclose(pFile);
			}
			else
			{
				std::cout << "Error! Cannot write file " << temp_name << std::endl;
				return 1;
			}

		}
	}
	return 0;
}

void help_info()
{
	std::cout << "Usage:" << std::endl;
	std::cout << "Standard output: BIRD_predict -b model_file.bin -i input_file.txt -o output_file.txt" << std::endl;
	std::cout << "Standard output will save a matrix contained all predited value in log scale (log2(x+1) transformed)." << std::endl;
	std::cout << "WIG output: BIRD_predict -b model_file.bin -i input_file.txt -o output_name -w" << std::endl;
	std::cout << "WIG output will save each sample as a WIG file." << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-b   Specify library file. If not sepecified,the program will search for model_file.bin in the current directory." << std::endl;
	std::cout << "-i   Specify input file (gene expression obtained from GeneBASE)." << std::endl;
	std::cout << "-o   Specify output file." << std::endl;
	std::cout << "-u   Set upper bound for predicted values (default:14)." << std::endl;
	std::cout << "-w   Output WIG file for each sample." << std::endl;
	std::cout << "-l   Use locus-level model for prediction." << std::endl;
	return;
}
