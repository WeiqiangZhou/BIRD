#include "predict_header.h"

int main(int argc, char *argv[])
{
	char infile[255];
	char outfile[255];
	char libfile[255]="./model_file.bin";
	int write_flag = 0;

	if (argc == 2 && strcmp(argv[1], "-h") == 0)
	{
		help_info();
		return 1;
	}
	else if (argc == 5 && strcmp(argv[1], "-i") == 0 && strcmp(argv[3], "-o") == 0 )
	{
		strcpy(infile, argv[2]);
		strcpy(outfile, argv[4]);
	}
	else if (argc == 6 && strcmp(argv[1], "-i") == 0 && strcmp(argv[3], "-o") == 0 && strcmp(argv[5], "-w") == 0)
	{
		strcpy(infile, argv[2]);
		strcpy(outfile, argv[4]);
		write_flag = 1;
	}
	else if (argc == 7 && strcmp(argv[1], "-b") == 0 && strcmp(argv[3], "-i") == 0 && strcmp(argv[5], "-o") == 0)
	{
		strcpy(libfile, argv[2]);
		strcpy(infile, argv[4]);
		strcpy(outfile, argv[6]);
	}
	else if (argc == 8 && strcmp(argv[1], "-b") == 0 && strcmp(argv[3], "-i") == 0 && strcmp(argv[5], "-o") == 0 && strcmp(argv[7], "-w") == 0)
	{
		strcpy(libfile, argv[2]);
		strcpy(infile, argv[4]);
		strcpy(outfile, argv[6]);
		write_flag = 1;
	}
	else
	{
		std::cout << "Please input the correct parameters." << std::endl;
		std::cout << "Standard output: BIRD_predict -b model_file.bin -i input_file.txt -o output_file.txt" << std::endl;
		std::cout << "WIG output: BIRD_predict -b model_file.bin -i input_file.txt -o output_name -w" << std::endl;
		return 1;
	}

	Exondata exonin;
	Exondata *indata = &exonin;
	int predictor_size, sample_size, cluster_size, var_size, loci_size, bin_size;

	if (ReadPar(libfile, loci_size, predictor_size, cluster_size, bin_size, var_size))
	{
		return 1;
	}

	double *quantile_in = new double[predictor_size];
	int *cluster_idx = new int[predictor_size];
	//int *TC_id = new int[predictor_size];
	double *exon_mean = new double[predictor_size];
	double *exon_sd = new double[predictor_size];
	double *DNase_mean = new double[loci_size];
	double *DNase_sd = new double[loci_size];

	double **coef = new double *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		coef[j] = new double[loci_size]();
	}

	int **pre_idx = new int *[var_size];
	for (int j = 0; j < var_size; j++)
	{
		pre_idx[j] = new int[loci_size]();
	}

	char **select_loci = new char *[loci_size];
	for (int j = 0; j < loci_size; j++)
	{
		select_loci[j] = new char[30];
	}

	char **TC_id = new char *[predictor_size];
	for (int j = 0; j < predictor_size; j++)
	{
		TC_id[j] = new char[30];
	}

	if (ReadinModel(libfile, quantile_in, exon_mean, exon_sd, coef, DNase_mean, DNase_sd, pre_idx, TC_id, cluster_idx, select_loci, predictor_size, var_size, loci_size))
	{
		return 1;
	}
	if (ReadinExon(infile, indata))
	{
		return 1;
	}
	std::cout << "Processing data..." << std::endl;
	if (CheckTCid(TC_id, exonin.TC_name,predictor_size))     //check if the input file contains the same predictor as library.
	{
		std::cout << "Exon array data format incorrect: Transcription cluster id is not the same as the library file." << std::endl;
		std::cout << "Please check sample exon file for reference." << std::endl;
		return 1;
	}

	sample_size = (int)exonin.sample_name.size();

	double ** data_norm = new double *[sample_size];
	double ** data_mean = new double *[sample_size];
	double ** output = new double *[sample_size];
	for (int i = 0; i < sample_size; i++)
	{
		data_norm[i] = new double[predictor_size]();
		data_mean[i] = new double[cluster_size]();
		output[i] = new double[loci_size]();
	}

	for (int i = 0; i < sample_size; i++)
	{
		for (int j = 0; j < predictor_size; j++)
		{
			data_norm[i][j] = log2(exonin.data[j][i] + 1);
		}
		QuantileNorm(data_norm[i], quantile_in, predictor_size);    //quantile normalization accroding to quantile from UW exon data.
		
	}

	StandardizeRow(data_norm, exon_mean, exon_sd, predictor_size, sample_size);
	ClusterMean(data_norm, data_mean, cluster_idx, predictor_size, cluster_size, sample_size); //get cluster mean;
	Regression(data_mean, output, coef, pre_idx, var_size, loci_size, sample_size);
	StandardizeRow_r(output, DNase_mean, DNase_sd, loci_size, sample_size);

	std::cout << "Writing output file..." << std::endl;
	if (WriteWIG(output, select_loci, exonin.sample_name, outfile, bin_size, loci_size, sample_size, write_flag))
	{
		return 1;
	}

	//release memory.
	ReleaseExondata(exonin);
	delete[] quantile_in;
	delete[] cluster_idx;
	for (int i = 0; i < sample_size; i++)
	{
		delete[] data_norm[i];
		delete[] data_mean[i];
		delete[] output[i];
	}
	for (int i = 0; i < var_size; i++)
	{
		delete[] coef[i];
		delete[] pre_idx[i];
	}
	for (int i = 0; i < predictor_size; i++)
	{
		delete[] TC_id[i];
	}
	for (int i = 0; i < loci_size; i++)
	{
		delete[] select_loci[i];
	}
	delete[] data_norm;
	delete[] data_mean;
	delete[] output;
	delete[] coef;
	delete[] pre_idx;
	delete[] TC_id;
	delete[] select_loci;
	return 0;
}
