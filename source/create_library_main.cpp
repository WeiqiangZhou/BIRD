#include "create_library_header.h" 

int main(int argc, char *argv[])
{
	char infile[255]="./par_file.txt";
	char outfile[255]="./model_file.bin";

	if (argc == 5 && strcmp(argv[1], "-i") == 0 && strcmp(argv[3], "-o") == 0)
	{
		strcpy(infile, argv[2]);
		strcpy(outfile, argv[4]);
	}
        else if (argc == 2 && strcmp(argv[1], "-h") == 0)
        {
                std::cout << "Usage: BIRD_build_library -i par_file.txt -o model_file.bin" << std::endl;
                std::cout << "par_file.txt: contains the required parameters and files list." << std::endl;
                std::cout << "model_file.bin: the model file to be used by BIRD." << std::endl;
                return 1;
        }
	else
	{
		std::cout << "Please input the correct parameters." << std::endl;
		std::cout << "Example: BIRD_build_library -i par_file.txt -o model_file.bin" << std::endl;
		return 1;
	}

	FILE * p_out;
	int param[5];
	int num_loci, num_gene, cluster_size, bin_size, num_var;
	char **filelist = new char *[10];
	for (int i = 0; i < 10; i++)
	{
		filelist[i] = new char[255];
	}

	if (ReadParam(infile, param, filelist)){ return 1; }

	num_loci = param[0];
	num_gene = param[1];
	cluster_size = param[2];
	bin_size = param[3];
	num_var = param[4];

	double *uw_quantile = new double[num_gene];
	double *exon_mean = new double[num_gene];
	double *exon_sd = new double[num_gene];
	double *DNase_mean = new double[num_loci];
	double *DNase_sd = new double[num_loci];
	int *cluster_idx = new int[num_gene];

	double **coef = new double *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		coef[j] = new double[num_loci]();
	}

	int **pre_idx = new int *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		pre_idx[j] = new int[num_loci]();
	}

	char **select_loci = new char *[num_loci];
	for (int j = 0; j < num_loci; j++)
	{
		select_loci[j] = new char[30];
	}

	char **row_name = new char *[num_gene];
	for (int j = 0; j < num_gene; j++)
	{
		row_name[j] = new char[30];
	}

	if (ReadVector_d(filelist[0], uw_quantile)){ return 1; }
	if (ReadVector_d(filelist[1], exon_mean)){ return 1; }
	if (ReadVector_d(filelist[2], exon_sd)){ return 1; }
	if (ReadVector_i(filelist[3], cluster_idx)){ return 1; }
	if (ReadVector_d(filelist[4], DNase_mean)){ return 1; }
	if (ReadVector_d(filelist[5], DNase_sd)){ return 1; }
	if (ReadMatrix_d(filelist[6], coef)){ return 1; }
	if (ReadMatrix_i(filelist[7], pre_idx)){ return 1; }
	if (ReadName(filelist[8], select_loci)){ return 1; }
	if (ReadName(filelist[9], row_name)){ return 1; }

	p_out = fopen(outfile, "wb");
	if (p_out != NULL)
	{
		std::cout << "Writing " << outfile << "..." << std::endl;
		fwrite(param, sizeof(int), 5, p_out);
		fwrite(uw_quantile, sizeof(double), num_gene, p_out);
		fwrite(exon_mean, sizeof(double), num_gene, p_out);
		fwrite(exon_sd, sizeof(double), num_gene, p_out);
		for (int k = 0; k < num_var; k++)
		{
			fwrite(coef[k], sizeof(double), num_loci, p_out);
		}
		fwrite(DNase_mean, sizeof(double), num_loci, p_out);
		fwrite(DNase_sd, sizeof(double), num_loci, p_out);
		for (int k = 0; k < num_var; k++)
		{
			fwrite(pre_idx[k], sizeof(int), num_loci, p_out);
		}
		fwrite(cluster_idx, sizeof(int), num_gene, p_out);
		for (int k = 0; k < num_loci; k++)
		{
			fwrite(select_loci[k], sizeof(char), 30, p_out);
		}
		for (int k = 0; k < num_gene; k++)
		{
			fwrite(row_name[k], sizeof(char), 30, p_out);
		}
		fclose(p_out);
		std::cout << "Successfully generated library file." << std::endl;
	}
	else
	{
		std::cout << "Error! Cannot write " << outfile << std::endl;
		return 1;
	}

	

	//Release memory
	delete[] uw_quantile;
	delete[] exon_mean;
	delete[] exon_sd;
	delete[] DNase_mean;
	delete[] DNase_sd;
	delete[] cluster_idx;

	for (int j = 0; j < num_var; j++)
	{
		delete[] coef[j];
		delete[] pre_idx[j];
	}
	delete[] coef;
	delete[] pre_idx;
	
	return 0;


}
