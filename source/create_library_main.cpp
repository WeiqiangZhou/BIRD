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
                std::cout << "-i Specify input file which contains the required parameters and files list." << std::endl;
                std::cout << "-o Specify output name of the model file to be used by BIRD." << std::endl;
                return 1;
    }
	else
	{
		std::cout << "Please input the correct parameters." << std::endl;
		std::cout << "Example: BIRD_build_library -i par_file.txt -o model_file.bin" << std::endl;
		return 1;
	}

	FILE * p_out;
	int param[8];
	int num_loci, num_gene, cluster_size, bin_size, num_var, DH_num1, DH_num2, DH_num3;
	char **filelist = new char *[18];
	for (int i = 0; i < 18; i++)
	{
		filelist[i] = new char[255];
	}

	if (ReadParam(infile, param, filelist)){ return 1; }

	num_loci = param[0];
	num_gene = param[1];
	cluster_size = param[2];
	bin_size = param[3];
	num_var = param[4];
	DH_num1 = param[5];
	DH_num2 = param[6];
	DH_num3 = param[7];

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

	double **dis_matrix = new double *[3];
	for (int j = 0; j < 3; j++)
	{
		dis_matrix[j] = new double[num_loci]();
	}

	double **DH_coef1 = new double *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		DH_coef1[j] = new double[DH_num1]();
	}

	double **DH_coef2 = new double *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		DH_coef2[j] = new double[DH_num2]();
	}

	double **DH_coef3 = new double *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		DH_coef3[j] = new double[DH_num3]();
	}

	int **pre_idx = new int *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		pre_idx[j] = new int[num_loci]();
	}

	int **DH_pre_idx1 = new int *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		DH_pre_idx1[j] = new int[DH_num1]();
	}

	int **DH_pre_idx2 = new int *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		DH_pre_idx2[j] = new int[DH_num2]();
	}

	int **DH_pre_idx3 = new int *[num_var];
	for (int j = 0; j < num_var; j++)
	{
		DH_pre_idx3[j] = new int[DH_num3]();
	}

	int **DH_cluster = new int *[3];
	for (int j = 0; j < 3; j++)
	{
		DH_cluster[j] = new int[num_loci]();
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
	if (ReadMatrix_d(filelist[10], dis_matrix)){ return 1; }
	if (ReadMatrix_i(filelist[11], DH_cluster)){ return 1; }
	if (ReadMatrix_d(filelist[12], DH_coef1)){ return 1; }
	if (ReadMatrix_d(filelist[13], DH_coef2)){ return 1; }
	if (ReadMatrix_d(filelist[14], DH_coef3)){ return 1; }
	if (ReadMatrix_i(filelist[15], DH_pre_idx1)){ return 1; }
	if (ReadMatrix_i(filelist[16], DH_pre_idx2)){ return 1; }
	if (ReadMatrix_i(filelist[17], DH_pre_idx3)){ return 1; }

	p_out = fopen(outfile, "wb");
	if (p_out != NULL)
	{
		std::cout << "Writing " << outfile << "..." << std::endl;
		fwrite(param, sizeof(int), 8, p_out);
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
		for (int k = 0; k < 3; k++)
		{
			fwrite(dis_matrix[k], sizeof(double), num_loci, p_out);
		}
		for (int k = 0; k < 3; k++)
		{
			fwrite(DH_cluster[k], sizeof(int), num_loci, p_out);
		}
		for (int k = 0; k < num_var; k++)
		{
			fwrite(DH_coef1[k], sizeof(double), DH_num1, p_out);
		}
		for (int k = 0; k < num_var; k++)
		{
			fwrite(DH_coef2[k], sizeof(double), DH_num2, p_out);
		}
		for (int k = 0; k < num_var; k++)
		{
			fwrite(DH_coef3[k], sizeof(double), DH_num3, p_out);
		}
		for (int k = 0; k < num_var; k++)
		{
			fwrite(DH_pre_idx1[k], sizeof(int), DH_num1, p_out);
		}
		for (int k = 0; k < num_var; k++)
		{
			fwrite(DH_pre_idx2[k], sizeof(int), DH_num2, p_out);
		}
		for (int k = 0; k < num_var; k++)
		{
			fwrite(DH_pre_idx3[k], sizeof(int), DH_num3, p_out);
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
	for (int j = 0; j < 3; j++)
	{
		delete[] dis_matrix[j];
		delete[] DH_cluster[j];
	}
	for (int j = 0; j < num_var; j++)
	{
		delete[] coef[j];
		delete[] pre_idx[j];
		delete[] DH_coef1[j];
		delete[] DH_coef2[j];
		delete[] DH_coef3[j];
		delete[] DH_pre_idx1[j];
		delete[] DH_pre_idx2[j];
		delete[] DH_pre_idx3[j];
	}
	delete[] coef;
	delete[] pre_idx;
	delete[] DH_coef1;
	delete[] DH_coef2;
	delete[] DH_coef3;
	delete[] DH_pre_idx1;
	delete[] DH_pre_idx2;
	delete[] DH_pre_idx3;
	delete[] dis_matrix;
	delete[] DH_cluster;

	return 0;


}
