#include"MatrixCalculator.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"malloc.h"
/*
create a matrix with 0
*/
Matrix_t create_mat(uint16_t row, uint16_t column)
{
	Matrix_t mat;
	if (row <= 0||column<=0)
	{
		printf("error, in create_mat: row <= 0||column<=0\n");
		exit(1);
	}
	if (row > 0 && column > 0)
	{
		mat.row = row;
		mat.column = column;
		mat.data = (float **)malloc(row*sizeof(float *));//先指针的指针
		if (mat.data == NULL)
		{
			printf("error, in create_mat: mat.data==NULL");
			exit(1);
		}
		uint16_t i;
		for (i = 0; i < row; i++)
		{
			*(mat.data + i) = (float *)malloc(column*sizeof(float));//再分配每行的指针
			if (mat.data[i] == NULL)
			{
			  printf("error, in create_mat: mat.data==NULL");
			  exit(1);
			}
		}
		clear_mat(&mat);
	}
	return mat;
}
/*
free a matrix
*/
void free_mat(Matrix_t *mat)
{
	uint16_t i;
	for (i = 0; i < mat->row; i++)
		free(mat->data[i]);/*释放行*/
	free(mat->data);/*释放头指针*/
}

/*
display the matrix
*/
void show_mat(const char* name,const Matrix_t *mat)
{
	uint16_t i, j;
	printf("%s=\n",name);
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			printf("%.6f\t", mat->data[i][j]);
		}
		printf("\n");
	}
}
/*
set datas to the matrix
*/
void set_mat_data(Matrix_t* mat,const float *data)
{
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat->data[i][j] = data[i*mat->column+j];
		}
	}
}
/*
mat=mat1+mat2
*/

Matrix_t add_mat(const Matrix_t* mat1, const Matrix_t* mat2)
{
	
	if (mat1->row != mat2->row)
	{
		printf("error, in add_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1->column != mat2->column)
	{
		printf("error, in add_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	Matrix_t mat;
	uint16_t i, j;
	mat = create_mat(mat1->row, mat1->column);
	for (i = 0; i < mat1->row; i++)
	{
		for (j = 0; j < mat1->column; j++)
			mat.data[i][j] = mat1->data[i][j] + mat2->data[i][j];
	}
	return mat;
}

void add_mat_2(Matrix_t* mat,const Matrix_t* mat1, const Matrix_t* mat2)
{

	if (mat1->row != mat2->row)
	{
		printf("error, in add_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1->column != mat2->column)
	{
		printf("error, in add_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	if (mat->column != mat1->column || mat->row != mat1->row)
	{
		printf("error, in add_mat: mat->column != mat1->column || mat->row != mat1->row\n");
		exit(1);
	}
	uint16_t i, j;
	for (i = 0; i < mat1->row; i++)
	{
		for (j = 0; j < mat1->column; j++)
			mat->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
	}
}

/*
mat=mat1-mat2
*/
void sub_mat_2(Matrix_t* mat,const Matrix_t* mat1, const Matrix_t* mat2)//mat1-mat2;
{
	if (mat1->row != mat2->row)
	{
		printf("error, in sub_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1->column != mat2->column)
	{
		printf("error, in sub_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	if (mat->column != mat1->column || mat->row != mat1->row)
	{
		printf("error, in sub_mat: mat->column != mat1->column || mat->row != mat1->row\n");
		exit(1);
	}
	
	uint16_t i, j;
	for (i = 0; i < mat1->row; i++)
	{
		for (j = 0; j < mat1->column; j++)
			mat->data[i][j] = mat1->data[i][j] - mat2->data[i][j];
	}
}
Matrix_t sub_mat(const Matrix_t* mat1, const Matrix_t* mat2)//mat1-mat2;
{
	if (mat1->row != mat2->row)
	{
		printf("error, in sub_mat: mat1->row != mat2->row\n");
		exit(1);
	}
	if (mat1->column != mat2->column)
	{
		printf("error, in sub_mat: mat1->column != mat2->column\n");
		exit(1);
	}
	Matrix_t mat;
	uint16_t i, j;
	mat = create_mat(mat1->row, mat1->column);
	for (i = 0; i < mat1->row; i++)
	{
		for (j = 0; j < mat1->column; j++)
			mat.data[i][j] = mat1->data[i][j] - mat2->data[i][j];
	}
	return mat;
}
/*
transpose the matrix, mat=mat'
*/
Matrix_t transpose_mat(const Matrix_t* mat)//mat'
{
	Matrix_t mat_T;
	mat_T = create_mat(mat->row, mat->column);
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat_T.data[i][j] = mat->data[j][i];
		}
	}
	return mat_T;
}
void transpose_mat_2(Matrix_t* mat_T,const Matrix_t* mat)//mat'
{
	if (mat_T->row != mat->column)
	{
		printf("error, in transpose_mat_2: mat_T->row != mat->column\n");
		exit(1);
	}
	if (mat_T->column != mat->row)
	{
		printf("error, in transpose_mat_2: mat_T->column != mat->row\n");
		exit(1);
	}
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat_T->data[i][j] = mat->data[j][i];
		}
	}
}
/*
mat=scaler*mat
every element in the matrix multiplys a scaler
*/
Matrix_t scale_mat(const Matrix_t* mat, const float scaler)//scaler*Mat
{
	Matrix_t mat1;
	mat1 = create_mat(mat->row, mat->column);
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat1.data[i][j] = mat->data[i][j]*scaler;
		}
	}
	return mat1;
}
void scale_mat_2(Matrix_t* mat, const float scaler)//Mat=scaler*Mat
{
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat->data[i][j] = mat->data[i][j] * scaler;
		}
	}
}
/*
set all datas in matrix to zero
*/
void clear_mat(Matrix_t* mat)
{
	uint16_t i, j;
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat->data[i][j] = 0;
		}
	}
}
/*
mat=mat1*mat2
*/
Matrix_t mult_mat(const Matrix_t *mat1, const Matrix_t* mat2)
{
	Matrix_t mat;
	if (mat1->column != mat2->row)
	{
		printf("error,In mult_mat: mat1->column != mat2->row\n");
		exit(1);
	}	
	else
	{
		mat = create_mat(mat1->row, mat2->column);
		clear_mat(&mat);
		uint16_t i, j;
		for (i = 0; i < mat1->row; i++)
		{
			for (j = 0; j < mat2->column; j++)
			{
				uint16_t m;
				for (m = 0; m < mat1->column; m++)
				{
					mat.data[i][j] += mat1->data[i][m] * mat2->data[m][j];
				}
			}
		}
	}
	return mat;
}
void mult_mat_2(Matrix_t *mat,const Matrix_t *mat1, const Matrix_t* mat2)
{
	if (mat1->column != mat2->row)
	{
		printf("error,In mult_mat: mat1->column != mat2->row\n");
		exit(1);
	}
	if (mat->column != mat1->row || mat->row != mat2->column)
	{
		printf("error,In mult_mat: mat->column != mat1->row || mat->row != mat2->column\n");
		exit(1);
	}
	else
	{
		
		clear_mat(&mat);
		uint16_t i, j;
		for (i = 0; i < mat1->row; i++)
		{
			for (j = 0; j < mat2->column; j++)
			{
				uint16_t m;
				for (m = 0; m < mat1->column; m++)
				{
					mat->data[i][j] += mat1->data[i][m] * mat2->data[m][j];
				}
			}
		}
	}
}
/*
result = mat1(vector) * mat2(vector)
*/
float dot_product(const Matrix_t *mat1, const Matrix_t* mat2)
{
	float result;
	if (mat1->column != mat2-> column || mat1->row != mat2->row )
	{
		printf("error,In dot_product: mat1->column != mat2-> column || mat1->row != mat2->row\n");
		exit(1);
	}
	else
	{
		result = 0;
		uint16_t i, j;
		for (i = 0; i < mat1->row; i++)
		{
			for (j = 0; j < mat1->column; j++)
			{
				//printf("i=%d ,j=%d \nm1=%f, m2=%f\n",i,j, mat1->data[i][j] , mat2->data[j][i]);
				result += mat1->data[i][j] * mat2->data[i][j];
				//printf("res=%f\n", result);
			}
		}
	}
	return result;

}
/*
generate a I(nxn) matrix
*/
Matrix_t eye(uint16_t n)
{
	if (n <= 0)
	{
		printf("error, in eye: n<0\n");
		exit(1);
	}
	Matrix_t mat = create_mat(n, n);
	uint16_t i,j;
	for (i = 0; i < n;i++)
	for (j = 0; j < n; j++)
	{
		if (i == j)
			mat.data[i][j] = 1;
		else
			mat.data[i][j] = 0;
	}
	return mat;
}
/*
generate a diagonal matrix with diag[n] as its main diagonal elements
*/
Matrix_t diag_mat(uint16_t n, float* diag)
{
	if (n <= 0){
		printf("error: in diag_mat(n<0)\n");
		exit(1);
	}
	Matrix_t mat = create_mat(n, n);
	uint16_t i, j;
	for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	{
		if (i == j)
			mat.data[i][j] = diag[i];
		else
			mat.data[i][j] = 0;
	}
	return mat;
}
/*
copy a matrix
*/
Matrix_t copy_mat(const Matrix_t *mat)
{
	uint16_t i, j;
	Matrix_t mat_copy = create_mat(mat->row, mat->column);
	for (i = 0; i < mat->row; i++)
	{
		for (j = 0; j < mat->column; j++)
		{
			mat_copy.data[i][j] = mat->data[i][j];
		}
	}
	return mat_copy;
}
/*
copy matrix's data to another matrix
*/
void copy_mat_data(Matrix_t* copy,const Matrix_t* data)
{
	if (data->row != copy->row || data->column != copy->column)
	{
		printf("error, in copy_mat_data: data->row != copy->row || data->column != copy->column\n");
		exit(1);
	}
	uint16_t i, j;
	for (i = 0; i < data->row;i++)
	for (j = 0; j < data->column; j++)
	{
		copy->data[i][j] = data->data[i][j];
	}
}
/*
get matrix's derterminent value
*/
float det_mat(Matrix_t *m)
{
	uint16_t i, j, n, max_row;
	char swap_f;
	float max, k;
	float det=1;
	swap_f = 0;
	if (m->column != m->row)
	{
		printf("error:In det_mat (m->column != m->row)\n");
		exit(1);
	}
	Matrix_t mat = copy_mat(m);
	for (i = 0; i < mat.row-1; i++)
	{
		max = abs(mat.data[i][i]);
		max_row = i;
		for (j = i + 1; j < mat.row; j++)
		{
			if (max < abs(mat.data[j][i]))
			{
				max = abs(mat.data[j][i]);
				max_row = j;
			}
		}
		if (i != max_row)
		{
			swap_row_mat(&mat, i, max_row);
			swap_f++;
		}
		for (j = i + 1; j < mat.row; j++)
		{
			k = -mat.data[j][i]/mat.data[i][i];
			for (n= 0; n < mat.column; n++)
			{
				mat.data[j][n] = mat.data[i][n] * k + mat.data[j][n];
			}
		}
	}
	if (swap_f%2==1)swap_f = -1;
	else swap_f = 1;
	det = 1;
	for (i = 0; i < mat.column; i++)
		det *= mat.data[i][i];
	det *= swap_f;
	//show_mat("in det_mat: mat_upper",&mat);
	free_mat(&mat);
	return det;
}
/*
get inverse matrix
use main column element of Gauss-Jordan algrithm: A|I  --->  I|A^(-1)
*/
Matrix_t inverse_mat(Matrix_t* m)
{
	if (det_mat(m) == 0)
	{
		printf("error: In inverse_mat(det_mat(mat) == 0)\n");
		exit(1);
	}

	Matrix_t mat = copy_mat(m);
	Matrix_t inv_mat = eye(m->row);
	
	int i, j, n,max_row;
	int swap_f = 0;
	float max,k;
	for (i = 0; i < mat.row - 1; i++)
	{
		max = abs(mat.data[i][i]);
		max_row = i;
		for (j = i + 1; j < mat.row; j++)
		{
			if (max < abs(mat.data[j][i]))
			{
				max = abs(mat.data[j][i]);
				max_row = j;
			}
		}
		if (i != max_row)
		{
			swap_row_mat(&mat, i, max_row);
			swap_row_mat(&inv_mat, i, max_row);
			swap_f++;
		}
		for (j = i + 1; j < mat.row; j++)
		{
			k = -mat.data[j][i] / mat.data[i][i];
			for (n = 0; n < mat.column; n++)
			{
				mat.data[j][n] = mat.data[i][n] * k + mat.data[j][n];
				inv_mat.data[j][n] = inv_mat.data[i][n] * k + inv_mat.data[j][n];
			}
		}
	}

	for (i = 0; i < mat.row; i++)
	{
		k = 1/mat.data[i][i];
		scale_row_mat(&mat,i, k);
		scale_row_mat(&inv_mat, i, k);
	}
	for (i = mat.row-1; i>0; i--)
	{
		for (j = i - 1; j >=0; j--)
		{
			k = -mat.data[j][i] / mat.data[i][i];
			for (n = 0; n < mat.column; n++)
			{
				mat.data[j][n] = mat.data[j][n] + k*mat.data[i][n];
				inv_mat.data[j][n] = inv_mat.data[j][n] + k*inv_mat.data[i][n];
			}
		}
		
	}
	//show_mat("in inverse_mat: mat", &mat);
	//show_mat("in inverse_mat: inv_mat", &inv_mat);
	free_mat(&mat);
	return inv_mat;
}
/*swap NO.m and NO.n row in mat*/
void swap_row_mat(Matrix_t* mat,uint16_t m,uint16_t n)
{
	float temp;
	uint16_t i;
	for (i = 0; i < mat->column; i++)
	{
		temp = mat->data[m][i];
		mat->data[m][i] = mat->data[n][i];
		mat->data[n][i] = temp;
	}
}
/*
NO.m row in matrix multiply a scaler
*/
void scale_row_mat(Matrix_t* mat, uint16_t m, float scaler)
{
	uint16_t i;
	for (i = 0; i < mat->column; i++)
		mat->data[m][i] *= scaler;
}

