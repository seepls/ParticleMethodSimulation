#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "define.h"
#include "inverse.h"


void INVERSE_matrixInverse(double matrix[10][10], int size){
	double det;

	det = INVERSE_determinant(matrix,size);

	if (fabs(det) < TOO_SMALL_VALUE){
		printf("Warning: Determinant is 0 in INVERSE_matrixInverse, setting identity matrix!\n");
		INVERSE_setIdendityMatrix(matrix,size);
		return;
	}

	INVERSE_cofactorSymmetric(matrix,size);

	INVERSE_transposeAndDivideDeterminant(matrix,size,det);

}

double INVERSE_determinant(double matrix[10][10], int size){
	double s = 1, det = 0, subMatrix[10][10];
	int i,j,m,n,c;

	if (size == 1)
		return matrix[0][0];
	else if (size == 2)
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	else if (size == 3){
		for(i = 0; i < 3; i++)
			det += (matrix[0][i] * (matrix[1][(i+1)%3] * matrix[2][(i+2)%3] - matrix[1][(i+2)%3] * matrix[2][(i+1)%3]));
		return det;
	}
	else{
		det = 0;
		for (c=0;c<size;c++){
			m = 0;
			n = 0;
			for (i=0; i<size; i++){
				for (j=0; j<size; j++){
					subMatrix[i][j] = 0;
					if ((i != 0) && (j != c)){
						subMatrix[m][n]=matrix[i][j];
						if (n<(size-2))
							n++;
						else{
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (matrix[0][c] * INVERSE_determinant(subMatrix, size-1));
			s = -1 * s;
		}
	}
	return det;
}

void INVERSE_setIdendityMatrix(double matrix[10][10],int size){
	int i,j;
	for (i=0; i<size; i++)
		for (j=0; j<size; j++){
			if(i == j)
				matrix[i][j] = 1;
			else
				matrix[i][j] = 0;
		}

}

void INVERSE_cofactor(double matrix[10][10], int size){
	double subMatrix[10][10], cofacMatrix[10][10];
	int p,q,m,n,i,j;
	
	for (q=0; q<size; q++){
		for(p=0; p<size; p++){
			m = 0;
			n = 0;
			for (i=0; i<size; i++){
				for (j=0; j<size; j++){
					if((i != q) && (j != p)){
						subMatrix[m][n] = matrix[i][j];
						if (n<(size-2))
							n++;
						else{
							n = 0;
							m++;
						}
					}
				}
			}
			cofacMatrix[q][p] = pow((double)-1,q+p) * INVERSE_determinant(subMatrix,size-1);
		}
	}

	for (i=0; i<size; i++)
		for (j=0; j<size; j++)
			matrix[i][j] = cofacMatrix[i][j];
}

void INVERSE_cofactorSymmetric(double matrix[10][10], int size){
	double subMatrix[10][10], cofacMatrix[10][10];
	int p,q,m,n,i,j;
	
	for (q=0; q<size; q++){
		for(p=0; p<=q; p++){
			m = 0;
			n = 0;
			for (i=0; i<size; i++){
				for (j=0; j<size; j++){
					if((i != q) && (j != p)){
						subMatrix[m][n] = matrix[i][j];
						if (n<(size-2))
							n++;
						else{
							n = 0;
							m++;
						}
					}
				}
			}
			cofacMatrix[q][p] = pow((double)-1,q+p) * INVERSE_determinant(subMatrix,size-1);
			if (p!=q)
				cofacMatrix[p][q] = cofacMatrix[q][p];
		}
	}

	for (i=0; i<size; i++)
		for (j=0; j<size; j++)
			matrix[i][j] = cofacMatrix[i][j];
}

void INVERSE_transposeAndDivideDeterminant(double matrix[10][10], int size, double det){
	int i,j;
	double transpose[10][10];

	for (i=0; i<size; i++)
		for (j=0; j<size; j++)
			transpose[i][j] = matrix[j][i];

	for (i=0; i<size; i++)
		for (j=0; j<size; j++)
			matrix[i][j] = transpose[i][j] / det;
	

}
