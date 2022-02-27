// Copyright 2021 Popa Bianca 311CA <biancapopa2002@gmail.com>

#include <stdio.h>
#include <stdlib.h>

//incarcarea in memorie a unei matrice

void cit_mat(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat)
{
	(*nr_mat)++;
	if (*nr_mat == 1) {
		(*vec_mat) = malloc(*nr_mat * sizeof(int **));
		(*vec_m) = malloc(*nr_mat * sizeof(int));
		(*vec_n) = malloc(*nr_mat * sizeof(int));
	} else {
		int ***aux_vec_mat, *aux_vec_m, *aux_vec_n;

		aux_vec_mat = (int ***)realloc(*vec_mat, *nr_mat * sizeof(int **));
		aux_vec_m = (int *)realloc(*vec_m, *nr_mat * sizeof(int));
		aux_vec_n = (int *)realloc(*vec_n, *nr_mat * sizeof(int));

		if (aux_vec_mat && aux_vec_m && aux_vec_n) {
			(*vec_mat) = aux_vec_mat;
			(*vec_m) = aux_vec_m;
			(*vec_n) = aux_vec_n;
		}
	}

	scanf("%d%d", &(*vec_m)[*nr_mat - 1], &(*vec_n)[*nr_mat - 1]);

	(*vec_mat)[*nr_mat - 1] =
	(int **)malloc((*vec_m)[*nr_mat - 1] * sizeof(int *));

	for (int i = 0; i < (*vec_m)[*nr_mat - 1]; i++)
		(*vec_mat)[*nr_mat - 1][i] =
		(int *)malloc((*vec_n)[*nr_mat - 1] * sizeof(int));
	for (int i = 0; i < (*vec_m)[*nr_mat - 1]; i++)
		for (int j = 0; j < (*vec_n)[*nr_mat - 1]; j++)
			scanf("%d", &(*vec_mat)[*nr_mat - 1][i][j]);
}

//determinarea dimensiunilor unei matrice

void dimensiune(int **vec_m, int **vec_n, int *nr_mat)
{
	int index;
	scanf("%d", &index);

	if (index >= *nr_mat || index < 0)
		printf("No matrix with the given index\n");
	else
		printf("%d %d\n", (*vec_m)[index], (*vec_n)[index]);
}

//afisarea unei matrice

void afisare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat)
{
	int index;
	scanf("%d", &index);

	if (index >= *nr_mat || index < 0)
		printf("No matrix with the given index\n");
	else
		for (int i = 0; i < (*vec_m)[index]; i++) {
			for (int j = 0; j < (*vec_n)[index]; j++)
				printf("%d ", (*vec_mat)[index][i][j]);
			printf("\n");
		}
}

//redimensionarea unei matrice

void redimensionare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat)
{
	int index;
	scanf("%d", &index);

	if (index >= *nr_mat || index < 0) {
		printf("No matrix with the given index\n");
	} else {
		int l, c, *vec_l, *vec_c;
		scanf("%d", &l);

		vec_l = malloc(l * sizeof(int));
		for (int i = 0; i < l; i++)
			scanf("%d", &vec_l[i]);

		scanf("%d", &c);

		vec_c = malloc(c * sizeof(int));
		for (int i = 0; i < c; i++)
			scanf("%d", &vec_c[i]);
		int **aux_mat;
		aux_mat = malloc(l * sizeof(int *));
		for (int i = 0; i < l; i++)
			aux_mat[i] = malloc(c * sizeof(int));

		for (int i = 0; i < l; i++)
			for (int j = 0; j < c; j++)
				aux_mat[i][j] = (*vec_mat)[index][vec_l[i]][vec_c[j]];

		for (int i = 0; i < (*vec_m)[index]; i++)
			free((*vec_mat)[index][i]);

		free((*vec_mat)[index]);

		(*vec_m)[index] = l;
		(*vec_n)[index] = c;
		(*vec_mat)[index] = aux_mat;
		free(vec_l), free(vec_c);
	}
}

//inmultirea a doua matrice

void inmultire(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat)
{
	int index1, index2;
	scanf("%d%d", &index1, &index2);

	if (index1 >= *nr_mat || index2 >= *nr_mat || index1 < 0 || index2 < 0) {
		printf("No matrix with the given index\n");
	} else if ((*vec_n)[index1] != (*vec_m)[index2]) {
		printf("Cannot perform matrix multiplication\n");
	} else {
		(*nr_mat)++;

		int ***aux_vec_mat, *aux_vec_m, *aux_vec_n;

		aux_vec_mat = (int ***)realloc(*vec_mat, *nr_mat * sizeof(int **));
		aux_vec_m = (int *)realloc(*vec_m, *nr_mat * sizeof(int));
		aux_vec_n = (int *)realloc(*vec_n, *nr_mat * sizeof(int));

		if (aux_vec_mat && aux_vec_m && aux_vec_n) {
			(*vec_mat) = aux_vec_mat;
			(*vec_m) = aux_vec_m;
			(*vec_n) = aux_vec_n;
		}

		(*vec_m)[*nr_mat - 1] = (*vec_m)[index1];
		(*vec_n)[*nr_mat - 1] = (*vec_n)[index2];

		(*vec_mat)[*nr_mat - 1] =
		(int **)malloc((*vec_m)[*nr_mat - 1] * sizeof(int *));
		for (int i = 0; i < (*vec_m)[*nr_mat - 1]; i++)
			(*vec_mat)[*nr_mat - 1][i] =
			(int *)malloc((*vec_n)[*nr_mat - 1] * sizeof(int));

		for (int i = 0; i < (*vec_m)[*nr_mat - 1]; i++) {
			for (int j = 0; j < (*vec_n)[*nr_mat - 1]; j++) {
				(*vec_mat)[*nr_mat - 1][i][j] = 0;
				for (int k = 0; k < (*vec_n)[index1]; k++) {
					(*vec_mat)[*nr_mat - 1][i][j] +=
					(*vec_mat)[index1][i][k] * (*vec_mat)[index2][k][j];
					(*vec_mat)[*nr_mat - 1][i][j] %= 10007;
					if ((*vec_mat)[*nr_mat - 1][i][j] < 0)
						(*vec_mat)[*nr_mat - 1][i][j] += 10007;
				}
			}
		}
	}
}

//sortarea matricelor

void sortare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat)
{
	int *sum;
	sum = (int *)calloc(*nr_mat, sizeof(int));

	for (int i = 0; i < *nr_mat; i++) {
		for (int j = 0; j < (*vec_m)[i]; j++) {
			for (int k = 0; k < (*vec_n)[i]; k++) {
				sum[i] += (*vec_mat)[i][j][k];
				sum[i] %= 10007;
				if (sum[i] < 0)
					sum[i] += 10007;
			}
		}
	}

	for (int i = 0; i < *nr_mat - 1; i++) {
		for (int j = i + 1; j < *nr_mat; j++) {
			if (sum[j] < sum[i]) {
				int aux;
				aux = sum[j];
				sum[j] = sum[i];
				sum[i] = aux;

				int **aux2;
				aux2 = (*vec_mat)[j];
				(*vec_mat)[j] = (*vec_mat)[i];
				(*vec_mat)[i] = aux2;

				int aux_3;
				aux_3 = (*vec_m)[i];
				(*vec_m)[i] = (*vec_m)[j];
				(*vec_m)[j] = aux_3;

				int aux_4;
				aux_4 = (*vec_n)[i];
				(*vec_n)[i] = (*vec_n)[j];
				(*vec_n)[j] = aux_4;
			}
		}
	}
	free(sum);
}

//transpunerea unei matrice

void transpusa(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat)
{
	int index;
	scanf("%d", &index);

	if (index >= *nr_mat || index < 0) {
		printf("No matrix with the given index\n");
	} else {
		int **aux_mat;
		aux_mat = malloc((*vec_n)[index] * sizeof(int *));
		for (int i = 0; i < (*vec_n)[index]; i++)
			aux_mat[i] = malloc((*vec_m)[index] * sizeof(int));

		for (int i = 0; i < (*vec_m)[index]; i++)
			for (int j = 0; j < (*vec_n)[index]; j++)
				aux_mat[j][i] = (*vec_mat)[index][i][j];

		for (int i = 0; i < (*vec_m)[index]; i++)
			free((*vec_mat)[index][i]);
		free((*vec_mat)[index]);

		int aux = (*vec_m)[index];
		(*vec_m)[index] = (*vec_n)[index];
		(*vec_n)[index] = aux;
		(*vec_mat)[index] = aux_mat;
	}
}

//eliberarea memoriei unei matrice

void eliberare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat)
{
	int index;
	scanf("%d", &index);

	if (index >= *nr_mat || index < 0) {
		printf("No matrix with the given index\n");
	} else {
		for (int i = 0; i < (*vec_m)[index]; i++)
			free((*vec_mat)[index][i]);
		free((*vec_mat)[index]);

		for (int i = index; i < *nr_mat - 1; i++) {
			(*vec_m)[i] = (*vec_m)[i + 1];
			(*vec_n)[i] = (*vec_n)[i + 1];
			(*vec_mat)[i] = (*vec_mat)[i + 1];
		}

		(*nr_mat)--;
		int ***aux_vec_mat, *aux_vec_m, *aux_vec_n;

		aux_vec_mat = (int ***)realloc(*vec_mat, *nr_mat * sizeof(int **));
		aux_vec_m = (int *)realloc(*vec_m, *nr_mat * sizeof(int));
		aux_vec_n = (int *)realloc(*vec_n, *nr_mat * sizeof(int));

		if (aux_vec_mat && aux_vec_m && aux_vec_n) {
			*vec_mat = aux_vec_mat;
			*vec_m = aux_vec_m;
			*vec_n = aux_vec_n;
		}
	}
}
