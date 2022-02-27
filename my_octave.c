// Copyright 2021 Popa Bianca 311CA <biancapopa2002@gmail.com>

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "citire.h"

int main(void)
{
	int ***vec_mat, nr_mat = 0, *vec_m, *vec_n;

	char ch;
	scanf("%c", &ch);
	while (ch != 'Q') {
		if (ch == 'L')
			//incarcarea in memorie a unei matrice
			cit_mat(&vec_mat, &vec_m, &vec_n, &nr_mat);

		else if (ch == 'D')
			//determinarea dimensiunilor unei matrice
			dimensiune(&vec_m, &vec_n, &nr_mat);

		else if (ch == 'P')
			//afisarea unei matrice
			afisare(&vec_mat, &vec_m, &vec_n, &nr_mat);

		else if (ch == 'C')
			//redimensionarea unei matrice
			redimensionare(&vec_mat, &vec_m, &vec_n, &nr_mat);

		else if (ch == 'M')
			//inmultirea a doua matrice
			inmultire(&vec_mat, &vec_m, &vec_n, &nr_mat);

		else if (ch == 'O')
			//sortarea matricelor
			sortare(&vec_mat, &vec_m, &vec_n, &nr_mat);

		else if (ch == 'T')
			//transpunerea unei matrice
			transpusa(&vec_mat, &vec_m, &vec_n, &nr_mat);

		else if (ch == 'F')
			//eliberarea memoriei unei matrice
			eliberare(&vec_mat, &vec_m, &vec_n, &nr_mat);

		else if (isalpha(ch))
			printf("Unrecognized command\n");

		scanf("%c", &ch);
	}
	//eliberarea tuturor resurselor
	for (int i = 0; i < nr_mat; i++)
		for (int j = 0; j < vec_m[i]; j++)
			free(vec_mat[i][j]);

	for (int i = 0; i < nr_mat; i++)
		free(vec_mat[i]);

	if (nr_mat)
		free(vec_mat), free(vec_m), free(vec_n);

	return 0;
}
