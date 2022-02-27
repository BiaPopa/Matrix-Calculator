// Copyright 2021 Popa Bianca 311CA <biancapopa2002@gmail.com>

#ifndef CITIRE
#define CITIRE

void cit_mat(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat);
void dimensiune(int **vec_m, int **vec_n, int *nr_mat);
void afisare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat);
void redimensionare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat);
void inmultire(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat);
void sortare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat);
void transpusa(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat);
void eliberare(int ****vec_mat, int **vec_m, int **vec_n, int *nr_mat);

#endif
