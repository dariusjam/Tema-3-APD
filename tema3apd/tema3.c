/* Mangea Liviu Darius 334CA */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

#define NUM_COLORS 256


int main(int argc, char **argv) {

	int numtasks, rank, tag = 1;  
	
	int tip_multime, MAX_STEPS;
	double x_min, x_max, y_min, y_max, resolution, cR = 0, cI = 0;
	int width, height;
	int i, j, k, step, contor, fin, part;
	double info[11], zr, zi, cr, ci, temp;;
	unsigned char *color;
	
	
	MPI_Status Stat;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	
	/* Procesul master (rank == 0) realizeaza citirea din fisier si pune
	informatiile intr-un vector info pe care il celorlalte procese. */
	if (rank == 0) {
		
		FILE *fr, *fw;

		fr = fopen(argv[1], "r");	
		fw = fopen(argv[2], "w");
		
		fscanf(fr, "%d", &tip_multime);
		
		fscanf(fr, "%lf", &x_min);
		fscanf(fr, "%lf", &x_max);
		fscanf(fr, "%lf", &y_min);
		fscanf(fr, "%lf", &y_max);
		fscanf(fr, "%lf", &resolution);
		fscanf(fr, "%d", &MAX_STEPS);
		
		if (tip_multime == 1) {
			fscanf(fr, "%lf", &cR);
			fscanf(fr, "%lf", &cI);
		}
		
		fclose(fr);
		
		double w = ((x_max - x_min) / resolution);
		double h = ((y_max - y_min) / resolution);
		width = (int)w;
		height = (int)h;
		
		info[0] = (double)tip_multime;
		info[1] = x_min;
		info[2] = x_max;
		info[3] = y_min;
		info[4] = y_max;
		info[5] = (double)width;
		info[6] = (double)height;
		info[7] = resolution;
		info[8] = (double)MAX_STEPS;
		info[9] = cR;
		info[10] = cI;
			
		/* Trimiterea vectorului cu informatiile necesare */
		for (i = 1; i < numtasks; i++) {
			MPI_Send(info, 11, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
		}
		
		fprintf(fw, "P2\n");
		fprintf(fw, "%d %d\n", width, height);
		fprintf(fw, "%d\n", NUM_COLORS - 1);
		
		color = (unsigned char*)malloc(width * sizeof(unsigned char));
		
		/* Pentru fiecare proces, de la cel cu rankul mai mare catre cel mai mic
		(deoarece matricea trebuie scrisa in fisier cu liniile in ordine inversa)
		se asteapta primirea liniilor din matrice, una cate una la fiecare
		receive, dupa care sunt scrise in fisier */
		for (i = numtasks - 1; i >= 1; i--) {
		
			/* Ultimul task executa calculeaza si restul de linii care raman din
			impartirea cu rest a liniilor la numarul de taskuri */
			if (i == numtasks - 1) {
				fin = (int)(height / numtasks) + height % numtasks;
			} else {
				fin = (int)(height / numtasks);
			}
			
			for(j = 0; j < fin; j++) {
				MPI_Recv(color, width, MPI_UNSIGNED_CHAR, i, tag, MPI_COMM_WORLD, &Stat);

				for (k = 0; k < width; k++) {
					fprintf(fw, "%d ", color[k]);
				}
				fprintf(fw, "\n");
			}
		}
		
		fclose(fw);
	
	}
	
	
	
	/* Daca procesul nu este master, se primesc informatiile in vectorul info */
	if(rank > 0) {
		MPI_Recv(info, 11, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &Stat);
	}
	
	tip_multime = info[0];
	x_min = info[1];
	x_max = info[2];
	y_min = info[3];
	y_max = info[4];
	width = (int)info[5];
	height = (int)info[6];
	resolution = info[7];
	MAX_STEPS = (int)info[8];
	cR = info[9];
	cI = info[10];
	

	color = (unsigned char*)malloc(width * sizeof(unsigned char));
	contor = 0;
	part = height / numtasks;
	
	/* In functie de rank se calculeaza numarul de linii pe care il va calcula
	fiecare proces */
	if(rank == numtasks - 1) {
		i = height - 1;
		fin = rank * part;
	} else {
		i = (rank + 1) * part-1;
		fin = rank * part;
	}
	
	/* Pentru fiecare linie din chunkul de linii al fiecarui proces, se executa
	in functie de tipul multimii, algoritmul prezentat in enunt */	
	for (; i >= fin; i--) {
		
		if (tip_multime == 0) {
			/* Mandelbrot */
			for(j = 0; j < width; j++) {
				cr = (double)(j * (x_max - x_min)) / width + x_min;
				ci = (double)(i * (y_max - y_min)) / height + y_min;
				zr = 0;
				zi = 0;
				step = 0;
				while(sqrt(zr*zr + zi*zi) < 2 && step < MAX_STEPS) {
					temp = zr;
					zr = zr * zr - zi * zi + cr;
					zi = 2 * zi * temp + ci;
					step++;
				}
				color[contor] = step % NUM_COLORS;
				contor++;			
			}
			
		} else {
			
			/* Julia */
			for(j = 0; j < width; j++) {
				zr = (double)(j * (x_max - x_min)) / width + x_min;
				zi = (double)(i * (y_max - y_min)) / height + y_min;

				step = 0;
				while(sqrt(zr*zr + zi*zi) < 2 && step < MAX_STEPS) {
					temp = zr;
					zr = zr * zr - zi * zi + cR;
					zi = 2 * zi * temp + cI;
					step++;
				}
				color[contor] = step % NUM_COLORS;
				contor++;	
			}
		}
		
		/* Daca rankul nu e 0 se trimite linia calculata la procesul master.
		Daca rank == 0 se scrie linia direct in fisierul de output */
		if(rank != 0) {
			MPI_Send(color, contor, MPI_UNSIGNED_CHAR, 0, tag, MPI_COMM_WORLD);
		} else {
			FILE *fw;
			fw = fopen(argv[2], "a");
	
			for (k = 0; k < width; k++) {
				fprintf(fw, "%d ", color[k]);
			}
			fprintf(fw, "\n");
		
			fclose(fw);
		}
		
		contor = 0;
	}
	
	
	MPI_Finalize();

	return 0;
}
