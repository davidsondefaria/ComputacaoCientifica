#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double *y, *fn, *xn, *V;
double Energia;

double sinal(double x, double y) {
	if (y >= 0)
		return fabs(x);
	else if (y < 0)
		return -fabs(x);
}

main(int argc, char *argv[]) {
	int i, j, n, conta0, cont = 0, pc;
	int tam;
	double eup, edn, ypc, norm, djump;
	double dx, xf;
	double maior=0, *p, mp=0;

	if (sscanf (argv[1], "%i", &n) != 1) {
    	fprintf(stderr, "error - not an integer");
	}

	if (sscanf (argv[2], "%i", &tam) != 1) {
    	fprintf(stderr, "error - not an integer");
	}

	if (sscanf (argv[3], "%lf", &xf) != 1) {
    	fprintf(stderr, "error - not an double");
	}
	printf("%lf\n", xf);

	y = (double *)malloc((tam) * sizeof(double));
	fn = (double *)malloc((tam) * sizeof(double));
	xn = (double *)malloc((tam) * sizeof(double));
	V = (double *)malloc((tam) * sizeof(double));
	p = (double *)malloc((tam) * sizeof(double));

	
	FILE *fnfile;
	fnfile = fopen("Energia.out", "a");
	
	dx = xf/tam;

	for (i = 0; i < tam; i++) {
		xn[i] = (double)i * dx;
		V[i] = 0.5 * xn[i] * xn[i];
	}

	eup = V[tam - 1];
	edn = eup;

	for (i = 0; i < tam; i++) {
		if (eup < V[i])
			eup = V[i];
		if (edn > V[i])
			edn = V[i];
	}

	Energia = (eup + edn) / 2.0;
	//printf("#%g\n", Energia);

	do {
		pc = -1;

		for (i = 0; i < tam; i++) {
			fn[i] = -2.0 * (Energia - V[i]) * dx * dx / 12.0;
			if (fn[i] == 0.)
				fn[i] = 1e-20;
			if (fn[i] != sinal(fn[i], fn[i + 1])) {		//verifica sinal
				pc = i;
			}
		}

		for (i = 0; i < tam; i++) {
			fn[i] = 1-fn[i];
			y[i] = 0;
		}

		printf("#Após y0: cont: %d\t\n", cont);

		//paridade
		if (n % 2 == 0) {
			y[0] = 1.0;
			y[1] = 0.5 * (12.0 - 10.0 * fn[0]) * y[0] / fn[1];
		} else {
			y[0] = 0.0;
			y[1] = dx;
		}

		printf("#Após paridade y: cont: %d\t\n", cont);

		conta0 = 0;
		printf("#pc: %d\n", pc);
		for (i = 1; i <=pc - 1; i++) {
			//y[i+1] = psi(i);
			y[i + 1] = ((12.0 - 10.0 * fn[i]) * y[i] - fn[i - 1] * y[i - 1]) / fn[i + 1];
			printf("#yi: %lf\tyi+1%lf\n", y[i], y[i+1]);
			if (y[i] != sinal(y[i], y[i + 1])) {		//verifica sinal
				//printf("#%lf\t%lf\n", y[i], y[i+1]);
				conta0++;		//conta zeros
			}
		}
		ypc = y[pc];

		printf("#Após yi+1: cont: %d\t\n", cont);

		if (n % 2 == 0)
			conta0 = 2 * conta0;
		else
			conta0 = 2 * conta0 + 1;

		//if (cont > 1) {
		if (conta0 != n) {
			if (conta0 > n) {
				eup = Energia;
			}
			else {
				edn = Energia;
			}
			Energia = 0.5 * (eup + edn);
		}
		//}

		printf("#Após Energia: cont: %d\t\n", cont);

		if (conta0 == n) {
			printf("#Entra no if cont ou conta0\n");
			y[tam - 1] = dx;
			y[tam - 2] = (12. - 10.*fn[tam - 1]) * y[tam - 1] / fn[tam - 2];

			printf("#\tApós condição inicial do ytam-1\n");

			for (i = tam - 2; i >= pc + 1; i--) {
				y[i - 1] = ((12. - 10.*fn[i]) * y[i] - fn[i + 1] * y[i + 1]) / fn[i - 1];
			}

			printf("#\tApós yi-1: cont: %d\t\n", cont);

			ypc = ypc / y[pc];
			for (i = pc; i < tam; i++) {
				y[i] = y[i] * ypc;
			}

			/*norm = 0;
			for (i = 1; i < tam; i++) {
				norm += y[i] * y[i];
			}
			norm = dx*(2.*norm + y[0] * y[0]);
			norm = sqrt(norm);
			printf("#norm: %g\n", norm);

			for (i = 0; i < tam; i++) {
				y[i] = y[i] / norm;
			}*/

			//if (cont > 1) {
			i = pc;
			djump = (y[i + 1] + y[i - 1] - (14 - 12 * fn[i]) * y[i]) / dx;
			if (djump * y[i] > 0.0)
				eup = Energia;
			else
				edn = Energia;
			Energia = 0.5 * (eup + edn);
			//}
		}
		printf("#Saiu ou não entrou no if cont ou conta0\n");	

		printf("#aqui: #%g\t%d\n", Energia, cont);
		cont++;
		printf("##########################\n\n");


	} while (fabs(eup - edn) > 1e-6 && cont<10000);


	for(i=0; i<tam; i++){
		if((xn[pc]*xn[pc] - xn[i]*xn[i]) <=0)
			p[i] = 0;
		else
			p[i] = 1.0/(M_PI*sqrt(xn[pc]*xn[pc] - xn[i]*xn[i]));
		if (fabs(p[i]) > fabs(mp))
			mp = p[i];	
		if (fabs(y[i]) > fabs(maior))
			maior = y[i];
	}
	
	for(i=0; i<tam; i++){
		y[i] = y[i]/maior;
		p[i] = p[i]/mp;
	}

	if(n%2==0){
		for (i = tam-1; i > 0; i--) {
			printf("%d\t%lf\t%g\t%g\t%g\t%g\t%g\n", i, -xn[i], y[i], y[i]*y[i], Energia, V[i], p[i]);
		}
	}else{
		for (i = tam-1; i > 0; i--) {
			printf("%d\t%lf\t%g\t%g\t%g\t%g\t%g\n", i, -xn[i], -y[i], y[i]*y[i], Energia, V[i], p[i]);
		}
	}
	for (i = 0; i < tam; i++) {
		printf("%d\t%lf\t%g\t%g\t%g\t%g\t%g\n", i, xn[i], y[i], y[i]*y[i], Energia, V[i], p[i]);
	}

	printf("#%d\t%d\t%d\n", n, conta0, cont);

	fprintf(fnfile, "%lf\t%lf\n", Energia, xn[pc]);
	fclose(fnfile);
	/*free(V);
	free(fn);
	free(y);
	free(xn)*/;
}