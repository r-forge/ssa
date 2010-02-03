
void
cconv (int *l, int *m, double *x, int *n, double *s)
{
double *y, *z, *u;
int i;
#pragma omp parallel for default(shared) private(i) schedule(dynamic)
for (i = 0; i < *m; i++)
    {
    y = x + (*n - l[i]);
    z = x + l[i];
    u = x;
    s[i] = 0.0;
    while (u < y)
    s[i] += *u++ * *z++;
    }
}

void
hankelC(double *x, double *c, int *l, int *t)
{
int i, j, h, K, L, T; double s;
L = *l; T=*t; K = T - L + 1; 
#pragma omp parallel for default(shared) private(i, j, s) schedule(dynamic)
for (j = 1; j <= L; j++) 
	for (h = j; h <= L; h++) {
		s = 0.0;
		for (i = 1; i <= K; i++) 
			s += x[i+j-2] * x[i+h-2];
		c[j + (h - 1) * L - 1] = c[h + (j - 1) * L - 1] = s;
		}
}

void
fromCrossC(double *x, double *s, double *a, int *l, int *t)
{
int i, j, k, h, ilw, iup, ni, K, L, T;
double sv, sw;
L = *l; T=*t; K = T - L + 1; 
#pragma omp parallel for default(shared) private(i, j, k, sv, sw) schedule(dynamic)
for (i = 1; i <= T; i++) {
	ilw = (i + 1) - L; if (ilw < 1) ilw = 1;
	iup = i; if (iup > K) iup = K;
	ni = iup - ilw + 1;
	for (j = 1; j <= L; j++) {
		sv = 0.0;
		for (k = ilw; k <= iup; k++) {
			sw = 0.0;
			for (h = 1; h <= L; h++)
				sw += x[k + h - 2] * s[h + (j - 1) * L - 1];
			sv += sw * s[(i + 1) - k + (j - 1) * L - 1];
		    }
	    a[i + (j - 1) * T - 1] = sv / ((double) ni);
	    }
    }
}
