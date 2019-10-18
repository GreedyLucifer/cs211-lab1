#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/


void dgemm0(const double* A, const double* B, double* C, const int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) 
		{
			register double r = C[i * n + j];
			for (int k = 0; k < n; k++)
				r += A[i * n + k] * B[k * n + j];
			C[i * n + j] = r;
		}
}

void dgemm2(const double *A, const double *B, double *C, const int n) 
{
	for (int i = 0; i < n; i += 2)
		for (int j = 0; j < n; j += 2)
		{
			register double c00 = C[i * n + j];
			register double c01 = C[i * n + (j + 1)];
			register double c10 = C[(i + 1) * n + j];
			register double c11 = C[(i + 1) * n + (j + 1)];
			for (int k = 0; k < n; k += 2)
			{
				register double a00 = A[i * n + k];
				register double a01 = A[i * n + (k + 1)];
				register double a10 = A[(i + 1) * n + k];
				register double a11 = A[(i + 1) * n + (k + 1)];
				register double b00 = B[k * n + j];
				register double b01 = B[k * n + (j + 1)];
				register double b10 = B[(k + 1) * n + j];
				register double b11 = B[(k + 1) * n + (j + 1)];
				c00 = a00 * b00 + a01 * b10 + c00;
				c01 = a00 * b01 + a01 * b11 + c01;
				c10 = a10 * b00 + a11 * b10 + c10;
				c11 = a10 * b01 + a11 * b11 + c11;
			}
			C[i * n + j] = c00;
			C[i * n + (j + 1)] = c01;
			C[(i + 1) * n + j] = c10;
			C[(i + 1) * n + (j + 1) = ]c11;
		}
}

void dgemm3(const double *A, const double *B, double *C, const int n) 
{
	for (int i = 0; i < n; i += 3)
		for (int j = 0; j < n; j += 3) 
		{
			int c0 = i * n + j; 
			int c1 = c0 + n;
			int c2 = c1 + n;
			register double c00 = C[c0]; 
			register double c01 = C[c0 + 1];
			register double c02 = C[c0 + 2];
			register double c10 = C[c1]; 
			register double c11 = C[c1 + 1];
			register double c12 = C[c1 + 2];
			register double c20 = C[c2];
			register double c21 = C[c2 + 1];
			register double c22 = C[c2 + 2];

			for (int k = 0; k < n; k += 3) 
			{
				int a0 = i * n + k; 
				int a1 = a0 + n; 
				int a2 = a1 + n;
				int b0 = k * n + j;
				int b1 = b0 + n;
				int b2 = b1 + n;
				register double a00 = A[a0];
				register double a10 = A[a1];
				register double a20 = A[a2];
				register double b00 = B[b0]; register double b01 = B[b0 + 1]; register double b02 = B[b0 + 2];

				c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
				c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
				c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

				a00 = A[a0 + 1];
				a10 = A[a1 + 1];
				a20 = A[a2 + 1];
				b00 = B[b1];b01 = B[b1 + 1]; b02 = B[b1 + 2];

				c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
				c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
				c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

				a00 = A[a0 + 2];
				a10 = A[a1 + 2];
				a20 = A[a2 + 2];
				b00 = B[b2]; b01 = B[b2 + 1]; b02 = B[b2 + 2];

				c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
				c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
				c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

			}
			C[c0] = c00;
			C[c0 + 1] = c01;
			C[c0 + 2] = c02;
			C[c1] = c10;
			C[c1 + 1] = c11;
			C[c1 + 2] = c12;
			C[c2] = c20;
			C[c2 + 1] = c21;
			C[c2 + 2] = c22;
		}
}

void ijk(const double *A, const double *B, double *C, const int n) 
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			register double r = C[i * n + j];
			for (int k = 0; k < n; k++)
				r += A[i * n + k] * B[k * n + j];
			C[i * n + j] = r;
		}
	}
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
	for (int i = 0; i < n; i += b)
		for (int j = 0; j < n; j += b)
			for (int k = 0; k < n; k += b)
				/* B x B mini matrix multiplications */
				for (int i1 = i; (i1 < i + b) && (i1 < n); i1++)
					for (int j1 = j; (j1 < j + b) && (j1 < n); j1++)
					{
						register double r = C[i1 * n + j1];
						for (int k1 = k; (k1 < k + b) && (k1 < n); k1++)
							r += A[i1 * n + k1] * B[k1 * n + j1];
						C[i1 * n + j1] = r;
					}
}

void jik(const double *A, const double *B, double *C, const int n) 
{
	for (int j = 0; j < n; j++) {
		 for (int i = 0; i < n; i++){
			register double r = C[i * n + j];
			for (int k = 0; k < n; k++)
				r += A[i * n + k] * B[k * n + j];
			C[i * n + j] = r;
		}
	}
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
	for (int j = 0; j < n; j += b)
		for (int i = 0; i < n; i += b)
			for (int k = 0; k < n; k += b)
				/* B x B mini matrix multiplications */
				for (int j1 = j; (j1 < j + b) && (j1 < n); j1++)
					for (int i1 = i; (i1 < i + b) && (i1 < n); i1++)
					{
						register double r = C[i1 * n + j1];
						for (int k1 = k; (k1 < k + b) && (k1 < n); k1++)
							r += A[i1 * n + k1] * B[k1 * n + j1];
						C[i1 * n + j1] = r;
					}
}

void kij(const double *A, const double *B, double *C, const int n) 
{
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			register double r = A[i * n + k];
			for (int j = 0; j < n; j++)
				C[i * n + j] += r * B[k * n + j];
		}
	}
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
	for (int k = 0; k < n; k += b)
		for (int i = 0; i < n; i += b)
			for (int j = 0; j < n; j += b)
				/* B x B mini matrix multiplications */
				for (int k1 = k; (k1 < k + b) && (k1 < n); k1++)
					for (int i1 = i; (i1 < i + b) && (i1 < n); i1++)
					{
						register double r = A[i1 * n + k1];
						for (int j1 = j; (j1 < j + b) && (j1 < n); j1++)
							C[i1 * n + j1] += r * B[k1 * n + j1];
					}
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			register double r = A[i * n + k];
			for (int j = 0; j < n; j++)
				C[i * n + j] += r * B[k * n + j];
		}
	}
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
	for (int i = 0; i < n; i += b)
		for (int k = 0; k < n; k += b)
			for (int j = 0; j < n; j += b)
				/* B x B mini matrix multiplications */
				for (int i1 = i; (i1 < i + b) && (i1 < n); i1++)
					for (int k1 = k; (k1 < k + b) && (k1 < n); k1++)
					{
						register double r = A[i1 * n + k1];
						for (int j1 = j; (j1 < j + b) && (j1 < n); j1++)
							C[i1 * n + j1] += r * B[k1 * n + j1];
					}
}

void jki(const double *A, const double *B, double *C, const int n) 
{
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			register double r = B[k * n + j];
			for (int j = 0; j < n; j++)
				C[i * n + j] += A[i * n + k] * r;
		}
	}
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
	for (int j = 0; j < n; j += b)
		for (int k = 0; k < n; k += b)
			for (int i = 0; i < n; i += b)
				/* B x B mini matrix multiplications */
				for (int j1 = j; (j1 < j + b) && (j1 < n); j1++)
					for (int k1 = k; (k1 < k + b) && (k1 < n); k1++)
					{
						register double r = B[k1 * n + j1];
						for (int i1 = i; (i1 < i + b) && (i1 < n); i1++)
							C[i1 * n + j1] += A[i1 * n + k1] * r;
					}
}

void kji(const double *A, const double *B, double *C, const int n) 
{
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			register double r = B[k * n + j];
			for (int j = 0; j < n; j++)
				C[i * n + j] += A[i * n + k] * r;
		}
	}
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
	for (int k = 0; k < n; k += b)
		for (int j = 0; j < n; j += b)
			for (int i = 0; i < n; i += b)
				/* B x B mini matrix multiplications */
				for (int k1 = k; (k1 < k + b) && (k1 < n); k1++)
					for (int j1 = j; (j1 < j + b) && (j1 < n); j1++)
					{
						register double r = B[k1 * n + j1];
						for (int i1 = i; (i1 < i + b) && (i1 < n); i1++)
							C[i1 * n + j1] += A[i1 * n + k1] * r;
					}
}

void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
	for (int i = 0; i < n; i += b)
		for (int j = 0; j < n; j += b)
			for (int k = 0; k < n; k += b)
				/* B x B mini matrix multiplications */
				for (int i1 = i; (i1 < i + b) && (i1 < n); i1+3)
					for (int j1 = j; (j1 < j + b) && (j1 < n); j1+3)
					{
						int c0 = i * n + j;
						int c1 = c0 + n;
						int c2 = c1 + n;
						register double c00 = C[c0];
						register double c01 = C[c0 + 1];
						register double c02 = C[c0 + 2];
						register double c10 = C[c1];
						register double c11 = C[c1 + 1];
						register double c12 = C[c1 + 2];
						register double c20 = C[c2];
						register double c21 = C[c2 + 1];
						register double c22 = C[c2 + 2];

						for (int k1 = k; (k1 < k + b) && (k1 < n); k1+3)
						{
							int a0 = i * n + k;
							int a1 = a0 + n;
							int a2 = a1 + n;
							int b0 = k * n + j;
							int b1 = b0 + n;
							int b2 = b1 + n;
							register double a00 = A[a0];
							register double a10 = A[a1];
							register double a20 = A[a2];
							register double b00 = B[b0]; register double b01 = B[b0 + 1]; register double b02 = B[b0 + 2];

							c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
							c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
							c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

							a00 = A[a0 + 1];
							a10 = A[a1 + 1];
							a20 = A[a2 + 1];
							b00 = B[b1]; b01 = B[b1 + 1]; b02 = B[b1 + 2];

							c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
							c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
							c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

							a00 = A[a0 + 2];
							a10 = A[a1 + 2];
							a20 = A[a2 + 2];
							b00 = B[b2]; b01 = B[b2 + 1]; b02 = B[b2 + 2];

							c00 += a00 * b00; c01 += a00 * b01; c02 += a00 * b02;
							c10 += a10 * b00; c11 += a10 * b01; c12 += a10 * b02;
							c20 += a20 * b00; c21 += a20 * b01; c22 += a20 * b02;

						}
						C[c0] = c00;
						C[c0 + 1] = c01;
						C[c0 + 2] = c02;
						C[c1] = c10;
						C[c1 + 1] = c11;
						C[c1 + 2] = c12;
						C[c2] = c20;
						C[c2 + 1] = c21;
						C[c2 + 2] = c22;
					}
}