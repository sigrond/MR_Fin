#include <complex>
#include <cmath>
#include <cstdio>
#include"globals.h"

const std::complex<float> I(0.0f,1.0f);

void calculateMieAB(int * nMaxTable, int nPiiTau, int sizeR, std::complex<float> m, float * x, float * afloat, float * aImag, float * bfloat, float * bImag) {
	std::complex<float> a;
	std::complex<float> b;
	std::complex<float> *Theta=new std::complex<float> [nPiiTau];
	std::complex<float> *Eta = new std::complex<float>[nPiiTau];
	std::complex<float> *Psi = new std::complex<float>[nPiiTau];
#pragma omp parallel for private(a, b, Theta, Eta, Psi) default(none) shared(sizeR, nMaxTable, x, m, afloat, nPiiTau, aImag, bfloat, bImag)
	for(int k=0;k<sizeR;k++) {
		const std::complex<float> invM = 1.0f/m;
		const std::complex<float> invX = 1.0f/x[k];
		int Nadn;
		const std::complex<float> r = x[k]*m;
		int const Nmax = nMaxTable[k];
		float j=abs(r);
		if (((int)j)>Nmax) Nadn=(int)ceil(j)+15;
		else Nadn=Nmax+15;

		std::complex<float> *D=new std::complex<float> [Nadn+1];

		/* Calculating D */
		D[Nadn]=0.0f;
		for (int i=Nadn;i>=1;--i) {
			const std::complex<float> aux = (std::complex<float>) i/r;
			D[i - 1] =  aux - 1.0f / (D[i] + aux); 
		}
		/*initial values */	
		Theta[0]=sin(x[k]);
		Theta[1]=Theta[0]/(x[k])-cos(x[k]);

		Eta[0]=cos(x[k]);
		Eta[1]=Eta[0]/(x[k])+sin(x[k]);

		Psi[0]=Theta[0]-I*Eta[0];
		Psi[1]=Theta[1]-I*Eta[1];

		a = ((D[1] * invM + (float)1 * invX) * Theta[1] - Theta[0])
			/ ((D[1] * invM + (float)1 * invX) * Psi[1] - Psi[0]);
		b = ((D[1] * m + (float)1 * invX) * Theta[1] - Theta[0])
			/ ((D[1] * m + (float)1 * invX) * Psi[1] - Psi[0]);
		afloat[nPiiTau*k] = a.real();
		aImag[nPiiTau*k] = imag(a);
		bfloat[nPiiTau*k] = b.real();
		bImag[nPiiTau*k] = imag(b);
		for (int i=2;i<=Nmax;++i) {
			const std::complex<float> aux = (2.0f*i - 1.0f)*invX;
			const std::complex<float> aux2 = (std::complex<float>)i*invX;
			const std::complex<float> aux3 = D[i] * invM + aux2;
			const std::complex<float> aux4 = D[i] * m + aux2;

			Theta[i] =  aux * Theta[i - 1] - Theta[i - 2];
			Eta[i]   =  aux * Eta[i - 1] - Eta[i - 2];
			Psi[i]   = Theta[i] - I*Eta[i];
			a=(aux3 * Theta[i] - Theta[i-1])
				/ (aux3 * Psi[i] - Psi[i-1]);
			b=(aux4 * Theta[i] - Theta[i-1])
				/ (aux4 * Psi[i] - Psi[i-1]);
			int index = i-1 + nPiiTau * k;
			afloat[index] = real(a);
			aImag[index] = imag(a);
			bfloat[index] = real(b);
			bImag[index] = imag(b);
		}
	}
}
