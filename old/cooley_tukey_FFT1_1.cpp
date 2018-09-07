/* ----------------------------------------------
	FFT�֐�
	2002/11/14  by Masahiro Nakayama 
	Ver. 1.1 (Cooley_Tukey_FFT1_1.cpp)
---------------------------------------------- */

#include <stdio.h>
#include <math.h>

#include "parameter1_1.h"
#include "Cooley_Tukey_FFT1_1.h"

/*-------------------------------------------------------------*/
/* ----- Main�֐� ----- */
int FFTCalculate(struct ComplexNumbers sif[], int bit_num){
	int bit_r[DATA];

	BitRev(bit_r, bit_num);	/*-- �r�b�g���]���s�� --*/
	BitRIndex(bit_r, sif, bit_num);	/*-- �r�b�g���]�C���f�b�N�V���O���s�� --*/
	CooleyTukeyFFT(sif, bit_num);	/*-- FFT�������s�� --*/

	return 0;
}

/*--------------------------------------------------------------*/
/*---------------- Cooley Tukey FFT �A���S���Y�� ---------------*/
/* ----- �r�b�g���] ----- */
void BitRev(int bit_r[], int bit_num){
	int n2, i, j, k, data_num;
	
	data_num = (int)pow(2, bit_num);

	n2 = data_num / 2;
	i = 0;	j = 0;
	
	for(;;){
		bit_r[i] = j;
		i = i + 1;
		if(i>=data_num){break;}
		else{
			k = n2;
			while(k<=j){
				j = j - k;
				k = k / 2;
			}
			j = j + k;
		}
	}
}

/* ----- �r�b�g���]�C���f�b�N�V���O ----- */
void BitRIndex(int bit_r[], struct ComplexNumbers x[], int bit_num){
	int i, j, data_num;
	struct ComplexNumbers t;

	data_num = (int)pow(2, bit_num);

	for(i=0; i<=data_num-1; i++){
		j=bit_r[i];
		if(j > i){
			t.re = x[i].re;
			t.im = x[i].im;
			x[i].re = x[j].re;
			x[i].im = x[j].im;
			x[j].re = t.re;
			x[j].im = t.im;
		}
	}
}

/* ----- FFT�{�� ----- */
void CooleyTukeyFFT(struct ComplexNumbers x[], int bit_num){
	int k, index0, index1, b, r, data_num;
	struct ComplexNumbers w0, w1, b0, b1;

	data_num = (int)pow(2, bit_num);

	for(k=0; k<=bit_num-1; k++){
		for(b=0; b<=pow(2,k)-1; b++){
			w0.re = PowReW(pow(2,bit_num-1-k) * b , data_num);
			w0.im = PowImW(pow(2,bit_num-1-k) * b , data_num );
			w1.re = PowReW(pow(2,bit_num-1-k) * (pow(2, k)+b) , data_num);
			w1.im = PowImW(pow(2,bit_num-1-k) * (pow(2, k)+b) , data_num);
			for(r=0; r<=pow(2,bit_num-1-k)-1; r++){
				index0 = r * (int)pow(2, k+1) + b;
				index1 = r * (int)pow(2, k+1) + (int)pow(2, k) + b;
				b0.re = x[index0].re;
				b0.im = x[index0].im;
				b1.re = x[index1].re;
				b1.im = x[index1].im;
				x[index0].re = b0.re + w0.re * b1.re - w0.im * b1.im;
				x[index0].im = b0.im + w0.im * b1.re + w0.re * b1.im;
				x[index1].re = b0.re + w1.re * b1.re - w1.im * b1.im;
				x[index1].im = b0.im + w1.im * b1.re + w1.re * b1.im;
			}
		}
	}
}

/*---------------------------------------------------------------*/
/* ---------------------- �e�� �����p �֐� ----------------------*/
/* ----- W(���f��)�ׂ̂�����v�Z���� ----- */
/* ������ */
double PowReW(double y, int data_num){
	y = cos(y*2*PAI/data_num);
	return y;
}
/* ������ */
double PowImW(double y, int data_num){
	y = sin(y*2*PAI/data_num);
	return y;
}
