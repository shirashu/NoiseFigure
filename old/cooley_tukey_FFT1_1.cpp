/* ----------------------------------------------
	FFT関数
	2002/11/14  by Masahiro Nakayama 
	Ver. 1.1 (Cooley_Tukey_FFT1_1.cpp)
---------------------------------------------- */

#include <stdio.h>
#include <math.h>

#include "parameter1_1.h"
#include "Cooley_Tukey_FFT1_1.h"

/*-------------------------------------------------------------*/
/* ----- Main関数 ----- */
int FFTCalculate(struct ComplexNumbers sif[], int bit_num){
	int bit_r[DATA];

	BitRev(bit_r, bit_num);	/*-- ビット反転を行う --*/
	BitRIndex(bit_r, sif, bit_num);	/*-- ビット反転インデックシングを行う --*/
	CooleyTukeyFFT(sif, bit_num);	/*-- FFT処理を行う --*/

	return 0;
}

/*--------------------------------------------------------------*/
/*---------------- Cooley Tukey FFT アルゴリズム ---------------*/
/* ----- ビット反転 ----- */
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

/* ----- ビット反転インデックシング ----- */
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

/* ----- FFT本体 ----- */
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
/* ---------------------- 各種 処理用 関数 ----------------------*/
/* ----- W(複素数)のべき乗を計算する ----- */
/* 整数部 */
double PowReW(double y, int data_num){
	y = cos(y*2*PAI/data_num);
	return y;
}
/* 虚数部 */
double PowImW(double y, int data_num){
	y = sin(y*2*PAI/data_num);
	return y;
}
