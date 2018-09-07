/* ----------------------------------------------------
	�m�C�Y�w���̌v�Z
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (noise_figure2_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "parameter2_0.h"
#include "common_function2_0.h" 
#include "new.h"
#include "y_parameter2_0.h"

/* ********************************** */
/* ----- �m�C�Y�w����Mian�֐� ----- */
void NoiseFigure(struct Spectrum si[], struct YPara y[], 
				 double fr[],struct CircuitPara circ[], struct SetPara para){

	/* ================================================================= */
	/*** --- �ϐ��錾 --- ***/
	double *p, *r, *c, *fmin, *gn,*gass;
	struct ComplexNumbers *zopt, *zc;

	/*** --- �v���g�^�C�v�錾 --- ***/
	void MinimumNF_Gonzalez(struct Spectrum si[], struct YPara y[], 
							double p[], double r[], double c[], double fmin[],
							double fr[], double gn[], struct ComplexNumbers zopt[],
							struct ComplexNumbers zc[],double gass[],struct CircuitPara circ[], struct SetPara para);
	void MinimumNF_Cappy(struct Spectrum si[], struct YPara y[], 
						 double p[], double r[], double c[], double fmin[],
						 double fr[], double gn[], struct ComplexNumbers zopt[], 
						 struct ComplexNumbers zc[],double gass[],struct CircuitPara circ[], struct SetPara para);
	void WriteNoise(double p[], double r[], double c[], double fmin[],
					double fr[], double gn[], struct ComplexNumbers zopt[], 
					struct ComplexNumbers zc[],double gass[],struct CircuitPara circ[], struct SetPara para);

	/*** ----- �ǂݍ��񂾃p�����[�^�ɏ]���A�z��̃��������m�� */
	_set_new_handler(error);
	p = new double[para.data];
	r = new double[para.data];
	c = new double[para.data];
	fmin = new double[para.data];
	gn = new double[para.data];
	zopt = new ComplexNumbers[para.data];
	zc = new ComplexNumbers[para.data];
	gass = new double[para.data];
	/* ================================================================= */


	/* --- �G���p�����[�^�̌v�Z --- */
	/* Cappy��Ă̎� */
	if(para.nway==1){
		MinimumNF_Cappy(si,	y,
						p, r, c, fmin, fr, gn, zopt, zc,gass,circ, para);
	}
	/* Gonzalez��Ă̎� */
  	else if(para.nway==0){
		MinimumNF_Gonzalez(si,	y, 
							p, r, c, fmin, fr, gn, zopt, zc,gass,circ, para);
	}

	/* --- �t�@�C���ɏ������� --- */
	WriteNoise(p, r, c, fmin, fr, gn, zopt, zc,gass,circ, para);

	/* ----- �������̊J�� */
	delete [] p;
	delete [] r;
	delete [] c;
	delete [] fmin;
	delete [] gn;
	delete [] zopt;

}



/* ******************************** */
/* ----- �m�C�Y�w������������ ----- */
void WriteNoise(double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[], 
				struct ComplexNumbers zc[],double gass[],struct CircuitPara circ[], struct SetPara para){
    FILE *fp;	//�ǉ�
	int jf;

    fp = OpenFile(para.noise, 'w');


	fprintf(fp, "f[Hz]	gm	fc	fmin	Gass	P	R	C	gn	ropt	xopt	rc	 \n");

	for(jf=0;jf<=para.data/2-1 ; jf++){
		if(fr[jf]>para.fmin_y){
			if(fr[jf]>para.fmax_y){break;}

			fprintf(fp, "%1.6e	", fr[jf]);
			fprintf(fp, "%f		%e ",  circ[jf].gm,  circ[jf].ft);
			fprintf(fp, "%lf	%lf	",fmin[jf],  gass[jf]);
			fprintf(fp, "%lf	%lf		%lf ",  p[jf], r[jf], c[jf]);
//			fprintf(fp, "%lf	",  fmin[jf]);
			fprintf(fp, "%lf	%lf		%lf	",  gn[jf], zopt[jf].re, zopt[jf].im);
			fprintf(fp, "%lf\n",  zc[jf].re);
			
		}
	}

	fclose(fp);
	printf("'%s' write\n", para.noise);

}
