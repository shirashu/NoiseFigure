/* --------------------------------------------------------------------
	Y�p�����[�^����ъe�함���萔�̌v�Z
	2002/11/13 by Masahiro Nakayama 
	Ver. 1.0 (y_parameter.cpp)
-------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_0.h"
#include "common_function1_0.h"
#include "y_parameter1_0.h"
#include "cooley_tukey_FFT1_0.h"

/* ************************************* 
     Y�p�����[�^�v�Z��Mian�֐�
************************************** */
void YParameter(double cur1[][JT][IN0], double ss[][IN0], double rex[][IN0],
				struct ComplexNumbers y11[], struct ComplexNumbers y21[],
				double fr[], char *fn_yparam, char *fn_out){
	int jf;
	double dv[2];
	double g[IN0][IN0][DATA], b[IN0][IN0][DATA];
	double cgd[50], cgs[50];
	double ri[50], gm0[50], tau[50], cds[50], gds[50];
	struct ComplexNumbers cur_f[IN0][DATA][IN0];

	// ���� --- �ϐ��̏�����
	dv[0] = DV0;	dv[1] = DV1;

	/* -----------------------------------------------------*/
	/* --- ���l�ϕ��@��p����Y�p�����[�^�v�Z*/
	/* -----------------------------------------------------*/
	if(WAY==0){
		// Y�p�����[�^�[�v�Z
		CalculateYPara_inte(cur1, g, b,	rex, ss, dv, fr);
		// �����萔�v�Z
		CalculateOtherPara(g, b, cgd, cgs, ri, gm0, tau, cds, gds, fr);
		//�O�p�`�E�B���h�ɂ�镽��
		FlatYPara_inte(g, b, fr); 
	}

	/* -----------------------------------------------------*/
	/* --- FFT�@��p����Y�p�����[�^�v�Z*/
	/* -----------------------------------------------------*/
	if(WAY==1){
		// �d���ӂ����FFT�ɂ�����
		CurrentFFT(cur1, cur_f);
		// Y�p�����[�^�[�v�Z
		CalculateYPara_FFT(cur_f, g, b,	rex, ss, dv, fr);
		// �����萔�v�Z
		CalculateOtherPara(g, b, cgd, cgs, ri, gm0, tau, cds, gds, fr);
	}

	// ���� --- Y�p�����[�^�[���֐��󂯓n���p�ϐ��ɑ��
	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} //�ő���g���Ōv�Z���X�g�b�v
		y11[jf].re = g[0][0][jf];
		y11[jf].im = b[0][0][jf];
		y21[jf].re = g[1][0][jf];
		y21[jf].im = b[1][0][jf];
	}

	// ���� --- �t�@�C���ɏ�������
	WriteParameter(g, b, cgd, cgs, ri, gm0, tau, cds, gds, fr, fn_yparam, fn_out);

}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  �ȉ��e��֐�
------------------------------------------------------------------------- */

/* ************************************** */
/* --- ���l�ϕ��@ ----------------------- */
/* --- Y�p�����[�^���v�Z ---------------- */
void CalculateYPara_inte(double cur1[][JT][IN0], double g[][IN0][DATA],
						 double b[][IN0][DATA], double rex[][IN0],
						 double ss[IN0][IN0], double dv[], double fr[]){
	int i, ii, jt, jf;
	double w, s1, s2, r1, r2, t;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} //�ő���g���Ōv�Z���X�g�b�v
		w = 2*PAI*fr[jf];
		for(ii=0; ii<=1; ii++){
		for(i=0 ; i<=IN0-1; i++){
			s1 = cur1[i][JTL0-1][ii]*sin(w*double(JTL0)*DT)/2;
			r1 = 0.0;
			s2 = cur1[i][JTL0-1][ii]*cos(w*double(JTL0)*DT)/2;
			r2 = 0.0;
			for(jt=1; jt<=JTL0-2; jt++){
		        t = DT*double(jt);
				Sigma(cur1[i][jt][ii]*sin(w*t), &s1, &r1);
				Sigma(cur1[i][jt][ii]*cos(w*t), &s2, &r2);
			}
			g[i][ii][jf] = (rex[i][ii] - ss[i][ii])/dv[ii] + w*s1*DT/dv[ii];
			b[i][ii][jf] = w*s2*DT/dv[ii];
		}
		}
	}
}

/* **************************************************** */
/* --- ���l�ϕ��@ ------------------------------------- */
/* --- Y�p�����[�^�𕽊����O�p�`�E�B���h�ɂ�镽���@--- */
void FlatYPara_inte(double g[][IN0][DATA], double b[][IN0][DATA], double fr[]){
	int jf, k, i, ii, l1[DATA], l2[DATA], fmax;
	double t1[DATA], t2[DATA], rs1, rs2;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){fmax=jf-1; break;} //�ő���g���Ōv�Z���X�g�b�v
	}
	
	for(ii=0; ii<=1; ii++){
	for(i=0 ; i<=IN0-1; i++){
		for(jf=0; jf<=fmax; jf++){
			rs1 = 0;	rs2 = 0;
			t1[jf] = 0; t2[jf] = 0;
			l1[jf] = YAVE; l2[jf] = YAVE;
			if(jf-l1[jf]+1<0)    {do{l1[jf]--;}while(jf-l1[jf]+1<0);}
			if(jf+l2[jf]-1>fmax)    {do{l2[jf]--;}while(jf+l2[jf]-1>fmax);}
			for(k=-l2[jf]+1; k<=l1[jf]-1; k++){
				Sigma((YAVE-abs(k))*g[i][ii][jf-k], &t1[jf], &rs1);
				Sigma((YAVE-abs(k))*b[i][ii][jf-k], &t2[jf], &rs2);
			}
		}
		for(jf=0; jf<=fmax; jf++){
			g[i][ii][jf] = t1[jf]/YAVE/(l1[jf]+l2[jf])*2;
			b[i][ii][jf] = t2[jf]/YAVE/(l1[jf]+l2[jf])*2;
		}
	}
	}
}



/* ************************* */
/* --- FFT�@ --------------- */
/* --- �d����FFT�ɂ����� --- */
void CurrentFFT(double cur1[][JT][IN0],
				struct ComplexNumbers cur_f[][DATA][IN0]){
	struct ComplexNumbers temp[DATA];
	int i, ii;

	/* --- �d���ӂ����FFT�ɂ����� --- */
	for(ii=0; ii<=1; ii++){
	for(i=0 ; i<=IN0-1; i++){
		ConverterToComp(temp, cur1, i, ii);
		FFTCalculate(&temp[0], JBIT);
		ConverterFromComp(temp, cur_f, i, ii);
	}
	}

}

/* ****************************************** */
/* --- FFT�@ -------------------------------- */
/* ----- �ӂ���d���𕡑f���̒l�ɕϊ� ----- */
void ConverterToComp(struct ComplexNumbers temp[], double cur1[][JT][IN0],
					 int i, int ii){
	int k;
	for(k=0; k<=DATA-1; k++){
		temp[k].re = cur1[i][k][ii];
		temp[k].im = 0;
	}
}
void ConverterFromComp(struct ComplexNumbers temp[], 
					   struct ComplexNumbers cur_f[][DATA][IN0], int i, int ii){
	int k;
	for(k=0; k<=DATA/2-1; k++){
		cur_f[i][k][ii].re = temp[k].re;
		cur_f[i][k][ii].im = temp[k].im;
	}
}

/* ********************************** */
/* --- FFT�@ ------------------------ */
/* --- Y�p�����[�^���v�Z ------------ */
void CalculateYPara_FFT(struct ComplexNumbers cur_f[][DATA][IN0],
						double g[][IN0][DATA], double b[][IN0][DATA],
						double rex[][IN0], double ss[][IN0],
						double dv[], double fr[]){
	int i, ii, jf;
	double w;

	for(jf=0; jf<=DATA/2-1; jf++){
		w = 2*PAI*fr[jf];
		for(ii=0; ii<=1; ii++){
		for(i=0 ; i<=IN0-1; i++){
			g[i][ii][jf] = (rex[i][ii] - ss[i][ii])/dv[ii] + w*cur_f[i][jf][ii].im*DT/dv[ii];
			b[i][ii][jf] = w*cur_f[i][jf][ii].re*DT/dv[ii];
		}
		}
	}
}



/* **************************** */
/* ���� ----------------------- */
/* --- �e�함���萔�����߂� --- */
void CalculateOtherPara(double g[][IN0][DATA], double b[][IN0][DATA],
						double cgd[], double cgs[],	double ri[],
						double gm0[], double tau[], double cds[],
						double gds[], double fr[]){
	int jf;
	double g11, b11, g12, b12, g21, b21, g22, b22;
	double d, sd, sg, w;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} //�ő���g���Ōv�Z���X�g�b�v
		w = 2*PAI*fr[jf];
		g11 = g[0][0][jf];		b11 = b[0][0][jf];
		g12 = g[0][1][jf];		b12 = b[0][1][jf];
		g21 = g[1][0][jf];		b21 = b[1][0][jf];
		g22 = g[1][1][jf];		b22 = b[1][1][jf];

		cgd[jf] = -b12/w;
		d   = b11+b12;
		sd  = d * d;
		sg  = g11*g11;
		cgs[jf] = d/w*(1+sg/sd);
		ri[jf]  = g11/(g11*g11 +b11*b11 +b12*b12 -2*b11*b12);
		gm0[jf] = sqrt((g21*g21+(b21-b12)*(b21-b12)) * (1+(w*cgs[jf]*ri[jf])*(w*cgs[jf]*ri[jf])));
		tau[jf] = 1/w * acos((g21 - w*b21*ri[jf]*cgs[jf] + w*b12*cgs[jf]*ri[jf]) / gm0[jf]);
		cds[jf] = (b22+b12)/w;
		gds[jf] = g22;
	}
}

/* ********************************** */
/* --- ���� ------------------------ */
/* ----- �p�����[�^�[���������� ----- */
void WriteParameter(double g[][IN0][DATA], double b[][IN0][DATA],
					double cgd[], double cgs[], double ri[],
					double gm0[], double tau[], double cds[],
					double gds[], double fr[], char *fn_yparam, char *fn_out){
	int jf;
	FILE *fp_yparam, *fp_out;

	// --- �t�@�C���I�[�v��
    fp_yparam = OpenFile(fn_yparam, 'w');
    fp_out    = OpenFile(fn_out, 'w');

	fprintf(fp_yparam, "f[Hz]	Re[Y11]	Im[Y11]	Re[Y12]	Im[Y12]	Re[Y21]	Im[Y21]	Re[Y22]	Im[Y22]\n");
	fprintf(fp_out, "f[Hz]	Cgd[F/m]	Cgs[F/m]	Cds[F/m]	Ri[ohm/m]	gm0[S/m]	gds[S/m]	tau[s]\n");

	for(jf=0; ; jf++){
		if(fr[jf]<=0){continue;}
		else if(fr[jf]>1.05e+11){break;}
		fprintf(fp_yparam, "%1.7e	" , fr[jf]);
		fprintf(fp_yparam, "%lf	%lf	" , g[0][0][jf], b[0][0][jf]);
		fprintf(fp_yparam, "%lf	%lf	" , g[0][1][jf], b[0][1][jf]);
		fprintf(fp_yparam, "%lf	%lf	" , g[1][0][jf], b[1][0][jf]);
		fprintf(fp_yparam, "%lf	%lf\n", g[1][1][jf], b[1][1][jf]);

		fprintf(fp_out, "%1.7e	%1.7e	%1.7e	%1.7e	", fr[jf], cgd[jf], cgs[jf], cds[jf]);
		fprintf(fp_out, "%1.7e	%lf	%lf	%1.7e\n", ri[jf], gm0[jf], gds[jf], tau[jf]);
	}

	// --- �t�@�C���N���[�Y
	fclose(fp_yparam);
	fclose(fp_out);

	printf("'%s' write\n", fn_yparam);
	printf("'%s' write\n", fn_out);
}
