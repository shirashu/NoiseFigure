/* --------------------------------------------------------------------
	Y�p�����[�^����ъe�함���萔�̌v�Z
	2003/4/11 by Masahiro Nakayama 
	Ver. 2.0 (y_parameter2_0.cpp)

	Y�p�����[�^�v�Z���̃t�[���G�ϊ���ύX
-------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_1.h"
#include "common_function1_1.h"
#include "y_parameter1_1.h"
#include "cooley_tukey_FFT1_1.h"

/* ************************************* 
     Y�p�����[�^�v�Z��Mian�֐�
************************************** */
void YParameter(double cur[][JT][IN0], double ss[][IN0], double rex[][IN0],
				struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				struct ComplexNumbers y21[], struct ComplexNumbers y22[],
				double fr[], char *fn_yparam, char *fn_out){
	int jf;
	double dv[2];
	double g[IN0][IN0][FIN], b[IN0][IN0][FIN];
	double cgd[FIN], cgs[FIN];
	double ri[FIN], gm0[FIN], tau[FIN], cds[FIN], gds[FIN];

	// �ϐ��̏�����
	dv[0] = DV0;	dv[1] = DV1;

	// Y�p�����[�^�[�v�Z
	CalculateYPara_inte(cur, g, b,	rex, ss, dv, fr);
	// �O�p�`�E�B���h�ɂ�镽��
	if(YAVE!=0){
		FlatYPara_inte(g, b, fr);
	}
	// �����萔�v�Z
	CalculateOtherPara(g, b, cgd, cgs, ri, gm0, tau, cds, gds, fr);

	// Y�p�����[�^�[���֐��󂯓n���p�ϐ��ɑ��
	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} //�ő���g���Ōv�Z���X�g�b�v
		y11[jf].re = g[0][0][jf];
		y21[jf].re = g[1][0][jf];
		y11[jf].im = b[0][0][jf];
		y21[jf].im = b[1][0][jf];

		y12[jf].re = g[0][1][jf];
		y22[jf].re = g[1][1][jf];
		y12[jf].im = b[0][1][jf];
		y22[jf].im = b[1][1][jf];
	}

	// �t�@�C���ɏ�������
	WriteParameter(g, b, cgd, cgs, ri, gm0, tau, cds, gds, fr, fn_yparam, fn_out);

}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  �ȉ��e��֐�
------------------------------------------------------------------------- */

/* ************************************** */
/* --- ���l�ϕ��@ ----------------------- */
/* --- Y�p�����[�^���v�Z ---------------- */
void CalculateYPara_inte(double cur[][JT][IN0], double g[][IN0][FIN],
						 double b[][IN0][FIN], double rex[][IN0],
						 double ss[IN0][IN0], double dv[], double fr[]){
	int i, ii, jt, jf;
	double w, s1, s2, r1, r2, t, r, x, a;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} //�ő���g���Ōv�Z���X�g�b�v
		w = 2*PAI*fr[jf];

		// --- 20030411�ǉ�v2.0(���R)
		r = (1 + w*DT*sin(w*DT) - cos(w*DT));
		x = (+w*DT*cos(w*DT) + sin(w*DT));
		a = r*r + x*x;
		r = -w*w * r;
		x =  w*w * x;
		// --- 20030411�ǉ�v2.0

		for(ii=0; ii<=1; ii++){
		for(i=0 ; i<=IN0-1; i++){

			/* --- 20030411�ύXv2.0(���R)
			s1 = cur[i][JTL0-1][ii]*sin(w*double(JTL0-1)*DT)/2;
			r1 = 0.0;
			s2 = cur[i][JTL0-1][ii]*cos(w*double(JTL0-1)*DT)/2;
			r2 = 0.0;*/
			s1 = cur[i][JTL0-1][ii] * (r*cos(w*double(JTL0-1)*DT) + x*sin(w*double(JTL0-1)*DT))/2;
			r1 = 0.0;
			s2 = cur[i][JTL0-1][ii] * (x*cos(w*double(JTL0-1)*DT) - r*sin(w*double(JTL0-1)*DT))/2;
			r2 = 0.0;
			// --- 20030411�ǉ�v2.0
			
			for(jt=1; jt<=JTL0-2; jt++){
		        t = DT*double(jt);
				/* --- 20030411�ύXv2.0(���R)
				Sigma(cur[i][jt][ii]*sin(w*t), &s1, &r1);
				Sigma(cur[i][jt][ii]*cos(w*t), &s2, &r2);*/
				Sigma(cur[i][jt][ii]*(r*cos(w*t) + x*sin(w*t)), &s1, &r1);
				Sigma(cur[i][jt][ii]*(x*cos(w*t) - r*sin(w*t)), &s2, &r2);
				// --- 20030411�ǉ�v2.0
			}
			/* --- 20030411�ύXv2.0(���R)
			g[i][ii][jf] = (rex[i][ii] - ss[i][ii])/dv[ii] + w*s1*DT/dv[ii];
			b[i][ii][jf] = w*s2*DT/dv[ii];*/
			g[i][ii][jf] = DT/dv[ii]/a * ( x/w*(rex[i][ii] - ss[i][ii]) + s1*DT);
			b[i][ii][jf] = DT/dv[ii]/a * (-r/w*(rex[i][ii] - ss[i][ii]) + s2*DT);
			// --- 20030411�ǉ�v2.0
		}
		}
	}
}


/* **************************************************** */
/* --- ���l�ϕ��@ ------------------------------------- */
/* --- Y�p�����[�^�𕽊����O�p�`�E�B���h�ɂ�镽���@--- */
void FlatYPara_inte(double g[][IN0][FIN], double b[][IN0][FIN], double fr[]){
	int jf, k, i, ii, l[FIN], fmax;
	double t1[FIN], t2[FIN], rs1, rs2;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){fmax=jf-1; break;} //�ő���g���Ōv�Z���X�g�b�v
	}

	for(ii=0; ii<=1; ii++){
	for(i=0 ; i<=IN0-1; i++){

		for(jf=0; jf<=fmax; jf++){
			t1[jf] = 0; t2[jf] = 0;
			rs1 = 0;	rs2 = 0;

			l[jf] = YAVE;
			if(jf-l[jf]<=0)    {do{l[jf]--;}while(jf-l[jf]<=0);}
			else if(jf+l[jf]>fmax)    {do{l[jf]--;}while(jf+l[jf]>fmax);}

			for(k=-l[jf]; k<=l[jf]; k++){
				Sigma((l[jf]-abs(k))*g[i][ii][jf-k], &t1[jf], &rs1);
				Sigma((l[jf]-abs(k))*b[i][ii][jf-k], &t2[jf], &rs2);
			}
		}
		for(jf=0; jf<=fmax; jf++){
			if(l[jf] == 0){continue;}
			g[i][ii][jf] = t1[jf]/l[jf]/l[jf];
			b[i][ii][jf] = t2[jf]/l[jf]/l[jf];
		}

	}
	}
}

/* **************************** */
/* --- �e�함���萔�����߂� --- */
void CalculateOtherPara(double g[][IN0][FIN], double b[][IN0][FIN],
						double cgd[], double cgs[],	double ri[],
						double gm0[], double tau[], double cds[],
						double gds[], double fr[]){
	int jf;
	double g11, b11, g12, b12, g21, b21, g22, b22;
	double d, sd, sg, w;

	for(jf=0; ; jf++){
		if(fr[jf]<=0){continue;} //0GHz�ȉ��̎��g���ł͌v�Z���Ȃ�
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
/* ----- �p�����[�^�[���������� ----- */
void WriteParameter(double g[][IN0][FIN], double b[][IN0][FIN],
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
		else if(fr[jf]>1.05e+11){break;} //�ő���g���Ōv�Z���X�g�b�v

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
