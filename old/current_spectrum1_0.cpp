/* ----------------------------------------------------
	�X�y�N�g���v�Z�p�w�b�_�t�@�C��
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (correlation_function1_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_0.h" 
#include "common_function1_0.h" 
#include "current_spectrum1_0.h" 
#include "cooley_tukey_FFT1_0.h"


/* ***************************+************** */
/* ----- �d���̃X�y�N�g���v�Z��Mian�֐� ----- */
void CurrentSpectrum(double cur1[][JT][IN0], double corr[][NUMS],
					 double sig[], double sid[], struct ComplexNumbers sigd[],
					 double fr[], char *fp_name_spec){

	struct ComplexNumbers tg[DATA], td[DATA];

	/* -----------------------------------------------------*/
	/* --- ���֊֐�����t�[���G�ϊ��ɂ��X�y�N�g�����v�Z   */
	/* -----------------------------------------------------*/
	if(WAY==0){
		// --- �d���̂ӂ�����瑊�֊֐����v�Z
		CalculateCorrelation_inte(cur1, corr);
		// --- ���֊֐�����X�y�N�g�����v�Z
		CalculateSpectrum_inte(corr, sig, sid, sigd, fr);
		// ----- FILE�փX�y�N�g���Ƒ��֊֐�����������
		WriteSpectrum_inte(corr, sig, sid, sigd, fr, fp_name_spec);
	}

	/* -----------------------------------------------------*/
	/* --- FFT��p�����X�y�N�g���v�Z*/
	/* -----------------------------------------------------*/
	if(WAY==1){
		// --- �ӂ���d���𕡑f���̒l�ɕϊ�
		ConvertCompNumber(tg, td, cur1);
		// --- �X�y�N�g���E�B���h�E��������
		if(SWIN==1){
			SpectrumWindows(tg);
			SpectrumWindows(td);
		}
		// --- �d���ӂ����FFT�ɂ�����
		FFTCalculate(&tg[0], JBIT);
		FFTCalculate(&td[0], JBIT);
		// --- �X�y�N�g�����v�Z����
		CalculateSpectrum_FFT(sig, sid, sigd, tg, td);
		// --- �X�y�N�g�����O�p�g�E�B���h�E�ɂ��␳
		FlatSpectrum_FFT(sig, sid, sigd);
		// --- �X�y�N�g�����瑊�֊֐����v�Z����
		CalculateCorrelation_FFT(sig, sid, sigd, corr);
		// ----- FILE�փX�y�N�g���Ƒ��֊֐�����������
		WriteSpectrum_FFT(corr, sig, sid, sigd, fp_name_spec);
	}

}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  �ȉ��e��֐�
------------------------------------------------------------------------- */


/* ****************************************************** */
/* ���l�ϕ��@(���֊֐������߂Čv�Z) --------------------- */
/* ****************************************************** */
/* ****************************************************** */
/* ----- �d���̂ӂ�����璼�ڂɑ��֊֐����v�Z���� ----- */
void CalculateCorrelation_inte(double cur1[][JT][IN0], double corr[][NUMS]){
	int dt, jt;
	double temp1, temp2, rs1, rs2;

	// --- ���ȑ��֊֐��̌v�Z
	for(jt=0; jt<=MSMALLS; jt++){
		rs1 = 0;	rs2 = 0;	temp1 = 0;	temp2 = 0;
		for(dt=(JTL0-NUMS)+1; dt<=(JTL0-NUMS)+MLARGE-MSMALLS; dt++){
			Sigma(cur1[0][dt][STEP]*cur1[0][dt+jt][STEP], &temp1, &rs1);
			Sigma(cur1[1][dt][STEP]*cur1[1][dt+jt][STEP], &temp2, &rs2);
		}
		corr[0][jt] = temp1/(double)(MLARGE-MSMALLS)/EPP;
		corr[1][jt] = temp2/(double)(MLARGE-MSMALLS)/EPP;
	}

	// --- ���ݑ��֊֐��̌v�Z
	for(jt=0; jt<=MSMALLC; jt++){
		rs1 = 0;	rs2 = 0;	temp1 = 0;	temp2 = 0;
		for(dt=(JTL0-NUMS)+1; dt<=(JTL0-NUMS)+MLARGE-MSMALLC; dt++){
			Sigma(cur1[0][dt][STEP]*cur1[1][dt+jt][STEP], &temp1, &rs1);
			Sigma(cur1[1][dt][STEP]*cur1[0][dt+jt][STEP], &temp2, &rs2);
		}
		corr[2][MSMALLC+1+jt] = temp1/(double)(MLARGE-MSMALLC)/EPP;
		corr[2][MSMALLC  -jt] = temp2/(double)(MLARGE-MSMALLC)/EPP;
	}

}

/* ********************************************************** */
/* --- ���֊֐�����t�[���G�ϊ��ɂ��X�y�N�g�����v�Z���� --- */
void CalculateSpectrum_inte(double corr[][NUMS], double sig[], double sid[],
							struct ComplexNumbers sigd[], double fr[]){
	int jt, fmax, fmin, jf;
	double rs1, rs2, rs3, rs4;

	// --- �X�y�N�g���̌v�Z
	for(jf=0; ; jf++){
		if(fr[jf]==0){fmin=jf;} //���g��0Hz�̃C���f�b�N�X���L�^
		else if(fr[jf]==32168){fmax=jf-1; break;} //�ő���g���Ōv�Z���X�g�b�v
		sig[jf] = 0.0;		sid[jf] = 0.0;
		sigd[jf].re = 0.0;	sigd[jf].im = 0.0;
		rs1 = 0;	rs2 = 0;	rs3 = 0;	rs4 = 0;

		// --- ���ȑ��֊֐�����̎��ȃX�y�N�g���v�Z
		for(jt=-MSMALLS; jt<=MSMALLS; jt++){
			Sigma(corr[0][abs(jt)]*cos(2*PAI*fr[jf]*jt*DT), &sig[jf], &rs1);
			Sigma(corr[1][abs(jt)]*cos(2*PAI*fr[jf]*jt*DT), &sid[jf], &rs2);
		}
		sig[jf] = sig[jf]*2*DT;
		sid[jf] = sid[jf]*2*DT;


		// --- ���ݑ��֊֐�����̑��݃X�y�N�g���v�Z
		for(jt=-MSMALLC; jt<=MSMALLC; jt++){
			Sigma(corr[2][abs(MSMALLC+jt)]*cos(2*PAI*fr[jf]*jt*DT), &sigd[jf].re, &rs3);
			Sigma(corr[2][abs(MSMALLC+jt)]*sin(2*PAI*fr[jf]*jt*DT), &sigd[jf].im, &rs4);
		}
		sigd[jf].re = sigd[jf].re*2*DT;
		sigd[jf].im = sigd[jf].im*2*DT;
	}

	// ---  ���g��0Hz�ŃX�y�N�g��sig�Asigd��0�Ƃ��āA�␳
	for(jf=fmax; jf>=0; jf--){
		sig[jf] = sig[jf] - sig[fmin];
		sigd[jf].re = sigd[jf].re - sigd[fmin].re;
		sigd[jf].im = sigd[jf].im - sigd[fmin].im;
	}
}

/* ************************************************ */
/* ----- FILE�փX�y�N�g���Ƒ��֊֐����������� ----- */
void WriteSpectrum_inte(double corr[][NUMS], double sig[], double sid[],
						struct ComplexNumbers sigd[], double fr[], char *fp_name){
	int k, fmax, fmin;
	char fn_spe[]=SPE;
    FILE *fp;

	// --- �t�@�C���I�[�v��
	fp = OpenFile(fn_spe, 'a');

	fprintf(fp, "��I_G	��I_D	");
	fprintf(fp, "NUMS	MSMALL_Self	MSMALL_Cross	MLARGE	DATA\n");
	fprintf(fp, "%lf	%lf	",sqrt(corr[0][0]), sqrt(corr[1][0]));
	fprintf(fp, "%d	%d	%d	%d	%d\n", NUMS, MSMALLS, MSMALLC, MLARGE, DATA);

	fprintf(fp, "t[sec]	Cor_GD(t)	");
	fprintf(fp, "t[sec]	Cor_G(t)	Cor_D(t)	");
	fprintf(fp, "f[Hz]	Si_G(f)	Si_D(f)	Si_gd(f)[re]	Si_gd(f)[im]\n");

	for(k=0; ; k++){
		if(fr[k]==0){fmin=k; continue;}
		else if(fr[k]==1.05e+11){fmax=k; break;}
	}

	for(k=0; k<=fmax-fmin; k++){
		fprintf(fp, "%1.4e	%1.4e	", (k-MSMALLC)*DT, corr[2][k]);
		fprintf(fp, "%1.4e	%1.4e	%1.4e	", k*DT, corr[0][k], corr[1][k]);
		fprintf(fp, "%1.4e	%1.4e	%1.4e	%1.4e	%1.4e\n", fr[k+fmin] , sig[k+fmin], sid[k+fmin], sigd[k+fmin].re, sigd[k+fmin].im);
	}

	for(k=fmax-fmin+1; k<=MSMALLS; k++){
		fprintf(fp, "%1.4e	%1.4e	", (k-MSMALLC)*DT, corr[2][k]);
		fprintf(fp, "%1.4e	%1.4e	%1.4e\n", k*DT , corr[0][k], corr[1][k]);
	}
	for(k=MSMALLS+1; k<=2*MSMALLC; k++){
		fprintf(fp, "%1.4e	%1.4e\n", (k-MSMALLC)*DT, corr[2][k]);
	}

	// --- �t�@�C���N���[�Y
	fclose(fp);
	printf("'%s' write\n", fn_spe);
}







/* ****************************************************** */
/* FFT�@(�X�y�N�g���𒼐ڌv�Z) -------------- */
/* ****************************************************** */
/* ****************************************** */
/* ----- �ӂ���d���𕡑f���̒l�ɕϊ� ----- */
void ConvertCompNumber(struct ComplexNumbers tg[],
					   struct ComplexNumbers td[], double cur1[][JT][IN0]){
	int k;
	for(k=0; k<=DATA-1; k++){
		tg[k].re = cur1[0][JTL0-NUMS+k][STEP];
		tg[k].im = 0.0;
		td[k].re = cur1[1][JTL0-NUMS+k][STEP];
		td[k].im = 0.0;
	}
}

/* **************************************** */
/* ----- �X�y�N�g���E�B���h�E���|���� ----- */
void SpectrumWindows(struct ComplexNumbers t[]){
	int i;
	// 0�`T/10 �܂ł̃E�B���h
	for(i=0; i<=(int)(DATA/10); i++){
		t[i].re = t[i].re * (1-cos(PAI*10*i/DATA))/2;
	}

	// 9T/10�`T �܂ł̃E�B���h
	for(i=(int)(DATA*9/10); i<=DATA-1; i++){
		t[i].re = t[i].re * (1+cos(PAI*10*((i-(int)(DATA*9/10))/DATA)))/2;
	}
}

/* ******************************** */
/* ----- �X�y�N�g�����v�Z���� ----- */
void CalculateSpectrum_FFT(double sig[], double sid[], struct ComplexNumbers sigd[],
						   struct ComplexNumbers tg[], struct ComplexNumbers td[]){
	int jf;
	double st;

	if(SWIN==0){st = 1.0;}
	else if(SWIN==1){st = 0.875;}//�X�y�N�g���E�B���h�E���|�����ۂ̋��x����

	for(jf=0; jf<=DATA/2-1; jf++){
		sig[jf]    = (tg[jf].re*tg[jf].re + tg[jf].im*tg[jf].im)/st/EPP;
		sid[jf]    = (td[jf].re*td[jf].re + td[jf].im*td[jf].im)/st/EPP;
		sigd[jf].re= (tg[jf].re*td[jf].re + tg[jf].im*td[jf].im)/st/EPP;
		sigd[jf].im= (tg[jf].im*td[jf].re - tg[jf].re*td[jf].im)/st/EPP;
	}
}

/* ************************************************************ */
/* ----- �O�p�`�E�B���h�ɂ��X�y�N�g���̎��g���������s�� ----- */
void FlatSpectrum_FFT(double sig[], double sid[], struct ComplexNumbers sigd[]){
	int dt, jf;
	double t1[DATA], t2[DATA], t3[DATA], t4[DATA];
	double rs1, rs2, rs3, rs4;

	for(jf=0; jf<=DATA/2-1; jf++){
		rs1 = 0.0;	rs2 = 0.0;	rs3 = 0.0;	rs4 = 0.0;
		t1[jf]= 0.0; t2[jf]= 0.0; t3[jf]= 0.0; t4[jf]= 0.0; 
		for(dt=-SAVE+1; dt<=SAVE-1; dt++){
			if(jf-dt>=0 && jf-dt<=DATA/2-1){
				Sigma((SAVE-abs(dt))*sig[jf-dt], &t1[jf], &rs1);
				Sigma((SAVE-abs(dt))*sid[jf-dt], &t2[jf], &rs2);
				Sigma((SAVE-abs(dt))*sigd[jf-dt].re, &t3[jf], &rs3);
				Sigma((SAVE-abs(dt))*sigd[jf-dt].im, &t4[jf], &rs4);
			}
		}
		sig[jf] = t1[jf]/SAVE/SAVE;
		sid[jf] = t2[jf]/SAVE/SAVE;
		sigd[jf].re = t3[jf]/SAVE/SAVE;
		sigd[jf].im = t4[jf]/SAVE/SAVE;
	}
}

/* ******************************************** */
/* ----- �X�y�N�g�����瑊�֊֐����v�Z���� ----- */
void CalculateCorrelation_FFT(double sig[], double sid[], 
								struct ComplexNumbers sigd[], double corr[][NUMS]){
	int i;
	struct ComplexNumbers temp1[DATA], temp2[DATA], temp3[DATA], temp4[DATA];

	// �X�y�N�g�����ꎞ�ϐ��z��Ɋi�[
	for(i=0; i<=DATA/2-1; i++){
		temp1[i].re= sig[i]; //���ȑ��֊֐� sig
		temp1[i].im= 0;
		temp2[i].re= sid[i]; //���ȑ��֊֐� sig
		temp2[i].im= 0;
		temp3[i].re= sigd[i].re; //���ݑ��֊֐� sigd
		temp3[i].im= 0;
		temp4[i].re= sigd[i].im; //���ݑ��֊֐� sigd
		temp4[i].im= 0;

	}

	// �X�y�N�g����FFT�ɂ�����
	FFTCalculate(&temp1[0], JBIT-1);
	FFTCalculate(&temp2[0], JBIT-1);
	FFTCalculate(&temp3[0], JBIT-1);
	FFTCalculate(&temp4[0], JBIT-1);

	// �X�y�N�g����FFT���瑊�֊֐����v�Z����
	for(i=0; i<=DATA/2-1; i++) {
		if(i<=DATA/4-1){
			corr[0][i] = 2*temp1[i].re/(DATA*DATA);
			corr[1][i] = 2*temp2[i].re/(DATA*DATA);
			corr[2][i] = 2*(temp3[DATA/4-1-i].re+temp4[DATA/4-1-i].im)/(DATA*DATA);
		}
		else{
			corr[2][i] = 2*(temp3[i-DATA/4+1].re-temp4[i-DATA/4+1].im)/(DATA*DATA);
		}
	}

}

/* ************************************************ */
/* ----- FILE�փX�y�N�g���Ƒ��֊֐����������� ----- */
void WriteSpectrum_FFT(double corr[][NUMS], double sig[], double sid[], 
				   struct ComplexNumbers sigd[], char *fp_name) {
	int k;
    FILE *fp;

	fp = OpenFile(fp_name, 'a');

	fprintf(fp, "��I_G	��I_D	DATA\n");
	fprintf(fp, "%lf	%lf	%d\n",sqrt(corr[0][0]), sqrt(corr[1][0]), DATA);

	fprintf(fp, "f[Hz]	Si_G(f)	Si_D(f)	Si_G-D(f).re	Si_G-D(f).im");
	fprintf(fp, "	t[sec]	CrossCor(t)");
	fprintf(fp, "	t[sec]	Cor_G(t)	Cor_D(t)");

	fprintf(fp, "\n");

	for(k=0; k<=DATA/2-1; k++){
		// --- �X�y�N�g���̏�������
		fprintf(fp, "%1.4e	", k/DT/DATA);
		fprintf(fp, "%1.4e	" , sig[k]*DT/DATA);
		fprintf(fp, "%1.4e	" , sid[k]*DT/DATA);
		fprintf(fp, "%1.4e	%1.4e	", sigd[k].re*DT/DATA, sigd[k].im*DT/DATA);

		// --- ���ݑ��֊֐��̏�������
		if(k<=DATA/4-1){
			fprintf(fp, "%1.4e	", (k-DATA/4)*DT);
			fprintf(fp, "%1.4e	", corr[2][k]);
		}

		else{
			fprintf(fp, "%1.4e	", (k-DATA/4)*DT);
			fprintf(fp, "%1.4e	", corr[2][k]);
		}

		// --- ���ȑ��֊֐��̏�������
		if(k<=DATA/4-1){
			fprintf(fp, "%1.4e	", k*DT);
			fprintf(fp, "%1.4e	%1.4e" , corr[0][k], corr[1][k]);
		}

		fprintf(fp, "\n");
	}

	fclose(fp);
	printf("'%s' write\n", fp_name);
}

