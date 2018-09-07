/* ----------------------------------------------------
	�X�y�N�g���v�Z --- ���ۂ̐��l�v�Z�֐�
	�i���֊֐�����t�[���G�ϊ��ɂ��X�y�N�g�����v�Z�j
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (spe_cor2_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "new.h"

#include "parameter2_0.h" 
#include "common_function2_0.h" 
#include "current_spectrum2_0.h" 

/* ***************************+********************** */
/* ----- ���֊֐�����̃X�y�N�g���v�Z��Mian�֐� ----- */
void Spectrum_Cor(struct Current cur[], struct Correlation corr[],
				  struct Spectrum si[], double fr[], struct SetPara para){

	/* ================================================================= */
	/*** --- �v���g�^�C�v�錾 --- ***/
	void CalculateCorrelation_inte(struct Current cur[], struct Correlation corr[], 
								   struct SetPara para);
	void CalculateSpectrum_inte(struct Correlation corr[], struct Spectrum si[], 
								double fr[], struct SetPara para);
	void SpectrumWindow_inte(struct Spectrum si[], double fr[], struct SetPara para);
	void WriteSpectrum_inte(struct Correlation corr[], struct Spectrum si[], 
							double fr[], struct SetPara para);
	/* ================================================================= */

	/* --- �d���̂ӂ�����瑊�֊֐����v�Z */
	CalculateCorrelation_inte(cur, corr, para);
	/* --- ���֊֐�����X�y�N�g�����v�Z */
	CalculateSpectrum_inte(corr, si, fr, para);
	/* --- �X�y�N�g�����O�p�g�E�B���h�E�ɂ��␳ */
	SpectrumWindow_inte(si, fr, para);
	/* ----- FILE�փX�y�N�g���Ƒ��֊֐����������� */
	WriteSpectrum_inte(corr, si, fr, para);

}



/* ********************************************************************** 
/* -------------------------------------------------------------------------
  �ȉ��e��֐�
------------------------------------------------------------------------- */

/* ****************************************************** */
/* ----- �d���̂ӂ�����璼�ڂɑ��֊֐����v�Z���� ----- */
void CalculateCorrelation_inte(struct Current cur[], struct Correlation corr[], 
							   struct SetPara para){
	int dt, jt;
	double temp1, temp2, temp3, temp4, rs1, rs2, rs3, rs4;

	for(jt=0; jt<=para.msmalls; jt++){
		if(jt%100 == 0){printf("SelfCorr:%d / %d\n", jt, para.msmalls);}
		temp1 = 0;	temp2 = 0;		rs1 = 0;	rs2 = 0;
		for(dt=(para.jtn0-para.numn)+1; dt<=(para.jtn0-para.numn)+para.mlarge-para.msmalls; dt++){
			Sigma(cur[dt]._[0]*cur[dt+jt]._[0], &temp1, &rs1);/* ���ȑ��� */
			Sigma(cur[dt]._[1]*cur[dt+jt]._[1], &temp2, &rs2);/* ���ȑ��� */
		}
		corr[jt]._[0] = temp1/(double)(para.mlarge-para.msmalls)/para.epp;
		corr[jt]._[1] = temp2/(double)(para.mlarge-para.msmalls)/para.epp;
	}

	for(jt=0; jt<=para.msmallc; jt++){
		if(jt%100 == 0){printf("CrossCorr:%d / %d\n", jt, para.msmallc);}
		temp3=0;	temp4=0;	rs3 = 0;	rs4 = 0;
		for(dt=(para.jtn0-para.numn)+1; dt<=(para.jtn0-para.numn)+para.mlarge-para.msmallc; dt++){
			Sigma(cur[dt]._[0]*cur[dt+jt]._[1], &temp3, &rs3);/* ���ݑ��� */
			Sigma(cur[dt]._[1]*cur[dt+jt]._[0], &temp4, &rs4);/* ���ݑ��� */
		}
		corr[jt]._[2] = temp3/(double)(para.mlarge-para.msmallc)/para.epp;
		corr[jt]._[3] = temp4/(double)(para.mlarge-para.msmallc)/para.epp;
	}
}

/* ********************************************************** */
/* --- ���֊֐�����t�[���G�ϊ��ɂ��X�y�N�g�����v�Z���� --- */
void CalculateSpectrum_inte(struct Correlation corr[], struct Spectrum si[], double fr[], struct SetPara para){

	/*** ----- �ϐ��錾 */
	int jt, fmax, jf, tau;
	double *tg, *td;
	struct ComplexNumbers *tgd;
	double rs1, rs2;
	double cmax, cgd, cdg, *a, *b;

	/* ----- �ǂݍ��񂾃p�����[�^�ɏ]���A�z��̃��������m�� */
	_set_new_handler(error);
	tg = new double [para.numn];
	td = new double [para.numn];
	tgd = new struct ComplexNumbers [para.numn];
	a = new double [para.numn];
	b = new double [para.numn];


	/* --- ���ȑ��֊֐�����̃p���[�X�y�N�g���v�Z   */
	for(jf=0; ; jf++){
		if(jf%100 == 0){printf("PowerSpe:%d (%1.2e Hz)\n", jf, fr[jf]);}
		if(fr[jf]==32168){fmax=jf-1; break;} //�ő���g���Ōv�Z���X�g�b�v
		rs1 = 0;	rs2 = 0;
		tg[jf] = corr[0]._[0] + corr[para.msmalls]._[0]*cos(2*para.pai*fr[jf]*para.msmalls*para.dt);
		td[jf] = corr[0]._[1] + corr[para.msmalls]._[1]*cos(2*para.pai*fr[jf]*para.msmalls*para.dt);
		for(jt=1; jt<=para.msmalls-1; jt++){
			Sigma(2*corr[jt]._[0]*cos(2*para.pai*fr[jf]*jt*para.dt), &tg[jf], &rs1);
			Sigma(2*corr[jt]._[1]*cos(2*para.pai*fr[jf]*jt*para.dt), &td[jf], &rs2);
		}
		tg[jf] = tg[jf]*para.dt;
		td[jf] = td[jf]*para.dt;
	}

	for(jf=0; jf<=fmax; jf++){
		si[jf].g = tg[jf];
		si[jf].d = td[jf];
	}


	/* --- ���ݑ��֊֐�����̃N���X�X�y�N�g���v�Z   */
	cmax = corr[0]._[2];	tau = 0;
	for(jt=1; jt<=para.msmallc; jt++){
		if(corr[jt]._[2]>cmax){tau=jt; cmax=corr[jt]._[2];}
	}

	for(jt=0; jt<=para.msmallc-tau; jt++){
		cgd=corr[jt+tau]._[2];
		if(jt-tau>=0){cdg=corr[jt-tau]._[3];}
		else {cdg=corr[-jt+tau]._[2];}
		a[jt] = (cgd + cdg)/2;
		b[jt] = (cgd - cdg)/2;
	}

	for(jf=0; jf<=fmax; jf++){
		if(jf%100 == 0){printf("PowerSpe:%d (%1.2e Hz)\n", jf, fr[jf]);}
		rs1 = 0;	rs2 = 0;
		tgd[jf].re = a[0] + a[para.msmallc-tau]*cos(2*para.pai*fr[jf]*(para.msmallc-tau)*para.dt);
		tgd[jf].im =        b[para.msmallc-tau]*sin(2*para.pai*fr[jf]*(para.msmallc-tau)*para.dt);
		for(jt=1; jt<=para.msmallc-tau-1; jt++){
			Sigma(2*a[jt]*cos(2*para.pai*fr[jf]*jt*para.dt), &tgd[jf].re, &rs1);
			Sigma(2*b[jt]*sin(2*para.pai*fr[jf]*jt*para.dt), &tgd[jf].im, &rs2);
		}
		tgd[jf].re = tgd[jf].re*para.dt;
		tgd[jf].im = tgd[jf].im*para.dt;
	}

	for(jf=0; jf<=fmax; jf++){
		a[jf] = tgd[jf].re;
		b[jf] = tgd[jf].im;
	}

	for(jf=0; jf<=fmax-1; jf++){
		si[jf].gd.re = (a[jf]*cos(2*para.pai*fr[jf]*tau*para.dt) -
					    b[jf]*sin(2*para.pai*fr[jf]*tau*para.dt));
		si[jf].gd.im = (a[jf]*sin(2*para.pai*fr[jf]*tau*para.dt) +
					    b[jf]*cos(2*para.pai*fr[jf]*tau*para.dt));
	}

	/* ----- �������̊J�� */
	delete [] tg;
	delete [] td;
	delete [] tgd;
	delete [] a;
	delete [] b;


}


/* ************************************************************ */
/* ----- �O�p�`�E�B���h�ɂ��X�y�N�g���̎��g���������s�� ----- */
void SpectrumWindow_inte(struct Spectrum si[], double fr[], struct SetPara para){
	int k, jf, fmax, fmin, *l1, *l2;
	double *t1, *t2, *t3, *t4;
	double rs1, rs2, rs3, rs4;

	/* ----- �ǂݍ��񂾃p�����[�^�ɏ]���A�z��̃��������m�� */
	_set_new_handler(error);
	l1 = new int [para.numn];
	l2 = new int [para.numn];
	t1 = new double [para.numn];
	t2 = new double [para.numn];
	t3 = new double [para.numn];
	t4 = new double [para.numn];


	for(jf=0; ; jf++){
		if(fr[jf]==0){fmin=jf;} /*���g��0Hz�̃C���f�b�N�X���L�^ */
		else if(fr[jf]==32168){fmax=jf-1; break;} /*�ő���g���Ōv�Z���X�g�b�v */
	}

	/* --- Sig�CSigd�����߂� */
	for(jf=0; jf<=fmax; jf++){
		t1[jf] = 0; t2[jf] = 0; t3[jf] = 0; t4[jf] = 0;
		rs1 = 0;	rs2 = 0;	rs3 = 0;	rs4 = 0;

		l1[jf] = para.saved; l2[jf] = para.saved;
		if(jf-l1[jf]+1<0)    {do{l1[jf]--;}while(jf-l1[jf]+1<0);}
		if(jf+l2[jf]-1>fmax)    {do{l2[jf]--;}while(jf+l2[jf]-1>fmax);}

		for(k=-l2[jf]+1; k<=l1[jf]-1; k++){
			Sigma((para.saveg-abs(k))*si[jf-k].g, &t1[jf], &rs1);
			Sigma((para.savegd-abs(k))*si[jf-k].gd.re, &t3[jf], &rs3);
			Sigma((para.savegd-abs(k))*si[jf-k].gd.im, &t4[jf], &rs4);
		}
	}

	for(jf=0; jf<=fmax; jf++){
		si[jf].g = t1[jf]/para.saveg/(l1[jf]+l2[jf])*2;
		si[jf].gd.re = t3[jf]/para.savegd/(l1[jf]+l2[jf])*2;
		si[jf].gd.im = t4[jf]/para.savegd/(l1[jf]+l2[jf])*2;
	}

	/* --- Sid�����߂� */
	for(jf=1; jf<=fmax; jf++){
		t2[jf] = 0;	rs2 = 0;
		l1[jf] = para.saved; l2[jf] = para.saved;
		if(jf-l1[jf]+1<1)    {do{l1[jf]--;}while(jf-l1[jf]+1<1);}
		if(jf+l2[jf]-1>fmax)    {do{l2[jf]--;}while(jf+l2[jf]-1>fmax);}

		for(k=-l2[jf]+1; k<=l1[jf]-1; k++){
			Sigma((para.saved-abs(k))*si[jf-k].d, &t2[jf], &rs2);
		}
	}

	for(jf=1; jf<=fmax; jf++){
		si[jf].d = t2[jf]/para.saved/(l1[jf]+l2[jf])*2;
	}


	/* ---  ���g��0Hz�ŃX�y�N�g��sig�Csigd��0�Ƃ��ĕ␳ */
	for(jf=fmax; jf>=0; jf--){
		si[jf].g = si[jf].g - si[fmin].g;
		si[jf].gd.im = si[jf].gd.im - si[fmin].gd.im;
		si[jf].gd.re = si[jf].gd.re - si[fmin].gd.re;

	}

	/* ----- �������̊J�� */
	delete [] l1;
	delete [] l2;
	delete [] t1;
	delete [] t2;
	delete [] t3;
	delete [] t4;

}

/* ************************************************ */
/* ----- FILE�փX�y�N�g���Ƒ��֊֐����������� ----- */
void WriteSpectrum_inte(struct Correlation corr[], struct Spectrum si[], 
						double fr[], struct SetPara para){
	int k;
    FILE *fp;

	/* --- ���ϓd���l�̏������� */
	if(para.current_out ==1){
		fp = OpenFile(para.current, 'a');
		fprintf(fp, "��I_G	��I_D\n");
		fprintf(fp, "%lf	%lf\n",sqrt(corr[0]._[0]), sqrt(corr[0]._[1]));
	fclose(fp);
	printf("'%s' write\n", para.current);
	}

	/* --- ���֊֐��̏������� */
	if(para.cor_out ==1){
		fp = OpenFile(para.cor, 'w');
		fprintf(fp, "t[sec]	Cor_d(t)	Cor_g(t)	CrossCor(t)\n");
		for(k=0; k<=para.data/4-1; k++){
			if(fr[k]==32168){break;} /*�ő���g���Ōv�Z���X�g�b�v */
			fprintf(fp, "%1.4e	", 2*k*para.dt);
			fprintf(fp, "%1.4e	%1.4e	%1.4e", corr[k]._[1], corr[k]._[0], corr[k]._[2]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("'%s' write\n", para.cor);
	}


	/* --- �X�y�N�g���̏������� */
	fp = OpenFile(para.spe, 'w');
	fprintf(fp, "f[Hz]	Si_D(f)	Si_G(f)	Si_G-D(f).re	Si_G-D(f).im\n");
	for(k=0; k<=para.data/4-1; k++){
		if(fr[k]==32168){break;} /*�ő���g���Ōv�Z���X�g�b�v */
		fprintf(fp, "%1.4e	", k/para.dt/para.data);
		fprintf(fp, "%1.4e	" , si[k].d);
		fprintf(fp, "%1.4e	" , si[k].g);
		fprintf(fp, "%1.4e	" , si[k].gd.re);
		fprintf(fp, "%1.4e\n"   , si[k].gd.im);
	}
	fclose(fp);
	printf("'%s' write\n", para.spe);

}





