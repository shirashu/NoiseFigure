/* ----------------------------------------------------
	�X�y�N�g���v�Z
	2015/ T.takahashi 
	Ver. 3.0 (current_spectrum2_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "new.h"

#include "parameter2_0.h" 
#include "common_function2_0.h" 
#include "current_spectrum2_0.h" 

/* ***************************+************** */
/* ----- �d���̃X�y�N�g���v�Z��Mian�֐� ----- */
void CurrentSpectrum(struct Spectrum si[], double fr[], struct SetPara para){


	/* ================================================================= */
	/* ================================================================= */
	/*** --- �ϐ��錾 --- ***/
	struct Current *cur;		//-- �m�C�Y������͗p�f�[�^
	struct Correlation *corr;	//-- ���֊֐�

	/*** --- �v���g�^�C�v�錾 --- ***/
	void ReadCurrent_Spe(struct Current cur[], struct SetPara para);
	void CurInit_Spe(struct Current cur[], struct SetPara para);
	void LowPassFilter_S(struct Current cur[], struct SetPara para);
	void Spectrum_Cor(struct Current cur[], struct Correlation corr[],
					  struct Spectrum si[], double fr[], struct SetPara para);
	void Spectrum_FFT(struct Current cur[], struct Correlation corr[],
					  struct Spectrum si[], double fr[], struct SetPara para);

	/* ----- �ǂݍ��񂾃p�����[�^�ɏ]���A�z��̃��������m�� */
	_set_new_handler(error);
	cur = new struct Current[para.data];
	corr = new struct Correlation[para.data];

	/* ================================================================= */
	/* ================================================================= */


	/* --- FILE����d���l��ǂݍ��� */
	ReadCurrent_Spe(cur, para);
	/* --- �d���̂ӂ���A����Ԃ̓d���l���v�Z */
	CurInit_Spe(cur, para);

	/* --- ���ʉ߃t�B���^ */
	if(para.lowpass==1){
		LowPassFilter_S(cur, para);
	}

	/* --- FFT�ɂ��X�y�N�g�����v�Z */
	if(para.way==1){
		Spectrum_FFT(cur, corr, si, fr, para);
	}
	/* --- ���֊֐�����X�y�N�g�����v�Z */
	else if(para.way==0){
		Spectrum_Cor(cur, corr, si, fr, para);
	}

	/* ----- �������̊J�� */
	delete [] cur;
	delete [] corr;


}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  �ȉ��e��֐�
------------------------------------------------------------------------- */

/* ****************************************************** */
/* ----- ���ʉ߃t�B���^ ----- */
void LowPassFilter_S(struct Current cur[], struct SetPara para){
	int jf;
	double *x;

	/* ----- �ǂݍ��񂾃p�����[�^�ɏ]���A�z��̃��������m�� */
	_set_new_handler(error);
	x = new double [(para.jtn0 * 4)];


	/* ----- ���[�p�X�t�B���^�[�������邽�߂̏��� */
	for(jf=0; jf<=para.jtl0-1; jf++){
		x[jf] = cur[jf]._[0];
		x[para.jtn0*1 + jf] = cur[jf]._[1];
		x[para.jtn0*2 + jf] = cur[jf]._[2];
		x[para.jtn0*3 + jf] = 0.0;
	}

	/* ----- ���[�p�X�t�B���^�[���� */
	LowPassFilter(x, para.jtn0, para.cutfreq_s, para.fwin_s, para);

	/* ----- ���[�p�X�t�B���^�[����������̏��� */
	for(jf=0; jf<=para.jtl0-1; jf++){
		cur[jf]._[0]= x[jf];
		cur[jf]._[1]= x[para.jtn0*1 + jf];
		cur[jf]._[2]= x[para.jtn0*2 + jf];
	}

	/* ----- �������̊J�� */
	delete [] x;




}
