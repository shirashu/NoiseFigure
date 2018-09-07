/* ----------------------------------------------
	���M���p�����[�^�A�G�������v�Z�v���O���� 
	2004/7/26 by Masahiro Nakayama 
	Ver.  2.0 (main_noise2_0.cpp)
---------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "new.h"

#include "parameter2_0.h"
#include "common_function2_0.h"
#include "y_parameter2_0.h"

/* ********************************************** */
/* ------------------ Main�֐� ------------------ */
int main() {

	/* ================================================================= */
	/*** --- �ϐ��錾 --- ***/
	struct SetPara para;	//-- �p�����[�^�ꗗ

	double *fr;				//-- ���g���C���f�b�N�X
	struct YPara *y;		//-- Y�p�����[�^
	struct Spectrum *si;	//-- �X�y�N�g�����x
	struct CircuitPara *circ;		/*--������H�p�����[�^ */

	/*** --- �v���g�^�C�v�錾 --- ***/
	void YParameter(struct YPara y[], double fr[],struct CircuitPara circ[], struct SetPara para);
	void CurrentSpectrum(struct Spectrum si[], double fr[], struct SetPara para);
	void NoiseFigure(struct Spectrum si[], struct YPara y[],
					 double fr[],struct CircuitPara circ[], struct SetPara para);
	/* ================================================================= */

	/* ----- �p�����[�^�̓ǂݍ��� */
	ReadParameter(&para);

	/* ----- �ǂݍ��񂾃p�����[�^�ɏ]���A�z��̃��������m�� */
	_set_new_handler(error);
	fr = new double[para.data];
	y = new struct YPara[para.data];
	si = new struct Spectrum[para.data];
	circ = new struct CircuitPara[para.data];

	/* ----- ���߂������g���̃C���f�b�N�X�𐶐� */
	FrequencyIndex(fr, para);


	/* ********* Y�p�����[�^�[�E������H�萔�v�Z ********* */
	/* --- Y�p�����[�^�[�̌v�Z */
	if(para.culc!=2){
		YParameter(y, fr,circ, para);
	}

	/* ********* �G���X�y�N�g�����x�̌v�Z ********* */
	/* --- �X�y�N�g���̌v�Z */
	if(para.culc!=1){
		CurrentSpectrum(si, fr, para);
	}
	
	/* ********* �G�������̌v�Z ********* */
	/* --- �m�C�Y�p�����[�^�E�G���w���̌v�Z */
	if(para.culc==0){
		NoiseFigure(si, y, fr,circ, para);
	}


	/* ----- �������̊J�� */
	delete [] fr;
	delete [] y;
	delete [] si;

	/* ----- �v���O�����̒�~���� */
	if(para.p_stop == 1){
		StopPG();
	}

    return 0;
}
