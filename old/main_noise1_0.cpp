/* ----------------------------------------------
	�m�C�Y�w���v�Z�v���O���� 
	2002/11/13 by Masahiro Nakayama 
	Ver.  1.0 (main_noise1_0.cpp)
	�p�����[�^�Ɋւ��ẮCparameter.h���Q�Ƃ̂��ƁD
---------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_0.h"
#include "current_spectrum1_0.h"
#include "cooley_tukey_FFT1_0.h"
#include "common_function1_0.h" 
#include "y_parameter1_0.h"
#include "noise_figure1_0.h"

/* ********************************************** */
/* ------------------ Main�֐� ------------------ */
int main() {
	char fn_vdc[]=VDC, fn_vgc[]=VGC;
	char fn_yparam[]=YPARAM, fn_out[]=OUTPUT;
	char fn_spe[]=SPE ,fn_noise[]=NIOSE;
	double cur[IN0][JT][IN0], cur1[IN0][JT][IN0], ss[IN0][IN0], rex[IN0][IN0];
		/*--cur:�d���l�Acur1:�ӂ���d���l*/
		/*--�d���l [0][*][*]:Gate,[1][*][*]:Drain,[2][*][*]:Source*/
		/*--       [*][*][0]:Gate Step�d�����,[*][*][1]:Drain Step�d�����*/
	struct ComplexNumbers y11[NUMS], y21[NUMS];//-- Y�p�����[�^
	double fr[NUMS];//���g���C���f�b�N�X
	double corr[IN0][NUMS]; //-- ���֊֐�
	double sig[NUMS], sid[NUMS]; //-- �Q�[�g�C�h���C���ɂ�����X�y�N�g��
	struct ComplexNumbers sigd[NUMS];//-- �Q�[�g�E�h���C�����݂ɂ�����X�y�N�g��


	/* ----- FILE����d���l��ǂݍ��݁A�ϐ��ɑ������ */
	ReadCurrent(cur, 0, &fn_vgc[0]);
	ReadCurrent(cur, 1, &fn_vdc[0]);
	/* ----- �d���̂ӂ���A����Ԃ̓d���l���v�Z */
	CurInit(cur, cur1, ss, rex, &fn_spe[0]);
	/* ----- ���߂������g���̃C���f�b�N�X�𐶐� */
	FrequencyIndex(fr);

	/* ----- Y�p�����[�^�[�E������H�萔�̌v�Z */
	YParameter(cur1, ss, rex, y11, y21, fr, &fn_yparam[0], &fn_out[0]);

	/* ----- �X�y�N�g������̃m�C�Y�v�Z ----- */
	CurrentSpectrum(cur1, corr, sig, sid, sigd, fr, &fn_spe[0]);
	NoiseFigure(sig, sid, sigd, y11, y21, fr, &fn_noise[0]);

    return 0;
}