/* ----------------------------------------------
	�m�C�Y�w���v�Z�v���O���� 
	2003/3/31 by Masahiro Nakayama 
	Ver.  1.1 (main_noise1_1.cpp)
---------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_1.h"
#include "current_spectrum1_1.h"
#include "cooley_tukey_FFT1_1.h"
#include "common_function1_1.h" 
#include "y_parameter1_1.h"
#include "noise_figure1_1.h"

/* ********************************************** */
/* ------------------ Main�֐� ------------------ */
int main() {
	char fn_vdc[]=VDC, fn_vgc[]=VGC, fn_vc[]=VC;
	char fn_yparam[]=YPARAM, fn_out[]=OUTPUT, fn_spe[]=SPE ,fn_noise[]=NIOSE;
	double fr[FIN];//���g���C���f�b�N�X
	double cur[IN0][JT][IN0], cur1[IN0][JT][IN0], ss[IN0][IN0], rex[IN0][IN0];
		/*--cur:�d���l�Acur1:�ӂ���d���l*/
		/*--�d���l [0][*][*]:Gate,[1][*][*]:Drain,[2][*][*]:Source*/
		/*--       [*][*][0]:Gate Step�d�����,[*][*][1]:Drain Step�d�����*/
	double cur2[IN0][JTN0]; //-- �m�C�Y������͗p�f�[�^
	struct ComplexNumbers y11[FIN], y21[FIN];//-- Y�p�����[�^
	struct ComplexNumbers y12[FIN], y22[FIN];//-- Y�p�����[�^
	double corr[IN0+1][JTN0]; //-- ���֊֐�
	double sig[NUMN], sid[NUMN]; //-- �Q�[�g�C�h���C���ɂ�����X�y�N�g��
	struct ComplexNumbers sigd[NUMN];//-- �Q�[�g�E�h���C�����݂ɂ�����X�y�N�g��


	// ----- ���߂������g���̃C���f�b�N�X�𐶐�
	FrequencyIndex(fr);

	// ********* Y�p�����[�^�[�E������H�萔�v�Z ********* 
	// --- FILE����d���l��ǂݍ��� 
	ReadCurrent(cur, 0, &fn_vgc[0]);
	ReadCurrent(cur, 1, &fn_vdc[0]);
	// --- �d���̂ӂ���A����Ԃ̓d���l���v�Z 
	CurInit(cur, cur1, ss, rex);
	// --- Y�p�����[�^�[�̌v�Z
	YParameter(cur1, ss, rex, y11, y12, y21, y22, fr, &fn_yparam[0], &fn_out[0]);

	// ********* �m�C�Y������͌v�Z ********* 
	// --- FILE����d���l��ǂݍ��� 
	ReadCurrentVC(cur2, &fn_vc[0]);
	// --- �d���̂ӂ���A����Ԃ̓d���l���v�Z 
	CurInitVC(cur2, &fn_spe[0]);
	// --- �X�y�N�g���̌v�Z
	CurrentSpectrum(cur2, corr, sig, sid, sigd, fr, &fn_spe[0]);
	// --- �m�C�Y�p�����[�^�E�w���̌v�Z
	NoiseFigure(sig, sid, sigd, y11, y12, y21, y22, fr, &fn_noise[0]);

    return 0;
}