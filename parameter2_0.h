/* ----------------------------------------------
	���ʃp�����[�^�ꗗ
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (parameter2_0.h)
---------------------------------------------- */

#define PARAMETER "__SetPara20.txt" /* �p�����[�^�ꗗ�t�@�C���� */

void ReadParameter(struct SetPara *para);

/* **************************************************** */
/* ----- ���f���\���̂��߂̍\���� -----*/
struct ComplexNumbers{
	double re;
	double im;
};


/* **************************************************** */
/* ----- Y�p�����[�^�̂��߂̍\���� -----*/
struct YPara{
	struct ComplexNumbers _11;
	struct ComplexNumbers _12;
	struct ComplexNumbers _21;
	struct ComplexNumbers _22;
};


/* **************************************************** */
/* ----- �G���X�y�N�g�����x�\���̂��߂̍\���� -----*/
struct Spectrum{
	double g;//-- �Q�[�g�C�h���C���ɂ�����X�y�N�g��
	double d;//-- �Q�[�g�C�h���C���ɂ�����X�y�N�g��
	struct ComplexNumbers gd;//-- �Q�[�g�E�h���C�����݂ɂ�����X�y�N�g��
};



/* ******************************************* */
/* ----- �ǂݍ��݃p�����[�^�̂��߂̍\���� -----*/
struct SetPara{

	int culc;	// �v�Z�ΏہA0:���ׂāA1:Y�p���̂�2:�X�y�N�g���̂�

	double dt;	// ��t�F1�X�e�b�v�̎���

	// --- Y�p�����[�^�v�Z
	int jtp0;	// �X�e�b�v�d����������O�̎��ԃX�e�b�v
	int jtl0;	// �X�e�b�v�d������������̎��ԃX�e�b�v
	int jt;		// JTL0+JTP0 --- Simulation���s���S���ԃX�e�b�v
	int nump;	// �v�Z�Ώۂ̎��ԃX�e�b�v���i�d����������O�̃X�e�b�v���j
	int numl;	// �v�Z�Ώۂ̎��ԃX�e�b�v���i�d������������̃X�e�b�v���j
	int cur0;	//������H�@�@0:Ig=0�AIs=Id�Ɖ���C1:���L����Ȃ��AIg=Is-Id
	int yave;	// Y�p�����[�^�[�̕��������s���Ƃ��́@�d�ݕt���̃X�e�b�v��

	double fmin_y;	// �o�̓t�@�C���̍ŏ����g��
	double fmax_y;	// �o�̓t�@�C���̍ŏ����g���̍ő���g��
	double fstep_y;	// �o�̓t�@�C���̍ŏ����g���̎��g�����iy�p���v�Z�ɂ̂ݗL���j

	int lowpass_y;	//0:�s��Ȃ��C1:�s��
	double cutfreq_y;	// �J�b�g�I�t���g�� */
	int fwin_y;		// ���g���E�B���h�E�̃t�B���^�W�� */

	double dv0;	// �Q�[�g�d�ɂɗ^����X�e�b�v�d��
	double dv1;	// �h���C���d�ɂɗ^����X�e�b�v�d��
	int fin;	// ���g���C���f�b�N�X�̔z��

	char vgc[128];		//   �Q�[�g�X�e�b�v�d�����͂̓d���l�t�@�C���� */
	char vdc[128];		// �h���C���X�e�b�v�d�����͓d���l�t�@�C���� */
	char yparam[128];	// Y�p�����[�^�[�̃t�@�C���� */
	char output[128];	// �A�E�g�v�b�g�̃t�@�C���� */


	// --- �X�y�N�g���̌v�Z ---
	int jtn0;	// �S���ԃX�e�b�v
	int numn;	// �v�Z�Ώۂ̎��ԃX�e�b�v���i�d������������̃X�e�b�v���j
	int fnum;
	int epp;	// �����q1������̓d�q��

	int way;	//0:���l�ϕ��@�C1:FFT�@
		// --- FFT�@�̂Ƃ��ɕK�v�Ȑݒ�
		int data;	// FFT�v�Z�Ƃ��ĎQ�Ƃ���f�[�^��
		int jbit;	// DATA = 2^JBIT�@�Ɛݒ肷�邱��
		int twin;	// �X�y�N�g���E�B���h�E�A0:���p�����C1:���p����

		int saved;	// �h���C���G���ɑ΂��Ă̕������̏d�ݕt����(�s��Ȃ��Ƃ���1)
		int appr_sd;	// �h���C���G���X�y�N�g�����x[Sid = ���]�A0:���p�����C1:���p����
		double fmin_sd;	// ���Ƃ�����g���͈̔� - �ŏ��l
		double fmax_sd;	// ���Ƃ�����g���͈̔� - �ő�l
		double fmax_sd2;	// ���Ƃ�����g���͈̔� - �ő�l

		int saveg;	// �Q�[�g�G���ɑ΂��Ă̕������̏d�ݕt����(�s��Ȃ��Ƃ���1)
		int appr_sg;	// �Q�[�g�G���X�y�N�g�����x[Sig �� f^2]�A0:���p�����C1:���p����
		double fmin_sg;	// ���Ƃ�����g���͈̔� - �ŏ��l
		double fmax_sg;	// ���Ƃ�����g���͈̔� - �ő�l
		double fmax_sg2;	// ���Ƃ�����g���͈̔� - �ő�l

		int savegd;	// �h���C���E�Q�[�g���݂ɑ΂��Ă̕�����(�s��Ȃ��Ƃ���1)
		int appr_sgd;	// ���݃X�y�N�g��[Sigid �� f]�A0:���p�����C1:���p����
		double fmin_sgd;	// ���Ƃ�����g���͈̔� - �ŏ��l
		double fmax_sgd;	// ���Ƃ�����g���͈̔� - �ő�l
		double fmax_sgd2;	// ���Ƃ�����g���͈̔� - �ő�l


		// --- ���l�ϕ��@�̂Ƃ��ɕK�v�Ȑݒ� --- 
		int msmalls;	// ���֊֐��̎��ԃX�e�b�vm ---���ȑ��֊֐�
		int msmallc;	// ���֊֐��̎��ԃX�e�b�vm ---���ݑ��֊֐�
		int mlarge;		// M (MLARGE > MSMALL)


	// --- ���ʉ߃t�B���^�̐ݒ� ---
	int lowpass;	// 0:�s��Ȃ��C1:�s��
	double cutfreq_s;	// �J�b�g�I�t���g�� */
	int fwin_s;		// ���g���E�B���h�E�̃t�B���^�W�� */

	char vc[256];		// �G����͗p�d���f�[�^�̓��̓t�@�C���� */
	char spe[256];		// �G���d���X�y�N�g�����x�o�̓t�@�C���� */
	int cor_out;
	int current_out;
	int fouricomp_out;
	char cor[256];		// ���֊֐��o�̓t�@�C���� */
	char current[256];		// �d�����ρE�G���d������ */
	char fouricomp[256];		// �m�C�Y�w���o�̓t�@�C���� */

	/*---------------------------------*/

	// --- �G����� --- 
	int nway;	// �m�C�Y�w���v�Z�� 0:Gonzalez�C1:Cappy
	char noise[256];		// �m�C�Y�w���o�̓t�@�C���� */


	double pai;	// �~������
	double bk;	// �{���c�}���萔
	int ta;		// ���x

	double wg;	// �G���w���~�̌v�Z���s���Q�[�g��


	int p_stop;
};


