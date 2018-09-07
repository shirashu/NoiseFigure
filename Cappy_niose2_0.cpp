/* ----------------------------------------------------
	ノイズ指数の計算
	2003/9/30 by Masahiro Nakayama 
	Ver. 2.0 (Cappy_niose2_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter2_0.h"
#include "common_function2_0.h" 
#include "y_parameter2_0.h"

/* ******************************************** */
/* --- ノイズ指数計算 Cappy提案の式による--- */
void MinimumNF_Cappy(struct Spectrum si[], 
				struct YPara y[], 
				double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[], 
				struct ComplexNumbers zc[],double gass[],struct CircuitPara circ[], struct SetPara para){
	int jf;
	double y11_abs, y21_abs, y12_abs, y22_abs;
	double rn, x, a, b;
	double gmat, bmat, d, e, f, g, h, i;
	double svg, svd;
	double re_yopt,im_yopt,re_ya,im_ya,re_yb,im_yb,re_yc,re_yd,re_ye; 			/*Gass計算用パラメータ*/

	struct ComplexNumbers svgd;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} /*最大周波数で計算をストップ */

		gmat = 	(y[jf]._11.re*y[jf]._22.re - y[jf]._11.im*y[jf]._22.im) - (y[jf]._12.re*y[jf]._21.re - y[jf]._12.im*y[jf]._21.im);
		bmat = 	(y[jf]._11.re*y[jf]._22.im + y[jf]._11.im*y[jf]._22.re) - (y[jf]._12.re*y[jf]._21.im + y[jf]._12.im*y[jf]._21.re);

		y11_abs = (y[jf]._11.re*y[jf]._11.re + y[jf]._11.im*y[jf]._11.im);
		y12_abs = (y[jf]._12.re*y[jf]._12.re + y[jf]._12.im*y[jf]._12.im);
		y21_abs = (y[jf]._21.re*y[jf]._21.re + y[jf]._21.im*y[jf]._21.im);
		y22_abs = (y[jf]._22.re*y[jf]._22.re + y[jf]._22.im*y[jf]._22.im);

		p[jf] = si[jf].d / (4*para.bk*para.ta*sqrt(y21_abs));
		r[jf] = si[jf].g*sqrt(y21_abs) / (4*para.bk*para.ta*y11_abs);
		c[jf] = sqrt(si[jf].gd.im*si[jf].gd.im+si[jf].gd.re*si[jf].gd.re) / sqrt(si[jf].g*si[jf].d);


		a = gmat*gmat+bmat*bmat;
		b = y22_abs * si[jf].g;
		x = y12_abs * si[jf].d;
		d = y[jf]._22.re*y[jf]._12.re + y[jf]._22.im*y[jf]._12.im;
		e =-y[jf]._22.re*y[jf]._12.im + y[jf]._22.im*y[jf]._12.re;
		f = 2* (d*si[jf].gd.re + e*si[jf].gd.im);
		svg = 1/a*(b+x-f);

		a = gmat*gmat+bmat*bmat;
		b = y21_abs * si[jf].g;
		x = y11_abs * si[jf].d;
		d = y[jf]._21.re*y[jf]._11.re + y[jf]._21.im*y[jf]._11.im;
		e =-y[jf]._21.re*y[jf]._11.im + y[jf]._21.im*y[jf]._11.re;
		f = 2* (d*si[jf].gd.re + e*si[jf].gd.im);
		svd = 1/(a)*(b+x-f);

		a = gmat*gmat+bmat*bmat;
		b = ( y[jf]._22.re*y[jf]._21.re + y[jf]._22.im*y[jf]._21.im) * si[jf].g;
		x = ( y[jf]._12.re*y[jf]._11.re + y[jf]._12.im*y[jf]._11.im) * si[jf].d;
		d = y[jf]._22.re*y[jf]._11.re + y[jf]._22.im*y[jf]._11.im;
		e = y[jf]._22.re*y[jf]._11.im - y[jf]._22.im*y[jf]._11.re;
		f = d*si[jf].gd.re - e*si[jf].gd.im;
		g = y[jf]._12.re*y[jf]._21.re + y[jf]._12.im*y[jf]._21.im;
		h = y[jf]._12.re*y[jf]._21.im - y[jf]._12.im*y[jf]._21.re;
		i = g*si[jf].gd.re - h*si[jf].gd.im;
		svgd.re = 1/(a)*(-b - x + f + i);

		a = gmat*gmat+bmat*bmat;
		b = ( y[jf]._22.re*y[jf]._21.im - y[jf]._22.im*y[jf]._21.re) * si[jf].g;
		x = ( y[jf]._12.re*y[jf]._11.im - y[jf]._12.im*y[jf]._11.re) * si[jf].d;
		d = y[jf]._22.re*y[jf]._11.re + y[jf]._22.im*y[jf]._11.im;
		e = y[jf]._22.re*y[jf]._11.im - y[jf]._22.im*y[jf]._11.re;
		f = d*si[jf].gd.im + e*si[jf].gd.re;
		g = y[jf]._12.re*y[jf]._21.re + y[jf]._12.im*y[jf]._21.im;
		h = y[jf]._12.re*y[jf]._21.im - y[jf]._12.im*y[jf]._21.re;
		i = g*si[jf].gd.im + h*si[jf].gd.re;
		svgd.im = 1/(a)*(-b - x + f + i);

		a = y[jf]._22.re + y[jf]._12.re*svgd.re / svd  - y[jf]._12.im*svgd.im / svd;
		b = y[jf]._22.im + y[jf]._12.re*svgd.im / svd  + y[jf]._12.im*svgd.re / svd;
		x = (gmat)/(gmat*gmat+bmat*bmat);
		d =-(bmat)/(gmat*gmat+bmat*bmat);
		zc[jf].re = (a*x - b*d);
		zc[jf].im = (a*d + b*x);

		gn[jf] = svd / (4*para.bk*para.ta)*(gmat*gmat+bmat*bmat)/y21_abs;

		a = svg / (4*para.bk*para.ta);
		b = svd / (4*para.bk*para.ta);
		x = y12_abs / y21_abs;
		d = (svgd.re*svgd.re + svgd.im*svgd.im)/svd/svd;
		rn = a-b*x*d;

		fmin[jf] = 1+2*gn[jf]*(zc[jf].re+sqrt(zc[jf].re*zc[jf].re+rn/gn[jf]));
		fmin[jf] = 10*log10(fmin[jf]);

		zopt[jf].re = sqrt(zc[jf].re*zc[jf].re+rn/gn[jf]);
		zopt[jf].im = -zc[jf].im;
//Gassの計算(Yopt=Gopt+jBopt:Monte Carlo analysis of the dynamic behavior of III-V MOSFETs for low-noise RF applicatons)
		y11_abs = sqrt(y[jf]._11.re*y[jf]._11.re + y[jf]._11.im*y[jf]._11.im);
		y12_abs = sqrt(y[jf]._12.re*y[jf]._12.re + y[jf]._12.im*y[jf]._12.im);
		y21_abs = sqrt(y[jf]._21.re*y[jf]._21.re + y[jf]._21.im*y[jf]._21.im);
		y22_abs = sqrt(y[jf]._22.re*y[jf]._22.re + y[jf]._22.im*y[jf]._22.im);

		re_yopt = 2*para.pai*fr[jf]*(circ[jf].cgs+circ[jf].cgd)*((c[jf]*sqrt(r[jf]/p[jf]))-1);
		im_yopt = 2*para.pai*fr[jf]*(circ[jf].cgs+circ[jf].cgd)*sqrt(r[jf]/p[jf])*sqrt(1-c[jf]*c[jf]);

		re_ya = y[jf]._21.re*y[jf]._12.re - y[jf]._21.im*y[jf]._12.im;		//Re{Y12Y21}
		im_ya = y[jf]._21.re*y[jf]._12.im + y[jf]._12.re*y[jf]._21.im;		//Im{Y12Y21}
		re_yb = re_yopt + y[jf]._11.re;							//Re{Yopt+Y11}
		im_yb = im_yopt + y[jf]._11.im;							//Im{Yopt+Y11}
		re_yc = (re_ya*re_yb + im_ya*im_yb) / (re_yb*re_yb +im_yb*im_yb);	// Y21Y12 / (Yopt+Y11)

	    re_yd = y[jf]._22.re - re_yc;			//Re[Y22-{Y21Y12/(Yopt+Y11)]

		re_ye= pow( fabs( y21_abs / sqrt(re_yb*re_yb +im_yb*im_yb) ),2.0 );	//|Y21/(Yopt+Y11)|**2

		gass[jf]= fabs((re_yopt / re_yd) * re_ye);
		gass[jf] = 10*log10(gass[jf]);
	}
}
