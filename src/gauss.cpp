#include <vector>
#include "gauss.h"

using namespace std;


gauss_iso::gauss_iso(vector<vector<vector<double>>>& keke, double& e,
	double& NU, int& gaussnum, vector<double>& xx, vector<double>& yy, vector<double>& zz)
	:ke(keke), E(e), nu(NU), gauss_num(gaussnum), x1(xx), y1(yy), z1(zz) {

	ga[0] = -0.577350269189626;
	ga[1] = 0.577350269189626;
	DD.resize(3, vector<double>(3, 0.0));
	double sum;
	//平面应力情况下的本构
	sum = E / (1 - nu * nu);
	DD[0][0] = sum; DD[0][1] = sum * nu;
	DD[1][0] = sum * nu; DD[1][1] = sum;
	DD[2][2] = sum * (1 - nu) / 2;
	/*平面应变的本构如下
	sum = E / ((1 + nu) * (1 - 2 * nu));
	DD[0][0] = sum * (1 - nu); DD[0][1] = sum * nu;
	DD[1][0] = sum * nu; DD[1][1] = sum * (1 - nu);
	DD[2][2] = sum * (1 - 2 * nu) / 2;
	*/
	ne.resize(4, 0.0);
	neds.resize(4, 0.0);
	nedt.resize(4, 0.0);
	nedn.resize(4, 0.0);
	J.resize(2, vector<double>(2, 0.0));
	dx.resize(2, 0.0);
	dy.resize(2, 0.0);
	dz.resize(2, 0.0);
	B.resize(8, vector<double>(3, 0.0));
}


void gauss_iso::compute(vector<vector<double>>& location) {
	int ii = 0;
	for (int i = 0; i < gauss_num; i++) {
		for (int j = 0; j < gauss_num; j++) {
			s = ga[i]; t = ga[j];
			Nene();
			Nedstn();
			Nxyz();
			location[ii][0] = x;
			location[ii][1] = y;

			xyzdsdtdn();
			jacobimatrix();
			Bs();
			//以下进行矩阵运算B'DB,使用其他的一些库最好
			for (int m = 0; m < 8; m++) {
				for (int n = 0; n < 8; n++) {
					ke[ii][m][n] = Jdet * BDB_line(B[m], B[n]);
				}
			}
			ii++;

		}
	}
}


//计算BDB结果某一个值2D
inline double gauss_iso::BDB_line(vector<double>& lineone, vector<double>& linetwo) {
	double sum = 0.0;
	double sum1 = 0.0;
	for (int i = 0; i < 3; i++) {
		sum1 = 0.0;
		for (int j = 0; j < 3; j++) {
			sum1 += DD[i][j] * linetwo[j];
		}
		sum += lineone[i] * sum1;
	}
	return sum;
}


//形状函数的计算2D

inline void gauss_iso::Nene() {
	ne[0] = (1.0 - s) * (1.0 - t) / 4;
	ne[1] = (1.0 - s) * (1.0 + t) / 4;
	ne[2] = (1.0 + s) * (1.0 + t) / 4;
	ne[3] = (1.0 + s) * (1.0 - t) / 4;
}


//计算坐标变换2D
inline void gauss_iso::Nxyz() {
	x = ne[0] * x1[0] + ne[1] * x1[1] + ne[2] * x1[2] + ne[3] * x1[3];
	y = ne[0] * y1[0] + ne[1] * y1[1] + ne[2] * y1[2] + ne[3] * y1[3];
}



//形状函数对s,t的导函数计算2D
inline void gauss_iso::Nedstn() {
	neds[0] = (-1) * (1.0 - t) / 4;
	neds[1] = (-1) * (1.0 + t) / 4;
	neds[2] = (1.0 + t) / 4;
	neds[3] = (1.0 - t) / 4;

	nedt[0] = (1.0 - s) * (-1) / 4;
	nedt[1] = (1.0 - s) / 4;
	nedt[2] = (1.0 + s) / 4;
	nedt[3] = (1.0 + s) * (-1) / 4;

}

//笛卡尔坐标系x,y对s,t的偏导数2D
inline void gauss_iso::xyzdsdtdn() {
	dx[0] = neds[0] * x1[0] + neds[1] * x1[1] + neds[2] * x1[2] + neds[3] * x1[3];
	dx[1] = nedt[0] * x1[0] + nedt[1] * x1[1] + nedt[2] * x1[2] + nedt[3] * x1[3];
	dy[0] = neds[0] * y1[0] + neds[1] * y1[1] + neds[2] * y1[2] + neds[3] * y1[3];
	dy[1] = nedt[0] * y1[0] + nedt[1] * y1[1] + nedt[2] * y1[2] + nedt[3] * y1[3];
}


//雅可比矩阵生成2D
inline void gauss_iso::jacobimatrix() {
	J[0][0] = dx[0]; J[0][1] = dy[0];
	J[1][0] = dx[1]; J[1][1] = dy[1];
	Jdet = dx[0] * dy[1] - dx[1] * dy[0];
}

//2D情况下B矩阵的生成，B：形函数偏导数矩阵   
inline void gauss_iso::Bs() {
	double a, b, c, d;  //B矩阵的部分系数计算
	a = dy[0];
	b = dy[1];
	c = dx[1];
	d = dx[0];
	//为了和matlab88行对应，这里变了，不是按照输入的来了，实际上是两者的正反手系不同
	//正确的如下
	//a = dy[1];
	//b = dy[0];
	//c = dx[0];
	//d = dx[1];
	for (int i = 0; i < 4; i++) {
		B[2 * i][0] = (a * neds[i] - b * nedt[i]) / Jdet;
		B[2 * i + 1][1] = (c * nedt[i] - d * neds[i]) / Jdet;
		B[2 * i][2] = (c * nedt[i] - d * neds[i]) / Jdet;
		B[2 * i + 1][2] = (a * neds[i] - b * nedt[i]) / Jdet;
	}
}




