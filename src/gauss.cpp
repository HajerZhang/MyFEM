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
	// plane stress
	sum = E / (1 - nu * nu);
	DD[0][0] = sum; DD[0][1] = sum * nu;
	DD[1][0] = sum * nu; DD[1][1] = sum;
	DD[2][2] = sum * (1 - nu) / 2;
	/* plane strain
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
			for (int m = 0; m < 8; m++) {
				for (int n = 0; n < 8; n++) {
					ke[ii][m][n] = Jdet * BDB_line(B[m], B[n]);
				}
			}
			ii++;

		}
	}
}

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

inline void gauss_iso::Nene() {
	ne[0] = (1.0 - s) * (1.0 - t) / 4;
	ne[1] = (1.0 - s) * (1.0 + t) / 4;
	ne[2] = (1.0 + s) * (1.0 + t) / 4;
	ne[3] = (1.0 + s) * (1.0 - t) / 4;
}

inline void gauss_iso::Nxyz() {
	x = ne[0] * x1[0] + ne[1] * x1[1] + ne[2] * x1[2] + ne[3] * x1[3];
	y = ne[0] * y1[0] + ne[1] * y1[1] + ne[2] * y1[2] + ne[3] * y1[3];
}

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

inline void gauss_iso::xyzdsdtdn() {
	dx[0] = neds[0] * x1[0] + neds[1] * x1[1] + neds[2] * x1[2] + neds[3] * x1[3];
	dx[1] = nedt[0] * x1[0] + nedt[1] * x1[1] + nedt[2] * x1[2] + nedt[3] * x1[3];
	dy[0] = neds[0] * y1[0] + neds[1] * y1[1] + neds[2] * y1[2] + neds[3] * y1[3];
	dy[1] = nedt[0] * y1[0] + nedt[1] * y1[1] + nedt[2] * y1[2] + nedt[3] * y1[3];
}

inline void gauss_iso::jacobimatrix() {
	J[0][0] = dx[0]; J[0][1] = dy[0];
	J[1][0] = dx[1]; J[1][1] = dy[1];
	Jdet = dx[0] * dy[1] - dx[1] * dy[0];
}

inline void gauss_iso::Bs() {
	double a, b, c, d;
	a = dy[0];
	b = dy[1];
	c = dx[1];
	d = dx[0];
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




