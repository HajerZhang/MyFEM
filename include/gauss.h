#pragma once
#include <vector>
using namespace std;

class gauss_iso
{
public:
    vector<vector<vector<double>>>& ke;
	
	double& E;
	double& nu;
	int& gauss_num;
	vector<double>& x1;
	vector<double>& y1;
	vector<double>& z1;
	gauss_iso(vector<vector<vector<double>>>& keke, double& e, double& NU, int& gaussnum, vector<double>& xx, vector<double>& yy, vector<double>& zz);

	void compute(vector<vector<double>>& location);
private:
    double ga[2];
	vector<vector<double>> DD;
	double s, t, n;
	vector<double> ne;
	vector<double> neds;
	vector<double> nedt;
	vector<double> nedn;
	vector<vector<double>> J;
	double Jdet;
	double x, y, z;
	vector<double> dx;
	vector<double> dy;
	vector<double> dz;
	vector<vector<double>> B;
	void Nene();
	void Nxyz();
	void Nedstn();
	void xyzdsdtdn();
	void jacobimatrix();
	void Bs();
	double BDB_line(vector<double>& lineone, vector<double>& linetwo);    
};