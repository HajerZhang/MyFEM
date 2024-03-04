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
	vector<vector<double>> DD;   //弹性本构
	//s,t,n代表高斯积分点坐标
	double s, t, n;
	vector<double> ne;  //形函数
	vector<double> neds; //形函数对s的偏导数
	vector<double> nedt; //形函数对t的偏导数
	vector<double> nedn; //形函数对n的偏导数
	vector<vector<double>> J;   //雅可比阵
	double Jdet;        //雅可比值
	double x, y, z;     //坐标变换,即笛卡尔坐标
	vector<double> dx;  //笛卡尔坐标系x,y,z对s,t,n的偏导数
	vector<double> dy;
	vector<double> dz;
	//B矩阵：6×3×8矩阵，这里写成6×24矩阵（3维）
	//B矩阵：3×2×4矩阵，这里写成3×8矩阵（2维）
	vector<vector<double>> B;
	void Nene();          //形状函数计算
	void Nxyz();          //坐标变换结果
	void Nedstn();        //形状函数对s,t,n的导函数计算
	void xyzdsdtdn();     //笛卡尔坐标系x,y,z对s,t,n的偏导数
	void jacobimatrix();  //雅可比矩阵生成
	void Bs();            //B矩阵的生成，B：形函数偏导数矩阵
	double BDB_line(vector<double>& lineone, vector<double>& linetwo);    //用于计算BDB结果矩阵中任意一项的值
};