#include <iostream>
#include "MPC2X2.h"
#include <Eigen/Dense>
#include "math.h"






float heaviside(float t, float dt) {
	bool bin;
	bin = t > dt;
	return (float)bin;
}

float y(float t, float deadTime, float coef, float toa) {
	float result;
	result = coef * (1 - exp(-(t - deadTime) / toa)) * heaviside(t, deadTime);
	return result;
}


int main() {
	int predH = 100;
	int manH = 2;
	float dt = 1;
	Eigen::Vector2f setPt;
	setPt(0) = 1;// simple setpoint change
	setPt(1) = 0;
	Eigen::Vector2f weights;
	weights(0) = .1;// weights
	weights(1) = .1;// weights
	// Constructor
	MPC2X2 MPC = MPC2X2::MPC2X2(predH, manH, dt, setPt, weights);
	float t, y11, y12, y21, y22;



	MPC.resizeAll();
	// Set the MV to CV vectors
	for (int i = 0; i < predH; i++) {
		
		t = (float)i * dt;

		y11 = y(t, 1, .5, 5);
		y12 = y(t, .2, .3, 3);
		y21 = y(t, 1, .2, 5);
		y22 = y(t, 0, .1, 1.5);

		MPC.y11Set(y11, i);
		MPC.y12Set(y12, i);
		MPC.y21Set(y21, i);
		MPC.y22Set(y22, i);
	}
	
	// make / resize matrices
	MPC.makeMatrices();

	Eigen::Vector2f y_test;
	y_test(0) = 1;
	y_test(1) = 3;
	MPC.y_hSet(y_test);

	Eigen::VectorXf y_h = MPC.get_y_h();
	
	Eigen::Vector2f du;
	for (int i = 0; i < 700; i++) {
		t = dt * (float)i;
		du = MPC.compute(y_h);
		y_h = MPC.yhAdd(du);
		std::cout << "dMV" << std::endl;
		std::cout << du << std::endl;
		std::cout << "outputs" << std::endl;
		std::cout << y_h(0) << std::endl;
		std::cout << y_h(predH) << std::endl;
		std::cout << "setpoints" << std::endl;
		std::cout << setPt(0) << std::endl;
		std::cout << setPt(1) << std::endl;
	}











	//Eigen::VectorXf setPoint = MPC.getSetPoint();
	//std::cout << setPoint << std::endl;
	//Eigen::MatrixXf W = MPC.getW();
	//std::cout << W << std::endl;


	// GENERAL GETTERS / SETTERS below


	//std::cout << "=============X=============" << std::endl;
	//Eigen::MatrixXf X = MPC.getX();
	//std::cout << X << std::endl;
	//std::cout << X.size() << std::endl;

	//std::cout << "=============H=============" << std::endl;
	//Eigen::MatrixXf H = MPC.getH_shift();
	//std::cout << H << std::endl;
	//std::cout << H.size() << std::endl;

	//std::cout << "=============y11=============" << std::endl;
	//eigen::vectorxf _y11 = mpc.get_y11();
	//std::cout << _y11 << std::endl;
	//std::cout << _y11.size() << std::endl;

	//eigen::vectorxf der;
	//der.setones(_y11.size(), 1);
	//der = der * 3;
	//std::cout <<der<< std::endl;
	//std::cout << der.transpose()*_y11  << std::endl;


	//std::cout << "=============Hy11=============" << std::endl;
	//std::cout << H*_y11 << std::endl;
	//


	//std::cout << "=============INV XTX==========" << std::endl;
	//Eigen::MatrixXf inv_XTX = MPC.getInv_XTX();
	//std::cout << X << std::endl;
}