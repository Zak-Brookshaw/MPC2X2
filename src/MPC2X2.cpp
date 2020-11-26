#include <iostream>
#include "MPC2X2.h"
#include "Eigen/Dense"



MPC2X2::MPC2X2(int predH, int manH, float dt, Eigen::Vector2f setPt, Eigen::Vector2f weights) :
	predH(predH), manH(manH), dt(dt), setPt(setPt), weights(weights) {

}	



Eigen::MatrixXf MPC2X2::makeInsertMatrix(Eigen::VectorXf y) {
	Eigen::MatrixXf insert(predH, manH);
	
	for (int k = 0; k < manH; k++) {
		insert.block(0, k, predH, 1) = y;
		y = H_shift * y;
	}
	
	return insert;
}

void MPC2X2::resizeAll() {

	// MV to CV vectors resized
	y11.resize(predH);
	y12.resize(predH);
	y21.resize(predH);
	y22.resize(predH);

	// Shift matrix
	H_shift.resize(predH, predH);

	// y_shift

	y_shift.resize(predH, predH);

	// Weights matrix
	W.resize(mvNum * manH, mvNum * manH);
	W.setZero();

	// X matrix
	X.resize(predH * cvNum, manH * mvNum);

	// X modified
	X_modified.resize(predH * cvNum, mvNum);

	// setPoint
	setPoint.resize(predH * cvNum);

	// y_h
	y_h.resize(predH * cvNum);
	y_h.setZero();

	// duFull
	duFull.resize(mvNum * manH);
}



void MPC2X2::makeMatrices() {

	

	// Shift matrix
	
	H_shift.setZero();
	H_shift.block(1, 0, predH - 1, predH - 1).setIdentity();
	H_shift(predH - 1, predH - 1) = 1;

	// y_shift
	y_shift.setZero();
	y_shift.block(0, 1, predH - 1, predH - 1).setIdentity();
	y_shift(predH - 1, predH - 1) = 1;
	
	// Fill the weights matrix
	int j;
	for (int i = 0; i < manH; i++) {
		j = i + manH;
		W(i, i) = weights(0);
		W(j, j) = weights(1);
	}

	// temporary variable for X matrix
	Eigen::MatrixXf _insert(predH, manH);

	// fill the X matrix
	//top left quadrant filled
	_insert = makeInsertMatrix(y11);
	X.block(0, 0, predH, manH) = _insert;

	// top right quadrant filled
	_insert = makeInsertMatrix(y12);
	X.block(0, manH, predH, manH) = _insert;

	// bottom left quadrant filled
	_insert = makeInsertMatrix(y21);
	X.block(predH, 0, predH, manH) = _insert;
	// bottom right quadrant filled
	_insert = makeInsertMatrix(y22);
	X.block(predH, manH, predH, manH) = _insert;

	X_modified.resize(predH*cvNum, mvNum);
	for (int col = 0; col < mvNum; col++) {
		X_modified.block(0, col, predH * cvNum, 1) = X.block(0, col*mvNum, predH * cvNum, 1);
	}

	// Transpose X matrix
	XT = X.transpose();
	// special inverse matrix
	inv_XTXW = (XT * X + W).inverse();

	// Fill the setpoint vector
	Eigen::VectorXf ones = createOnes(predH);
	
	// create setPoint vector same shape as y_h, so they can be added ... 
	for (int i = 0; i < cvNum; i++) {
		setPoint.block(i * predH, 0, predH, 1) = ones * setPt(i);
	}

}



Eigen::Vector2f MPC2X2::compute(Eigen::VectorXf y_ti) {

	Eigen::VectorXf ones = createOnes(predH);
	// get error
	Eigen::VectorXf error(predH*cvNum); // wont work
	// transform e_ti into proper form for matrix mult
	//Eigen::Vector2f e_ti  // To Do
	for (int i = 0; i < cvNum; i++) {
		error.block(i * predH, 0, predH, 1) = setPoint.block(i * predH, 0, predH, 1) - y_ti.block(i * predH, 0, predH, 1);
	}
	duFull = inv_XTXW * XT * error;
	for (int i = 0; i < mvNum; i++) {
		du(i) = duFull(i * mvNum);
	}
	return du;
}

// The Vector2f prevents the extension from 2X2 
Eigen::VectorXf MPC2X2::yhAdd(Eigen::Vector2f du_) {

	for (int i = 0; i < cvNum; i++) {
		y_h.block(i * predH, 0, predH, 1) = y_shift * y_h.block(i * predH, 0, predH, 1);
	}

	y_h = y_h + X_modified * du_;
	return y_h;
}


// Setters

// The Vector2f prevents the extension from 2X2 
void MPC2X2::modelError(Eigen::Vector2f y_measure) {
	Eigen::Vector2f de;
	de.setZero();

	Eigen::VectorXf ones = createOnes(predH);

	for (int i = 0; i < cvNum; i++) {
		de(i) = y_measure(i) - y_h(i*predH);
		y_h.block(i * predH, 0, predH, 1) = y_h.block(i * predH, 0, predH, 1) + de(i) * ones;
	}



}

void MPC2X2::y_hSet(Eigen::VectorXf y) {
	Eigen::VectorXf ones = createOnes(predH);

	for (int i = 0; i < cvNum; i++) {
		y_h.block(i*predH, 0, predH, 1) = y(i)*ones;
	}

}

void MPC2X2::y11Set(float y, int i) {
	y11(i) = y;
}

void MPC2X2::y12Set(float y, int i) {
	y12(i) = y;
}
void MPC2X2::y21Set(float y, int i) {
	y21(i) = y;
}
void MPC2X2::y22Set(float y, int i) {
	y22(i) = y;
}
// Getters
Eigen::MatrixXf MPC2X2::getX() {
	return X;
}

Eigen::MatrixXf MPC2X2::getXT() {
	return XT;
}

Eigen::MatrixXf MPC2X2::getInv_XTXW() {
	return inv_XTXW;

}

Eigen::MatrixXf MPC2X2::getH_shift() {
	return H_shift;
}
Eigen::VectorXf MPC2X2::get_y11() {
	return y11;
}

Eigen::MatrixXf MPC2X2::getX_modified() {
	return X_modified;
}

Eigen::VectorXf MPC2X2::getSetPoint() {
	return setPoint;
}

Eigen::MatrixXf MPC2X2::getW() {
	return W;
}

Eigen::VectorXf MPC2X2::get_y_h() {
	return y_h;
}




Eigen::VectorXf MPC2X2::createOnes(int length) {
	Eigen::VectorXf result(length);
	result.setOnes();
	return result;
}

