#pragma once
#include "Eigen/Dense"
#include <iostream>


class MPC2X2 {
private:
	const int predH = 100;
	const int manH = 5;
	const int mvNum = 2;
	const int cvNum = 2;
	const int Xrow = predH * cvNum;
	const int Xcol = manH * mvNum;
	float dt;

	Eigen::Vector2f setPt;
	Eigen::VectorXf setPoint;
	
	Eigen::VectorXf y11;
	Eigen::VectorXf y12;
	Eigen::VectorXf y21;
	Eigen::VectorXf y22;
	Eigen::MatrixXf H_shift;  // used to shift the values down
	Eigen::MatrixXf y_shift;  // used to shift y_h vector down a timestep

	Eigen::MatrixXf X;
	Eigen::MatrixXf XT;  // store this so I compute less

	Eigen::MatrixXf W;  // weight matrix
	Eigen::Vector2f weights; // weights for individual MVs

	
	Eigen::MatrixXf inv_XTXW; // manipulated matrix of X and W
	Eigen::VectorXf duFull;
	Eigen::Vector2f du; // manipulated variables
	
	Eigen::MatrixXf X_modified;  // modified X for du -> multiplier
	Eigen::VectorXf y_h; // homogeneous y

	/// <summary>
	/// Creates a predH X manH matrix, used to insert into the X matrix
	/// </summary>
	/// 
	Eigen::MatrixXf makeInsertMatrix(Eigen::VectorXf y);




	/// <summary>
	/// Creates a ones Vector on size length
	/// </summary>
	/// <param name="length"></param>
	/// <returns></returns>
	Eigen::VectorXf createOnes(int length);



public:
	/// <summary>
	/// Constructor for this MPC class
	/// </summary>
	/// <param name="predH"></param> Prediction Horizon
	/// <param name="manH"></param> Manipulation Horizon
	/// <param name="dt"></param> the time value of each discrete increment
	MPC2X2(int predH, int manH, float dt, Eigen::Vector2f setPoint, Eigen::Vector2f weights);

	/// <summary>
	/// This function sets the y11 vector to a particular value at the index i;
	/// therefore create a loop to add all the values to the y11 vector
	/// </summary>
	/// <param name="y"></param> The value
	/// <param name="i"></param> The index
	void y11Set(float y, int i);

	/// <summary>
	/// This function sets the y11 vector to a particular value at the index i;
	/// therefore create a loop to add all the values to the y11 vector
	/// </summary>
	/// <param name="y"></param> The value
	/// <param name="i"></param> The index
	void y12Set(float y, int i);

	/// <summary>
	/// This function sets the y11 vector to a particular value at the index i;
	/// therefore create a loop to add all the values to the y11 vector
	/// </summary>
	/// <param name="y"></param> The value
	/// <param name="i"></param> The index
	void y21Set(float y, int i);

	/// <summary>
	/// This function sets the y11 vector to a particular value at the index i;
	/// therefore create a loop to add all the values to the y11 vector
	/// </summary>
	/// <param name="y"></param> The value
	/// <param name="i"></param> The index
	void y22Set(float y, int i);


	/// <summary>
	/// This sets an initial steady-state value for the control variables
	/// </summary>
	/// <param name="y"></param>
	void y_hSet(Eigen::VectorXf y);


	/// <summary>
	/// This function makes the X and inv_XTX matrices from the provided information
	/// </summary>
	/// <param name=""></param>
	void makeMatrices(void);


	/// <summary>
	/// This method resizes all the matrices and vectors based on
	/// Constructor inputs
	/// </summary>
	void resizeAll();

	/// <summary>
	/// Adds constant model error to the y_h vector
	/// based on the difference between measured and 
	/// predicted values
	/// </summary>
	/// <param name="y_measure"></param>
	void modelError(Eigen::Vector2f y_measure);

	/// <summary>
	/// This adds the new MV changes to the homogeneous solution
	/// </summary>
	/// <param name="du"></param>
	Eigen::VectorXf yhAdd(Eigen::Vector2f du);

	Eigen::MatrixXf getX();
	Eigen::MatrixXf getXT();
	Eigen::MatrixXf getInv_XTXW();
	Eigen::MatrixXf getH_shift();
	Eigen::VectorXf get_y11();
	Eigen::MatrixXf getX_modified();
	Eigen::VectorXf getSetPoint();
	
	Eigen::MatrixXf getW();
	Eigen::VectorXf get_y_h();


	/// <summary>
	/// This is the work horse of this program
	/// input is the homogeneous solution y_h
	/// User is forced to insert it, in order to be 
	/// aware of its values
	/// </summary>
	/// <param name="error"></param>
	/// <returns></returns>
	Eigen::Vector2f compute(Eigen::VectorXf y_h);


};