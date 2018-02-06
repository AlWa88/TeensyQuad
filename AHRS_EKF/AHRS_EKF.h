/*
AHRS_EKF.h
Alexander Waldejer
alexander.waldejer@hotmail.com

Orientation and altitude estimation using the MPU9250.
Sensor fusion for orientation estimation is based on the (nonlinear) Extended Kalman Filter.
Altitude measurement frequency and accuracy is enhanced based on the (linear) Kalman Filter,
combining the barometer with linear acceleration (achieved by rotating the raw accelerometer (sensor frame)
with the resulting orientation estimation.
*/

#ifndef AHRS_EKF_h
#define AHRS_EKF_h

#include "Eigen30.h"
#include "Eigen/LU"
#include "Arduino.h"
#include <math.h>

//using namespace Eigen;

class AHRS_EKF{
	public:
	// Variables
	
	// Functions
	AHRS_EKF();
	void updateEKF(float* yAccMeas, float* yGyrMeas, float* yMagMeas);
	void getQuaternions(float* result);
	void getEulers(float* result);
	void getEulers2(float* result);
	void updateTime();
	
	private:
	// Variables
	Eigen::MatrixXf Qq;
	Eigen::MatrixXf P;
	Eigen::MatrixXf RAcc;
	Eigen::MatrixXf RGyr;
	Eigen::MatrixXf RMag;	
	Eigen::MatrixXf S;
	Eigen::MatrixXf K;	
	Eigen::MatrixXf h1;
	Eigen::MatrixXf h2;
	Eigen::MatrixXf h3;
	Eigen::MatrixXf h4;
	Eigen::MatrixXf hd;
	Eigen::MatrixXf Sq;
	Eigen::MatrixXf Sw;
	Eigen::MatrixXf G;
	Eigen::MatrixXf F;
	Eigen::MatrixXf vAcc;
	Eigen::MatrixXf vGyr;
	Eigen::MatrixXf vMag;
	
	Eigen::Vector4f q;
	Eigen::Vector3f g0;	
	Eigen::Vector3f m0;	
	Eigen::Vector3f fka;
	Eigen::Vector3f fkm;
	Eigen::Vector3f yka;
	Eigen::Vector3f yAcc;
	Eigen::Vector3f yGyr;
	Eigen::Vector3f yMag;
	Eigen::Vector3f mAcc;
	Eigen::Vector3f mGyr;
	Eigen::Vector3f mMag;
	Eigen::Vector3f vTemp;
	
	float ts; // delta t since last measurement
	float eulerTemp[5]; // Temp variable used for q2euler convertion
	const int limitCalibration = 100; // desired filter calibration limit of data samples to reach
	int counterCalibration; // counter for calibration
	uint32_t timerStart = 0, timerStop = 0;
	float beta = 0.5, alpha = 0.1;
	
	// Functions
	void updateAccelerometer();
	void updateGyroscope();
	void updateMagnetometer();
	void qNorm();
	void initVariables();
	void print_mtxf(const Eigen::MatrixXf& X);
};
#endif
