/*
AHRS_EKF.cpp
Alexander Waldejer
alexander.waldejer@hotmail.com

Orientation and altitude estimation using the MPU9250.
Sensor fusion for orientation estimation is based on the (nonlinear) Extended Kalman Filter.
Altitude measurement frequency and accuracy is enhanced based on the (linear) Kalman Filter,
combining the barometer with linear acceleration (achieved by rotating the raw accelerometer (sensor frame)
with the resulting orientation estimation.
*/

#include "AHRS_EKF.h"

/* AHRS_EKF object */
AHRS_EKF::AHRS_EKF(){
	// Initialize all class variables with values
	initVariables();
}

/* Public funtions */
void AHRS_EKF::updateEKF(float* yAccMeas, float* yGyrMeas, float* yMagMeas){
	// Move over measurements to private variables
	yAcc << yAccMeas[0], yAccMeas[1], yAccMeas[2];
	yGyr << yGyrMeas[0], yGyrMeas[1], yGyrMeas[2]; 
	yMag << yMagMeas[0], yMagMeas[1], yMagMeas[2];

	if (counterCalibration > limitCalibration){
		// Calibration complete: Call each individual step of EKF
		updateAccelerometer();
		qNorm();
		updateGyroscope();
		qNorm();
		//updateMagnetometer();
		//qNorm();
	}
	else if (counterCalibration == limitCalibration){
		Serial.printf("counterCalibration: %i\n", counterCalibration);
		
		// Finish data collection
		// Find mean of raw sensor data		
		mAcc /= (float)limitCalibration;
		mGyr /= (float)limitCalibration;
		mMag /= (float)limitCalibration;
		
		// Assign filter variables
		// g0 = mAcc;
		
		// Find sum of all variances of raw accelerometer sensor data
		for (int i = 0; i < limitCalibration; i++){
			vTemp = vAcc.col(i) - mAcc;
			vTemp << pow(vTemp(0),2), pow(vTemp(1),2), pow(vTemp(2),2);
			RAcc += vTemp.asDiagonal();
		}
		
		// Find sum of all variances of raw gyroscope sensor data
		for (int i = 0; i < limitCalibration; i++){
			vTemp = vGyr.col(i) - mGyr;
			vTemp << pow(vTemp(0),2), pow(vTemp(1),2), pow(vTemp(2),2);
			RGyr += vTemp.asDiagonal();
		}
		
		// Find sum of all variances of raw magnetometer sensor data
		for (int i = 0; i < limitCalibration; i++){
			vTemp = vMag.col(i) - mMag;
			vTemp << pow(vTemp(0),2), pow(vTemp(1),2), pow(vTemp(2),2);
			RMag += vTemp.asDiagonal();
		}
		
		// Divide by total number of samples to get the variance of each sensor
		RAcc /= ((float)limitCalibration - 1.0);
		RGyr /= ((float)limitCalibration - 1.0);
		RMag /= ((float)limitCalibration - 1.0);
		
		// Print resulting calibration data
		Serial.printf("mAcc = % f % f % f\n", mAcc(0), mAcc(1), mAcc(2));
		Serial.printf("mGyr = % f % f % f\n", mGyr(0), mGyr(1), mGyr(2));
		Serial.printf("mMag = % f % f % f\n", mMag(0), mMag(1), mMag(2));
		Serial.println();
		Serial.printf("RAcc:\n");
		print_mtxf(RAcc);
		Serial.printf("RGyr:\n");
		print_mtxf(RGyr);
		Serial.printf("RMag:\n");
		print_mtxf(RMag);
		
		// Increment calibration counter one last time to end calibration
		counterCalibration++;
	}
	else{
		Serial.printf("counterCalibration: %i\n", counterCalibration);
		
		// Collect data for calibration mean
		mAcc += yAcc;
		mGyr += yGyr;
		mMag += yMag;
		
		// Collect data for calibration variance
		vAcc.col(counterCalibration) = yAcc;
		vGyr.col(counterCalibration) = yGyr;
		vMag.col(counterCalibration) = yMag;
									
		// Increment calibration counter
		counterCalibration++;
	}
}

/* Return Euler angles */
void AHRS_EKF::getEulers(float* result){
	// Quaternions to Eulers according to ZYX rotation
	eulerTemp[0] = 2.*pow(q(0),2)-1+2.*pow(q(1),2);
    eulerTemp[1] = 2.*(q(1)*q(2)-q(0)*q(3));
    eulerTemp[2] = 2.*(q(1)*q(3)+q(0)*q(2));
    eulerTemp[3] = 2.*(q(2)*q(3)-q(0)*q(1));
    eulerTemp[4] = 2.*pow(q(0),2)-1+2.*pow(q(3),2);

    result[2] = atan2(eulerTemp[3],eulerTemp[4]); // phi
    result[1] = -atan(eulerTemp[2]/(sqrt(1-pow(eulerTemp[2],2)))); // theta
    result[0] = atan2(eulerTemp[1],eulerTemp[0]); // psi
}

/* Return Euler angles */
void AHRS_EKF::getEulers2(float* result){
	float test = q(0)*q(1) + q(2)*q(3);
	if (test > 0.499) { // singularity at north pole
		result[2] = 2 * atan2(q(0),q(3));
		result[1] = M_PI/2;
		result[0] = 0;
		return;
	}
	if (test < -0.499) { // singularity at south pole
		result[2] = -2 * atan2(q(0),q(3));
		result[1] = - M_PI/2;
		result[0] = 0;
		return;
	}
    double sqx = q(0)*q(0);
    double sqy = q(1)*q(1);
    double sqz = q(2)*q(2);
    result[2] = atan2(2*q(1)*q(3)-2*q(0)*q(2) , 1 - 2*sqy - 2*sqz);
	result[1] = asin(2*test);
	result[0] = atan2(2*q(0)*q(3)-2*q(1)*q(2) , 1 - 2*sqx - 2*sqz);
	
}


/* Return Quaternions */
void AHRS_EKF::getQuaternions(float* result){
	result[0] = q(0);
	result[1] = q(1);
	result[2] = q(2);
	result[3] = q(3);
}

/* Calculate sampling time (delta t sinze last update)*/ 
void AHRS_EKF::updateTime(){
	// Set integration time by time elapsed since last filter update
	timerStop = micros();
	ts = ((timerStop - timerStart) / 1000000.0f);
	timerStart = timerStop;
	//ts=0.5;
	//Serial.printf("Execution time %1.15f\n",ts);
	//updateAccelerometer();
	//updateGyroscope();
	//updateMagnetometer();
}

/* Private functions */
/* Accelerometer step of EKF */
void AHRS_EKF::updateAccelerometer(){
	// Check if measurement is valid
	//  && vAcc1(1) < ( vAcc2(1) * beta ) && vAcc1(2) < ( vAcc2(2) * beta )) && ( vAcc1(0) > ( vAcc2(0) / beta ) && vAcc1(1) > ( vAcc2(1) / beta ) && vAcc1(2) > ( vAcc2(2) / beta )
	if((sqrt(pow(yAcc(0),2)) > beta * sqrt(pow(mAcc(0),2))|| sqrt(pow(yAcc(1),2)) > beta * sqrt(pow(mAcc(1),2)) || sqrt(pow(yAcc(2),2)) > beta * sqrt(pow(mAcc(2),2))) && (sqrt(pow(yAcc(0),2)) < sqrt(pow(mAcc(0),2)) / beta || sqrt(pow(yAcc(1),2)) < sqrt(pow(mAcc(1),2)) / beta || sqrt(pow(yAcc(2),2)) < sqrt(pow(mAcc(2),2)) / beta)){
		//Serial.printf("Acc outlier!\n");
	}
	else{		
		// fka=0.0; 
		// g0=[0;0;9.8];
		// yacc=[0;0;0];
		//yka=Qq(x)'*(g0+fka);
		Qq << 2.0*(pow(q(0),2.0)+pow(q(1),2.0)) - 1.0,  2.0*(q(1)*q(2)-q(0)*q(3)),    2.0*(q(1)*q(3)+q(0)*q(2)),
			2.0*(q(1)*q(2)+q(0)*q(3)),    2.0*(pow(q(0),2)+pow(q(2),2)) - 1.0,  2.0*(q(2)*q(3)-q(0)*q(1)),
			2.0*(q(1)*q(3)-q(0)*q(2)),    2.0*(q(2)*q(3)+q(0)*q(1)),    2.0*(pow(q(0),2.0)+pow(q(3),2.0)) - 1.0;
		yka = Qq.transpose() * (mAcc + fka);
		
		// [h1 h2 h3 h4]=dQqdq(x);
		// hd=[h1'*g0 h2'*g0 h3'*g0 h4'*g0];
		// S=hd*P*hd'+RAcc;
		// K=P*hd'*S^-1;
		
		// The derivative of Qq wrt q(i), i={0,1,2,3}
		h1 << 2.0*q(0),   -q(3),    q(2),
			  q(3),  2.0*q(0),   -q(1),
			 -q(2),    q(1),  2.0*q(0);
		h2 << 2.0*q(1),    q(2),    q(3),
			  q(2),     0.0,   -q(0),
			  q(3),    q(0),     0.0;
		h3 << 0.0,    q(1),    q(0),
			  q(1),  2.0*q(2),    q(3),
			 -q(0),    q(3),     0.0;
		h4 << 0.0,   -q(0),    q(1),
			  q(0),     0.0,    q(2),
			  q(1),    q(2),  2.0*q(3);
		h1 *= 2.0;	
		h2 *= 2.0;	
		h3 *= 2.0;	
		h4 *= 2.0;	
		hd << h1.transpose() * g0,  h2.transpose() * g0,  h3.transpose() * g0,  h4.transpose() * g0;
		S = hd * P * hd.transpose() + RAcc;
		K = P * hd.transpose() * S.inverse();

		// update
		//x=x+K*(yAcc-yka);
		//P=P-K*S*K';
		q += K * (yAcc - yka);
		P -= K * S * K.transpose();
		
		//q = q.normalized();
		//Serial.printf("q=%1.5f, %1.5f, %1.5f, %1.5f\n", q[0], q[1], q[2], q[3]);
		//print_mtxf(P);
	}
}

/* Gyroscope step of EKF */
void AHRS_EKF::updateGyroscope(){
	//G=Sq(x)*0.5*T;
    //F=eye(4)+Somega(omega)*0.5*T;
    //x=F*x; % omega measured contains both signal and noise
    //P=F*P*F'+G*Rw*G';
	Sq << -q(1), -q(2), -q(3),
		q(0), -q(3),  q(2),
		q(3),  q(0), -q(1),
		-q(2),  q(1),  q(0);	
	G = Sq * 0.5 * ts;
	//Serial.printf("G\n");
	//print_mtxf(G);
	
	Sw << 0.0,  -yGyr(0),  -yGyr(1),  -yGyr(2),
		yGyr(0),    0.0,   yGyr(2),  -yGyr(1),
		yGyr(1),  -yGyr(2),    0.0,   yGyr(0),
		yGyr(2),   yGyr(1),  -yGyr(0),    0.0;
	
	F = Eigen::MatrixXf::Identity(4,4) + Sw * 0.5 * ts;
	
	//Serial.printf("F\n");
	//print_mtxf(F);
	
	q = F * q;
	P = F * P * F.transpose() + G * RGyr * G.transpose();
	
	//Serial.printf("q=%1.5f, %1.5f, %1.5f, %1.5f\n", q[0], q[1], q[2], q[3]);
	//qNorm();
	//Serial.printf("qnorm=%1.5f, %1.5f, %1.5f, %1.5f\n", q[0], q[1], q[2], q[3]);
	//print_mtxf(P);
}

/* Magnetometer step of EKF */
void AHRS_EKF::updateMagnetometer(){
	// fka=0.0; 
	// g0=[0;0;9.8];
	// yacc=[0;0;0];
	//yka=Qq(x)'*(g0+fka);
	m0 << 0, sqrt(pow(yMag(0),2)+pow(yMag(1),2)), yMag(2);
	
	//Serial.printf("m0: %f, %f, %f\n", m0(0), m0(1), m0(2));
	
	Qq << 2.0*(pow(q(0),2.0)+pow(q(1),2.0)) - 1.0,  2.0*(q(1)*q(2)-q(0)*q(3)),    2.0*(q(1)*q(3)+q(0)*q(2)),
		2.0*(q(1)*q(2)+q(0)*q(3)),    2.0*(pow(q(0),2)+pow(q(2),2)) - 1.0,  2.0*(q(2)*q(3)-q(0)*q(1)),
		2.0*(q(1)*q(3)-q(0)*q(2)),    2.0*(q(2)*q(3)+q(0)*q(1)),    2.0*(pow(q(0),2.0)+pow(q(3),2.0)) - 1.0;
	
	//Serial.printf("Qq\n");
	//print_mtxf(Qq);
	
	yka = Qq.transpose() * (m0 + fkm);
	
	//Serial.printf("yka: %f, %f, %f\n", yka(0), yka(1), yka(2));
	
	// [h1 h2 h3 h4]=dQqdq(x);
	// hd=[h1'*m0 h2'*m0 h3'*m0 h4'*m0];
	// S=hd*P*hd'+RMag;
	// K=P*hd'*S^-1;
	
	// The derivative of Qq wrt q(i), i={0,1,2,3}
	h1 << 2.0*q(0),   -q(3),    q(2),
		  q(3),  2.0*q(0),   -q(1),
		 -q(2),    q(1),  2.0*q(0);
	h2 << 2.0*q(1),    q(2),    q(3),
		  q(2),     0.0,   -q(0),
		  q(3),    q(0),     0.0;
	h3 << 0.0,    q(1),    q(0),
		  q(1),  2.0*q(2),    q(3),
		 -q(0),    q(3),     0.0;
	h4 << 0.0,   -q(0),    q(1),
		  q(0),     0.0,    q(2),
		  q(1),    q(2),  2.0*q(3);
	h1 *= 2.0;	
	h2 *= 2.0;	
	h3 *= 2.0;	
	h4 *= 2.0;	
	hd << h1.transpose() * m0,  h2.transpose() * m0,  h3.transpose() * m0,  h4.transpose() * m0;
	S = hd * P * hd.transpose() + RMag;
	K = P * hd.transpose() * S.inverse();

	// update
	//x=x+K*(yAcc-yka);
	//P=P-K*S*K';
	q += K * (yMag - yka);
	P -= K * S * K.transpose();
	
	//Serial.printf("q=%1.5f, %1.5f, %1.5f, %1.5f\n", q[0], q[1], q[2], q[3]);
	//qNorm();
	//Serial.printf("qnorm=%1.5f, %1.5f, %1.5f, %1.5f\n", q[0], q[1], q[2], q[3]);
	//print_mtxf(P);
}

/* Initialize variables in class constructor */
void AHRS_EKF::initVariables(){
	// Matrices
	Qq = Eigen::MatrixXf::Zero(3,3);
	P = Eigen::MatrixXf::Identity(4,4);
	RAcc = Eigen::MatrixXf::Identity(3,3);
	RGyr = Eigen::MatrixXf::Identity(3,3);
	RMag = Eigen::MatrixXf::Identity(3,3);
	S = Eigen::MatrixXf::Zero(3,3);
	K = Eigen::MatrixXf::Zero(4,3);
	h1 = Eigen::MatrixXf::Zero(3,3);
	h2 = Eigen::MatrixXf::Zero(3,3);
	h3 = Eigen::MatrixXf::Zero(3,3);
	h4 = Eigen::MatrixXf::Zero(3,3);
	hd = Eigen::MatrixXf::Zero(3,4);
	Sq = Eigen::MatrixXf::Zero(4,3);
	Sw = Eigen::MatrixXf::Zero(4,4);
	G = Eigen::MatrixXf::Zero(4,3);
	F = Eigen::MatrixXf::Zero(4,4);
	vAcc = Eigen::MatrixXf::Zero(3,limitCalibration);
	vGyr = Eigen::MatrixXf::Zero(3,limitCalibration);
	vMag = Eigen::MatrixXf::Zero(3,limitCalibration);
	
	// Vectors
	q << 1.0, 0.0, 0.0, 0.0;
	g0 << 0.0, 0.0, 9.8;
	m0 << 0.0, 0.0, 0.0;	
	fka << 0.0, 0.0, 0.0;
	fkm << 0.0, 0.0, 0.0;
	yka << 0.0, 0.0, 0.0;
	yAcc << 0.0, 0.0, 0.0;
	yGyr << 0.0, 0.0, 0.0;
	yMag << 0.0, 0.0, 0.0;
	mAcc << 0.0, 0.0, 0.0;
	mGyr << 0.0, 0.0, 0.0;
	mMag << 0.0, 0.0, 0.0;
	vTemp << 0.0, 0.0, 0.0;
	
	// Variables
	counterCalibration = 0;
}

/* Print matrix function */
void AHRS_EKF::print_mtxf(const Eigen::MatrixXf& X){
   int i, j, nrow, ncol;
   
   nrow = X.rows();
   ncol = X.cols();

   //Serial.print("nrow: "); Serial.println(nrow);
   //Serial.print("ncol: "); Serial.println(ncol);       
   //Serial.println();
   
   for (i=0; i<nrow; i++)
   {
       for (j=0; j<ncol; j++)
       {
           Serial.print(X(i,j), 6);   // print 6 decimal places
           Serial.print(", ");
       }
       Serial.println();
   }
   Serial.println();
}

/* Quaternions normalized */
void AHRS_EKF::qNorm(){
	q = q.normalized();
	if(q(0) < 0.0){
		q *= -1;
	}
}




















