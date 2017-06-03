#include <iostream>
#include <math.h> 
#include "kalman_filter.h"
#include "tools.h"

using namespace std;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    // Applying state transition matrix
    x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
    
    // Recomputing state covariance matrix
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    // Getting prediction in measurement space
    VectorXd z_pred = H_ * x_;
    cout << "z_pred = " << z_pred << endl;
	
    // Computing difference between prediction and measurement
    cout << "z = " << z << endl;
    VectorXd y = z - z_pred;
	
    // Adding state transition covariance and measurement covariance matrices
    MatrixXd Ht = H_.transpose();
    MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
    
    // Computing Kalman gain
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	// Updating estimate
	x_ = x_ + (K * y);
	
    // Recomputing state covariance matrix
    long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // Getting state parameters
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	// Computing radar range
    float rho = sqrt(px * px + py * py);
	
    // Checking denominator for division
    cout << "rho = " << rho << endl;
    if (rho < 0.0001) {
	    cout << "Error in KalmanFilter::UpdateEKF() - No division by zero." << endl;
	}
    
    else {
        // Computing radar angle and range rate from Cartesion coordinates
        float phi = atan2(py, px);
        float rho_dot = (px * vx + py * vy) / rho;
        
        // Instantiating measurement vector in polar coordinates
        VectorXd z_pred = VectorXd(3);
        z_pred << rho, phi, rho_dot;
        cout << "z_pred size = " << z_pred.size() << endl;
        
        // Computing difference between prediction and measurement
        VectorXd y = z - z_pred;
        
        // Normalizing angle in [-pi,pi]
        y(1) = remainder(y(1), 2.0 * M_PI);
        
        // Computing measurement jacobian matrix Hjs
        cout << ".x_ = " << x_ << endl;
        Tools tools;
        MatrixXd Hj_ = tools.CalculateJacobian(x_);
        
        // Adding state transition covariance and measurement covariance matrices
        MatrixXd Hjt = Hj_.transpose();
        MatrixXd PHt = P_ * Hjt;
        MatrixXd S = Hj_ * PHt + R_;
        
        // Computing Kalman gain
        MatrixXd Si = S.inverse();
        MatrixXd K = PHt * Si;

        // Updating estimate
        x_ = x_ + (K * y);
        
        // Recomputing state covariance matrix
        long x_size = x_.size();
        MatrixXd I = MatrixXd::Identity(x_size, x_size);
        P_ = (I - K * Hj_) * P_;
    }
}
