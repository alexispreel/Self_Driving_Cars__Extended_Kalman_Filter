#include <iostream>
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    previous_timestamp_ = 0;
    
    // Creating state vector
	ekf_.x_ = VectorXd(4);

	// Initializing state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;
    
    // Setting measurement extraction matrix (laser only)
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;
    
    // Initializing measurement jacobian
    Hj_ = MatrixXd(3, 4);

    // Initializing measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
            0, 0.0225;
    
    // Initializing measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;
    
    // Initializing transition matrix F
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;
    
	// Setting acceleration noise components
	noise_ax = 9;
	noise_ay = 9;
    
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            // Getting initial radar measurement: range, angle, range rate
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            float rho_dot = measurement_pack.raw_measurements_[2];
            
            // Converting to cartesian coordinates
            ekf_.x_ << rho * cos(phi), rho * sin(phi), rho_dot * cos(phi), rho_dot * sin(phi);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            // Getting initial lidar measurement: position
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }

        // Done initializing, no need to predict or update
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/
    // Computing time elapsed between measurements in seconds
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;
    
    // Modifying transition F matrix with time elapsed
	ekf_.F_(0,2) = dt;
	ekf_.F_(1,3) = dt;
    
    // Computing differential factors for Q matrix
    float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;
    
    // Setting process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4 * noise_ax / 4, 0, dt_3 * noise_ax / 2, 0,
			  0, dt_4 * noise_ay / 4, 0, dt_3 * noise_ay / 2,
			  dt_3 * noise_ax / 2, 0, dt_2 * noise_ax, 0,
			  0, dt_3 * noise_ay / 2, 0, dt_2 * noise_ay;
    
    // Predicting
    ekf_.Predict();

    /*****************************************************************************
    *  Update
    ****************************************************************************/
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Choosing radar measurement covariance matrix
        ekf_.R_ = R_radar_;
        
        // Updating state with radar measurement
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        
    } else {
        // Choosing Lidar measurement extraction matrix
        ekf_.H_ = H_laser_;
        
        // Choosing Lidar measurement covariance matrix
        ekf_.R_ = R_laser_;
        
        // Updating state with lidar measurement
        cout << "measurement_pack.raw_measurements_ = " << measurement_pack.raw_measurements_ << endl;
        ekf_.Update(measurement_pack.raw_measurements_);
        
    }
    
    // Printing output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
