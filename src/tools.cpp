#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
    TODO:
    * Calculate the RMSE here.
    */
    // Initializing RMSE
    VectorXd rmse(4);
	rmse << 0,0,0,0;
    
    // Checking inputs
	if(estimations.size() != ground_truth.size() || estimations.size() == 0){
		cout << "Invalid inputs for RMSE" << endl;
		return rmse;
	}
    
	// Cumulating squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
        
		VectorXd r = estimations[i] - ground_truth[i];
        
		//coefficient-wise multiplication
		r = r.array() * r.array();
		rmse += r;
	}
    
	// Calculating mean
	rmse = rmse / estimations.size();
    
	// Calculating squared root
	rmse = rmse.array().sqrt();
    
	// Return result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    // Initializing jacobian
    MatrixXd Hj(3,4);
    
	// Getting state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// Computing terms
    float d1 = px * px + py * py;
	float d2 = sqrt(d1);
	float d3 = (d1 * d2);
	
    // Checking denominator for division
    if (d1 < 0.0001) {
	    cout << "Error in CalculateJacobian() - No division by zero." << endl;
	}
    
    else {
	// Computing Jacobian
        Hj <<   px / d2, py / d2, 0, 0,
                -py / d1, px / d1, 0, 0,
                py * (vx * py - vy * px) / d3, px * (vy * px - vx * py) / d3, px / d2, py / d2;
	}
	return Hj;
}
