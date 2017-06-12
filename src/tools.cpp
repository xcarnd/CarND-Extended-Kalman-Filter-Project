#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    assert(estimations.size() == ground_truth.size() && estimations.size() > 0);

    VectorXd v(4);
    v << 0, 0, 0, 0;

    for (std::size_t i = 0; i < estimations.size(); ++i) {
        VectorXd diff = estimations[i] - ground_truth[i];
        v += diff.array() * diff.array();
    }
    v /= estimations.size();
    return v.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    double px = x_state[0];
    double py = x_state[1];
    double vx = x_state[2];
    double vy = x_state[3];
    double b = px * px + py * py;
    double a = sqrt(b);
    double c = b * a;
    MatrixXd H(3, 4);
    H << px / a, py / a, 0, 0,
         -py / b,-px / b, 0, 0,
         py * (vx * py - vy * px) / c, px * (vy * px - vx * py) / c, px / a, py / a;
    return H;
}
