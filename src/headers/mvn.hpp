#ifndef MVN_HPP
#define MVN_HPP


#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <vector>
#include <random>

std::random_device dev;
std::mt19937 engine(dev());
typedef boost::mt19937 PRNG_s;
PRNG_s rng(engine()); //Generate non-static random numbers (pick different numbers from prio distribution each run)

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]


class Mvn
{
public:
    
    Eigen::VectorXd mean;
    Eigen::MatrixXd sigma;
    
    Mvn(const Eigen::VectorXd& mu, const Eigen::MatrixXd& s) :mean(mu), sigma(s) {
        
    };
    
    Eigen::VectorXd sample(int nr_iterations)
    {
        int n = mean.rows();
        
        // Generate x from the N(0, I) distribution
        Eigen::VectorXd x(n);
        Eigen::VectorXd sum(n);
        sum.setZero();
        
        for (unsigned int i = 0; i < nr_iterations; i++)
        {
            x.setRandom();
            x = 0.5 * (x + Eigen::VectorXd::Ones(n));
            sum = sum + x;
        }
        sum = sum - (static_cast<double>(nr_iterations) / 2) * Eigen::VectorXd::Ones(n);
        x = sum / (std::sqrt(static_cast<double>(nr_iterations) / 12));
    
        
        // Find the eigen vectors of the covariance matrix
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd > eigen_solver(sigma);  //computes Eigenvalues and vectors of self adjoint matrices
        Eigen::MatrixXd eigenvectors = eigen_solver.eigenvectors().real();   // get Eigenvalues
        Eigen::MatrixXd eigenvalues = eigen_solver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();  //find max then square root of eigen values x times vectors and take diganoal
        
        
        // Find the transformation matrix
        Eigen::MatrixXd Q = eigenvectors * eigenvalues;
        
        return Q * x + mean;
    }
    
    void checkTrunc(Eigen::VectorXd sampleCheck, Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound, bool& parOOS)
    {
        for (int i = 0; i < sampleCheck.size(); i++){
            if ((sampleCheck(i) < lowerBound(i)) || (sampleCheck(i) > upperBound(i))){
                parOOS = true;
                return;
            }
        }
        parOOS = false;
    }
    
    Eigen::VectorXd sampleTrunc(Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound, int nr_iterations)
    {
        bool parOOS = FALSE;
        Eigen::VectorXd sampleCheck;
        sampleCheck = sample(nr_iterations);
        checkTrunc(sampleCheck, lowerBound, upperBound, parOOS);
        
        while (parOOS){
            sampleCheck = sample(nr_iterations);
            checkTrunc(sampleCheck, lowerBound, upperBound, parOOS);
        }
        return sampleCheck;
    }
};


#endif
