#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <Eigen/Sparse>

#define BAND_SIZE 50
#define EIGEN_INTERPOLANT 0.4


#define MAXNCV_FACTOR 3.5 // taken from GRAPHITE
#define EIG_EPSILON 1E-10
#define MAX_ITERATIONS 10000
#define FI(x) x-1

class SparseEigenSolver{
public:
    SparseEigenSolver(Eigen::SparseMatrix <double > _matrix, int _spectrumSize);
    ~SparseEigenSolver();
    void computeSpectrum();

    std::vector<double> getEigenvalues();
    std::vector<Eigen::VectorXd> getEigenvectors();
    double getEigenvalue(int);
    Eigen::VectorXd getEigenvector(int);

private:
    SparseEigenSolver();
    class EigenPair{
    public:
        EigenPair(double _eigenvalue, double * _eigenvector):eigenvalue(_eigenvalue),eigenvector(_eigenvector){}
        inline bool operator<(const EigenPair& pair2) const  {
            return eigenvalue < pair2.eigenvalue;
        }
        double eigenvalue;
        double * eigenvector;
    };
    double getEigenvalueformBand(int);
    double * getEigenvectorformBand(int);
    void solverOneBand();
    void shiftMatrix(double);


    Eigen::SparseMatrix <double >matrix;
    int bandSize;
    int spectrumSize;
    double * workspace;
    double * workspaced;
    std::vector<double> eigenvalues;
    std::vector<Eigen::VectorXd> eigenvectors;
};

#endif
