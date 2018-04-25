#include "eigensolver.h"
#include <iostream>

#ifdef USING_CHOLMOD_MODULE
#include <cholmod.h>
#include <SuiteSparse_config.h>
typedef SuiteSparse_long UF_long;
#include <Eigen/CholmodSupport>
#endif

extern "C"  {
void dsaupd_(
        int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid,
        int *ncv, double *V, int *ldv, int *iparam, int *ipntr, double *workd,
        double *workl, int *lworkl, int *info
);

void dseupd_(
        int *rvec, char *howmany, int *select, double *d, double *Z, int *ldz, double *sigma, char *bmat,
        int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv, int *iparam,
        int *ipntr, double *workd, double *workl, int *lworkl, int *info
);
}

SparseEigenSolver::~SparseEigenSolver() {
    delete [] workspace;
    delete [] workspaced;
}

SparseEigenSolver::SparseEigenSolver(Eigen::SparseMatrix<double> _matrix, int _spectrumSize):
        matrix(_matrix),spectrumSize(_spectrumSize),bandSize(BAND_SIZE){
    eigenvalues.reserve(spectrumSize);
    eigenvectors.reserve(spectrumSize);
    int n_a = matrix.rows();
    int ncv_a = int(bandSize*MAXNCV_FACTOR);
    if(ncv_a > n_a) ncv_a = n_a;
    workspace = new double [n_a*ncv_a];
    workspaced = new double [bandSize];
#ifdef USING_CHOLMOD_MODULE
    std::cout << "Using SUITESPARSE Cholmod module for solving linear system " << std::endl << std::flush;
#endif
}

void SparseEigenSolver::computeSpectrum() {
    int dim = matrix.rows();
    int eigen_index = 0;
    double current_eigenvalue = -1e-8, last_eigenvalue = -1e-8;
    shiftMatrix(current_eigenvalue);
    while (eigen_index < spectrumSize){
        solverOneBand();
        double min_eigenvalue = 1e20, max_eigenvalue = -1e20;
        int last_index = eigen_index;
        std::vector<EigenPair> eigenpairs;
        for (int i = 0 ; i < bandSize; ++i){
            double next_eigenvalue =  current_eigenvalue + getEigenvalueformBand(i);
            if (next_eigenvalue > last_eigenvalue){
                min_eigenvalue = next_eigenvalue < min_eigenvalue ? next_eigenvalue : min_eigenvalue;
                max_eigenvalue = next_eigenvalue > max_eigenvalue ? next_eigenvalue : max_eigenvalue;
                eigenpairs.push_back(EigenPair(next_eigenvalue,getEigenvectorformBand(i)));
                eigen_index ++;
            }
        }
        std::sort(eigenpairs.begin(),eigenpairs.end());
        for (int i = 0; i < eigenpairs.size(); ++i){
            int next_index = last_index + i;
            if (next_index == spectrumSize) break;
            EigenPair pair = eigenpairs[i];
            eigenvalues.push_back(pair.eigenvalue);
            Eigen::VectorXd eigenvector(dim);
            for (int j = 0; j < dim; ++j) eigenvector(j) = pair.eigenvector[j];
            eigenvectors.push_back(eigenvector);
        }
        double increase_value = max_eigenvalue + EIGEN_INTERPOLANT*(max_eigenvalue-min_eigenvalue) - current_eigenvalue;
        shiftMatrix(increase_value);
        current_eigenvalue += increase_value;
        last_eigenvalue = max_eigenvalue;
    }
}

void SparseEigenSolver::solverOneBand() {

#ifdef USING_CHOLMOD_MODULE
    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
#else
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
#endif

    solver.compute(matrix);
    if(solver.info()!=Eigen::Success) {
        // decomposition failed
        std::cout << "decomposition failed" << std::endl;
        exit(0);
    }

    int ido_a = 0;
    char* bmat_a = new char[1];
    bmat_a[0] = 'I';
    int n_a = matrix.rows();
    char* which_a = new char[2];
    which_a[0] = 'L';which_a[1] = 'M';
    int nev_a = bandSize;
    double tol_a = EIG_EPSILON;
    double * resid_a = new double[n_a];
    int ncv_a = int(nev_a*MAXNCV_FACTOR);
    if(ncv_a > n_a) ncv_a = n_a;
    int ldv_a = n_a;
    double * v_a = workspace; //new double[ldv_a * ncv_a];
    int *iparam_a = new int[12];

    int ishfts_a = 1;
    int maxitr_a = MAX_ITERATIONS;
    int mode_a = 1;

    iparam_a[FI(1)] = ishfts_a;
    iparam_a[FI(3)] = maxitr_a;
    iparam_a[FI(7)] = mode_a;

    int *ipntr_a = new int [12];
    double *workd_a = new double [3*n_a];
    for (int i = 0; i < 3*n_a; i++) workd_a[i] = (double)rand() / (double)RAND_MAX;

    int lworkl_a = ncv_a*(ncv_a+8);
    double *workl_a = new double [lworkl_a];

    int info_a = 0;

    bool iterationDone = false;

    while(!iterationDone){
        dsaupd_(&ido_a,bmat_a,&n_a,which_a,&nev_a,&tol_a,resid_a,&ncv_a,v_a,&ldv_a,iparam_a,ipntr_a,workd_a,workl_a,&lworkl_a,&info_a);
        if (ido_a != 1 && ido_a != -1 ) iterationDone = true;
        if (!iterationDone){
            double* rhs_data = workd_a + FI(ipntr_a[FI(1)]);
            double* unknown_data = workd_a + FI(ipntr_a[FI(2)]);
            Eigen::VectorXd rhs(n_a);
            Eigen::VectorXd unknown(n_a);
            for (int i = 0; i < n_a; ++i) rhs(i) = rhs_data[i];
            unknown = solver.solve(rhs);
            for (int i = 0; i < n_a; ++i) unknown_data[i] = unknown(i);
        }
    }

    int rvecs_a = true;
    char *howmany_a = new char[1];
    howmany_a[0] = 'A';
    int * select_a = new int[ncv_a];
    double * d_a =  workspaced; //new double [nev_a];
    double sigma_a;
    int ierr_a;

    dseupd_(&rvecs_a, howmany_a, select_a, d_a, v_a, &ldv_a, &sigma_a,
            bmat_a, &n_a, which_a, &nev_a, &tol_a, resid_a, &ncv_a, v_a,
            &ldv_a, iparam_a, ipntr_a, workd_a, workl_a, &lworkl_a,
            &ierr_a);

    delete [] bmat_a;
    delete [] which_a;
    delete [] resid_a;
    //delete [] v_a;
    delete [] iparam_a;
    delete [] ipntr_a;
    delete [] workd_a;
    delete [] workl_a;
    delete [] howmany_a;
    delete [] select_a;
    //delete [] d_a;
}

void SparseEigenSolver::shiftMatrix(double shift) {
    int dim = matrix.rows();
    Eigen::SparseMatrix <double > identity(dim,dim);
    identity.setIdentity();
    matrix = matrix - shift * identity;
}

double SparseEigenSolver::getEigenvalueformBand(int index){
    int proper_ind = bandSize - 1 - index;
    return 1.0/ (workspaced[proper_ind]);
}

double * SparseEigenSolver::getEigenvectorformBand(int index){
    int proper_ind = bandSize - 1 - index;
    int n_a = matrix.rows();
    return workspace + n_a*proper_ind;
}

Eigen::VectorXd SparseEigenSolver::getEigenvector(int index){
    return eigenvectors[index];
}
double SparseEigenSolver::getEigenvalue(int index) {
    return eigenvalues[index];
}
std::vector<double> SparseEigenSolver::getEigenvalues(){
    return eigenvalues;
}
std::vector<Eigen::VectorXd> SparseEigenSolver::getEigenvectors(){
    return eigenvectors;
}