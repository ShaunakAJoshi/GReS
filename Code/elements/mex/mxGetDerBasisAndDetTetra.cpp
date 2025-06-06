#include "mex.h"
#include <cmath>

// Compute determinant of 3x3 matrix (column-major)
double det3x3(const double* M) {
    return M[0]*(M[4]*M[8] - M[5]*M[7])
         - M[1]*(M[3]*M[8] - M[5]*M[6])
         + M[2]*(M[3]*M[7] - M[4]*M[6]);
}

// Compute inverse of 3x3 matrix (column-major)
void inv3x3(const double* M, double* invM) {
    double det = det3x3(M);
    if (fabs(det) < 1e-15) {
        mexErrMsgIdAndTxt("mexFunction:singularMatrix", "Jacobian matrix is singular.");
    }
    double invDet = 1.0 / det;

    invM[0] =  (M[4]*M[8] - M[5]*M[7]) * invDet;
    invM[1] = -(M[1]*M[8] - M[2]*M[7]) * invDet;
    invM[2] =  (M[1]*M[5] - M[2]*M[4]) * invDet;
    invM[3] = -(M[3]*M[8] - M[5]*M[6]) * invDet;
    invM[4] =  (M[0]*M[8] - M[2]*M[6]) * invDet;
    invM[5] = -(M[0]*M[5] - M[2]*M[3]) * invDet;
    invM[6] =  (M[3]*M[7] - M[4]*M[6]) * invDet;
    invM[7] = -(M[0]*M[7] - M[1]*M[6]) * invDet;
    invM[8] =  (M[0]*M[4] - M[1]*M[3]) * invDet;
}

// Multiply 3x4 by 4x3 to get 3x3: out = A*B
void mul3x4_4x3(const double* A, const double* B, double* out) {
    for (int col = 0; col < 3; ++col) {
        for (int row = 0; row < 3; ++row) {
            double sum = 0.0;
            for (int k = 0; k < 4; ++k) {
                sum += A[row + 3*k] * B[k + 4*col];
            }
            out[row + 3*col] = sum;
        }
    }
}

// Multiply 3x3 by 3x4 to get 3x4: out = A*B
void mul3x3_3x4(const double* A, const double* B, double* out) {
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 3; ++row) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) {
                sum += A[row + 3*k] * B[k + 3*col];
            }
            out[row + 3*col] = sum;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("mexFunction:invalidNumInputs", "3 inputs required: J1(3x4), coords(4x3), weights(nGauss)");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("mexFunction:invalidNumOutputs", "2 outputs required: N (3x4xN), detJ (scalar)");
    }

    // Inputs
    const mxArray* J1_mx = prhs[0];
    const mxArray* coords_mx = prhs[1];
    const mxArray* weights_mx = prhs[2];

    // Validate sizes
    const mwSize* dimsJ1 = mxGetDimensions(J1_mx);
    if (mxGetNumberOfDimensions(J1_mx) != 2 || dimsJ1[0] != 3 || dimsJ1[1] != 4) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "J1 must be 3 x 4");
    }
    const mwSize* dimsCoords = mxGetDimensions(coords_mx);
    if (mxGetNumberOfDimensions(coords_mx) != 2 || dimsCoords[0] != 4 || dimsCoords[1] != 3) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "coords must be 4 x 3");
    }

    const mwSize* dimsW = mxGetDimensions(weights_mx);
    mwSize nGauss = (dimsW[0] > dimsW[1]) ? dimsW[0] : dimsW[1];

    const double* J1 = mxGetPr(J1_mx);
    const double* coords = mxGetPr(coords_mx);
    const double* weights = mxGetPr(weights_mx);

    // Allocate outputs
    mwSize dimsN[3] = {3, 4, nGauss};
    plhs[0] = mxCreateNumericArray(3, dimsN, mxDOUBLE_CLASS, mxREAL);
    double* N = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleScalar(0.0); // scalar detJ
    double* detJ = mxGetPr(plhs[1]);

    // Compute Jacobian J = J1 * coords (3x4 * 4x3 = 3x3)
    double J[9], invJ[9], Ntemp[12];  // J: 3x3, invJ: 3x3, Ntemp: 3x4
    mul3x4_4x3(J1, coords, J);
    *detJ = det3x3(J);

    // Apply weight sum (since constant detJ and constant integrand)
    double weightSum = 0.0;
    for (mwSize i = 0; i < nGauss; ++i) {
        weightSum += weights[i];
    }
    *detJ *= weightSum;

    // Compute N = invJ * J1 (3x3 * 3x4 = 3x4)
    inv3x3(J, invJ);
    mul3x3_3x4(invJ, J1, Ntemp);

    // Fill all slices of output N(:,:,i) = Ntemp (same for each Gauss point)
    for (mwSize i = 0; i < nGauss; ++i) {
        for (int j = 0; j < 3*4; ++j) {
            N[i*3*4 + j] = Ntemp[j];
        }
    }
}
