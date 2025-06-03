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

// Multiply 3x8 by 8x3 to get 3x3: out = A*B
void mul3x8_8x3(const double* A, const double* B, double* out) {
    for (int col=0; col<3; ++col) {
        for (int row=0; row<3; ++row) {
            double sum = 0.0;
            for (int k=0; k<8; ++k) {
                sum += A[row + 3*k] * B[k + 8*col];
            }
            out[row + 3*col] = sum;
        }
    }
}

// Multiply 3x3 by 3x8 to get 3x8: out = A*B
void mul3x3_3x8(const double* A, const double* B, double* out) {
    for (int col=0; col<8; ++col) {
        for (int row=0; row<3; ++row) {
            double sum = 0.0;
            for (int k=0; k<3; ++k) {
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
        mexErrMsgIdAndTxt("mexFunction:invalidNumInputs", "3 inputs required: J1(3x8xN), coords(8x3), weights(N)");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("mexFunction:invalidNumOutputs", "2 outputs required: N (3x8xN), detJ (Nx1)");
    }

    // Inputs
    const mxArray* J1_mx = prhs[0];
    const mxArray* coords_mx = prhs[1];
    const mxArray* weights_mx = prhs[2];

    // Validate sizes
    const mwSize* dimsJ1 = mxGetDimensions(J1_mx);
    mwSize nGauss = dimsJ1[2];
    if (mxGetNumberOfDimensions(J1_mx) != 3 || dimsJ1[0] != 3 || dimsJ1[1] != 8) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "J1 must be 3 x 8 x nGauss");
    }
    const mwSize* dimsCoords = mxGetDimensions(coords_mx);
    if (mxGetNumberOfDimensions(coords_mx) != 2 || dimsCoords[0] != 8 || dimsCoords[1] != 3) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "coords must be 8 x 3");
    }
    const mwSize* dimsW = mxGetDimensions(weights_mx);
    if (!( (dimsW[0] == nGauss && dimsW[1] == 1) || (dimsW[1] == nGauss && dimsW[0] == 1) )) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "weights must be a vector of length nGauss");
    }

    const double* J1 = mxGetPr(J1_mx);
    const double* coords = mxGetPr(coords_mx);
    const double* weights = mxGetPr(weights_mx);

    // Allocate outputs
    mwSize dimsN[3] = {3, 8, nGauss};
    plhs[0] = mxCreateNumericArray(3, dimsN, mxDOUBLE_CLASS, mxREAL);
    double* N = mxGetPr(plhs[0]);

    mwSize dimsDet[2] = {nGauss, 1};
    plhs[1] = mxCreateNumericArray(2, dimsDet, mxDOUBLE_CLASS, mxREAL);
    double* detJ = mxGetPr(plhs[1]);

    // Temporary matrices
    double J[9];       // Jacobian 3x3
    double invJ[9];    // Inverse Jacobian 3x3
    double Ntemp[24];  // 3x8 matrix temporary

    for (mwSize gp = 0; gp < nGauss; ++gp) {
        // Compute J(:,:,gp) = J1(:,:,gp) * coords (3x8 * 8x3 = 3x3)
        mul3x8_8x3(J1 + gp*3*8, coords, J);

        // Compute determinant times weight
        detJ[gp] = det3x3(J) * weights[gp];

        // Compute inverse Jacobian
        inv3x3(J, invJ);

        // Compute N(:,:,gp) = invJ * J1(:,:,gp) (3x3 * 3x8 = 3x8)
        mul3x3_3x8(invJ, J1 + gp*3*8, Ntemp);

        // Copy to output N
        for (int i = 0; i < 3*8; ++i) {
            N[gp*3*8 + i] = Ntemp[i];
        }
    }
}