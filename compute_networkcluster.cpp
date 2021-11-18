// compute_networkcluster.cpp
#include <stdlib.h>
#include "mex.h"
#include <forward_list>
#include <iostream>

using namespace std;

/* ************************************************************************
                         Clustering Algorithm
 *
 * We are assuming that the indices for row (irow) and column (icol) as passed
 * in ivert2irow and ivert2icol always satisfy the condition irow < icol  as
 * they correspond to an upper triangular matrix.
 *
 * NOTE:
 * ivert2irow and ivert2icol are considered arrays of uint16_T because the
 *      number of dipoles we have to dealt with in this program is about
 *      16403, under 65535 which is the max int of the type uint16.
 *4294967295
 *
***************************************************************************/
void compute_cluster(const mxArray *E, const mxArray *ivert2irow, const mxArray *ivert2icol,
        const mxArray *ivertNeighBeg, const mxArray *indsel, mxArray *cluster, mxArray *cluster_ext) {
    int ic, k;
    int *ptr_E, *ptr_ivertNeighBeg, *ptr_cluster, *ptr_cluster_ext;
    uint16_T * ptr_ivert2irow, *ptr_ivert2icol;
    uint32_T * ptr_indsel;
    uint16_T i1r, i1c, i2r, i2c, inddip;
    uint32_T iter, indconn, indconn_block, nnzC, nE;
    uint16_T Ndip;
    uint32_T nFC;
    
    Ndip = mxGetNumberOfElements(ivertNeighBeg);
    nFC = mxGetNumberOfElements(ivert2irow);
    
    int * exist_connection = new int[nFC];
    for (uint32_T it = 0; it < nFC; it++) exist_connection[it] = -1;
    
    nE = (uint32_T)mxGetM(E);
    nnzC = (uint32_T)mxGetN(indsel);
    
    ptr_E = (int *)mxGetData(E);
    ptr_ivert2irow = (uint16_T *)mxGetData(ivert2irow);
    ptr_ivert2icol = (uint16_T *)mxGetData(ivert2icol);
    ptr_ivertNeighBeg = (int *)mxGetData(ivertNeighBeg);
    ptr_indsel = (uint32_T *)mxGetData(indsel);
    ptr_cluster = (int *)mxGetData(cluster);
    ptr_cluster_ext = (int *)mxGetData(cluster_ext);
       
    /* --- Main --- */
    // mark the connections in <exist_connection> array
    for (iter = 0; iter < nnzC; iter++) {
        indconn_block = ptr_indsel[iter]; // this indices are following the index order according to the blockwise FC analysis
        i1r = ptr_ivert2irow[indconn_block];
        i1c = ptr_ivert2icol[indconn_block];
        indconn = i1r*Ndip - (i1r+1)*(i1r+2)/2 + i1c; // recodify for the FC indices correspoding to the upper triangular matrix (column-first order)
        exist_connection[indconn] = iter;
    }
    
    // run along all the connections until there is none unassigned to any cluster
    ic = -1;
    forward_list<uint32_T> slink;
    for (uint32_T it = 0; it < nnzC; it++) {
        indconn_block = ptr_indsel[it];
        i1r = ptr_ivert2irow[indconn_block];
        i1c = ptr_ivert2icol[indconn_block];
        indconn = i1r*Ndip - (i1r+1)*(i1r+2)/2 + i1c;
        if (exist_connection[indconn] < 0) continue;
        exist_connection[indconn] = -1; // remove this connection
        // init a new cluster
        ic++;
        // extend this cluster to all neighbour links starting from link <it>
        slink.push_front(it);
        while (!slink.empty()) {
            // peek at the top and remove the node
            iter = slink.front();
            slink.pop_front();
            // assign connection <iter> to current cluster
            ptr_cluster[iter] = ic;
            ptr_cluster_ext[ic]++;
            // read corresponding indices
            indconn_block = ptr_indsel[iter];
            i1r = ptr_ivert2irow[indconn_block];
            i1c = ptr_ivert2icol[indconn_block];
            // extend this cluster by following neighbour connections that have in common the same index i1r
            k = ptr_ivertNeighBeg[i1c];
            while ((k < nE) && (ptr_E[k] == i1c)) { // the dissimilar indices have to be neighbours in the cortical surface
                inddip = ptr_E[nE+k];
                k++;
                if (inddip == i1r) continue;
                // enforce that i2r < i2c
                if (inddip < i1r) {
                    i2r = inddip;
                    i2c = i1r;
                }
                else {
                    i2r = i1r;
                    i2c = inddip;
                }
                // connections (i1r,i1c) and (i2r,i2c) are neighbours IF (i2r,i2c) is found in exist_connection
                indconn = i2r*Ndip - (i2r+1)*(i2r+2)/2 + i2c;                
                if (exist_connection[indconn] > -1) {
                    slink.push_front(exist_connection[indconn]);
                    exist_connection[indconn] = -1; // remove this connection
                }
            }            
            // continue extending this cluster by following neighbour connections that have in common the same index i1c
            k = ptr_ivertNeighBeg[i1r];
            while ((k < nE) && (ptr_E[k] == i1r)) { // the dissimilar indices have to be neighbours in the cortical surface
                inddip = ptr_E[nE+k];
                k++;
                if (inddip == i1c) continue;
                // enforce that i2r < i2c
                if (inddip < i1c) {
                    i2r = inddip;
                    i2c = i1c;
                }
                else {
                    i2r = i1c;
                    i2c = inddip;
                }
                // connections (i1r,i1c) and (i2r,i2c) are neighbours IF (i2r,i2c) is found in exist_connection
                indconn = i2r*Ndip - (i2r+1)*(i2r+2)/2 + i2c;
                if (exist_connection[indconn] > -1) {
                    slink.push_front(exist_connection[indconn]);
                    exist_connection[indconn] = -1; // remove this connection
                }
            }
        }
    }
    // release memory
    delete[] exist_connection;
}

/* ************************************************************************
                         Gateway MEX Function
***************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize nE, N;
    uint64_T Ndip, nFC;
        
    /* check for proper number of arguments */
    if (nrhs != 5) mexErrMsgIdAndTxt("ClusterLinks:nrhs", "Five inputs are required.");
    if (nlhs>2) mexErrMsgIdAndTxt("ClusterLinks:nlhs", "No more than two outputs.");
    /* make sure the 1st argument is a nEx2 matrix, type int32 */
    if ((mxGetClassID(prhs[0]) != mxINT32_CLASS) ||
         mxIsComplex(prhs[0]) ||
        (mxGetN(prhs[0]) != 2)) {
        mexErrMsgIdAndTxt("ClusterLinks:notInt32Matrix", "First input is a nEx2 matrix, type integer (int32).");
    }
    nE = mxGetM(prhs[0]);
    /* make sure the 2nd argument is a vector with nFC elements, type uint16 */
    nFC = mxGetNumberOfElements(prhs[1]);
    if (nFC > 0.5*65535*(65535-1))
        mexErrMsgIdAndTxt("LinkedList:UInt32_overflow", "Number of connection greater than 4294967295, causing an uint32 overflow");
    if ((mxGetClassID(prhs[1]) != mxUINT16_CLASS) ||
         mxIsComplex(prhs[1]) ||
        (mxGetM(prhs[1]) != nFC) ||
        (mxGetN(prhs[1]) != 1)) {
        mexErrMsgIdAndTxt("ClusterLinks:notUInt16Matrix", "Second input is a vector with <nFC> elements, type integer (uint16).");
    }
    /* make sure the 3rd argument is a vector with nFC elements, type uint16 */
    if ((mxGetClassID(prhs[2]) != mxUINT16_CLASS) ||
         mxIsComplex(prhs[2]) ||
        (mxGetM(prhs[2]) != nFC) ||
        (mxGetN(prhs[2]) != 1)) {
        mexErrMsgIdAndTxt("LinkedList:notUInt16Matrix", "Third input is a column vector with <nFC> elements, type integer (uint16).");
    }
    /* make sure the 4th argument is a vector, type int32 */
    Ndip = mxGetNumberOfElements(prhs[3]);
    if (Ndip > 65535)
        mexErrMsgIdAndTxt("LinkedList:UInt16_overflow", "Number of dipoles greater than 65535, causing an uint16 overflow");
    if ((mxGetClassID(prhs[3]) != mxINT32_CLASS) ||
         mxIsComplex(prhs[3]) ||
        (mxGetM(prhs[3]) != 1) ||
        (mxGetN(prhs[3]) != Ndip)) {
        mexErrMsgIdAndTxt("LinkedList:notInt32Matrix", "Fourth input is a column vector with <Ndip> elements, type integer (int32).");
    }
    /* make sure the 5th argument is a vector, type uint32 */
    if ((mxGetClassID(prhs[4]) != mxUINT32_CLASS) ||
         mxIsComplex(prhs[4]) ||
       ((mxGetM(prhs[4]) != 1) &&
        (mxGetN(prhs[4]) != 1))) {
        mexErrMsgIdAndTxt("ClusterLinks:notUInt32Matrix", "Fifth input is a vector, type integer (uint32).");
    }
    N = mxGetN(prhs[4]);
    
    /* create output arguments */
    mwSize dims[] = {1,N};
    plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    
    /* call the computational routine */
    compute_cluster(prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], plhs[0], plhs[1]);
}