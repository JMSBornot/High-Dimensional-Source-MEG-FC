// compute_networkcluster_hashKdeep
#include <stdlib.h>
#include "mex.h"
#include <list>
#include <forward_list>
#include <iostream>

using namespace std;

int K = 5000*2000;
int Ndip;
list<uint32_T> *hash;
uint32_T * ptr_indsel;
uint16_T * ptr_ivert2irow;
uint16_T * ptr_ivert2icol;

int64_T remove_connection(uint16_T irow, uint16_T icol) {
    bool found = false;
    uint32_T iter, indconn;
    int loc = (ir*Ndip - (ir+1)*(ir+2)/2 + ic) % K;
    
    // if present, it has to be listed in hash[loc]
    for (auto ind = hash[loc].begin(); ind != hash[loc].end(); ind++) {
        iter = *ind;
        indconn = ptr_indsel[iter];
        if ((irow == ptr_ivert2irow[indconn]) && (icol == ptr_ivert2icol[indconn])) {
            found = true;
            hash[loc].erase(ind); // remove ind from the hash table
            break;
        }
    }
    
    if (found) return iter;
    else return -1;
}

/* ************************************************************************
                         Clustering Algorithm
 *
 * We are assuming that the indices for row (irow) and column (icol) as passed
 * in ivert2irow and ivert2icol always satisfy the condition icol > irow as
 * they correspond to an upper triangular matrix.
 *
***************************************************************************/
void compute_cluster(const mxArray *E, const mxArray *ivert2irow, const mxArray *ivert2icol,
        const mxArray *ivertNeighBeg, const mxArray *indsel, int Ndip, mxArray *cluster, mxArray *cluster_ext) {
    mwSize nE, nnzC;
    uint32_T iter, indconn;
    int ic, k, loc;
    int *ptr_E, *ptr_ivertNeighBeg, *ptr_cluster, *ptr_cluster_ext;
    uint16_T i1r, i1c, i2r, i2c, inddip;
    int64_T iter_next;
    
    nE = mxGetM(E);
    nnzC = mxGetN(indsel);
    
    ptr_E = (int *)mxGetData(E);
    ptr_ivert2irow = (uint16_T *)mxGetData(ivert2irow);
    ptr_ivert2icol = (uint16_T *)mxGetData(ivert2icol);
    ptr_ivertNeighBeg = (int *)mxGetData(ivertNeighBeg);
    ptr_indsel = (uint32_T *)mxGetData(indsel);
    ptr_cluster = (int *)mxGetData(cluster);
    ptr_cluster_ext = (int *)mxGetData(cluster_ext);
       
    /* --- Main --- */
    // Initialize the Hash Table (pointers to doubly linked list) and accompanying
    // fast localizer doubly linked list
    hash = new list<uint32_T>[Ndip];
    for (iter = 0; iter < nnzC; iter++) {
        indconn = ptr_indsel[iter];
        i1r = ptr_ivert2irow[indconn];
        i1c = ptr_ivert2icol[indconn];
        
        int loc = (ir*Ndip - (ir+1)*(ir+2)/2 + ic) % K;
        hash[loc].push_back(iter);
    }
    
    // initialize the double linked list
    bool * visited = new bool[nnzC];
    for (iter = 0; iter < nnzC; iter++) visited[iter] = false;
    
    // run along all the elements until there is none remaining
    ic = -1;
    forward_list<uint32_T> slink;
    for (uint16_T it = 0; it < nnzC; it++)
        if (visited[it]) continue;
        // remove this connection from the hash table
        indconn = ptr_indsel[it];
        i1r = ptr_ivert2irow[indconn];
        i1c = ptr_ivert2icol[indconn];
        remove_connection(i1r, i1c);
        // init a new cluster
        ic++;
        // extend this cluster to all neighboring links starting from link <iter>
        slink.push_front(it);
        while (!slink.empty()) {
            iter = slink.front(); // peek at the top
            slink.pop_front(); // detach node at the top
            indconn = ptr_indsel[iter];
            // assign connection <indconn> to current cluster
            ptr_cluster[iter] = ic;
            ptr_cluster_ext[ic]++;
            // read its corresponding row and column subindices
            i1r = ptr_ivert2irow[indconn];
            i1c = ptr_ivert2icol[indconn];
            
            //cout << "First half search:" << endl;
            
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
                
                //if (ic == 0) mexPrintf("(%d, %d) and (%d, %d)\n", i1r, i1c, i2r, i2c);
                
                // connections (i1r,i1c) and (i2r,i2c) are neighbours IF (i2r,i2c) is found in the hash table
                iter_next = remove_connection(i2r, i2c, hash, ptr_indsel, ptr_ivert2irow, ptr_ivert2icol);
                if (iter_next != -1) {
                    //mexPrintf("(%d, %d) and (%d, %d)\n", i1r, i1c, i2r, i2c);
                    //cout << "Removed: " << iter_next << endl;
                    visited[iter_next] = true;
                    slink.push_front(iter_next);
                }
            }
            
            //cout << "Second half search:" << endl;
            
            // extend this cluster by following neighbour connections that have in common the same index i1c
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
                
                //if (ic == 0) mexPrintf("(%d, %d) and (%d, %d)\n", i1r, i1c, i2r, i2c);
                
                // connections (i1r,i1c) and (i2r,i2c) are neighbours IF (i2r,i2c) is found in the hash table
                iter_next = remove_connection(i2r, i2c, hash, ptr_indsel, ptr_ivert2irow, ptr_ivert2icol);
                if (iter_next != -1) {
                    //mexPrintf("(%d, %d) and (%d, %d)\n", i1r, i1c, i2r, i2c);
                    //cout << "Removed: " << iter_next << endl;
                    visited[iter_next] = true;
                    slink.push_front(iter_next);
                }
            }
        }
    }
    // release memory
    delete[] visited;
    delete[] hash;
}

/* ************************************************************************
                         Gateway MEX Function
***************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize nE, nFC, Ndip, N;
        
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
    compute_cluster(prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], Ndip, plhs[0], plhs[1]);
}