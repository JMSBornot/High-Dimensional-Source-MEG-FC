#include <stdlib.h>
#include "mex.h"
#include <list>
#include <iterator>

using namespace std;

/* ************************************************************************
                         Clustering Algorithm
***************************************************************************/
void compute_cluster(const mxArray *E, const mxArray *ivert2irow, const mxArray *ivert2icol,
        const mxArray *ivertNeighBeg, const mxArray *indsel, mxArray *cluster, mxArray *cluster_ext) {
    mwSize nE, nnzE;
    uint32_T it, it1, i1, i2;
    uint32_T *ptr_indsel;
    int idip1, idip2, itmp, ic, k, cnt;
    int *ptr_E, *ptr_ivertNeighBeg, *ptr_cluster, *ptr_cluster_ext;
    uint16_T *ptr_ivert2irow, *ptr_ivert2icol;
    uint16_T i1r, i1c, i2r, i2c;
    
    nE = mxGetM(E);
    nnzE = mxGetN(indsel);
    
    list<uint32_T> dlink, slink;
    
    ptr_E = (int *)mxGetData(E);
    ptr_ivert2irow = (uint16_T *)mxGetData(ivert2irow);
    ptr_ivert2icol = (uint16_T *)mxGetData(ivert2icol);
    ptr_ivertNeighBeg = (int *)mxGetData(ivertNeighBeg);
    ptr_indsel = (uint32_T *)mxGetData(indsel);
    ptr_cluster = (int *)mxGetData(cluster);
    ptr_cluster_ext = (int *)mxGetData(cluster_ext);
       
    /* --- Main --- */
    // initialize the double linked list
    for (uint32_T iter = 0; iter < nnzE; iter++) dlink.push_back(iter);
    // run along dlink elements until they are all assigned to any cluster
    ic = -1;
    while (!dlink.empty()) {
        ic++;
        it = dlink.front(); // peek at the top
        dlink.pop_front(); // detach node at the top
        ptr_cluster[it] = ic;
        ptr_cluster_ext[ic] = 1;
        // extend this cluster to all neighboring links starting from link <it>
        slink.push_back(it);
        while (!slink.empty()) {
            it1 = slink.front(); // peek at the top
            slink.pop_front(); // detach node at the top
            i1 = ptr_indsel[it1];
            i1r = ptr_ivert2irow[i1];
            i1c = ptr_ivert2icol[i1];
            
            // visit all the elements of dlink as looking for neighbor links
            for (auto it2 = dlink.begin(); it2 != dlink.end(); ) {
                i2 = ptr_indsel[*it2];
                i2r = ptr_ivert2irow[i2];
                i2c = ptr_ivert2icol[i2];
                // These edges are neighbors if they share one (tied) point while
                // the other (loose) points are neighbors in the cortical surface
                // or in the brain volume.
                idip1 = idip2 = -1;
                if (i1r == i2r) {
                    idip1 = i1c; idip2 = i2c;
                }
                else if (i1r == i2c) {
                    idip1 = i1c; idip2 = i2r;
                }
                else if (i1c == i2r) {
                    idip1 = i1r; idip2 = i2c;
                }
                else if (i1c == i2c) {
                    idip1 = i1r; idip2 = i2r;
                }
                bool flagNeigh = false; // checking if the loose points are neighbors
                if (idip1 != -1) {
                    // ensure that idip1 is the greatest index
                    if (idip2 > idip1) {
                        itmp = idip1;
                        idip1 = idip2;
                        idip2 = itmp;
                    }
                    k = ptr_ivertNeighBeg[idip1];
                    while ((k<nE) && (ptr_E[k] == idip1)) {
                        if (ptr_E[nE+k] == idip2) { // sharing an Edge in the surface?
                            flagNeigh = true;
                            break;
                        }
                        k++;
                    }
                }
                if (flagNeigh) { // assign it2 to the current (expanding) cluster
                    ptr_cluster[*it2] = ic;
                    ptr_cluster_ext[ic]++;
                    slink.push_back(*it2); // adding the node at the end of the single linked list
                    it2 = dlink.erase(it2); // remove this node from the double linked list
                }
                else ++it2;
            }
        }
    }
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
    /* make sure the 2nd argument is a vector with nFC elements, type int32 */
    nFC = mxGetNumberOfElements(prhs[1]);
    if ((mxGetClassID(prhs[1]) != mxUINT16_CLASS) ||
         mxIsComplex(prhs[1]) ||
        (mxGetM(prhs[1]) != nFC) ||
        (mxGetN(prhs[1]) != 1)) {
        mexErrMsgIdAndTxt("ClusterLinks:notUInt16Matrix", "Second input is a vector with <nFC> elements, type integer (uint16).");
    }
    /* make sure the 3rd argument is a vector with nFC elements, type int32 */
    if ((mxGetClassID(prhs[2]) != mxUINT16_CLASS) ||
         mxIsComplex(prhs[2]) ||
        (mxGetM(prhs[2]) != nFC) ||
        (mxGetN(prhs[2]) != 1)) {
        mexErrMsgIdAndTxt("LinkedList:notUInt16Matrix", "Third input is a vector with <nFC> elements, type integer (uint16).");
    }
    /* make sure the 4th argument is a vector, type int32 */
    Ndip = mxGetNumberOfElements(prhs[3]);
    if ((mxGetClassID(prhs[3]) != mxINT32_CLASS) ||
         mxIsComplex(prhs[3]) ||
        (mxGetM(prhs[3]) != 1) ||
        (mxGetN(prhs[3]) != Ndip)) {
        mexErrMsgIdAndTxt("LinkedList:notInt32Matrix", "Fourth input is a vector with <Ndip> elements, type integer (int32).");
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