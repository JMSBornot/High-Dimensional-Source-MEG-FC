#include "mex.h"

/* ************************************************************************
                          Data Structures 
***************************************************************************/
struct Node {
    uint64_T ind;
    Node *prev;
    Node *next;
};
/* Called at the begginning for starting the list */
void start_list(Node **head, Node **tail) {
    *head = *tail = NULL;
}
/* insert at top */
void insert_front(Node **head, Node **tail, uint64_T ind) {
    Node *p = new Node;
    p->ind = ind;
    p->next = *head;
    p->prev = NULL;
    *head = p;
    if (*tail == NULL) *tail = p;
    if (p->next != NULL) p->next->prev = p;
}
/* insert at end */
void insert_back(Node **head, Node **tail, uint64_T ind) {
    Node *p = new Node;
    p->ind = ind;
    p->next = NULL;
    p->prev = *tail;
    *tail = p;
    if (*head == NULL) *head = p;
    if (p->prev != NULL) p->prev->next = p;
}
/* detach a node */
void detach_node(Node **head, Node **tail, Node *tmp) {
    Node *p = tmp;
    tmp = tmp->next;
    if (p->prev == NULL)
        *head = tmp;
    else p->prev->next = tmp;
    if (tmp == NULL)
        *tail = p->prev;
    else tmp->prev = p->prev;
    p->prev = NULL;
    p->next = NULL;
    delete(p);
}
/* Destroy and free all memory resources kept by the list */
void destroy_list(Node **head, Node **tail) {
    Node *p, *pnext;
    int cnt = 0;
    for (p = *head; p != NULL; p = pnext) {
        cnt++;
        pnext = p->next;
        delete(p);
    }
    *head = *tail = NULL; // prevent dangling pointer
    mexPrintf("Removed %d elements\n", cnt);
}

/* ************************************************************************
           Computational Routine Based on Ranksum Statistics
***************************************************************************/
void compute_indices_ranksum(const mxArray *xs, const mxArray *isort, const mxArray *iperm, const mxArray *Thr, mwSize mA, mxArray * cellind) {    
    mxArray *rank; // to compute the rank for each permutation resampling
    mxArray *up; // to keep position for unique elements and its number of repetitions
    mxArray *uc; // to keep number of repetitions (counter) for unique elements
    mxArray *ui; // to keep initial position for unique elements
    mxArray *tmp, *aux, *item;
    Node * dlhead, *dltail, *dlink;
    double *ptr_xs, *ptr_rank, *ptr_Thr, *ptr_tmp, *ptr_aux;
    int *ptr_up, *ptr_uc, *ptr_ui, *ptr_isort, *ptr_iperm;
    int indcurr;
    double r, cutval;
    uint64_T *ptr_item;
    mwSize M, N, Nr, Nth;
    
    M = mxGetM(xs);
    N = mxGetN(xs);
    Nr = mxGetN(iperm);
    Nth = mxGetNumberOfElements(Thr);
    
    start_list(&dlhead, &dltail);
    
    /* create 1-by-N arrays of int-32 integers */
    mwSize dims[] = {1,M};
    up = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    uc = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    ui = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    ptr_up = (int *)mxGetData(up);
    ptr_uc = (int *)mxGetData(uc);
    ptr_ui = (int *)mxGetData(ui);
    ptr_ui[0] = 1;
    tmp = mxCreateDoubleMatrix(1, M, mxREAL);
    aux = mxCreateDoubleMatrix(1, M, mxREAL);
    rank = mxCreateDoubleMatrix(1, N, mxREAL);
    ptr_Thr = mxGetPr(Thr);
    for (int ir = 0; ir < Nr; ir++) {
        ptr_rank = mxGetPr(rank);
        for (int iter = 0; iter < N; iter++) {
            ptr_xs = mxGetPr(xs) + iter*M;
            indcurr = 0;
            for (int it = 0; it < M; it++) {
                ptr_uc[it] = 0; // clean counter array
            }
            for (int it = 1; it < M; it++) {
                if (ptr_xs[it] == ptr_xs[it-1]) {
                    ptr_uc[indcurr] = ptr_uc[indcurr] + 1;
                    ptr_ui[it] = ptr_ui[it-1];
                }
                else {
                    indcurr = indcurr + 1;
                    ptr_ui[it] = it+1;
                }
                ptr_up[it] = indcurr;
            }
            // MATLAB: tmp = 0.5*uc(1:indcurr);
            ptr_tmp = mxGetPr(tmp);
            for (int it = 0; it <= indcurr; it++) {
                *ptr_tmp++ = 0.5*ptr_uc[it];
            }
            // MATLAB: rs = tmp(up) + ui; (compute the ranks for the sorted array)
            ptr_tmp = mxGetPr(tmp);
            ptr_aux = mxGetPr(aux);
            for (int it = 0; it < M; it++) {
                *ptr_aux++ = ptr_tmp[ptr_up[it]] + ptr_ui[it];
            }
            // MATLAB: R(isort) = Rs; (set the ranks for the original array positions)
            ptr_aux = mxGetPr(aux);
            ptr_isort = (int *)mxGetData(isort) + iter*M;
            for (int it = 0; it < M; it++) {
                ptr_tmp[*ptr_isort++] = *ptr_aux++;
            }
            // compute the ranksum stats
            r = 0;
            ptr_iperm = (int *)mxGetData(iperm) + ir*M;
            for (int it = 0; it < mA; it++) {
                r += ptr_tmp[*ptr_iperm++];
            }
            *ptr_rank++ = r;
        }
        for (int k = 0; k < Nth; k++) {
            ptr_rank = mxGetPr(rank);
            cutval = ptr_Thr[k];
            int cnt = 0;
            if (k % 2 == 0) // check for rank values in the lower tail
                for (uint64_T icol = 0; icol < N; icol++) {
                    r = *ptr_rank++;
                    if (r < cutval) {
                        insert_back(&dlhead, &dltail, icol);
                        cnt++;
                    }
                }
            else // check for rank values in the upper tail
                for (uint64_T icol = 0; icol < N; icol++) {
                    r = *ptr_rank++;
                    if (r > cutval) {
                        insert_back(&dlhead, &dltail, icol);
                        cnt++;
                    }
                }
            // Create an auxiliar mxArray and fill up with the elements in the linked list
            mwSize dims[] = {1,cnt};
            item = mxCreateNumericArray(2, dims, mxUINT64_CLASS, mxREAL);
            ptr_item = (uint64_T *)mxGetData(item);
            for (dlink = dlhead; dlink != NULL; dlink = dlhead) {
                *ptr_item++ = dlink->ind; // peek at the top
                detach_node(&dlhead, &dltail, dlink); // detach node at the top
            }
            mxSetCell(cellind, k*Nr+ir, mxDuplicateArray(item));
            mxDestroyArray(item);
        }
    }

    // release memory (mxArray)
    mxDestroyArray(up);
    mxDestroyArray(uc);
    mxDestroyArray(ui);
    mxDestroyArray(tmp);
    mxDestroyArray(aux);
    mxDestroyArray(rank);
    destroy_list(&dlhead, &dltail);
}

/* ************************************************************************
                         Gateway MEX Function
***************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize M, N, mA, Nr, Nth;
    mwSize *tmp;
    
    /* check for proper number of arguments */
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("Ranksum:nrhs", "Five inputs are required.");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("Ranksum:nlhs", "No more than one output is required.");
    }
    /* make sure the 1st argument is an MxN matrix, type double */
    if (!mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("Ranksum:notRealMatrix", "First input is an MxN matrix, type double.");
    }
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    /* make sure the 2nd argument is an MxN matrix, type int32 */
    if ((mxGetClassID(prhs[1]) != mxINT32_CLASS) ||
         mxIsComplex(prhs[1]) ||
        (mxGetM(prhs[1]) != M) ||
        (mxGetN(prhs[1]) != N)) {
        mexErrMsgIdAndTxt("Ranksum:notInt32Matrix", "Second input is a MxN matrix, type integer (int32).");
    }
    
    /* make sure the 3rd argument is an MxNr matrix, type int32 */
    if ((mxGetClassID(prhs[2]) != mxINT32_CLASS) ||
         mxIsComplex(prhs[2]) ||
        (mxGetM(prhs[2]) != M)) {
        mexErrMsgIdAndTxt("Ranksum:notInt32Matrix", "Third input is a MxNr matrix, type integer (int32).");
    }
    Nr = mxGetN(prhs[2]);
    /* make sure the 4th argument is a scalar, type int32 */
    if ( !(mxIsScalar(prhs[3])) ||
         mxIsComplex(prhs[3]) ||
        (mxGetClassID(prhs[3]) != mxINT32_CLASS) ) {
        mexErrMsgIdAndTxt("Ranksum:notInt32Scalar", "Four input is a scalar, type integer (int32).");
    }
    tmp = (mwSize *)mxGetData(prhs[3]);
    mA = tmp[0];
    /* make sure the 5th argument is a 2xNth matrix, type double */
    if (!mxIsDouble(prhs[4]) || 
         mxIsComplex(prhs[4]) ||
        (mxGetM(prhs[4]) != 2)) {
        mexErrMsgIdAndTxt("Ranksum:notRealMatrix", "Fifth input is a 2xNth matrix, type double.");
    }
    Nth = mxGetN(prhs[4]);
        
    /* create the unique output argument */
    plhs[0] = mxCreateCellMatrix(Nr,2*Nth);
    
    /* call the computational routine */
    compute_indices_ranksum(prhs[0], prhs[1], prhs[2], prhs[4], mA, plhs[0]);
}