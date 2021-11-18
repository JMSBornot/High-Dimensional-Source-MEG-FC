#include <math.h>
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

/**************************************************************************
                          Compute Ranking Order
***************************************************************************/
void compute_rankorder(const mxArray *xs, const mxArray *isort, mxArray * rank) {
    mxArray *up; // to keep position for unique elements and its number of repetitions
    mxArray *uc; // to keep number of repetitions (counter) for unique elements
    mxArray *ui; // to keep initial position for unique elements
    mxArray *tmp, *aux;
    double *ptr_xs, *ptr_rank, *ptr_tmp, *ptr_aux;
    int *ptr_up, *ptr_uc, *ptr_ui, *ptr_isort;
    int indcurr;
    mwSize M, N;
    
    M = mxGetM(xs);
    N = mxGetN(xs);
    
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
    for (int iter = 0; iter < N; iter++) {
        ptr_xs = mxGetPr(xs) + iter*M;
        ptr_rank = mxGetPr(rank) + iter*M;
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
            ptr_rank[*ptr_isort++] = *ptr_aux++;
        }
    }
    // release memory (mxArray)
    mxDestroyArray(up);
    mxDestroyArray(uc);
    mxDestroyArray(ui);
    mxDestroyArray(tmp);
    mxDestroyArray(aux);
}

/**************************************************************************
           Computational Routine Based on Correlation Statistics
***************************************************************************/
void compute_indices_correlation(const mxArray *X, const mxArray *y, const mxArray *iperm, const mxArray *Thr, mxArray * cellind) {
    double *ptrh, *ptlh, *ptrh2, *ptlh2, *ptrh3, *ptrh4, *ptr_Thr;
    double x, x2, tmp, xy, sum_y, sum_y2, r, cutval;
    int ind;
    int *ptr_iperm;
    mxArray *X2, *sum_x, *sum_x2, *ysurr, *corr, *item;
    mwSize M, N, Nr, Nth;
    Node * dlhead, *dltail, *dlink;
    uint64_T *ptr_item;
    
    M = mxGetM(X);
    N = mxGetN(X);
    Nr = mxGetN(iperm);
    Nth = mxGetNumberOfElements(Thr);
    
    // auxiliar arrays
    X2 = mxCreateDoubleMatrix(M, N, mxREAL);
    sum_x = mxCreateDoubleMatrix(1, N, mxREAL);
    sum_x2 = mxCreateDoubleMatrix(1, N, mxREAL);
    ysurr = mxCreateDoubleMatrix(M, 1, mxREAL);
    corr = mxCreateDoubleMatrix(1, N, mxREAL);
    
    // square X elements for efficient computations
    ptrh = mxGetPr(X);
    ptlh = mxGetPr(X2);
    for (int it = 0; it < M*N; it++) {
        x = *ptrh++;
        *ptlh++ = x*x;
    }
    // compute sum_y and sum_y2
    ptrh = mxGetPr(y);
    tmp = *ptrh++;
    sum_y = tmp;
    sum_y2 = tmp*tmp;
    for (int it = 1; it < M; it++) {
        tmp = *ptrh++;
        sum_y += tmp;
        sum_y2 += tmp*tmp;
    }    
    // compute sum_x and sum_x2 per column
    ptrh = mxGetPr(X);
    ptrh2 = mxGetPr(X2);
    ptlh = mxGetPr(sum_x);
    ptlh2 = mxGetPr(sum_x2);
    for (int icol = 0; icol < N; icol++) {
        x = *ptrh++;
        tmp = *ptrh2++;
        for (int irow = 1; irow < M; irow++) {
            x += *ptrh++;
            tmp += *ptrh2++;
        }
        *ptlh++ = x;
        *ptlh2++ = tmp;
    } 
    // Initialize double linked list
    start_list(&dlhead, &dltail);
    // pointer to threshold values
    ptr_Thr = mxGetPr(Thr);
    /* --- compute correlation stat for surrogate data --- */
    tmp = M*sum_y2 - sum_y*sum_y;
    for (int iter=0; iter < Nr; iter++) {
        // create surrogate data
        ptrh = mxGetPr(y);
        ptlh = mxGetPr(ysurr);
        ptr_iperm = (int *)mxGetData(iperm) + iter*M;
        for (int irow = 0; irow < M; irow++) {
            ind = *ptr_iperm++;
            *ptlh++ = ptrh[ind];
        }
        // compute correlation per column
        ptrh = mxGetPr(X);
        ptrh3 = mxGetPr(sum_x);
        ptrh4 = mxGetPr(sum_x2);
        ptlh = mxGetPr(corr);
        for (int icol = 0; icol < N; icol++) {
            ptrh2 = mxGetPr(ysurr);
            xy = (*ptrh++)*(*ptrh2++);
            for (int irow = 1; irow < M; irow++) {
                xy += (*ptrh++)*(*ptrh2++);
            }
            x = *ptrh3++;
            x2 = *ptrh4++;
            *ptlh++ = (M*xy - x*sum_y)/sqrt((M*x2 - x*x)*tmp);
        }
        // threshold correlation values and select indices for values in the tails
        for (int k = 0; k < Nth; k++) {
            ptrh = mxGetPr(corr);
            cutval = ptr_Thr[k];
            int cnt = 0;
            if (k % 2 == 0) // check for rank values in the lower tail
                for (uint64_T icol = 0; icol < N; icol++) {
                    r = *ptrh++;
                    if (r < cutval) {
                        insert_back(&dlhead, &dltail, icol);
                        cnt++;
                    }
                }
            else // check for rank values in the upper tail
                for (uint64_T icol = 0; icol < N; icol++) {
                    r = *ptrh++;
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
            mxSetCell(cellind, k*Nr+iter, mxDuplicateArray(item));
            mxDestroyArray(item);
        }
    }
    
    // release memory (mxArray)
    mxDestroyArray(X2);
    mxDestroyArray(sum_x);
    mxDestroyArray(sum_x2);
    mxDestroyArray(ysurr);
    mxDestroyArray(corr);
    destroy_list(&dlhead, &dltail);
}

/* ************************************************************************
                         Gateway MEX Function
***************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mxArray *rx, *ry; // rank of X and y
    mwSize M, N, Nr, Nth;
    int itype;
    int *tmp;
    
    /* check for proper number of arguments */
    if ((nrhs != 5) && (nrhs != 7)) {
        mexErrMsgIdAndTxt("ClustPermCorrelation:nrhs", "Five or seven inputs are required.");
    }
    if (nlhs>1) {
        mexErrMsgIdAndTxt("ClustPermCorrelation:nlhs", "No more than one output is required.");
    }
    /* make sure the 1st argument is an MxN matrix, type double */
    if (!mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("ClustPermCorrelation:notRealMatrix", "First input is an MxN matrix, type double.");
    }
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    /* make sure the 2nd argument is an Mx1 column vector, type double */
    if (!mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1]) ||
        (mxGetM(prhs[1]) != M) ||
        (mxGetN(prhs[1]) != 1)) {
        mexErrMsgIdAndTxt("ClustPermCorrelation:notRealMatrix", "Second input is an Mx1 column vector, type double.");
    }    
    /* make sure the 3rd argument is an MxNr matrix, type int32 */
    if ((mxGetClassID(prhs[2]) != mxINT32_CLASS) ||
         mxIsComplex(prhs[2]) ||
        (mxGetM(prhs[2]) != M)) {
        mexErrMsgIdAndTxt("ClustPermCorrelation:notInt32Matrix", "Third input is a MxNr matrix, type integer (int32).");
    }
    Nr = mxGetN(prhs[2]);
    /* make sure the 4th argument is a scalar, type int32 */
    if ( !(mxIsScalar(prhs[3])) ||
         mxIsComplex(prhs[3]) ||
        (mxGetClassID(prhs[3]) != mxINT32_CLASS) ) {
        mexErrMsgIdAndTxt("ClustPermCorrelation:notInt32Scalar", "Four input is a scalar, type integer (int32).");
    }
    tmp = (int *)mxGetData(prhs[3]);
    itype = tmp[0];
    /* make sure the 5th argument is a 2xNth matrix, type double */
    if (!mxIsDouble(prhs[4]) || 
         mxIsComplex(prhs[4]) ||
        (mxGetM(prhs[4]) != 2)) {
        mexErrMsgIdAndTxt("ClustPermCorrelation:notRealMatrix", "Fifth input is a 2xNth matrix, type double.");
    }
    Nth = mxGetN(prhs[4]);
    
    if (nrhs == 7) {
        /* make sure the 6th argument is an MxN matrix, type int32 */
        if ((mxGetClassID(prhs[5]) != mxINT32_CLASS) ||
                mxIsComplex(prhs[5]) ||
                (mxGetM(prhs[5]) != M) ||
                (mxGetN(prhs[5]) != N)) {
            mexErrMsgIdAndTxt("ClustPermCorrelation:notInt32Matrix", "Sixth input is a MxN matrix, type integer (int32).");
        }
        /* make sure the 7th argument is an Mx1 column vector, type int32 */
        if ((mxGetClassID(prhs[6]) != mxINT32_CLASS) ||
                mxIsComplex(prhs[6]) ||
                (mxGetM(prhs[6]) != M) ||
                (mxGetN(prhs[6]) != 1)) {
            mexErrMsgIdAndTxt("ClustPermCorrelation:notInt32Matrix", "Seventh input is a Mx1 column vector, type integer (int32).");
        }
    }
        
    /* create the unique output argument */
    plhs[0] = mxCreateCellMatrix(Nr,2*Nth);
    
    /* call the computational routine */
    switch (itype) {
        case 1: // classical correlation analysis
            // compute_indices_correlation(const mxArray *X, const mxArray *y, const mxArray *iperm, const mxArray *Thr, mwSize mA, mxArray * cellind)
            compute_indices_correlation(prhs[0], prhs[1], prhs[2], prhs[4], plhs[0]);
            break;
        case 2: // Spearman's correlation analysis
            // compute rank of X and y
            rx = mxCreateDoubleMatrix(M, N, mxREAL);
            ry = mxCreateDoubleMatrix(M, 1, mxREAL);
            compute_rankorder(prhs[0], prhs[5], rx);
            compute_rankorder(prhs[1], prhs[6], ry);
            // compute indices
            compute_indices_correlation(rx, ry, prhs[2], prhs[4], plhs[0]);
            // release memory (mxArray)
            mxDestroyArray(rx);
            mxDestroyArray(ry);
            break;
        default:
            mexErrMsgIdAndTxt("ClustPermCorrelation:unexpected", "This option is not implemented. Please check the help.");
            break;
    }
}