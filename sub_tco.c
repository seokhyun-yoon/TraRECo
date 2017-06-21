#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *pseq1, *pseq2;
    register char *ps1, *ps2, *pst;
    size_t len1, len2;                   /* size of matrix */
    double Norm_dist_threshold;
    int Threshold;
    int pos = 0, dst = 10000, ref = 0;
    int K, dist, ss_ind;
    register int k, m;

    /* create a pointer to the real data in the input matrix  */
    pseq1 = (char*)mxGetPr(prhs[0]);
    pseq2 = (char*)mxGetPr(prhs[1]);
    len1 = mxGetN(prhs[0]);
    len2 = mxGetN(prhs[1]);
    ref = 0;
    if( len1 < len2 )
    {
        len2 = mxGetN(prhs[0]);
        len1 = mxGetN(prhs[1]);
        pseq2 = (char*)mxGetPr(prhs[0]);
        pseq1 = (char*)mxGetPr(prhs[1]);
        ref = 1;
    }
    Norm_dist_threshold = mxGetScalar(prhs[2]);
    ss_ind = (int)mxGetScalar(prhs[3]);
    Threshold = (int)( Norm_dist_threshold*(int)len2 + 0.5 );

    pos = 0;
    dst = 10000;
    K = (int)len1 - (int)len2 + 1;
    for( k = 0; k < K; k++ )
    {
        ps2 = pseq2;
        ps1 = pseq1 + k;
        for( m = 0, dist = 0; m<len2; m++ )
        {
            if( *ps1++ != *ps2++ )
            {
                dist ++;
                if( dist > Threshold ) break;
            }
        }
        if( dist <= Threshold )
        {
           pos = k + 1;
           dst = dist;
           break;
        }
    }
    if( pos == 0 && ss_ind == 0 )
    {
        for( k = 0; k < K; k++ )
        {
            ps2 = pseq2+len2-1;
            ps1 = pseq1 + k;
            for( m = 0, dist = 0; m<len2; m++ )
            {
                if( *ps1++ != (3 - (*ps2--)) )
                {
                    dist ++;
                    if( dist > Threshold ) break;
                }
            }
            if( dist <= Threshold )
            {
               pos = -(k + 1);
               dst = dist;
               break;
            }
        }
    }
    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar((double)pos);
    plhs[1] = mxCreateDoubleScalar((double)ref);
    plhs[2] = mxCreateDoubleScalar((double)dst);
}
