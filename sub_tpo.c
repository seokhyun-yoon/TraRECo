#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *pseq1, *pseq2;
    register char *ps1, *ps2, *pst;
    size_t len1, len2;                   /* size of matrix */
    double *p_param;
    static double Norm_dist_threshold = 0;
    static int Threshold[1000], Max_overlap = 0, Min_overlap = 0;
    int min_overlap, max_overlap;
    double norm_dist_threshold;
    int pos = 0, dst = 10000, dep = 0;
    int K, dist, min_len, len, ss_ind;
    register int k, m, th;

    /* create a pointer to the real data in the input matrix  */
    pseq1 = (char*)mxGetPr(prhs[0]);
    pseq2 = (char*)mxGetPr(prhs[1]);
    len1 = mxGetN(prhs[0]);
    len2 = mxGetN(prhs[1]);
    p_param = mxGetPr(prhs[2]);
    norm_dist_threshold = p_param[0];
    min_overlap = (int)p_param[1];
    max_overlap = (int)p_param[2];
    ss_ind = (int)p_param[3];
    
    if( len1 > len2 ) min_len = (int)len2;
    else              min_len = (int)len1;
    if( max_overlap > min_len ) max_overlap = min_len;
    
    K = (int)max_overlap - (int)min_overlap + 1;
    if( Max_overlap < max_overlap || Min_overlap != min_overlap || Norm_dist_threshold != norm_dist_threshold  )
    {
        Norm_dist_threshold = norm_dist_threshold;
        Max_overlap = max_overlap;
        Min_overlap = min_overlap;
        for( k = 0; k < K; k++ ) Threshold[k] = (int)( (min_overlap+k)*norm_dist_threshold + 0.5 );
    }
    
    if( max_overlap >= min_overlap )
    {
        K = (int)max_overlap - (int)min_overlap + 1;
        for( k = 0; k < K; k++ )
        {
            ps1 = pseq1;
            ps2 = pseq2 + (int)(len2-min_overlap-k);
            len = (int)(min_overlap+k);
            dist = 0;
            th = Threshold[k];
            for( m = 0; m<len; m++ )
            {
                if( *ps1++ != *ps2++ )
                {
                    dist ++;
                    if( dist > th ) break;
                }
            }
            if( dist <= th )
            {
               dep = min_overlap+k;
               pos = len2 - dep;
               dst = dist;
               break;
            }
        }
        if( pos == 0 )
        {
            for( k = 0; k < K; k++ )
            {
                ps1 = pseq1 + (int)(len1-min_overlap-k);
                ps2 = pseq2;
                len = (int)(min_overlap+k);
                dist = 0;
                th = Threshold[k];
                for( m = 0; m<len; m++ )
                {
                    if( *ps1++ != *ps2++ )
                    {
                        dist ++;
                        if( dist > th ) break;
                    }
                }
                if( dist <= th )
                {
                   dep = -(min_overlap+k);
                   pos = len2 + dep;
                   dst = dist;
                   break;
                }
            }
        }
        if( pos == 0 && ss_ind == 0 )
        {
            pst = (char*)mxMalloc(len2);
            ps1 = pst;
            ps2 = pseq2+len2-1;
            for( m = 0; m<len2; m++ ) *ps1++ = 3-(*ps2--); // s2 = 3 - s2(end:-1:1);

            //ps2 = pseq2;
            for( k = 0; k < K; k++ )
            {
                ps1 = pseq1;
                ps2 = pst + (int)(len2-min_overlap-k);
                len = (int)(min_overlap+k);
                dist = 0;
                th = Threshold[k];
                for( m = 0; m<len; m++ )
                {
                    if( *ps1++ != *ps2++ )
                    {
                        dist ++;
                        if( dist > th ) break;
                    }
                }
                if( dist <= th )
                {
                   dep = (min_overlap+k);
                   pos = -(len2 - dep);
                   dst = dist;
                   break;
                }
            }
            if( pos != 0 )
            {
                mxFree(pst);
            }
            else
            {
                for( k = 0; k < K; k++ )
                {
                    ps1 = pseq1 + (int)(len1-min_overlap-k);
                    ps2 = pst;
                    len = (int)(min_overlap+k);
                    dist = 0;
                    th = Threshold[k];
                    for( m = 0; m<len; m++ )
                    {
                        if( *ps1++ != *ps2++ )
                        {
                            dist ++;
                            if( dist > th ) break;
                        }
                    }
                    if( dist <= th )
                    {
                       dep = -(min_overlap+k);
                       pos = -(len2 + dep);
                       dst = dist;
                       break;
                    }
                }
                mxFree(pst);
            }
        }
    }
    
    plhs[0] = mxCreateDoubleScalar((double)pos);
    plhs[1] = mxCreateDoubleScalar((double)dep);
    plhs[2] = mxCreateDoubleScalar((double)dst);
}
