#include "mex.h"

void find_max_run_length( char* pseq1, char *pseq2, int len1, int *dst_threshold_len, int flag, int *run, int *dst )
{
    register int last_0_pos = 0;  
    register int last_0_dst = 10000;  
    register int dst_acc = 0;
    register int k;
    register char *seq1, *seq2;
    register int *pdst;
    pdst = dst_threshold_len;
    seq1 = pseq1;
    seq2 = pseq2;
    if( flag == 0 )
    {
        for(k=0;k<len1;k++)
        {
            if( *seq1++ == *seq2++ )
            {
                last_0_pos = k;
                last_0_dst = dst_acc;
            }
            else
            {
                if( ++dst_acc > *pdst ) break;
            }
            pdst ++;
        }
    }
    else
    {
        for(k=0;k<len1;k++)
        {
            if( *seq1-- == *seq2-- )
            {
                last_0_pos = k;
                last_0_dst = dst_acc;
            }
            else
            {
                if( ++dst_acc > *pdst ) break;
            }
            pdst ++;
        }
    }
    *run = last_0_pos;
    *dst = last_0_dst;
}

int find_max( int *run, int len )
{
    register int k, tmp, idx = 0, max = *run++;
    for(k=1;k<len;k++)
    {
        if( max < (tmp = *run++) ){ max = tmp; idx = k; } 
    }
    return idx;
}

void set_zero( int *run, int len )
{
    register int k;
    for(k=0;k<len;k++) *run++ = 0;
}
void set_val( int *run, int len, int val )
{
    register int k;
    for(k=0;k<len;k++) *run++ = val;
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *pseq1, *pseq2, *psdvt;
    static char psdv1[1000], psdv2[1000];
    register char *ps1, *ps2, *pst;
    size_t len1, len2;                   /* size of matrix */
    double *p_param;
    static double Norm_dist_threshold = 0;
    static int Threshold[1000], Max_overlap = 0, Min_overlap = 0;
    int min_overlap, max_overlap;
    double norm_dist_threshold;
    int j_mode = 0, dst = 10000, pos = 0, dep = 0, n_junctions = 0;
    static int n_run[1000], dist[1000];
    int *pitmp;
    int K, min_len, len, sel, run_threshold, lent, idx, ss_ind;
    register int k, m, ln;

    /* create a pointer to the real data in the input matrix  */
    pseq1 = (char*)mxGetPr(prhs[0]);
    pseq2 = (char*)mxGetPr(prhs[1]);
    len1 = mxGetN(prhs[0]);
    len2 = mxGetN(prhs[1]);
    p_param = (double*)mxGetPr(prhs[2]);
    norm_dist_threshold = p_param[0];
    min_overlap = (int)p_param[3];
    max_overlap = (int)p_param[2];
    ss_ind = (int)p_param[4];
    sel = (int)mxGetScalar(prhs[3]);
    run_threshold = (int)mxGetScalar(prhs[4]);
    
    if( len1 > len2 ) min_len = (int)len2;
    else              min_len = (int)len1;
    if( max_overlap < min_len ) len = max_overlap;
    else len = min_len;
    
    K = (int)max_overlap; //  - (int)min_overlap + 1;
    if( Max_overlap < max_overlap || Min_overlap != min_overlap || Norm_dist_threshold != norm_dist_threshold  )
    {
        Norm_dist_threshold = norm_dist_threshold;
        Max_overlap = max_overlap;
        Min_overlap = min_overlap;
        for( k = 0; k < K; k++ )
        {
            if( k < min_overlap ) Threshold[k] = (int)( (min_overlap)*norm_dist_threshold + 0.5 );
            else  Threshold[k] = (int)( k*norm_dist_threshold + 0.5 );
        }
    }

    if(sel == 0)
    {
        for(k=0;k<len;k++)
        {
            psdv1[k] = pseq2[k];
            psdv2[k] = 3 - pseq2[len-k-1];
        }
    }
    else
    {
        for(k=0;k<len;k++)
        {
            psdv1[k] = 3 - pseq2[len2-1-k];
            psdv2[k] = pseq2[len2-len+k];
        }
    }
    
    if( len1 >= len )
    {
        if( sel == 0 || ss_ind == 0 )
        {
            K = (int)len1 - (int)len + 1;
            set_zero( n_run, K );
            set_val( dist, K, 10000 );
            // n = 1;
            ps1 = pseq1;
            ps2 = psdv1;
            for( k = 0; k < K; k++, ps1 ++ )
            {
                find_max_run_length( ps1, ps2, len, Threshold, 0, n_run+k, dist+k );            
                if( n_run[k] > run_threshold ) break;
            }

            idx = find_max( n_run, K );
            lent = n_run[idx];
            if( lent >= min_overlap )
            {
                pos = idx + 1;
                dep = lent;
                if( sel == 0 ) j_mode = 1; // n*2-1
                else           j_mode = 4; // (3-n)*2
                dst = dist[idx];
            }
        }
        if( pos == 0 && (sel != 0 || ss_ind == 0) )
        {
            K = (int)len1 - (int)len + 1;
            set_zero( n_run, K );
            set_val( dist, K, 10000 );
            // n = 2;
            ps1 = pseq1 + len - 1;
            ps2 = psdv2 + len - 1;
            for( k = 0; k < K; k++, ps1 ++ )
            {
                find_max_run_length( ps1, ps2, len, Threshold, 1, n_run+k, dist+k );            
                if( n_run[k] > run_threshold ) break;
            }
            idx = find_max( n_run, K );
            lent = n_run[idx];
            if( lent >= min_overlap )
            {
                pos = idx + len - lent + 1;
                dep = lent;
                if( sel == 0 ) j_mode = 3; // n*2-1
                else           j_mode = 2; // (3-n)*2
                dst = dist[idx];
            }
        }
         
    }

    if(pos == 0 && len >= min_overlap) 
    {
        if( sel == 0 || ss_ind == 0 )
        {
            K = len - min_overlap + 1;
            set_zero( n_run, K );
            set_val( dist, K, 10000 );
            // n = 1;
            ps1 = pseq1 + len1 - min_overlap;
            ps2 = psdv1;
            ln = min_overlap;
            for( k = 0; k < K; k++, ps1 --, ln ++ )
            {
                find_max_run_length( ps1, ps2, ln, Threshold, 0, n_run+k, dist+k );            
                if( n_run[k] > run_threshold ) break;
            }
            idx = find_max( n_run, K );
            lent = n_run[idx];
            if( lent >= min_overlap )
            {
                pos = len1 - min_overlap - idx + 1;
                dep = lent;
                if( sel == 0 ) j_mode = 1; // n*2-1
                else           j_mode = 4; // (3-n)*2
                dst = dist[idx];
            } 
        }
        if( pos == 0 && ( sel != 0 || ss_ind == 0 ) )
        {
            K = len - min_overlap + 1;
            set_zero( n_run, K );
            set_val( dist, K, 10000 );
            // n = 2;
            ln = min_overlap;
            ps1 = pseq1 + ln - 1;
            ps2 = psdv2 + len - 1;
            for( k = 0; k < K; k++, ps1 ++, ln ++ )
            {
                find_max_run_length( ps1, ps2, ln, Threshold, 1, n_run+k, dist+k );            
                if( n_run[k] > run_threshold ) break;
            }
            idx = find_max( n_run, K );
            lent = n_run[idx];
            if( lent >= min_overlap )
            {
                pos = min_overlap - lent + idx + 1;
                dep = lent;
                if( sel == 0 ) j_mode = 3; // n*2-1
                else           j_mode = 2; // (3-n)*2
                dst = dist[idx];
            }
        }
    }
    
    plhs[0] = mxCreateDoubleScalar((double)pos);
    plhs[1] = mxCreateDoubleScalar((double)dep);
    plhs[2] = mxCreateDoubleScalar((double)j_mode);
    plhs[3] = mxCreateDoubleScalar((double)dst);
    plhs[4] = mxCreateDoubleScalar((double)n_junctions);
}
