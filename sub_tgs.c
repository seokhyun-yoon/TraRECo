#include "mex.h"

/* The computational routine */
void MVProduct(register double *mat, register double *ivec, register double *ovec, mwSize nr, mwSize nc)
{
    register mwSize i, j;
    register double *z, *t;
    /* multiply each element y by x */
    for (i=0; i<nr; i++) 
    {
        z = ivec;
        t = mat+i;
        *ovec = 0;
        for (j=0; j<nc; j++){ *ovec += *t * *z++; t += nr; } 
        ovec ++;
    }
}
double ssum(register double *ivec, mwSize len )
{
    register int i;
    register double ss = 0;
    /* multiply each element y by x */
    for (i=0; i<len; i++){ ss += (*ivec)*(*ivec); ivec++; }
    return ss;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    register double *y, *T, *x, *z, *zt, *yt, sm;
    double er = 0, er_prev = 0, swgt, lambda, step, loop;
    register int nrow, ncol, leny, k, m, n;

    /* create a pointer to the real data in the input matrix  */
    y = (double *)mxGetPr(prhs[0]);
    T = (double *)mxGetPr(prhs[1]);
    leny = (int)mxGetN(prhs[0]);
    nrow = (int)mxGetM(prhs[1]);
    ncol = (int)mxGetN(prhs[1]);
    
    swgt = (double)mxGetScalar(prhs[2]);
    lambda = (double)mxGetScalar(prhs[3]);
    step = (double)mxGetScalar(prhs[4]);
    loop = (double)mxGetScalar(prhs[5]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)ncol,1,mxREAL);
    x = (double *)mxGetPr(plhs[0]);
/*
    MVProduct( T, y, x, nrow, ncol );
*/
    z = (double*)mxMalloc(nrow*sizeof(double));
    
    for( k = 0; k<1; k ++)
    {
        MVProduct( T, x, z, nrow, ncol );
        zt = z;
        yt = y;
        for( m=0; m<nrow; m++ ){ *zt = (*yt++ - *zt)*swgt; zt++; }
        er = ssum( z, nrow );
        for( m=0; m<ncol; m++ )
        {
            sm = 0;
            zt = z;
            yt = T+m*nrow;
            for( n=0; n<nrow; n++ ) sm += *zt++ * *yt++;
            sm = 2*sm - lambda;
            x[m] = x[m] + step*sm;
            if( x[m] < 0 ) x[m] = 0;
        }
        er_prev = er;
    }
    for( k = 1; k<loop; k ++)
    {
        MVProduct( T, x, z, nrow, ncol );
        zt = z;
        yt = y;
        for( m=0; m<nrow; m++ ) *zt++ =  (*yt++ - *zt)*swgt;
        er = ssum( z, nrow );
        if( er > er_prev ) step = step*0.5;
        for( m=0; m<ncol; m++ )
        {
            sm = 0;
            zt = z;
            yt = T+m*nrow;
            for( n=0; n<nrow; n++ ) sm += *zt++ * *yt++;
            sm = 2*sm - lambda;
            x[m] = x[m] + step*sm;
            if( x[m] < 0 ) x[m] = 0;
        }
        er_prev = er;
    }
    mxFree(z);
    
/*
    function [x, er] = f04_Lasso_fw2( y, T, wgt, Lambda, Step_size, N_loop )

    [nr,nc] = size(T);
    x = zeros(nc,1);
    er = zeros(N_loop,1);
    for k = 1:1:1
        ev = (y-T*x).*sqrt(wgt);
        er(k,1) = sum( abs(ev).^2 );
        for m = 1:1:nc
            sm = 2*ev'*T(:,m) - Lambda;
            x(m) = x(m) + Step_size*sm;
            if x(m) < 0
                x(m) = 0;
            end
        end
    end
    
    for k = 2:1:N_loop
        ev = (y-T*x).*sqrt(wgt);
        er(k,1) = sum( abs(ev).^2 );
        if er(k,1) > er(k-1,1)
            Step_size = Step_size/2;
        end
        for m = 1:1:nc
            sm = 2*ev'*T(:,m) - Lambda;
            x(m) = x(m) + Step_size*sm;
            if x(m) < 0
                x(m) = 0;
            end
        end
    end
*/
    /* create the output matrix */
    plhs[1] = mxCreateDoubleScalar((double)er);
}
