// [[Rcpp::interfaces(r, cpp)]]
//'

#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
//#include <complex.h>
#include <fftw3.h>

#define REAL 0
#define IMAG 1
#define PI2 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651e+00

//using namespace Rcpp;
//using namespace std;
//using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]
int mod(int a, int n)
{
    return a - floor(a/n)*n;
}  


void l_vec_compute_arma(int k, arma::vec& l_vec, arma::vec& cn_vec, int m)
{
  int i, aa, bb;
  
  for(i=0; i<m-1; i++)
  {
    aa=mod(k, cn_vec(i));
    bb=(k-aa)/cn_vec(i);
    l_vec(i)=bb;
    k=aa;
  }

  return;
}


//[[Rcpp::export]]  
arma::vec pmn_mdfft_arma(int nnt, arma::mat pp, arma::vec nn_vec, arma::vec l_vec, arma::vec cn_vec)
{  
  arma::vec res(nnt, arma::fill::zeros);

  int mm=pp.n_cols;  
  int nn=pp.n_rows;
  
  //int nn_vec_a[mm-1] = {0};
  int* nn_vec_a = new int[mm-1];
  
  fftw_complex *in, *out;
  int i, j, k;
  int n, m, nt;
  fftw_plan p;
  double tmp, con, pij, pim, ww;
  arma::cx_double ctmp, ctmp1, ctmp2, qval, a1, a2;
  
  nt=nnt;
  n=nn;
  m=mm;
  
  
  for(i=0; i<mm-1; i++)
  {
    nn_vec_a[i]=nn_vec(i);
  }
  
  //Rprintf("nt %u, n %u, m %u \n", nt, n, m);
 
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt); 
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);  
  
  ww=2*PI2/(n+1);

  //Rprintf("ww %lf \n", ww);
  
  for(k=0; k<nt; k++)
  {
    //qval=0.0 + 0.0 * I;
    qval=arma::cx_double(0.0, 0.0); 


    l_vec_compute_arma(k, l_vec, cn_vec, m);
    
    //for(ii=0; ii<m-1; ii++)
    //{
    //  Rprintf("l_vec %u \n", l_vec[ii]);
    //}
      
    for(i=0; i<n; i++)
    {
      //ctmp=0.0 + 0.0 * I;
      
      ctmp=arma::cx_double(0.0, 0.0);
       
      
      for(j=0; j<m-1; j++)
      {
        pij=pp(i,j);

        //printf("pij: %lf, l_vec[j], %u, ww, %lf, \n", pij, l_vec[j], ww);
        
        //ctmp1=0.0+l_vec(j)*ww*I;
        
        ctmp1=arma::cx_double(0.0, l_vec(j)*ww);
        
        //printf("ctmp1: %lf +%lf*i\n", creal(ctmp1), cimag(ctmp1));

        //a1=pij+0.0*I;
        
        a1=arma::cx_double(pij, 0.0);
        
        a2=exp(ctmp1);
        //printf("a1: %lf +%lf*i\n", creal(a1), cimag(a1));
        //printf("a2: %lf +%lf*i\n", creal(a2), cimag(a2));
                
        ctmp2=a1*a2;

        //printf("ctmp2: %lf +%lf*i\n\n", creal(ctmp2), cimag(ctmp2));
        
        ctmp+=ctmp2;
      }
      
      pim=pp(i, m-1);
      ctmp+=pim;
      
      ctmp=log(ctmp);
      qval+=ctmp;
    }
    
    qval=exp(qval);
    
    in[k][REAL]= real(qval);
    in[k][IMAG]= imag(qval);
    
    //printf("qval: %lf +%lf*i\n", creal(qval), cimag(qval));
    //printf("in[k]: %lf +%lf*i\n", creal(in[k]), cimag(in[k]));

  
  }
  
  p=fftw_plan_dft(m-1, nn_vec_a, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  fftw_execute(p);
  
  con=pow(n+1, m-1);
  
  for(k=0; k<nt; k++)
  {
    tmp=out[k][REAL];
    res[k]=tmp/con;
    //printf("out[k]: %lf +%lf*i\n", creal(out[k]), cimag(out[k]));
  }
  
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);  
  
  delete[] nn_vec_a;
  
  return res;
}


Rcpp::IntegerVector rmultinom_1(unsigned int &size,Rcpp:: NumericVector &probs, unsigned int &N) {
    Rcpp::IntegerVector outcome(N);
    rmultinom(size, probs.begin(), N, outcome.begin());
    return outcome;
}


Rcpp::IntegerMatrix rmultinom_rcpp(unsigned int &n, unsigned int &size, Rcpp::NumericVector &probs) {
    unsigned int N = probs.length();
    Rcpp::IntegerMatrix sim(N, n);
    for (unsigned int i = 0; i < n; i++) {
        sim(Rcpp::_,i) = rmultinom_1(size, probs, N);
    }
    return sim;
}


//[[Rcpp::export]]
arma::vec rpmd_arma(arma::mat pp)
{
  int mm=pp.n_cols;
  int nn=pp.n_rows;
  int i;
  
  unsigned int n=1;
  unsigned int size=1;
  
  arma::mat tmp(nn, mm, arma::fill::zeros);
  
  for(i=0;i<nn;i++){
    Rcpp::NumericVector prob=Rcpp::wrap(pp.row(i));
    
    Rcpp::IntegerMatrix res=rmultinom_rcpp(n, size, prob);
    /*tmp.row(i)=arma::conv_to<arma::rowvec>::from(res);*/

    arma::mat xx = Rcpp::as<arma::vec>(res);
    tmp.row(i)=xx.as_row();
  
  }
  
  arma::vec finalres(mm, arma::fill::zeros);
  finalres=sum(tmp).as_col();
  
  return finalres;
}


//[[Rcpp::export]]
double pmd_simulation_singlepoint(arma::mat pp, arma::vec x_vec, int t)
{
    /*arma::vec res(nnt, arma::fill::zeros);*/
    double res=0;
    int mm=pp.n_cols;
    int k;
    double count=0;
  
  arma::vec sim(mm,arma::fill::zeros);
  for(k=0;k<t;k++){
    //arma::mat tmp(nn, mm, arma::fill::zeros);
    //rpmd(pp).as_row().print();
    sim=rpmd_arma(pp);
    //sim.print();
    if(all(sim==x_vec)){
        count++;
    }
  }
  res=count/t;


    return res;
}



//[[Rcpp::export]]
arma::vec pmd_simulation_allpoints(arma::mat pp, int nnt, arma::vec l_vec, arma::vec cn_vec, int t){
    arma::vec res(nnt, arma::fill::zeros);

    int mm=pp.n_cols;
    int nn=pp.n_rows;

    int j, k;
    int n, m, nt;
    
    nt=nnt;
    n=nn;
    m=mm;
    
    arma::vec x_vec(m, arma::fill::zeros);
    
    for(k=0; k<nt; k++)
    {
        l_vec_compute_arma(k, l_vec, cn_vec, m);
       /* for (int i=0; i<m-1; i++) {
            Rcpp::Rcout << l_vec(i)<<arma::endl;
       }
        Rcpp::Rcout <<"print x_vec" <<arma::endl; */
        for (j=0; j<m-1; j++) {
            x_vec(j) = l_vec(j);
        }
        x_vec(m-1)=n-sum(l_vec);
        /* for(int i=0;i<m;i++){
            Rcpp::Rcout << x_vec(i)<<arma::endl;
        } */
        res(k) = pmd_simulation_singlepoint(pp, x_vec, t);
        
    }
    
    return res;
}











