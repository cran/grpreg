#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

static double *vector(int n);
static void free_vector(double *v);
static double **matrix(int nr, int nc);
static void free_matrix(double **M, int nr);
static double **as_matrix(double *v, int nr, int nc);
static void as_vector(double *v, double **M, int nr, int nc);

static double *vector(int n)
{
  double *v;
  v = Calloc(n, double);
  return v;
}

static void free_vector(double *v)
{
  Free(v);
}

static double **matrix(int nr, int nc)
{
  int   i;
  double **M;
  M = Calloc(nr, double *);
  for (i = 0; i < nr; i++) M[i] = Calloc(nc, double);
  return M;
}

static void free_matrix(double **M, int nr)
{
  int   i;
  for (i = 0; i < nr; i++) Free(M[i]);
  Free(M);
}

static double **as_matrix(double *v, int nr, int nc)
{
  int i,j;
  double **M;

  M = Calloc(nr, double *);

  for (i = 0; i < nr; i++) M[i] = Calloc(nc, double);
  for (j = 0; j < nc; j++)
    {
      for (i = 0; i < nr; i++) M[i][j] = v[j*nr+i];
    }
  return M;
}

static double **t_as_matrix(double *v, int nr, int nc)
{
  int i,j;
  double **M;

  M = Calloc(nr, double *);

  for (i = 0; i < nr; i++) M[i] = Calloc(nc, double);
  for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++) M[i][j] = v[j*nr+i];
    }
  return M;
}

static void as_vector(double *v,double **M, int nr, int nc)
{
  int i,j;

  for (j = 0; j < nc; j++)
    {
      for (i = 0; i < nr; i++) v[j*nr+i] = M[i][j];
    }
}

static int checkConvergence(double *beta_new, double *beta_old, double eps, int p)
{
  int j;
  int converged = 1;
  for (j=0; j < p; j++)
    {
      if (beta_new[j]!=0 & beta_old[j]!=0)
	{
	  if (fabs((beta_new[j]-beta_old[j])/beta_old[j]) > eps)
	    {
	      converged = 0;
	      break;
	    }
	}
      else if (beta_new[j]==0 & beta_old[j]!=0)
	{
	  converged = 0;
	  break;
	}
      else if (beta_new[j]!=0 & beta_old[j]==0)
	{
	  converged = 0;
	  break;
	}
    }
  return(converged);
}

static double S(double x,double y)
{
  if (fabs(x) <= y) return(0);
  else
    {
      if (x > 0) return(x-y);
      else return(x+y);
    }
}

static void gLasso(double *beta, double *x, double *w, double *r, int K0, int Kj, int n, double lambda, double *penpars, double lambda2, double *df)
{
  int i, j, k, j1, j2, k1, k2, K;
  K = Kj - K0;
  double sxr, sxx, oldbeta, gradient_norm, sbb, ljk, delta, s;
  double *u;
  u = vector(K);
  delta = .00000001 + penpars[0]*lambda;

  /* Calculate gradient_norm */
  gradient_norm = 0;
  for (j1=K0; j1<Kj; j1++)
    {
      u[j1-K0] = 0;
      for (i=0; i<n; i++)
	{
	  u[j1-K0] = u[j1-K0] + x[n*j1+i]*w[i]*r[i];
	  for (j2=K0; j2<Kj; j2++)
	    {
	      u[j1-K0] = u[j1-K0] + x[n*j1+i]*w[i]*x[n*j2+i]*beta[j2];
	    }
	}
      gradient_norm = gradient_norm + pow(u[j1-K0],2);
    }
  gradient_norm = sqrt(gradient_norm);
  /*gradient_norm = sqrt(gradient_norm);*/
  /*Rprintf("%d %f %f \n",K0,gradient_norm/n,sqrt(K)*lambda);*/
  if (gradient_norm/n > sqrt(K)*lambda + delta)
    {
      sbb = delta;
      for (j=K0; j<Kj; j++) sbb = sbb + pow(beta[j],2);
      for (j=K0; j<Kj; j++)
	{
	  oldbeta = beta[j];
	  sxr=0;
	  for (i=0; i<n; i++) sxr = sxr + w[i]*x[n*j+i]*r[i];
	  if (w[0]==1) sxx=n;
	  else
	    {
	      sxx=0;
	      for (i=0; i<n; i++) sxx = sxx + w[i]*pow(x[n*j+i],2);
	    }
	  ljk = sqrt(K)*lambda/sqrt(sbb);
	  beta[j] = (sxr+sxx*oldbeta)/(sxx+n*(ljk + lambda2));
	  if (fabs(beta[j]-oldbeta) > .000001)
	    {
	      for (i=0; i<n; i++) r[i] = r[i] - (beta[j]-oldbeta)*x[n*j+i];
	      sbb = sbb + pow(beta[j],2) - pow(oldbeta,2);
	    }
	  df[0] = df[0] + fabs(beta[j])/fabs(sxr/sxx+beta[j]);
	}
    }
  else
    {
      for (j=K0; j<Kj; j++)
	{
	  if (beta[j]==0) continue;
	  else
	    {
	      oldbeta = beta[j];
	      beta[j] = 0;
	      for (i=0; i<n; i++) r[i] = r[i] + oldbeta*x[n*j+i];
	    }
	}
    }
  free_vector(u);
}

static void gBridge(double *beta, double *x, double *w, double *r, int K0, int Kj, int n, double lambda, double *penpars, double lambda2, double *df)
{
  int i, j, K;
  K = Kj - K0;
  double sxr, sxx, oldbeta, sab, gamma, delta, ljk;
  gamma = penpars[0];
  delta = .0000001;

  sab = 0;
  for (j=K0; j<Kj; j++) sab = sab + fabs(beta[j]);
  if (sab==0) return;

  for (j=K0; j<Kj; j++)
    {
      if (sab > delta)
	{
	  oldbeta = beta[j];
	  sxr=0;
	  for (i=0; i<n; i++) sxr = sxr + w[i]*x[n*j+i]*r[i];
	  if (w[0]==1) sxx=n;
	  else
	    {
	      sxx=0;
	      for (i=0; i<n; i++) sxx = sxx + w[i]*pow(x[n*j+i],2);
	    }
	  ljk = sqrt(K)*lambda*gamma*pow(sab,gamma-1);
	  beta[j] = S(sxr+sxx*oldbeta,n*ljk)/(sxx+n*lambda2);
	}
      else beta[j] = 0;
      if (fabs(beta[j]-oldbeta) > delta)
	{
	  for (i=0; i<n; i++) r[i] = r[i] - (beta[j]-oldbeta)*x[n*j+i];
	  sab = sab + fabs(beta[j]) - fabs(oldbeta);
	}
      df[0] = df[0] + fabs(beta[j])/fabs(sxr/sxx+beta[j]);
    }
}

static double MCP(double theta, double l, double a)
{
  theta = fabs(theta);
  if (theta <= a*l) return(l*theta - pow(theta,2)/(2*a));
  else return(a*pow(l,2)/2);
}

static double dMCP(double theta, double l, double a)
{
  theta = fabs(theta);
  if (theta < a*l) return(l-theta/a);
  else return(0);
}

static void gMCP(double *beta, double *x, double *w, double *r, int K0, int Kj, int n, double lambda, double lambda2, double *penpars, double *df)
{
  int i, j, K;
  K = Kj - K0;
  double sxr, sxx, oldbeta, sMCP, ljk;
  double phi,a,l1,l2,a1,a2;
  a = penpars[0];
  l1 = lambda;
  l2 = lambda;
  a1 = (K*a*pow(lambda,2))/(2*l1);
  a2 = a;

  sMCP = 0;
  for (j=K0; j<Kj; j++) sMCP = sMCP + MCP(beta[j],l2,a2);

  for (j=K0; j<Kj; j++)
    {
      oldbeta = beta[j];
      sxr=0;
      for (i=0; i<n; i++) sxr = sxr + w[i]*x[j*n+i]*r[i];
      if (w[0]==1) sxx=n;
      else
	{
	  sxx=0;
	  for (i=0; i<n; i++) sxx = sxx + w[i]*pow(x[j*n+i],2);
	}
      if (lambda==0) ljk = 0;
      else ljk = dMCP(sMCP,l1,a1)*dMCP(beta[j],l2,a2);

      beta[j] = S(sxr+sxx*oldbeta,n*ljk)/(sxx+n*lambda2);
      /*Rprintf("%d %f %f %f %f %f %f\n",j,lambda,dMCP(sMCP,l1,a1),dMCP(beta[j],l2,a2),ljk,oldbeta,beta[j]);*/
      df[0] = df[0] + fabs(beta[j])/fabs(sxr/sxx+beta[j]);

      if (fabs(beta[j]-oldbeta) > .00001)
	{
	  for (i=0; i<n; i++) r[i] = r[i] - (beta[j]-oldbeta)*x[j*n+i];
	  sMCP = sMCP + MCP(beta[j],l2,a2) - MCP(oldbeta,l2,a2);
	}
    }
}

static void gpPathFit(double *beta, int *counter, double *df, double *x, double *y, int *group, char **family, int *n, int *p, int *J, int *K, char **penalty, double *lambda, int *nlambda, double *eps, int *max_iter, int *verbose, int *monitor, int *n_monitor, double *penpars, double *lambda2, int *warn_conv)
{
  int l, i, j, g, K0, Kj, status, converged;
  converged=0;
  double sxr, sxx, oldbeta;
  double ybar, yp, yy;
  double eta, pi;
  double **Beta, *r, *w, *beta_old;
  Beta = matrix(*nlambda,*p);
  r = vector(*n);
  beta_old = vector(*p);
  w = vector(*n);

  /* Initial setup */
  if (strcmp(penalty[0],"gBridge")==0)
    {
      if (strcmp(family[0],"gaussian")==0)
	{
	  if (group[0]==0)
	    {
	      sxr=0;
	      for (i=0; i<*n; i++) sxr = sxr + y[i];
	      beta_old[0] = sxr / *n;
	      for (i=0;i<*n;i++)
		{
		  w[i] = 1;
		  r[i] = y[i] - beta_old[0];
		}
	      j=1;
	    }
	  else for (i=0;i<*n;i++)
	    {
	      w[i] = 1;
	      r[i] = y[i];
	      j=0;
	    }
	  for (;j<*p;j++)
	    {
	      sxr=0;
	      for (i=0; i<*n; i++) sxr = sxr + x[j*n[0]+i]*r[i];
	      beta_old[j] = sxr / *n;
	    }
	  for (i=0;i<n[0];i++)
	    {
	      r[i] = y[i];
	      for (j=0;j<*p;j++) r[i] = r[i] - x[j*n[0] + i]*beta_old[j];
	    }
	}
      else
	{
	  if (group[0]==0)
	    {
	      sxr=0;
	      for (i=0; i<*n; i++) sxr = sxr + y[i];
	      pi = sxr/n[0];
	      beta_old[0] = log((pi)/(1-pi));
	      for (i=0;i<*n;i++)
		{
		  w[i] = sqrt(pi*(1-pi));
		  r[i] = (y[i] - pi)/w[i];
		}
	      j=1;
	    }
	  else for (i=0;i<*n;i++)
	    {
	      w[i] = 1;
	      r[i] = (y[i] - 0.5)/0.5;
	      j=0;
	    }
	  for (;j<*p;j++)
	    {
	      sxr=0;
	      for (i=0; i<*n; i++) sxr = sxr + w[i]*x[j*n[0]+i]*r[i];
	      beta_old[j] = sxr / *n;
	    }
	}
    }
  else
    {
      for (j=0;j<*p;j++) beta_old[j] = 0;
      if (strcmp(family[0],"gaussian")==0)
	{
	  for (i=0;i<*n;i++)
	    {
	      w[i] = 1;
	      r[i] = y[i];
	    }
	}
    }

  /* Path */
  for (l=0;l<*nlambda;l++)
    {
      if (l==0) for (j=0;j<*p;j++) Beta[0][j] = beta_old[j];
      else for (j=0;j<*p;j++) Beta[l][j] = beta_old[j] = Beta[l-1][j];
      if (*verbose) Rprintf("Starting new fit: lambda = %f\n",lambda[l]);
      while (counter[l] < *max_iter)
	{
	  converged = 0;
	  status = 0;
	  counter[l] = counter[l] + 1;
	  if (*verbose) Rprintf("Iteration: %d\n",counter[l]);

	  /* Approximate L                 */
	  if (strcmp(family[0],"binomial")==0)
	    {
	      ybar = 0;
	      yp = 0;
	      yy = 0; for (i=0;i<*n;i++) ybar = ybar + y[i];
	      ybar = ybar / *n;
	      for (i=0;i<*n;i++)
		{
		  eta = 0;
		  for (j=0;j<*p;j++) eta = eta + x[j * *n + i]*Beta[l][j];
		  pi = exp(eta)/(1+exp(eta));
		  if (pi > .999)
		    {
		      pi = 1;
		      w[i] = .001;
		    }
		  else if (pi < .001)
		    {
		      pi = 0;
		      w[i] = .001;
		    }
		  else w[i] = sqrt(pi*(1-pi));
		  r[i] = (y[i] - pi)/w[i];
		  yp = yp + pow(y[i]-pi,2);
		  yy = yy + pow(y[i]-ybar,2);
		}
	    }
	  df[l] = 0;
	  /*Rprintf("%f %f\n",r[1],X_[1][1]);*/
  
	  /* Update unpenalized covariates */
	  for (j=0;; j++)
	    {
	      if (group[j]!=0) break;
	      sxr=0;
	      for (i=0; i<*n; i++) sxr = sxr + w[i]*x[j * *n + i]*r[i];
	      if (w[0]==1) sxx=*n;
	      else
		{
		  sxx=0;
		  for (i=0; i<*n; i++) sxx = sxx + w[i]*pow(x[j * *n + i],2);
		}

	      oldbeta = Beta[l][j];
	      Beta[l][j] = sxr/sxx + Beta[l][j];
	      /*Rprintf("%d %f %f %f\n",j,oldbeta,Beta[l][j],sxr);*/
	      for (i=0; i<*n; i++) r[i] = r[i] - (Beta[l][j]-oldbeta)*x[j * *n + i];
	      df[l] = df[l] + 1;
	    }
	  for (g=0; g<*J; g++)
	    {
	      K0 = j;
	      Kj = j + K[g];
	      /*Rprintf("%d %d\n",K0,Kj);*/
	      /*Rprintf("%f %f %f %f %f\n",beta[0],beta[1],beta[2],beta[3],beta[4]);*/
	      if (strcmp(penalty[0],"gLasso")==0) gLasso(Beta[l],x,w,r,K0,Kj,*n,lambda[l],penpars,lambda2[l],&df[l]);
	      if (strcmp(penalty[0],"gBridge")==0) gBridge(Beta[l],x,w,r,K0,Kj,*n,lambda[l],penpars,lambda2[l],&df[l]);
	      if (strcmp(penalty[0],"gMCP")==0) gMCP(Beta[l],x,w,r,K0,Kj,*n,lambda[l],lambda2[l],penpars,&df[l]);
	      j = Kj;
	    }
	  if (strcmp(family[0],"binomial")==0)
	    {
	      if (yp/yy < .01)
		{
		  warning("Model saturated; exiting...");
		  status = 1;
		}
	    }
	  if (status==1) break;
	  /*Rprintf("%d %d\n",counter[0],max_iter[0]);
	    Rprintf("%f %f %f\n",beta_new[0],beta_new[1],beta_new[2]);
	    Rprintf("%f %f %f\n",beta_old[0],beta_old[1],beta_old[2]);*/
	  if (*n_monitor != 0)
	    {
	      for (i=0; i<*n_monitor;i++) Rprintf("%.3f ",Beta[l][monitor[i]]);
	      Rprintf("\n");
	    }
	  if (checkConvergence(Beta[l],beta_old,*eps,*p))
	    {
	      converged  = 1;
	      break;
	    }
	  for (j=0;j<*p;j++) beta_old[j] = Beta[l][j];
	}
      if (converged==0 & warn_conv[0]) warning("Failed to converge");
    }
  as_vector(beta,Beta,*nlambda,*p);

  free_matrix(Beta,*nlambda);
  free_vector(beta_old);
  free_vector(w);
  free_vector(r);
}

static void gLassoMax(double *max, double *x, double *r, double *w, int K0, int Kj, int n)
{
  int i, j, K;
  K = Kj - K0;
  double gradient_norm, u, max_g;

  /* Calculate gradient_norm */
  gradient_norm = 0;
  for (j=K0; j<Kj; j++)
    {
      u = 0;
      for (i=0; i<n; i++) u = u + w[i]*x[n*j+i]*r[i];
      gradient_norm = gradient_norm + pow(u,2);
    }
  gradient_norm = sqrt(gradient_norm);
  max_g = gradient_norm/(sqrt(K)*n);
  if (*max < max_g) *max = max_g;
}

static void gBridgeMax(double *max, double *x, double *r, double *w, int K0, int Kj, int n, double *penpars)
{
  int i, j, K;
  K = Kj - K0;
  double sxr, max_g, gamma, delta;
  gamma = penpars[0];
  if (w[0]==1) delta = .35;
  else delta = .2;

  for (j=K0; j<Kj; j++)
    {
      sxr=0;
      for (i=0; i<n; i++) sxr = sxr + w[i]*x[n*j+i]*r[i];
      max_g = fabs(sxr)*pow(delta,1-gamma)/(n*gamma*pow(K,gamma));
      if (*max < max_g) *max = max_g;
    }
}

static void gMCPMax(double *max, double *x, double *r, double *w, int K0, int Kj, int n, double *penpars)
{
  int i, j, K;
  K = Kj - K0;
  double sxr, max_g, a, phi;
  a = penpars[0];

  for (j=K0; j<Kj; j++)
    {
      sxr=0;
      for (i=0; i<n; i++) sxr = sxr + w[i]*x[n*j+i]*r[i];
      max_g = sqrt(fabs(sxr)/n);
      if (*max < max_g) *max = max_g;
    }
}

static void determineMax(double *max, double *x, double *r, double *w, int *group, char **family, int *n, int *p, int *J, int *K, char **penalty, double *penpars)
{
  int j, g, K0, Kj;

  j=0;
  for (g=0; g<*J; g++)
    {
      K0 = j;
      Kj = j + K[g];
      /*Rprintf("%d %d\n",K0,Kj);*/
      /*Rprintf("%f %f %f %f %f\n",beta[0],beta[1],beta[2],beta[3],beta[4]);*/
      if (strcmp(penalty[0],"gLasso")==0) gLassoMax(max,x,r,w,K0,Kj,*n);
      if (strcmp(penalty[0],"gBridge")==0) gBridgeMax(max,x,r,w,K0,Kj,*n,penpars);
      if (strcmp(penalty[0],"gMCP")==0) gMCPMax(max,x,r,w,K0,Kj,*n,penpars);
      j = Kj;
    }
}

static const R_CMethodDef cMethods[] = {
  {"gpPathFit", (DL_FUNC) &gpPathFit, 22},
  {"determineMax", (DL_FUNC) &determineMax, 12},
  NULL
};

void R_init_grpreg(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
