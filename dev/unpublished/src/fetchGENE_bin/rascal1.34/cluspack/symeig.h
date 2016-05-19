/**********************************************************************/
/*                                                                    */
/*        Symmetric Householder reduction to tridiagonal form.        */
/*                                                                    */
/**********************************************************************/

void tred2(double **V, int n, double *d, double *e);

/***********************************************************************************/
/*      ---------------------------------------------------------------------------*/
/*      Symmetric tridiagonal QL algorithm.                                        */
/*      ---------------------------------------------------------------------------*/
/*                                                                                 */
/***********************************************************************************/
void tql2(double **V, int n, double *d, double *e);

/*     ---------------------------------------------------------------------------*/
/*      Eigenvalues and eigenvectors of a real matrix.                            */
/*                                                                                */
/*    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is         */
/*    diagonal and the eigenvector matrix V is orthogonal.                        */
/*    I.e. A = V.times(D.times(V.transpose())) and                                */
/*    V.times(V.transpose()) equals the identity matrix.                          */
/*     ---------------------------------------------------------------------------*/
/*                                                                                */
/**********************************************************************************/

void symeig(double **A, int n, double **V, double *D);
