/*
This code is a free software. 
It can be redistributed and/or modified under the terms of the MIT License.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Copyright (c) 2021 by Jun Liu (juliu@siue.edu).
Version date 4/13/2021, created by Jun Liu. 
*/
static char help[] = "Solves a linear 2D wave equation with leapfrog scheme and PinT solver.\n\n";
#include <petscksp.h> //based on MPI 

/* define solution function */
PetscReal y_sol(PetscReal t, PetscReal x, PetscReal y){
    return sin(2 * PETSC_PI * t) * x * (x - 1.0) * y * (y - 1.0);
}
/* define RHS function */
PetscReal f(PetscReal t, PetscReal x, PetscReal y){
    return -4 * PetscPowReal(PETSC_PI, 2) * sin(2 * PETSC_PI * t) * x * (x - 1.0) * y * (y - 1.0) - 2 * sin(2 * PETSC_PI * t) * (x * (x - 1.0) + y * (y - 1.0));
}
/* define initial condition */
PetscReal y0ic(PetscReal x, PetscReal y){
    return 0.0;
}
/* define velocity initial condition */
PetscReal y1ic(PetscReal x, PetscReal y){
    return 2 * PETSC_PI * x * (x - 1.0) * y * (y - 1.0);
}
/*define complex Chebyshev function of first kind*/
PetscScalar ChebyT(PetscInt n, PetscScalar x){
  return PetscCosComplex(n * PetscAcosComplex(x));  
}
/*define complex Chebyshev function of second kind*/
PetscScalar ChebyU(PetscInt n, PetscScalar x){
  return PetscSinComplex((n + 1) * PetscAcosComplex(x))/PetscSinComplex(PetscAcosComplex(x));
}

// Thomas algorithm for solve S*x=b with S=pentadiagonal(e,0,d,0,e) by L'*D*L factorization of S
void ThomasPendiag(int n, PetscScalar *d, PetscScalar *e, PetscScalar *b, PetscScalar *x){
  PetscScalar alpha[n];
  PetscScalar delta[n - 2];
  // Factor A=LDL'
  alpha[0] = d[0];
  delta[0] = e[0] / alpha[0];
  alpha[1] = d[1];
  delta[1] = e[1] / alpha[1];
 
  for (int k = 2; k <= n - 3; k++) {
    alpha[k] = d[k] - e[k - 2] * delta[k - 2];
    delta[k] = e[k] / alpha[k];
  }
  alpha[n - 2] = d[n - 2] - e[n - 4] * delta[n - 4];
  alpha[n - 1] = d[n - 1] - e[n - 3] * delta[n - 3];

  // Update Lz=b
  PetscScalar z[n];
  z[0] = b[0];
  z[1] = b[1];
  for (int k = 2; k <= n - 1; k++)  {
    z[k] = b[k] - delta[k - 2] * z[k - 2];
  }
  //Dc=z combined into next step
  //Backsubstitution L'x=c
  x[n - 1] = z[n - 1] / alpha[n - 1];
  x[n - 2] = z[n - 2] / alpha[n - 2];
  for(int k = n - 3; k >= 0; k--)  {
    x[k] = z[k] / alpha[k] - delta[k] * x[k + 2];
  }
}
//Thomas algorithm for tridiag(b,a,c)*y=f system
void ThomasTridiag(int n, PetscScalar *b, PetscScalar *a, PetscScalar *c, PetscScalar *f, PetscScalar *y){
  PetscScalar v[n]; 
  PetscScalar w = a[0];
  y[0] = f[0] / w;
  for (int i = 1; i <= n - 1; i++)  {
    v[i - 1] = c[i - 1] / w;
    w = a[i] - b[i] * v[i - 1];
    y[i] = (f[i] - b[i] * y[i - 1]) / w;
  }
  for (int j = n - 2; j >= 0; j--)  {
    y[j] = y[j] - v[j] * y[j + 1];
  }
}

/*define functions for finding eigenvalues by Newton iterations:
p=@(n,theta) sin(n*theta)-1i*cos(n*theta).*sin(theta); % polynomial
pp=@(n,theta) n*cos(n*theta)-1i*cos(n*theta).*cos(theta)+1i*n*sin(n*theta).*sin(theta); 
*/
PetscScalar pfun(int n, PetscScalar theta){
  return PetscSinComplex(n * theta) - I * PetscCosComplex(n * theta) * PetscSinComplex(theta);
}
PetscScalar ppfun(int n, PetscScalar theta){
  return n * PetscCosComplex(n * theta) - I * PetscCosComplex(n * theta) * PetscCosComplex(theta) + I * n * PetscSinComplex(n * theta) * PetscSinComplex(theta);
}
PetscInt fasteigB(int nt, double dt, PetscScalar *eig, PetscScalar *iV){
  //first compute eigenvalues by Newton iterations
  int it, itmax = 0;
  double tol = 1e-10; //Newton iteration stopping tolerance
  //Newton iteration for each eigenvalue seperately
  PetscScalar eig0;
  for (int k = 0; k < nt; k++)  {
    if (k < (nt / 2))    {
      eig0 = 0.5 * (k + 1) * PETSC_PI * (1.0 / nt + 1.0 / (nt + 1)) + PETSC_i / nt; 
      it = 0;
      while (PetscAbsComplex(pfun(nt, eig0)) > tol)      {
        eig0 = eig0 - pfun(nt, eig0) / ppfun(nt, eig0);
        it++;
      }
      if (it > itmax){ itmax = it; }
      eig[k] = eig0;
      eig[nt - 1 - k] = PETSC_PI - conj(eig[k]); //assuming nt even
    }
    eig[k] = PetscCosComplex(eig[k]) * PETSC_i / dt;
  }
  //construct inverse of V using fast O(n^2) algorithm
  //construct vectors for algorithm use
  PetscScalar r[nt],b[nt],pn[nt],e[nt],d[nt],xj;
  for(int k=0;k<nt;k++){
    if(k==(nt-1)) {r[k]=2.0;}
    else if(k==(nt-2)) {r[k]=PETSC_i;}
    else {r[k]=0.0;}
    e[k]=-1.0;  
    if(k==0||k==(nt-1)){d[k]=3.0;}
    else {d[k]=2.0;}
  } 
  //STEP-1, get b
  ThomasPendiag(nt,d,e,r,b); 
  //STEP-2, get a_j=the j-th row of iV, notice iV is  row-major stored
  PetscScalar av[nt],bv[nt],ipowk;
  PetscScalar *A = malloc(nt * nt*sizeof(PetscScalar));//for computation use
  for(int j=0;j<nt;j++){ 
    xj=-PETSC_i*eig[j]*dt; //back to xj
    pn[j]=(-nt*ChebyT(nt,xj)+xj*ChebyU(nt-1,xj))/(1-xj*xj)-PETSC_i*nt*ChebyU(nt-1,xj);
    //solve a tridiagonal system for each xj
    for(int k=0;k<nt;k++){ av[k]=2.0*xj; bv[k]=-b[k]/pn[j];} //drop 2* since /2 in STEP-3
    ThomasTridiag(nt,e,av,e,bv,&(A[j*nt])); 
  }
  //STEP-3: W=A*S/2; since S is sparse matrix, we can get O(n^2) implementation
  for(int k=0;k<nt;k++){ //column of A
    ipowk=PetscPowComplex(-PETSC_i, k);
    for(int j=0;j<nt;j++){ //row of A
      if(k==0) {iV[j*nt+k]=3.0*A[j*nt+k]-A[j*nt+k+2];}
      else if(k==1) {iV[j*nt+k]=2.0*A[j*nt+k]-A[j*nt+k+2];}
      else if(k==(nt-2)) {iV[j*nt+k]=2.0*A[j*nt+k]-A[j*nt+k-2];}
      else if(k==(nt-1)) {iV[j*nt+k]=3.0*A[j*nt+k]-A[j*nt+k-2];}
      else {iV[j*nt+k]=2.0*A[j*nt+k]-A[j*nt+k+2]-A[j*nt+k-2];}
      iV[j*nt+k]=ipowk*iV[j*nt+k];  
    }
  } 
  free(A);
  return itmax;
}

int main(int argc, char **args) //Start of program
{
  Vec x, b, max;                                       
  Mat Ax, V,  U, X, VX, Xloc, Vloc, Floc, iV;  
  KSP kspb;                                          
  PC pc;                                      
  MPI_Comm comm;
  PetscReal   err, xi, yk, tj, dt, dx, dx2, T = 2.0, xa = 0.0, xb = 1.0; 
  PetscInt i, j, Ii, J, lev, maxL = 7, maxN = 7, nx, nt, *idxm, nx0, k, *idxn, eigiter,l;
  PetscMPIInt size, rank;
  PetscScalar value[3], v, *aX, *aXlocal, *eg, *aVlocal, *aVX, *aVXlocal,  *aFlocal;
  PetscLogDouble t1, t2, tb1, tb2, ta1, ta2, tc1, tc2, teig1, teig2;
  
//initilize MPI
   PetscInitialize(&argc, &args, (char *)0, help); 
   MPI_Comm_size(PETSC_COMM_WORLD, &size); 
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  comm = MPI_COMM_SELF; 

  //read inputed mesh level if provided
  PetscOptionsGetInt(NULL, NULL, "-maxL", &maxL, NULL); //mesh level in time
   PetscOptionsGetInt(NULL, NULL, "-maxN", &maxN, NULL); //mesh level in space   
   PetscPrintf(PETSC_COMM_WORLD, "Proc#\t (nx,nx,nt) \t& Error \t& Iter \t& CPU \t (eig(B), \t Step-a, \t Step-b, \t Step-c)\n");
  
  lev = maxL; nt = T * PetscPowReal(2, lev); dt = T / nt; //time step size
  nx0 = PetscPowReal(2, maxN); nx = nx0 * nx0;
  dx = (xb - xa) / (nx0 + 1); dx2 = PetscPowReal(dx, 2); //space step size
    
  if (rank == 0) { //setup solution matrix U on root process
    MatCreateSeqDense(comm, nt, nx, NULL, &U);
    //true solution
    for (j = 0; j < nt; j++)  { //j-th time step, row
      tj = (j + 1) * dt;
      for (i = 0; i < nx0; i++)   {
        xi = xa + (i + 1) * dx;
        for (k = 0; k < nx0; k++)   {
          yk = xa + (k + 1) * dx;
          MatSetValue(U, j, i * nx0 + k, y_sol(tj, xi, yk), INSERT_VALUES); 
        }
      }
    }
     MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY); 
     MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY); 
  }
  int NPROWS=1, NPCOLS=size, BLOCKROWS3 = nt, BLOCKCOLS3 = nx / NPCOLS;  
  //construct partial columns of F locally, avoid communication 
  MatCreateSeqDense(comm, BLOCKROWS3, BLOCKCOLS3, NULL, &Floc);
 for (j = 0; j < nt; j++) { //j-th time step, row
      tj = (j + 1) * dt;
      for(l=0;l<BLOCKCOLS3;l++)   {
        i=(l+rank*BLOCKCOLS3)/nx0; // global column index =i * nx0 + k 
        k=(l+rank*BLOCKCOLS3)%nx0;
        xi = xa + (i + 1) * dx; 
          yk = xa + (k + 1) * dx; 
          if (j == 0)        {
                MatSetValue(Floc, j, l, f(tj, xi, yk) + y1ic(xi, yk) / (2 * dt), INSERT_VALUES);
            } else if (j == 1) {
                MatSetValue(Floc, j, l, f(tj, xi, yk) - y0ic(xi, yk) / (4 * dt * dt), INSERT_VALUES);
            } else  {
                MatSetValue(Floc, j, l, f(tj, xi, yk), INSERT_VALUES);
            } 
      }
    }
    MatAssemblyBegin(Floc, MAT_FINAL_ASSEMBLY); 
    MatAssemblyEnd(Floc, MAT_FINAL_ASSEMBLY); 

  PetscTime(&t1);
  MatCreateSeqDense(comm, nt, nt, NULL, &V);
  MatCreateSeqDense(comm, nt, nx, NULL, &X); 
  PetscScalar *mem1 = malloc(nt * nt*sizeof(PetscScalar));
  PetscScalar (*iV2)[nt] = (PetscScalar(*)[nt])mem1; //eigenvector matrix 
  PetscScalar *mem2 = malloc(nt * nt*sizeof(PetscScalar));
  PetscScalar (*V2)[nt] = (PetscScalar(*)[nt])mem2; //eigenvector matrix
  
  PetscTime(&teig1);
  eg = malloc(nt* sizeof(PetscScalar));
  eigiter = fasteigB(nt, dt, eg, mem1); //O(n^2) algorithm for eigen decomposition
  
  //Construct V from paper's explicit formula, not normalized
  for (i = 0; i < nt; i++){
    for (j = 0; j < nt; j++){ 
      value[0] = PetscPowComplex(PETSC_i, i) * ChebyU(i, -dt * eg[j] * PETSC_i);
      V2[i][j] = value[0]; //save for compute inverse of V
       MatSetValues(V, 1, &i, 1, &j, value, INSERT_VALUES); 
    }
  }
   MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY); 
   MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY); 
  for (i = 0; i < nt; i++){
        eg[i] = eg[i] * eg[i]; //square eigenvalue for wave equation, but keep same eigenvector matrix
    }  

  MatCreateSeqDense(comm, nt, nt, NULL, &iV);
  for (i = 0; i < nt; i++){
    for (j = 0; j < nt; j++){
      value[0] = iV2[i][j];
       MatSetValues(iV, 1, &i, 1, &j, value, INSERT_VALUES); 
    }
  }
   MatAssemblyBegin(iV, MAT_FINAL_ASSEMBLY); 
   MatAssemblyEnd(iV, MAT_FINAL_ASSEMBLY);  
  PetscTime(&teig2); 

  //Direct PinT solver 3-Steps  ==========================================
  PetscTime(&ta1);
  if (rank == 0){ aX = malloc(nt * nx*sizeof(PetscScalar));}

  aFlocal = malloc(BLOCKROWS3 * BLOCKCOLS3* sizeof(PetscScalar));
  MPI_Datatype blocktypeF, blocktype2F;
  MPI_Type_vector(BLOCKROWS3, BLOCKCOLS3, nx, MPI_C_DOUBLE_COMPLEX, &blocktype2F);
  MPI_Type_create_resized(blocktype2F, 0, sizeof(PetscScalar), &blocktypeF);
  MPI_Type_commit(&blocktypeF);
  int disps3[NPROWS * NPCOLS];
  int counts3[NPROWS * NPCOLS];
  for (int ii = 0; ii < NPROWS; ii++){
    for (int jj = 0; jj < NPCOLS; jj++){
      disps3[ii * NPCOLS + jj] = ii * nx * BLOCKROWS3 + jj * BLOCKCOLS3;
      counts3[ii * NPCOLS + jj] = 1;
    }
  } 
  MatCreateSeqDense(comm, BLOCKROWS3, BLOCKCOLS3, NULL, &Xloc);
  //Since V^{-1}=iV, Xloc=iV*Floc, no need to solve anymore
  MatMatMult(iV, Floc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Xloc);
  MatTranspose(Xloc, MAT_INPLACE_MATRIX, &Xloc); //transpose to get row major format
  MatDenseGetArray(Xloc, &aFlocal);
  //gather results from all processes by columns
  MPI_Gatherv(aFlocal, BLOCKROWS3 * BLOCKCOLS3, MPI_C_DOUBLE_COMPLEX, aX, counts3, disps3, blocktypeF, 0, PETSC_COMM_WORLD);
  PetscTime(&ta2);

  //step (b)=================================================================
  PetscTime(&tb1);
  NPROWS = size;  NPCOLS = 1;   
  int BLOCKROWS = nt / NPROWS, BLOCKCOLS = nx / NPCOLS;  
  aXlocal =malloc(BLOCKROWS * BLOCKCOLS*sizeof(PetscScalar));

  MPI_Datatype blocktype, blocktype2;
  MPI_Type_vector(BLOCKROWS, BLOCKCOLS, nx, MPI_C_DOUBLE_COMPLEX, &blocktype2);
  MPI_Type_create_resized(blocktype2, 0, sizeof(PetscScalar), &blocktype);
  MPI_Type_commit(&blocktype);
  int disps[NPROWS * NPCOLS], counts[NPROWS * NPCOLS];
  for (int ii = 0; ii < NPROWS; ii++) {
    for (int jj = 0; jj < NPCOLS; jj++) {
      disps[ii * NPCOLS + jj] = ii * nx * BLOCKROWS + jj * BLOCKCOLS;
      counts[ii * NPCOLS + jj] = 1;
    }
  }
   //scatter data aX to each process as aXlocal
  MPI_Scatterv(aX, counts, disps, blocktype, aXlocal, BLOCKROWS * BLOCKCOLS, MPI_C_DOUBLE_COMPLEX, 0, PETSC_COMM_WORLD);
  /*
     Assemble space matrix Ax on each process to avoid communication cost
  */
  MatCreateSeqAIJ(comm, nx, nx, 5, 0, &Ax);
  for (Ii = 0; Ii < nx; Ii++){
    v = -1.0 / dx2;
    i = Ii / nx0; j = Ii - i * nx0;
    if (i > 0) {
      J = Ii - nx0;
       MatSetValues(Ax, 1, &Ii, 1, &J, &v, ADD_VALUES); 
    }
    if (i < nx0 - 1) {
      J = Ii + nx0;
       MatSetValues(Ax, 1, &Ii, 1, &J, &v, ADD_VALUES); 
    }
    if (j > 0) {
      J = Ii - 1;
       MatSetValues(Ax, 1, &Ii, 1, &J, &v, ADD_VALUES); 
    }
    if (j < nx0 - 1) {
      J = Ii + 1;
       MatSetValues(Ax, 1, &Ii, 1, &J, &v, ADD_VALUES); 
    }
    v = 4.0 / dx2;
     MatSetValues(Ax, 1, &Ii, 1, &Ii, &v, ADD_VALUES); 
  }
   MatAssemblyBegin(Ax, MAT_FINAL_ASSEMBLY); 
   MatAssemblyEnd(Ax, MAT_FINAL_ASSEMBLY); 
  
   VecCreate(comm, &x); 
   VecSetSizes(x, BLOCKCOLS, BLOCKCOLS); 
   VecSetFromOptions(x); 
   VecDuplicate(x, &b); 

  PetscCalloc1(BLOCKCOLS, &idxm);
  for (j = 0; j < BLOCKCOLS; j++) { idxm[j] = j;} //for index

   KSPCreate(comm, &kspb); 
  //each process solve its own portion local systems 
  for (int ii = 0; ii < BLOCKROWS; ii++) { 
    MatShift(Ax,eg[ii + rank * BLOCKROWS]); //shift
     MatSetOption(Ax, MAT_SYMMETRIC, PETSC_TRUE); 
     KSPSetOperators(kspb, Ax, Ax); 
     KSPGetPC(kspb, &pc); 
     KSPSetType(kspb, KSPPREONLY); 
     PCSetType(pc, PCLU); //use LU-based direct solver 
    // PCSetType(pc, PCMG);  //use iterative solver if needed
     KSPSetFromOptions(kspb); 
    //exact row of aXlocal as Vec rhs
    VecSetValues(b, BLOCKCOLS, idxm, &(aXlocal[ii * BLOCKCOLS]), INSERT_VALUES);
    VecAssemblyBegin(b); VecAssemblyEnd(b);
     KSPSolve(kspb, b, x); //solve step-b systems
    VecGetValues(x, BLOCKCOLS, idxm, &(aXlocal[ii * BLOCKCOLS]));
    MatShift(Ax, -eg[ii + rank * BLOCKROWS]); //shift back for next solve
  }  
  PetscTime(&tb2);
//step (c)===================================================

  PetscTime(&tc1);
  if (rank == 0){ aVX = malloc(nt * nx*sizeof(PetscScalar)); }  
  //take part of V as local V for local matrix multiplication, no communication
  NPROWS = 1; NPCOLS = size;
  int BLOCKROWS2 = nt / NPROWS, BLOCKCOLS2 = nt / NPCOLS; /* number of cols in _block_ */
  aVlocal = malloc(BLOCKROWS2 * BLOCKCOLS2*sizeof(PetscScalar));
  PetscCalloc1(nt, &idxm);
  for (j = 0; j < BLOCKROWS2; j++){idxm[j] = j;}
  PetscCalloc1(BLOCKCOLS2, &idxn);
  for (j = 0; j < BLOCKCOLS2; j++){
    idxn[j] = j + rank * BLOCKCOLS2;
  }
  MatGetValues(V, BLOCKROWS2, idxm, BLOCKCOLS2, idxn, aVlocal); //extract local V by columns
  MatCreateSeqDense(comm, BLOCKCOLS2, BLOCKROWS2, aVlocal, &Vloc); //column-wise
  MatTranspose(Vloc, MAT_INPLACE_MATRIX, &Vloc);
  MatCreateSeqDense(comm, BLOCKCOLS, BLOCKROWS, aXlocal, &Xloc); //column-wise
  MatTranspose(Xloc, MAT_INPLACE_MATRIX, &Xloc);
  MatMatMult(Vloc, Xloc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Xloc);
  MatDenseGetArray(Xloc, &aVXlocal);
  //sum the submatrices of different processes in to aVX
  MPI_Reduce(aVXlocal, aVX, BLOCKROWS2 * BLOCKCOLS, MPI_C_DOUBLE_COMPLEX, MPI_SUM, 0, PETSC_COMM_WORLD);
  PetscTime(&tc2);  

  PetscTime(&t2); //finish all computation
  //clear up memory
  free(aVlocal);  free(aXlocal); free(mem1);free(mem2); free(eg);
   MatDestroy(&Xloc); MatDestroy(&Vloc); MatDestroy(&Floc); MatDestroy(&Ax); 
   VecDestroy(&x); VecDestroy(&b);  KSPDestroy(&kspb); 
  //compute errors in root process
  if (rank == 0){
    MatCreateSeqDense(comm, nt, nx, aVX, &VX); //column-wise
     MatAssemblyBegin(VX, MAT_FINAL_ASSEMBLY); 
     MatAssemblyEnd(VX, MAT_FINAL_ASSEMBLY); 
     MatAXPY(VX, -1.0, U, SAME_NONZERO_PATTERN); //Check the solution error
    VecCreateSeq(comm, nt, &max);
    MatGetRowMaxAbs(VX, max, NULL);
    VecMax(max, NULL, &err); //compute maximum error
    VecDestroy(&max); MatDestroy(&V);  MatDestroy(&X);   MatDestroy(&VX);  MatDestroy(&U);     
  }
   
  PetscPrintf(PETSC_COMM_WORLD, "%2d\t&(%3d,%3d,%3d) \t& %1.2e \t&%d \t& %1.2f\t&(%1.2f,& \t %1.2f,& \t %1.2f,& \t %1.2f)\\\\\n",
                     size, nx0, nx0, nt, (double)err,  eigiter, t2 - t1, teig2 - teig1, ta2 - ta1, tb2 - tb1, tc2 - tc1);
  PetscFinalize();
  return 0;
}
