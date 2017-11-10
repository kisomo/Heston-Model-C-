// Nelder Mead minimization algorithm

// Function to calculate the mean value of
// a set of n vectors each of dimension n
// namely a (n x n) matrix

vector<double> VMean(vector<vector<double> > X, int n)
{
    vector<double> meanX(n);
	for (int i=0; i<=n-1; i++)
	{
		meanX[i]=0.0;
		for (int j=0; j<=n-1; j++)
            {
				meanX[i] += X[i][j];
			}
		meanX[i] = meanX[i] / n;
	}
return meanX;
}

// Function to add two vectors together
vector<double> VAdd(vector<double> x, vector<double> y)
{
    int n = x.size();
    vector<double> z(n, 0.0);
	for (int j=0; j<=n-1; j++)
		z[j] = x[j] + y[j];
return z;
}

// Function to subtract two vectors
vector<double> VSub(vector<double> x, vector<double> y)
{
    int n = x.size();
    vector<double> z(n, 0.0);
	for (int j=0; j<=n-1; j++)
		z[j] = x[j] - y[j];
return z;
}

// Function to multiply a vector by a constant
vector<double> VMult(vector<double> x, double a)
{
    int n = x.size();
    vector<double> z(n, 0.0);
	for (int j=0; j<=n-1; j++)
		z[j] = a*x[j];
return z;
}
// Nelder Mead Algorithm
vector<double> NelderMead(double (*f)(vector<double>, double r), 
						  int N, double NumIters, double MaxIters, 
						  double Tolerance, vector<vector<double> > x, double r)
{
int i,j;

// Value of the function at the vertices
vector<vector<double> > F(N+1, vector<double>(2));

// Step 0.  Ordering and Best and Worst points
// Order according to the functional values, compute the best and worst points
step0:
NumIters = NumIters + 1;
for (j=0; j<=N; j++){
    vector<double> z(N, 0.0);								// Create vector to contain
	for (i=0; i<=N-1; i++)
		z[i] = x[i][j];
	F[j][0] = f(z, r);				       					// Function values
	F[j][1] = j;											// Original index positions
}
sort(F.begin(), F.end());

// New vertices order first N best initial vectors and
// last (N+1)st vertice is the worst vector
// y is the matrix of vertices, ordered so that the worst vertice is last
vector<vector<double> > y(N, vector<double>(N+1));
for (j=0; j<=N; j++) {
	for (i=0; i<=N-1; i++) {
		y[i][j] = x[i][F[j][1]];
	}
}

//  First best vector y(1) and function value f1
vector<double> x1(N, 0.0); for (i=0; i<=N-1; i++) x1[i] = y[i][0];
double f1 = f(x1, r);

// Last best vector y(N) and function value fn
vector<double> xn(N, 0.0); for (i=0; i<=N-1; i++) xn[i] = y[i][N-1];
double fn = f(xn, r);

// Worst vector y(N+1) and function value fn1
vector<double> xn1(N, 0.0); for (i=0; i<=N-1; i++) xn1[i] = y[i][N];
double fn1 = f(xn1, r);

// z is the first N vectors from y, excludes the worst y(N+1)
vector<vector<double> > z(N, vector<double>(N));
for (j=0; j<=N-1; j++) {
	for (i=0; i<=N-1; i++)
		z[i][j] = y[i][j];
}

// Mean of best N values and function value fm
vector<double> xm(N, 0.0); xm = VMean(z,N);
double fm = f(xm, r);

// Reflection point xr and function fr
vector<double> xr(N, 0.0); xr = VSub(VAdd(xm, xm), xn1);
double fr = f(xr, r);

// Expansion point xe and function fe
vector<double> xe(N, 0.0); xe = VSub(VAdd(xr, xr), xm);
double fe = f(xe, r);

// Outside contraction point and function foc
vector<double> xoc(N, 0.0);	xoc = VAdd(VMult(xr, 0.5), VMult(xm, 0.5));
double foc = f(xoc, r);

// Inside contraction point and function foc
vector<double> xic(N, 0.0);	xic = VAdd(VMult(xm, 0.5), VMult(xn1, 0.5));
double fic = f(xic, r);

while ((NumIters <= MaxIters) && (abs(f1-fn1) >= Tolerance))
{
// Step 1. Reflection Rule
if ((f1<=fr) && (fr<fn)) {
    for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
		for (i=0; i<=N-1; i++)	x[i][N] = xr[i];
	goto step0;
}

// Step 2.  Expansion Rule
if (fr<f1) {
	for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
	if (fe<fr)
		for (i=0; i<=N-1; i++)	x[i][N] = xe[i];
	else
		for (i=0; i<=N-1; i++)	x[i][N] = xr[i];
	goto step0;
}

// Step 3.  Outside contraction Rule
if ((fn<=fr) && (fr<fn1) && (foc<=fr)) {
	for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
		for (i=0; i<=N-1; i++)	x[i][N] = xoc[i];
	goto step0;
}

// Step 4.  Inside contraction Rule
if ((fr>=fn1) && (fic<fn1)) {
	for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
		for (i=0; i<=N-1; i++)	x[i][N] = xic[i];
	goto step0;
}

// Step 5. Shrink Step
for (i=0; i<=N-1; i++)
	x[i][0] = y[i][0];
for (i=0; i<=N-1; i++) {
    for (j=1; j<=N; j++)
		x[i][j] = 0.5*(y[i][j] + x[i][0]);
}
	goto step0;
}

// Output component
// Return N parameter values, value of objective function, and number of iterations
vector<double> out(N+2);
for (i=0; i<=N-1; i++)
	out[i] = x1[i];
	out[N] = f1;
	out[N+1] = NumIters;
return out;
}
