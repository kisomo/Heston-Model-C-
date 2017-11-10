// Heston and Nandi Integration Files

// Trapezoidal Rule passing two vectors
double trapz(vector<double> X, vector<double> Y)
{
    int n = X.size();
	double sum = 0;
	for (int i=1; i<=n-1; i++)
		sum += 0.5*(X[i] - X[i-1])*(Y[i-1] + Y[i]);
	return sum;
}

// HNC_f returns the real part of the Heston & Nandi integral
double HNC_f(complex<double> phi, double alpha, double beta, double gamma, double omega, 
			 double lambda, double V, double S, double K, double r, int T, int FuncNum)
{
	int t;
	vector<complex<double> > A(T+1);
	vector<complex<double> > B(T+1);
	complex<double> zero(0.0, 0.0);
	complex<double>    i(0.0, 1.0);
	A[T] = zero;
	B[T] = zero;
	complex<double> one(1.0, 0.0);
	for (t=T-1; t>=0; t--) {
		if (FuncNum==1) {
			A[t] = A[t+1] + (phi+one)*r + B[t+1]*omega - 0.5*log(1.0 - 2.0*alpha*B[t+1]);
			B[t] = (phi+one)*(lambda+gamma) - 0.5*gamma*gamma + beta*B[t+1] 
				 + 0.5*(phi+one-gamma)*(phi+one-gamma)/(1.0 - 2.0*alpha*B[t+1]);
			}
		else  {
			A[t] = A[t+1] + (phi)*r + B[t+1]*omega - 0.5*log(1.0 - 2.0*alpha*B[t+1]);
			B[t] = phi*(lambda+gamma) - 0.5*gamma*gamma + beta*B[t+1] 
				 + 0.5*(phi-gamma)*(phi-gamma)/(1.0 - 2.0*alpha*B[t+1]);
		}
	}
	if (FuncNum==1)
        return real(pow(K, -phi)*pow(S, phi+one)*exp(A[0] + B[0]*V) / phi);
	else
		return real(pow(K, -phi)*pow(S, phi    )*exp(A[0] + B[0]*V) / phi);
}

// Returns the Heston and Nandi option price
double HNC(double alpha, double beta, double gamma, double omega, double lambda, 
		   double V, double S, double K, double r, int T, int PutCall)
{
	const double pi = 4.0*atan(1.0);
	double High   = 100;
	double Increment = 0.25;
	int NumPoints = int(High/Increment);
	vector<double> X(NumPoints);
	vector<double> Y1(NumPoints);
	vector<double> Y2(NumPoints);
	complex<double>   i(0.0, 1.0);
	complex<double> phi;
	for (int j=0; j<=NumPoints-1; j++) {
		if (j==0) X[j] = 0.0000001; else
		X[j]  = double(j)*Increment;
		phi   = X[j]*i;
		Y1[j] = HNC_f(phi, alpha, beta, gamma, omega, lambda, V, S, K, r, T, 1);
		Y2[j] = HNC_f(phi, alpha, beta, gamma, omega, lambda, V, S, K, r, T, 2);
	}
	double int1 = trapz(X, Y1);
	double int2 = trapz(X, Y2);
	double P1 = 1/2 + exp(-r*T)*int1/S/pi;
	double P2 = 1/2 + int2/pi;
	if (P1<0) P1=0;
	if (P1>1) P1=1;
	if (P2<0) P2=0;
	if (P2>1) P2=1;
	double Call = S/2 + exp(-r*T)*int1/pi - K*exp(-r*T)*(0.5 + int2/pi);
    double Put  = Call + K*exp(-r*T) - S;
	if (PutCall==1) return Call;
	else return Put;
	return 0;
}
