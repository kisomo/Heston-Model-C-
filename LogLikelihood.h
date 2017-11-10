// Log Likelihood function for Heston and Nandi Model

// Returns the sum of a vector's elements
double VecSum(vector<double> x) {
	int n = x.size();
	double Sum = 0.0;
	for (int i=0; i<=n-1; i++)
		Sum += x[i];
	return Sum;
}

// Returns the sample variance
double VecVar(vector<double> x) {
	int n = x.size();
	double mean = VecSum(x) / n;
	double sumV = 0.0;
	for (int i=0; i<=n-1; i++)
		sumV += (x[i] - mean)*(x[i] - mean);
	return sumV / (n-1);
}

// Returns the log-likelihood based on S&P500 returns
double LogLike(vector<double> B, double r) {
	// Open the S&P500 Levels
	ifstream inPrices;
	inPrices.open("SP500.txt");
	vector<double> Price(0, 0.0);
	double n1;
	int i=0;
	while (!inPrices.eof()) {
		inPrices >> n1;
		Price.push_back(n1);
	i++ ;
	}
	int N = Price.size();

	// Calculate S&P500 returns
	vector<double> ret(N-1);
	for (i=0; i<=N-2; i++) {
		ret[i] = log(Price[i]/Price[i+1]);
	}
	double Variance = VecVar(ret);
	vector<double> h(N-1);
	vector<double> Z(N-1);
	vector<double> L(N-1);

	// Construct GARCH(1,1) process by working back in time
	h[N-2] = Variance;
    Z[N-2] = (ret[N-2] - r - B[4]*h[N-2]) / sqrt(h[N-2]);
    L[N-2] = -log(h[N-2]) - pow(ret[N-2], 2) / h[N-2];
	for (i=N-3; i>=0; i--) {
		h[i] = B[0] + B[2]*h[i+1] + B[1]*pow(Z[i+1]-B[3]*sqrt(h[i+1]), 2);
		Z[i] = (ret[i] - r - B[4]*h[i]) / sqrt(h[i]);
		L[i] = - log(h[i]) - pow(ret[i], 2) / h[i];
	}
	double LogL = VecSum(L);
	if ((B[0]<0) | (B[1]<0) |(B[2]<0) |(B[3]<0) |(B[4]<0)) // | (B[2]+B[1]*pow(B[3],2)>=1))
	   return 1e50;
	else
	return -LogL;						// Minimize -Log-Like(Beta)
}

