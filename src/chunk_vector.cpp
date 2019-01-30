#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector c_chunk_sums(NumericVector vec, IntegerVector sizes) {
		int nchunks = sizes.length();

		NumericVector val(nchunks);

		double* pvec = vec.begin();
		int*    psizes = sizes.begin();
		double* pval = val.begin();

		for (int i = 0; i < nchunks; i++) {
			double sum = 0;
			for (int j = 0; j < psizes[i]; j++) {
				sum += *pvec;
				++pvec;
			}
			pval[i] = sum;
		}
		return val;
}

// [[Rcpp::export]]
NumericVector c_chunk_maxs(NumericVector vec,IntegerVector sizes) {
		int nchunks = sizes.length();
		NumericVector val(nchunks);

		double* pvec = vec.begin();
		int*    psizes = sizes.begin();
		double* pval = val.begin();

		for (int i = 0; i < nchunks; i++) {
			if (psizes[i]>0) {
				double max = *pvec;
				++pvec;
				for (int j = 1; j < psizes[i]; j++) {
					if (*pvec>max) {
						max = *pvec;
					}
					++pvec;
				}
				pval[i] = max;
			} else {
				pval[i] = 0;
			}
		}
		return val;
}

// [[Rcpp::export]]
IntegerVector c_which_chunk_maxs(NumericVector vec,IntegerVector sizes) {
		int nchunks = sizes.length();

		IntegerVector val(nchunks);

		double* pvec = vec.begin();
		int*    psizes = sizes.begin();
		int*    pval = val.begin();

		int counter = 0;
		for (int i = 0; i < nchunks; i++) {

			if (psizes[i]>0) {
				counter++;
				double max = *pvec;
				int maxind = counter;
				pvec++;

				for (int j = 1; j < psizes[i]; j++) {
					counter++;
					if (*pvec > max) {
						max = *pvec;
						maxind = counter;
					}
					pvec++;
				}
				pval[i] = maxind;
			} else {
				pval[i] = 0;
			}

		}

		return val;
}



