#include <getPoly.h>

using namespace std;


// Gets all potential combinations of the values in mymat
// [[Rcpp::export]]
arma::umat get_all_combins(arma::umat mymat, unsigned int polyorder){
	unsigned int d = mymat.n_cols;
	unsigned int total_combins = choose(d+polyorder,d)-1;
	unsigned int num_rows = mymat.n_rows;
	arma::uvec OneToD = arma::linspace<arma::uvec>(1, d, d);
	
	arma::umat final_mat(d,total_combins);
	
	arma::uvec uniqs;
	unsigned int total_uniqs;
	
	unsigned int index_final = 0;
	
	for (unsigned int j = 0; j<num_rows; j++){
		uniqs = arma::sort(arma::unique(mymat.row(j))).t();
		total_uniqs = uniqs.n_rows;
		
		arma::uvec num_per_uniqs(total_uniqs);
		for (unsigned int k = 0; k<total_uniqs; k++){
			arma::uvec temp_found = arma::find(mymat.row(j)==uniqs[k]);
			num_per_uniqs(k) = temp_found.n_rows;
		}
		
		unsigned int curr_ncol_double = static_cast<double>(choose(d,num_per_uniqs(0)));
		if (total_uniqs>2){
			unsigned int total_uniqs_removed = num_per_uniqs(0);
			for (unsigned int k = 1; k < (total_uniqs - 1); k++) {
				curr_ncol_double = curr_ncol_double*static_cast<double>(choose(d-total_uniqs_removed,num_per_uniqs(k)));
				total_uniqs_removed += num_per_uniqs(k);
			}
		}
		unsigned int curr_ncol = static_cast<unsigned int>(curr_ncol_double);
		
		unsigned int num_done = num_per_uniqs(0);
		unsigned int row_curr_start = 0; //index after this next line
		unsigned int row_curr_end = num_per_uniqs(0) - 1; //index after this next line
		unsigned int col_curr_start = 0; //index after this next line
		unsigned int col_curr_end = choose(d,num_per_uniqs(0)) - 1;
		arma::umat curr_mat_inds(d,curr_ncol);
		curr_mat_inds(arma::linspace<arma::uvec>(row_curr_start, row_curr_end,row_curr_end - row_curr_start + 1),arma::linspace<arma::uvec>(col_curr_start, col_curr_end, col_curr_end - col_curr_start + 1)) = combn(OneToD, num_per_uniqs(0));
		
		if (total_uniqs>1){
			for (unsigned int k = 2; k<=total_uniqs; k++){
				arma::umat diff_to_comb(d-num_done,col_curr_end + 1);
				arma::umat y = curr_mat_inds.rows(arma::linspace<arma::uvec>(0, row_curr_end,row_curr_end + 1));
				
				for (unsigned int jj = 0; jj < (col_curr_end + 1); jj++){
					diff_to_comb.col(jj) = setdiff(OneToD,y.col(jj));
				}
				
				if (k==total_uniqs){
					row_curr_start = row_curr_end + 1;
					row_curr_end = row_curr_start + num_per_uniqs(k-1) - 1;
					curr_mat_inds(arma::linspace<arma::uvec>(row_curr_start, row_curr_end,row_curr_end-row_curr_start + 1), arma::linspace<arma::uvec>(0, diff_to_comb.n_cols - 1, diff_to_comb.n_cols)) = diff_to_comb;
				} else{
					arma::uvec OneToDiffRows = arma::linspace<arma::uvec>(1, diff_to_comb.n_rows, diff_to_comb.n_rows);
					arma::umat combins_for_smaller_set = combn(OneToDiffRows, num_per_uniqs(k-1));
					
					arma::uvec rows_old = arma::linspace<arma::uvec>(row_curr_start, row_curr_end,row_curr_end-row_curr_start + 1);
					arma::uvec rows_new = arma::linspace<arma::uvec>(row_curr_end + 1, row_curr_end + num_per_uniqs(k-1),num_per_uniqs(k-1));
					arma::uvec cols_old = arma::linspace<arma::uvec>(col_curr_start, col_curr_end, col_curr_end - col_curr_start + 1);
					
					curr_mat_inds.submat(rows_new, arma::linspace<arma::uvec>(col_curr_start, col_curr_end, col_curr_end - col_curr_start + 1) ) = diff_to_comb.rows(combins_for_smaller_set.col(0)-1);
					
					for (unsigned int zz = 1; zz < combins_for_smaller_set.n_cols; zz++) {
						col_curr_start = col_curr_end + 1;
						col_curr_end = col_curr_start + choose(d,num_per_uniqs(k-2)) - 1;
						curr_mat_inds(rows_old, arma::linspace<arma::uvec>(col_curr_start, col_curr_end, col_curr_end - col_curr_start + 1) ) = curr_mat_inds.submat(rows_old,cols_old);
						curr_mat_inds(rows_new, arma::linspace<arma::uvec>(col_curr_start, col_curr_end, col_curr_end - col_curr_start + 1) ) = diff_to_comb.rows(combins_for_smaller_set.col(zz)-1);
					}
					row_curr_start = row_curr_end + 1;
					row_curr_end = row_curr_start + num_per_uniqs(k-1) - 1;
					
				}
				num_done += num_per_uniqs(k-1);
			}
		}
		
		arma::uvec curr_sorted = arma::sort(mymat.row(j)).t();
		for (unsigned int z = 0; z < curr_mat_inds.n_cols; z++) {
			for (unsigned int k = 0; k < total_uniqs; k++) {
				arma::uvec temp_inds = curr_mat_inds.col(z);
				temp_inds = temp_inds.rows(find(curr_sorted==uniqs(k)));
				for (unsigned int bb = 0; bb < temp_inds.n_rows; bb++) {
					final_mat(temp_inds(bb)-1,index_final) = uniqs(k);
				}
			}
			index_final = index_final + 1;
		}
	}
	return(final_mat);
	
}

// Turns the matrix of all combinations of polynomial orders into the matrix of covariates using the samples and derivatives. Does this WITHOUT using the R package partitions which cannot be used on Linux.
// [[Rcpp::export]]
arma::mat getPoly_withoutpackage(arma::mat samples, arma::mat derivatives, arma::mat combinations){
	
	unsigned long N = samples.n_rows;
	unsigned long d = samples.n_cols;
	
	derivatives = -0.5*derivatives;
	
	unsigned long total_monos = combinations.n_cols;
	
	unsigned int k, j, jj;
	arma::vec term;
	
	arma::mat X(N,total_monos);
	
	for (k = 0; k<total_monos; k++){
		term = arma::zeros<arma::vec>(N);
		// First order derivatives wrt parameter j
		for (j = 0; j<d; j++){
			arma::vec myprod = combinations(j,k)*pow(samples.col(j),combinations(j,k) - 1.0) % derivatives.col(j);
			for (jj=0; jj<d; jj++){
				if (jj!=j){myprod = myprod % pow(samples.col(jj),combinations(jj,k));}
			}
			term = term + myprod;
		}
		// Second order derivatives wrt parameter j
		for (j = 0; j<d; j++){
			arma::vec myprod = std::max(combinations(j,k)*(combinations(j,k)-1.0),0.0) * pow(samples.col(j),combinations(j,k) - 2.0);
			for (jj=0; jj<d; jj++){
				if (jj!=j){myprod = myprod % pow(samples.col(jj),combinations(jj,k));}
			}
			term = term + -0.5*myprod;
		}
		X.col(k) = term;
	}
	
	return X;
}

// Turns the matrix of all combinations of polynomial orders into the matrix of covariates using the samples and derivatives. Does this WITH the R package partitions which cannot be used on Linux.
// [[Rcpp::export]]
arma::mat getPoly_withpackage(arma::mat samples, arma::mat derivatives, Rcpp::List combinations){ 
	
	unsigned long N = samples.n_rows;
	unsigned long d = samples.n_cols;
	
	derivatives = -0.5*derivatives;
	
	unsigned long polyorder = combinations.size();
	
	unsigned int mypol, k, j, jj;
	
	
	unsigned long num_monos, total_monos;
	arma::vec term;
	
	total_monos = 0;
	for (mypol = 0; mypol<polyorder; mypol++){
		arma::Mat<int> alpha = combinations[mypol];
		total_monos += alpha.n_cols;
	}
	
	arma::mat X(N,total_monos);
	unsigned long ii = 0;
	
	for (mypol = 0; mypol<polyorder; mypol++){
		arma::Mat<int> alpha = combinations[mypol];
		num_monos = alpha.n_cols;
		for (k = 0; k<num_monos; k++){
			term = arma::zeros<arma::vec>(N);
			// First order derivatives wrt parameter j
			for (j = 0; j<d; j++){
				arma::vec myprod = std::max(alpha(j,k),0)*pow(samples.col(j),alpha(j,k) - 1) % derivatives.col(j);
				for (jj=0; jj<d; jj++){
					if (jj!=j){myprod = myprod % pow(samples.col(jj),alpha(jj,k));}
				}
				term = term + myprod;
			}
			// Second order derivatives wrt parameter j
			for (j = 0; j<d; j++){
				arma::vec myprod = std::max(alpha(j,k)*(alpha(j,k)-1),0) * pow(samples.col(j),alpha(j,k) - 2);
				for (jj=0; jj<d; jj++){
					if (jj!=j){myprod = myprod % pow(samples.col(jj),alpha(jj,k));}
				}
				term = term + -0.5*myprod;
			}
			X.col(ii) = term;
			ii++;
		}
	}
	
	return X;
}


// The number of ways you can choose k values from a list of n without taking order into account (i.e. the number of combinations)
unsigned int choose(unsigned int n, unsigned int k){
	
	double factorial_n_over_factorial_k = 1.0;
	double curr = static_cast<double>(k);
	for (unsigned int z = (k+1); z<=n; z++){
		curr += 1.0;
		factorial_n_over_factorial_k *= curr;
	}
	
	double factorial_n_minus_k = 1.0;
	curr = 0.0;
	for (unsigned int z = 0; z<(n-k); z++){
		curr += 1.0;
		factorial_n_minus_k *= curr;
	}
	
	return (static_cast<unsigned int>(factorial_n_over_factorial_k/factorial_n_minus_k));
}

// Gets the matrix of combinations of m values from the vector x (which is length > m)
arma::umat combn(arma::uvec x, unsigned int m){ // Converted from base R's combn function
	unsigned int n = x.n_rows;
	unsigned int e = 0;
	unsigned int h = m;
	arma::uvec a = arma::linspace<arma::uvec>(1, m, m);
	arma::uvec r = x.rows(a-1);
	unsigned int count = choose(n, m);
	arma::umat out(m,count);
	for (unsigned int z = 0; z<count; z++){
		out.col(z) = r;
	}
	unsigned int i = 2;
	unsigned int nmmp1 = n - m + 1;
	while (a(0) != nmmp1) {
		if (e < n - h) {
			h = 1;
			e = a(m - 1);
			unsigned int j = 1;
			a(m - h + j - 1) = e + j;
			r = a; // Just to set the dimension
			for (unsigned int kk = 0; kk<r.n_rows; kk++){
				r(kk) = x(a(kk)-1);
			}
		}
		else {
			e = a(m - h - 1);
			h = h + 1;
			arma::uvec j = arma::linspace<arma::uvec>(1, h, h); // 1:h
			
			for (unsigned int kk = 0; kk<h; kk++){
				a(m - h + j(kk) - 1) = e + j(kk);
			}
			
			r = a; // Just to set the dimension
			for (unsigned int kk = 0; kk<r.n_rows; kk++){
				r(kk) = x(a(kk)-1);
			}
		}
		out.col(i-1) = r;
		i = i + 1;
	}
	return(out);
}

// Gets the vector of elements that are in x but not y
arma::uvec setdiff(arma::uvec x, arma::uvec y){
	for (unsigned int j = 0; j < y.n_rows; j++) {
		arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y(j)));
		x.shed_row(q1);
	}
	return x;
}