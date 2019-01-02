#include <Rcpp.h>
using namespace Rcpp;

// Assume vec would be a matrix with nrow rows and ncol cols
// Compute for each column the max value and set all values
// of the column to this max

int which_max_with_cond(NumericVector vec, LogicalVector cond) {
  int res=-1;
  double max = -DBL_MAX;
  for (int i=0;i<vec.size();i++) {
    if (cond[i]) {
      if (vec[i]>max) {
        max = vec[i];
        res = i;
      }
    }
  }
  return res;
}


double max_with_cond(NumericVector vec, LogicalVector cond) {
  double max = -DBL_MAX;
  for (int i=0;i<vec.size();i++) {
    if (cond[i]) {
      if (vec[i]>max) max = vec[i];
    }
  }
  return max;
}


double min_with_cond(NumericVector vec, LogicalVector cond) {
  double min = DBL_MAX;
  for (int i=0;i<vec.size();i++) {
    if (cond[i]) {
      if (vec[i]<min) min = vec[i];
    }
  }
  return min;
}


// [[Rcpp::export]]
NumericVector vecLongColMaxs(NumericVector vec,int nrow) {
  int n = vec.size();
  int ncol = n / nrow;

  NumericVector res(n);

  int index = 0;
  for (int col =0;col < ncol; col++) {
    double max_val = vec[index];
    int start_index = index;
    for (int row=1; row < nrow; row++) {
      index++;
      if (vec[index] > max_val) max_val =  vec[index];
    }
    index = start_index;
    for (int row=0; row < nrow; row++) {
      res[index] = max_val;
      index++;
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericVector vecLongRowMaxs(NumericVector vec,int nrow) {
  int n = vec.size();
  int ncol = n / nrow;

  NumericVector res(n);

  int index = 0;
  for (int row =0;row < nrow; row++) {
    index = row;
    double max_val = vec[index];
    for (int col=1; col < ncol; col++) {
      index=index+nrow;
      if (vec[index] > max_val) max_val =  vec[index];
    }
    index = row;
    for (int col=0; col< ncol; col++) {
      res[index] = max_val;
      index=index+nrow;
    }
  }
  return res;
}


// [[Rcpp::export]]
NumericVector mat_times_vec_rows(NumericMatrix mat, NumericVector vec, IntegerVector vec_rows) {
  int ncol = mat.ncol();
  int c = 0;
  NumericMatrix::Column mcol = mat( _ , c);
  NumericVector res = mcol *vec[vec_rows[c]];

  for (c=1; c<ncol;c++) {
    mcol = mat( _ , c);
    res += mcol * vec[vec_rows[c]];
  }
  return res;
}

// [[Rcpp::export]]
NumericVector mat_times_vec(NumericMatrix mat, NumericVector vec) {
  int ncol = mat.ncol();
  int c = 0;
  NumericMatrix::Column mcol = mat( _ , c);
  NumericVector res = mcol *vec;

  for (c=1; c<ncol;c++) {
    mcol = mat( _ , c);
    res += mcol * vec;
  }
  return res;
}


// [[Rcpp::export]]
IntegerVector cpp_capped_rne_find_actions(
    double U,double v1, double v2,
    NumericVector U_hat, NumericVector v1_hat, NumericVector v2_hat,
    LogicalVector IC_holds,
    NumericVector next_r1, NumericVector next_r2,
    NumericMatrix trans_mat, IntegerVector dest_rows,
    String tie_breaking,
    double tol=1e-12)
{
  // Pick dynamic equilibrium actions
  // using the specified tie.breaking rule
  double tb_const = 1;
  NumericVector tb;
  if (tie_breaking=="slack") {
    NumericVector slack = U_hat - (v1_hat + v2_hat);
    tb = slack;
  } else if (tie_breaking=="last") {
    tb = seq_along(U_hat);
  } else if (tie_breaking=="first") {
    tb = seq(U_hat.size(),0);
  } else if (tie_breaking=="equal_r") {
    NumericVector next_r_diff = abs(next_r1-next_r2);
    tb = mat_times_vec_rows(trans_mat,next_r_diff,dest_rows);
    tb_const = min(tb) + max(tb)-min(tb);
  } else {
    tb = runif(U_hat.size());
  }

  // Base Indexes on 1
  int ae = 1L+which_max_with_cond(tb_const+tb,(abs(U_hat-U)<tol) & IC_holds);
  int a1 = 1L+which_max_with_cond(tb_const+tb,(abs(v1_hat-v1)<tol) & IC_holds);
  int a2 = 1L+which_max_with_cond(tb_const+tb,(abs(v2_hat-v2)<tol) & IC_holds);
  return IntegerVector::create(ae,a1,a2);
}


// [[Rcpp::export]]
DataFrame cpp_capped_rne_iterations(int T,
    DataFrame sdf, DataFrame rne,List transmats,
    double delta, double rho, double beta1,
    String tie_breaking,double tol=1e-10, int debug_row=-1)
{

  NumericVector next_U = rne["U"];
  NumericVector next_r1 = rne["r1"];
  NumericVector next_r2 = rne["r2"];
  NumericVector next_v1 = rne["v1"];
  NumericVector next_v2 = rne["v2"];

  NumericVector next_v1_r1 = (1-rho)*next_v1 + rho*next_r1;
  NumericVector next_v2_r2 = (1-rho)*next_v2 + rho*next_r2;


  int nx = next_U.size();



  IntegerVector vec_na1 = sdf["na1"];
  IntegerVector vec_na2 = sdf["na2"];
  List li_pi1 = sdf["pi1"];
  List li_pi2 = sdf["pi2"];

  StringVector vec_x = sdf["x"];
  DataFrame res;


  NumericVector res_U(nx);
  NumericVector res_r1(nx);
  NumericVector res_r2(nx);
  NumericVector res_v1(nx);
  NumericVector res_v2(nx);

  IntegerVector res_ae(nx);
  IntegerVector res_a1(nx);
  IntegerVector res_a2(nx);


  // Obtaining function to get action profiles
  //Environment pkg = Environment::namespace_env("RelationalContractsCpp");
  //Function find_a_fun = pkg["capped.rne.find.actions"];


  for (int iter=0; iter < T; iter++) {
    bool last_iter = (iter == T-1);
    // Loop over states
    for (int row=0; row < nx; row++) {
      String x = vec_x[row];
      //int na1 = vec_na1[row];
      int na2 = vec_na2[row];

      NumericMatrix trans_mat = transmats[row];
      NumericVector pi1 = li_pi1[row];
      NumericVector pi2 = li_pi2[row];

      StringVector xd = colnames(trans_mat);

      // Note: match returns 1-based index
      // but we need 0-based index
      IntegerVector dest_rows = match(xd, vec_x)-1;


      //# Include code to compute U v and r for the current state
      //U.hat = (1-delta)*(sdf$pi1[[srow]] + sdf$pi2[[srow]]) +
      //  delta * (trans.mat %*% rne$U[dest.rows])


      NumericVector U_hat = (1-delta)*(pi1+pi2) +
        delta* mat_times_vec_rows(trans_mat,next_U, dest_rows);


      //# "q-value" of punishment payoffs
      //q1.hat = (1-delta)*sdf$pi1[[srow]] +
      //  delta * (trans.mat %*% ( (1-rho)*rne$v1[dest.rows] + rho*rne$r1[dest.rows] ))
      NumericVector q1_hat = (1-delta)*pi1 +
        delta*mat_times_vec_rows(trans_mat,next_v1_r1, dest_rows);

      NumericVector q2_hat = (1-delta)*pi2 +
        delta*mat_times_vec_rows(trans_mat,next_v2_r2, dest_rows);

      // # v1.hat is best reply q for player 1
      // # Note player 1 is col player
      //NumericMatrix q1_hat = matrix(q1.hat,na2, na1)
      //v1.hat.short = rowMaxs(q1.hat)
      //v1.hat = rep(v1.hat.short, times=na1)

      NumericVector v1_hat = vecLongRowMaxs(q1_hat,na2);

      // v2_hat is best reply for player 2 (row player)
      NumericVector v2_hat = vecLongColMaxs(q2_hat,na2);


      // Compute which action profiles are implementable
      LogicalVector IC_holds = (U_hat+tol >= v1_hat + v2_hat);

      if (row==debug_row) {
        res=DataFrame::create(
          Named("row_minus1")= debug_row,

          Named("pi1")=pi1,
          Named("pi2")=pi2,
          Named("U_cont")=mat_times_vec_rows(trans_mat,next_U, dest_rows),
          Named("next_U")=next_U,
	        Named("U_hat")=U_hat,
	        Named("v1_hat")=v1_hat,
	        Named("v2_hat")=v2_hat,
	        Named("IC_holds")=IC_holds,
	        Named("q1_hat")=q1_hat,
	        Named("q2_hat")=q2_hat
        );
  		  return res;
      }

      if (is_false(any(IC_holds))) {
        //std::string msg =  Failed in state "+x;
        stop("No pure strategy RNE exists.");
      }

      double U = max_with_cond(U_hat,IC_holds);
      double v1 = min_with_cond(v1_hat, IC_holds);
      double v2 = min_with_cond(v2_hat, IC_holds);

      double r1 = v1+beta1*(U-v1-v2);
      double r2 = U-r1;

      res_U[row] = U;
      res_v1[row] = v1;
      res_v2[row] = v2;
      res_r1[row] = r1;
      res_r2[row] = r2;

      if (last_iter) {
        IntegerVector a_row =
          cpp_capped_rne_find_actions(
          //find_a_fun(
            U,v1,v2,
            U_hat,v1_hat,v2_hat,IC_holds,next_r1,next_r2,
            trans_mat,dest_rows,
            tie_breaking, tol
          );
        res_ae[row] = a_row[0];
        res_a1[row] = a_row[1];
        res_a2[row] = a_row[2];
      }
    }
    if (last_iter) break;

    next_U = clone(res_U);
    next_v1 = clone(res_v1);
    next_v2 = clone(res_v2);
    next_r1 = clone(res_r1);
    next_r2 = clone(res_r2);

    next_v1_r1 = (1-rho)*next_v1 + rho*next_r1;
    next_v2_r2 = (1-rho)*next_v2 + rho*next_r2;

  }

  res =
	  DataFrame::create(
	    Named("x")= vec_x,
	    Named("r1")=res_r1,
	    Named("r2")=res_r2,
	    Named("U")=res_U,
	    Named("v1")=res_v1,
	    Named("v2")=res_v2,
	    Named("ae")=res_ae,
	    Named("a1")=res_a1,
	    Named("a2")=res_a2
		);
  return res;
}




/*
capped.rne.find.actions = function(U,v1,v2,U.hat,v1.hat,v2.hat, IC.holds, tie.breaking=c("slack","random","first","last")[1], tol=1e-12, next.r1=NULL, next.r2=NULL, trans.mat=NULL) {
  restore.point("capped.rne.find.actions")
  # Pick dynamic equilibrium actions
  # using the specified tie.breaking rule
  slack = U.hat - (v1.hat + v2.hat)
  const = 1
  if (tie.breaking=="slack") {
    tb = slack
  } else if (tie.breaking=="last") {
    tb = seq_len(NROW(U.hat))
  } else if (tie.breaking=="first") {
    tb = rev(seq_len(NROW(U.hat)))
  } else if (tie.breaking=="equal_r") {
    next_r_diff = abs(next.r1-next.r2)
    tb = as.vector(trans.mat %*% next_r_diff)
    const = min(tb) + max(tb)-min(tb)
  } else {
    tb = runif(NROW(U.hat))
  }

  ae = which.max((const+tb) * (abs(U.hat-U)<tol & IC.holds))
  a1 = which.max((const+tb) * (abs(v1.hat-v1)<tol & IC.holds))
  a2 = which.max((const+tb) * (abs(v2.hat-v2)<tol & IC.holds))
  c(ae,a1,a2)
}
*/

