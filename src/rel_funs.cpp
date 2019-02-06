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
NumericVector mat_times_vec_rows(NumericMatrix mat, NumericVector vec, IntegerVector vec_rows, int target_nrow) {
  int ncol = mat.ncol();
  int nrow = mat.nrow();

  // Special case: A zero row matrix encodes a terminal state
  // in which all actions lead back to itself
  if ( (ncol == 1) & (nrow==0)) {
    NumericVector res0(target_nrow,vec[vec_rows[0]]);
    return res0;
  }

  NumericVector res(nrow);
  for (int c=0; c<ncol;c++) {
    NumericMatrix::Column mcol = mat( _ , c);
    res += mcol * vec[vec_rows[c]];
  }
  return res;
}

// [[Rcpp::export]]
NumericVector mat_times_vec(NumericMatrix mat, NumericVector vec, int target_nrow=0) {
  int ncol = mat.ncol();
  int nrow = mat.nrow();
  NumericVector res(nrow);
  for (int c=0; c<ncol;c++) {
    NumericMatrix::Column mcol = mat( _ , c);
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
  if (tie_breaking=="equal_r") {
    NumericVector next_r_diff = -abs(next_r1-next_r2);
    tb = mat_times_vec_rows(trans_mat,next_r_diff,dest_rows, U_hat.size());
    tb_const = -min(tb) + max(tb)-min(tb);
  } else if (tie_breaking=="random") {
    tb = runif(U_hat.size());
  } else if (tie_breaking=="slack") {
    NumericVector slack = U_hat - (v1_hat + v2_hat);
    tb = slack;
  } else if (tie_breaking=="last") {
    tb = seq_along(U_hat);
  } else if (tie_breaking=="first") {
    tb = seq(U_hat.size(),0);
  } else if (tie_breaking=="max_r1") {
    tb = mat_times_vec_rows(trans_mat,next_r1,dest_rows,U_hat.size());
    tb_const = -min(tb) + max(tb)-min(tb);
  } else if (tie_breaking=="max_r2") {
    tb = mat_times_vec_rows(trans_mat,next_r2,dest_rows,U_hat.size());
    tb_const = -min(tb) + max(tb)-min(tb);
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
    String tie_breaking,double tol=1e-12, int debug_row=-1)
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
      int na = pi1.size();

      StringVector xd = colnames(trans_mat);

      // Note: match returns 1-based index
      // but we need 0-based index
      IntegerVector dest_rows = match(xd, vec_x)-1;


      //# Include code to compute U v and r for the current state
      //U.hat = (1-delta)*(sdf$pi1[[srow]] + sdf$pi2[[srow]]) +
      //  delta * (trans.mat %*% rne$U[dest.rows])


      NumericVector U_hat = (1-delta)*(pi1+pi2) +
        delta* mat_times_vec_rows(trans_mat,next_U, dest_rows,na);


      //# "q-value" of punishment payoffs
      //q1.hat = (1-delta)*sdf$pi1[[srow]] +
      //  delta * (trans.mat %*% ( (1-rho)*rne$v1[dest.rows] + rho*rne$r1[dest.rows] ))
      NumericVector q1_hat = (1-delta)*pi1 +
        delta*mat_times_vec_rows(trans_mat,next_v1_r1, dest_rows,na);

      NumericVector q2_hat = (1-delta)*pi2 +
        delta*mat_times_vec_rows(trans_mat,next_v2_r2, dest_rows,na);

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
          Named("U_cont")=mat_times_vec_rows(trans_mat,next_U, dest_rows,na),
          Named("next_U")=next_U,
	        Named("U_hat")=U_hat,
	        Named("U") = max_with_cond(U_hat,IC_holds),
	        Named("v1_hat")=v1_hat,
	        Named("v2_hat")=v2_hat,
	        Named("IC_holds")=IC_holds,
	        Named("q1_hat")=q1_hat,
	        Named("q2_hat")=q2_hat,
	        Named("stringsAsFactors")=false
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
	    Named("a2")=res_a2,
	    Named("stringsAsFactors")=false
		);
  return res;
}


// [[Rcpp::export]]
DataFrame cpp_capped_rne_multistage_iterations(int T,
    DataFrame sdf, DataFrame rne,List transmats,
    List static_rep_li,
    double delta, double rho, double beta1,
    String tie_breaking,double tol=1e-12, int debug_row=-1)
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

  NumericVector res_static_Pi(nx);
  NumericVector res_static_c1(nx);
  NumericVector res_static_c2(nx);

  IntegerVector res_s_ae(nx);
  IntegerVector res_s_a1(nx);
  IntegerVector res_s_a2(nx);

  IntegerVector res_d_ae(nx);
  IntegerVector res_d_a1(nx);
  IntegerVector res_d_a2(nx);



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

      int na = pi1.size();

      StringVector xd = colnames(trans_mat);

      // Note: match returns 1-based index
      // but we need 0-based index
      IntegerVector dest_rows = match(xd, vec_x)-1;


      //# Include code to compute U v and r for the current state
      //U.hat = (1-delta)*(sdf$pi1[[srow]] + sdf$pi2[[srow]]) +
      //  delta * (trans.mat %*% rne$U[dest.rows])


      NumericVector U_hat = (1-delta)*(pi1+pi2) +
        delta* mat_times_vec_rows(trans_mat,next_U, dest_rows,na);


      //# "q-value" of punishment payoffs
      //q1.hat = (1-delta)*sdf$pi1[[srow]] +
      //  delta * (trans.mat %*% ( (1-rho)*rne$v1[dest.rows] + rho*rne$r1[dest.rows] ))
      NumericVector q1_hat = (1-delta)*pi1 +
        delta*mat_times_vec_rows(trans_mat,next_v1_r1, dest_rows,na);

      NumericVector q2_hat = (1-delta)*pi2 +
        delta*mat_times_vec_rows(trans_mat,next_v2_r2, dest_rows,na);

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
          Named("U_cont")=mat_times_vec_rows(trans_mat,next_U, dest_rows,na),
          Named("next_U")=next_U,
	        Named("U_hat")=U_hat,
	        Named("v1_hat")=v1_hat,
	        Named("v2_hat")=v2_hat,
	        Named("IC_holds")=IC_holds,
	        Named("q1_hat")=q1_hat,
	        Named("q2_hat")=q2_hat,
	        Named("stringsAsFactors")=false
        );
  		  return res;
      }

      if (is_false(any(IC_holds))) {
        //std::string msg =  Failed in state "+x;
        stop("No pure strategy RNE exists. No incentive constraint is satisfied in the dynamic stage.");
      }

      double dU = max_with_cond(U_hat,IC_holds);
      double dv1 = min_with_cond(v1_hat, IC_holds);
      double dv2 = min_with_cond(v2_hat, IC_holds);


      // 2. Solve static stage

      // Action lists of the repeated game
      List s_li = static_rep_li[row];
      NumericMatrix ae_df = s_li["ae.df"];
      NumericMatrix a1_df = s_li["a1.df"];
      NumericMatrix a2_df = s_li["a2.df"];

      /*
      NumericVector ae_L = ae_df[,"L"];
      NumericVector a1_L = a1_df[,"L"];
      NumericVector a2_L = a2_df[,"L"];

      NumericVector ae_G = ae_df["G"];
      NumericVector a1_c1 = a1_df["c1"];
      NumericVector a2_c2 = a2_df["c2"];
      */
      NumericMatrix::Column ae_L = ae_df(_,4);
      NumericMatrix::Column a1_L = a1_df(_,4);
      NumericMatrix::Column a2_L = a2_df(_,4);

      NumericMatrix::Column ae_G = ae_df(_,1);
      NumericMatrix::Column a1_c1 = a1_df(_,2);
      NumericMatrix::Column a2_c2 = a2_df(_,3);

      // Available liquidity after static stage
      double L_av = 1/(1-delta)*(dU-dv1-dv2);

      // Use previously computed list
      // of candidates for optimal profiles

      int ae_row = -1;
      for (int i=0; i < ae_L.size(); i++) {
        if (L_av -ae_L[i] >= -tol) {
          ae_row = i; break;
        }
      }
      if (ae_row==-1) {
        stop("No pure strategy RNE exists. In static stage no incentive constraint is satisfied.");
      }

      int a1_row = -1;
      for (int i=0; i < a1_L.size(); i++) {
        if (L_av -a1_L[i] >= -tol) {
          a1_row = i; break;
        }
      }
      if (a1_row==-1) {
        stop("No pure strategy RNE exists. In static stage no incentive constraint is satisfied.");
      }

      int a2_row = -1;
      for (int i=0; i < a2_L.size(); i++) {
        if (L_av -a2_L[i] >= -tol) {
          a2_row = i; break;
        }
      }
      if (a2_row==-1) {
        stop("No pure strategy RNE exists. In static stage no incentive constraint is satisfied.");
      }

      // 3. Compute total payoffs by combining static and dynamic stage
      double U = (1-delta)*ae_G[ae_row]+dU;
      double v1 = (1-delta)*a1_c1[a1_row] + dv1;
      double v2 = (1-delta)*a2_c2[a2_row] + dv2;

      double r1 = v1+beta1*(U-v1-v2);
      double r2 = U-r1;

      res_U[row] = U;
      res_v1[row] = v1;
      res_v2[row] = v2;
      res_r1[row] = r1;
      res_r2[row] = r2;

      // Save action profiles only in the last iteration
      // Otherwise only payoffs are relevant
      if (last_iter) {
        IntegerVector a_row =
          cpp_capped_rne_find_actions(
          //find_a_fun(
            dU,dv1,dv2,
            U_hat,v1_hat,v2_hat,IC_holds,next_r1,next_r2,
            trans_mat,dest_rows,
            tie_breaking, tol
          );
        res_d_ae[row] = a_row[0];
        res_d_a1[row] = a_row[1];
        res_d_a2[row] = a_row[2];

        /*
        IntegerVector ae_a = ae_df[".a"];
        IntegerVector a1_a = a1_df[".a"];
        IntegerVector a2_a = a2_df[".a"];
        */
        int ae_a =  (int)(ae_df(ae_row,0)+0.1);
        int a1_a =  (int)(a1_df(a1_row,0)+0.1);
        int a2_a =  (int)(a2_df(a2_row,0)+0.1);


        res_s_ae[row] = ae_a;
        res_s_a1[row] = a1_a;
        res_s_a2[row] = a2_a;

        res_static_Pi[row] = ae_G[ae_row];
        res_static_c1[row] = a1_c1[a1_row];
        res_static_c2[row] = a2_c2[a2_row];

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
	    Named("s.ae")=res_s_ae,
	    Named("s.a1")=res_s_a1,
	    Named("s.a2")=res_s_a2,
	    Named("ae")=res_d_ae,
	    Named("a1")=res_d_a1,
	    Named("a2")=res_d_a2,
	    Named("static.Pi") = res_static_Pi,
	    Named("static.c1") = res_static_c1,
	    Named("static.c2") = res_static_c2,
	    Named("stringsAsFactors")=false
		);
  return res;
}
