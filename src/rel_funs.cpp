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
NumericVector mat_times_vec_rows(NumericMatrix mat, NumericVector vec, IntegerVector vec_rows) {
  int ncol = mat.ncol();
  int nrow = mat.nrow();
  NumericVector res(nrow);
  for (int c=0; c<ncol;c++) {
    NumericMatrix::Column mcol = mat( _ , c);
    res += mcol * vec[vec_rows[c]];
  }
  return res;
}

// [[Rcpp::export]]
NumericVector mat_times_vec(NumericMatrix mat, NumericVector vec) {
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
    tb = mat_times_vec_rows(trans_mat,next_r_diff,dest_rows);
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
    tb = mat_times_vec_rows(trans_mat,next_r1,dest_rows);
    tb_const = -min(tb) + max(tb)-min(tb);
  } else if (tie_breaking=="max_r2") {
    tb = mat_times_vec_rows(trans_mat,next_r2,dest_rows);
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
        stop("No pure strategy RNE exists. No incentive constraint is satisfied in the dynamic stage.");
      }

      double dU = max_with_cond(U_hat,IC_holds);
      double dv1 = min_with_cond(v1_hat, IC_holds);
      double dv2 = min_with_cond(v2_hat, IC_holds);


      // 2. Solve static stage

      // Action lists of the repeated game
      List s_li = static_rep_li[row];
      DataFrame ae_df = s_li["ae.df"];
      DataFrame a1_df = s_li["a1.df"];
      DataFrame a2_df = s_li["a2.df"];

      NumericVector ae_L = ae_df["L"];
      NumericVector a1_L = a1_df["L"];
      NumericVector a2_L = a2_df["L"];

      NumericVector ae_G = ae_df["G"];
      NumericVector a1_c1 = a1_df["c1"];
      NumericVector a2_c2 = a2_df["c2"];

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

        IntegerVector ae_a = ae_df[".a"];
        IntegerVector a1_a = a1_df[".a"];
        IntegerVector a2_a = a2_df[".a"];

        res_s_ae[row] = ae_a[ae_row];
        res_s_a1[row] = a1_a[a2_row];
        res_s_a2[row] = a2_a[a2_row];

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
	    Named("d.ae")=res_d_ae,
	    Named("d.a1")=res_d_a1,
	    Named("d.a2")=res_d_a2
		);
  return res;
}

/*
#' Solve an RNE for a capped version of a multistage game
rel.capped.rne.multistage = function(g,T, save.details=FALSE, tol=1e-10,  delta=g$param$delta, rho=g$param$rho, res.field="rne", tie.breaking=c("slack","random","first","last")[1], add=TRUE, keep.all.t=FALSE) {
  restore.point("rel.capped.rne.multistage")
  if (!g$is_compiled) g = rel_compile(g)

  if (add) {
    pinfo = g$prev.capped.rne.info
    if (!is.null(pinfo)) {
      #if (pinfo$delta)
    }
  }

  g$param$delta = delta
  g$param$rho = rho

  sdf = g$sdf
  adj_delta = (1-rho)*delta
  beta1 = g$param$beta1
  beta2 = 1-beta1

  if (save.details) {
    x.df = non.null(g$x.df, quick_df(x=sdf$x))
  }

  # Use vectors for higher speed
  rne.x = rep(sdf$x,times=T)
  rne.t=rep(T:1,each=NROW(sdf))
  n = length(rne.x)
  rne.r1 = rep(NA_real_,n)
  rne.r2 = rep(NA_real_,n)
  rne.U = rep(NA_real_,n)
  rne.v1 = rep(NA_real_,n)
  rne.v2 = rep(NA_real_,n)

  rne.s.ae = rep(NA_integer_,n)
  rne.s.a1 = rep(NA_integer_,n)
  rne.s.a2 = rep(NA_integer_,n)

  rne.d.ae = rep(NA_integer_,n)
  rne.d.a1 = rep(NA_integer_,n)
  rne.d.a2 = rep(NA_integer_,n)


  rne.details = NULL
  if (save.details)
    rne.details = vector("list",NROW(rne))

  # First solve repeated games for all states
  # These are the continuation payoffs in state T
  rows = 1:NROW(sdf)
  for (row in rows) {
    if (is.null(sdf$rep[[row]])) {
      sdf$rep[[row]] = solve.x.rep.multistage(g,row=row)
    }

    # Compute U, v, r
    rep = sdf$rep[[row]] %>%
      filter(adj_delta >= delta_min, adj_delta < delta_max)

    rne.U[row] = rep[1,"U"]
    rne.r1[row] = rep[1,"r1"]
    rne.r2[row] = rep[1,"r2"]
    rne.s.ae[row] = rep[1,"s.ae"]
    rne.s.a1[row] = rep[1,"s.a1"]
    rne.s.a2[row] = rep[1,"s.a2"]
    rne.d.ae[row] = rep[1,"d.ae"]
    rne.d.a1[row] = rep[1,"d.a1"]
    rne.d.a2[row] = rep[1,"d.a2"]

    w = ((1-delta) / (1-adj_delta))
    v1 = w*rep$v1_rep + (1-w)*rep$r1
    v2 = w*rep$v2_rep + (1-w)*rep$r2
    rne.v1[row] = v1
    rne.v2[row] = v2
  }

  g$sdf = sdf


  t = T-1
  srow = 1
  # Compute all remaining periods
  for (t in rev(seq_len(T-1))) {
    is.final = (!is.null(g$final.tdf)) & (t == T-1)
    for (srow in 1:NROW(sdf)) {
      row = srow + (T-t)*NROW(sdf)
      x = sdf$x[srow]

      # 1. Solve dynamic stage

      na1 = sdf$na1[srow]
      na2 = sdf$na2[srow]
      if (!is.final) {
        trans.mat = sdf$trans.mat[[srow]]
      } else {
        trans.mat = sdf$final.trans.mat[[srow]]
      }
        #rownames(trans.mat) = make.state.lab.a(sdf[srow,])

      if (is.null(trans.mat)) {
        trans.mat = matrix(1,na1*na2,1)
        colnames(trans.mat) = x
      }

      xd = colnames(trans.mat)

      dest.rows = match(xd, sdf$x) + (T-(t+1))*NROW(sdf)


      # Include code to compute U v and r for the current state
      U.hat = (1-delta)*(sdf$pi1[[srow]] + sdf$pi2[[srow]]) +
        delta * (trans.mat %*% rne.U[dest.rows])
      U.hat = as.vector(U.hat)

      # "q-value" of punishment payoffs
      q1.hat = (1-delta)*sdf$pi1[[srow]] +
        delta * (trans.mat %*% ( (1-rho)*rne.v1[dest.rows] + rho*rne.r1[dest.rows] ))

      q2.hat = (1-delta)*sdf$pi2[[srow]] +
        delta * (trans.mat %*% ( (1-rho)*rne.v2[dest.rows] + rho*rne.r2[dest.rows] ))


      # v1.hat is best reply q for player 1
      # Note player 1 is col player
      q1.hat = matrix(q1.hat,na2, na1)
      v1.hat.short = rowMaxs(q1.hat)
      v1.hat = rep(v1.hat.short, times=na1)


      # v2.hat is best reply q for player 2
      q2.hat = matrix(q2.hat,na2, na1)
      v2.hat.short = colMaxs(q2.hat)
      v2.hat = rep(v2.hat.short, each=na2)

      # Compute which action profiles are implementable
      IC.holds = U.hat+tol >= v1.hat + v2.hat

      # Can at least one action profile be implemented?

      if (sum(IC.holds)==0) {
        # Maybe just return empty RNE
        # instead
        stop(paste0("In state ", x, " period ", t," no pure action profile can satisfy the incentive constraint. Thus no pure RNE exists in the capped game."))
      }

      U = max(U.hat[IC.holds])
      v1 = min(v1.hat[IC.holds])
      v2 = min(v2.hat[IC.holds])

      # Pick dynamic equilibrium actions
      # using the specified tie.breaking rule
      slack = U.hat - (v1.hat + v2.hat)
      if (tie.breaking=="slack") {
        tb = slack
        const = 1
      } else if (tie.breaking=="last") {
        tb = seq_len(NROW(U.hat))
        const = 1
      } else if (tie.breaking=="first") {
        tb = rev(seq_len(NROW(U.hat)))
        const=1
      } else {
        const = 1
        tb = runif(NROW(U.hat))
        #restore.point("hdfhdf")
        #if (t==1 & x=="0 0") stop()
        const = 1
      }

      rne.d.ae[row] = which.max((const+tb) * (abs(U.hat-U)<tol & IC.holds))
      rne.d.a1[row] = which.max((const+tb) * (abs(v1.hat-v1)<tol & IC.holds))
      rne.d.a2[row] = which.max((const+tb) * (abs(v2.hat-v2)<tol & IC.holds))

      # 2. Solve static stage

      # Action lists of the repeated game
      s.li = g$static.rep.li[[srow]]

      dU = U; dv1=v1; dv2=v2
      # Available liquidity after static stage
      # TO DO: Check formula
      L.av = 1/(1-delta)*(dU-dv1-dv2)

      # Use previously computed list
      # of candidates for optimal profiles
      #
      # Filter takes too long
      #s.e = filter(s.li$ae.df,L.av-L >= -tol)[1,]
      #s.1 = filter(s.li$a1.df,L.av-L >= -tol)[1,]
      #s.2 = filter(s.li$a2.df,L.av-L >= -tol)[1,]

      rows = which(L.av - s.li$ae.df$L >= -tol)
      if (length(rows)==0) stop(paste0("No incentive compatible pure static action profile exists in period ",t))
      s.e = s.li$ae.df[rows[1],]

      rows = which(L.av - s.li$a1.df$L >= -tol)
      if (length(rows)==0) stop(paste0("No incentive compatible pure static action profile exists in period ",t))
      s.1 = s.li$a1.df[rows[1],]

      rows = which(L.av - s.li$a2.df$L >= -tol)
      if (length(rows)==0) stop(paste0("No incentive compatible pure static action profile exists in period ",t))
      s.2 = s.li$a2.df[rows[1],]



      U = (1-delta)*s.e$G+dU
      v1 = (1-delta)*s.1$c1 + dv1
      v2 = (1-delta)*s.2$c2 + dv2

      r1 = v1 + beta1*(U-v1-v2)
      r2 = v2 + beta2*(U-v1-v2)


      rne.U[row] = U;
      rne.v1[row] = v1; rne.v2[row] = v2
      rne.r1[row] = r1; rne.r2[row] = r2
      rne.s.ae[row] = s.e$.a
      rne.s.a1[row] = s.1$.a
      rne.s.a2[row] = s.2$.a




      # Save only details about dynamic stage
      if (save.details) {
        pi1 = sdf$pi1[[srow]]
        Er1 = as.vector(trans.mat %*% (rne.r1[dest.rows]))
        # Continuation payoff if new negotiation in next period
        u1_neg = (1-delta)*pi1 + delta*Er1

        pi2 = sdf$pi2[[srow]]
        Er2 = as.vector(trans.mat %*% (rne.r2[dest.rows]))
        # Continuation payoff if new negotiation in next period
        u2_neg = (1-delta)*pi2 + delta*Er2


        arows = seq_along(IC.holds)
        a.info = cbind(
          quick_df(t=t),
          x.df[x.df$x==x,],
          sdf$a.grid[[srow]],
          quick_df(
            d.can.ae = (abs(U.hat-dU)<tol & IC.holds)*1 + (arows==rne.d.ae[row]),
            d.can.a1 = (abs(v1.hat-dv1)<tol & IC.holds)*1 + (arows==rne.d.a1[row]),
            d.can.a2 = (abs(v2.hat-dv2)<tol & IC.holds)*1 + (arows==rne.d.a2[row]),
            IC.holds=IC.holds,
            slack=slack,

            pi1 = pi1,
            Er1 = Er1,
            u1_neg = u1_neg,

            pi2 = pi2,
            Er2 = Er2,
            u2_neg = u2_neg,

            r1=r1,
            r2=r2,

            U.hat = U.hat,
            v1.hat=v1.hat,
            v2.hat=v2.hat,
            U=U,
            v1=v1,
            v2=v2
          )
        )
        rne.details[[row]] = a.info
      }

    }
  }
  rne = quick_df(
    x = rne.x,
    t = rne.t,
    r1 = rne.r1,
    r2 = rne.r2,
    U=rne.U,
    v1=rne.v1,
    v2=rne.v2,

    s.ae=rne.s.ae,
    s.a1=rne.s.a1,
    s.a2=rne.s.a2,

    d.ae=rne.d.ae,
    d.a1=rne.d.a1,
    d.a2=rne.d.a2
  )

  # Add some additional info

  rows = match.by.cols(rne,g$a.labs.df, cols1=c("x","d.ae"), cols2=c("x","a"))
  d.lab = g$a.labs.df$lab[rows]
  rows = match.by.cols(rne,g$gs$a.labs.df, cols1=c("x","s.ae"), cols2=c("x","a"))
  s.lab = g$gs$a.labs.df$lab[rows]
  rne$ae.lab = paste0(s.lab," | ", d.lab)

  rows = match.by.cols(rne,g$a.labs.df, cols1=c("x","d.a1"), cols2=c("x","a"))
  d.lab = g$a.labs.df$lab[rows]
  rows = match.by.cols(rne,g$gs$a.labs.df, cols1=c("x","s.a1"), cols2=c("x","a"))
  s.lab = g$gs$a.labs.df$lab[rows]
  rne$a1.lab = paste0(s.lab," | ", d.lab)

  rows = match.by.cols(rne,g$a.labs.df, cols1=c("x","d.a2"), cols2=c("x","a"))
  d.lab = g$a.labs.df$lab[rows]
  rows = match.by.cols(rne,g$gs$a.labs.df, cols1=c("x","s.a2"), cols2=c("x","a"))
  s.lab = g$gs$a.labs.df$lab[rows]
  rne$a2.lab = paste0(s.lab," | ", d.lab)



  if (!is.null(g$x.df))
    rne = left_join(rne, g$x.df, by="x")

  g$sdf = sdf
  g[[res.field]] = rne
  g[[paste0(res.field,".details")]] = rne.details
  g
}
*/
