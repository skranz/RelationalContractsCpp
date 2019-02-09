#include <Rcpp.h>
using namespace Rcpp;

// ax specifies on action profile per state
// returns a vector of ax indices that contains
// all possible replies for player 1 in all states
// Note we assume ax is indexed starting with 1
// [[Rcpp::export]]
IntegerVector c_pl1_ax_replies(
  IntegerVector ax, IntegerVector na1, IntegerVector na2
) {
  int nx = na1.size();
  int nres = sum(na1);
  IntegerVector replies(nres);

  int offset = 0;
  int ind = 0;
  // Loop through all states
  for (int xrow=0; xrow < nx; xrow++) {
    int a = ax[xrow] - offset-1; // -1 because ax is 1-based
    int a2 = a % na2[xrow];
    // Insert all replies of player 1
    int rep_ind = offset+a2;
    for (int a1=0; a1<na1[xrow]; a1++) {
      replies[ind] = rep_ind+1; // +1 because replies are 1-based
      //Rcout << "xrow " << xrow << " a " << a << " a1 "<<a1 << " a2 " << a2 << std::endl <<
      //  " offset" << offset << " ind " << ind << " reply "<< rep_ind+1 << std::endl;

      ind++;
      rep_ind += na2[xrow];
    }
    offset += na1[xrow]*na2[xrow];
  }
  return replies;
}

// [[Rcpp::export]]
IntegerVector c_pl2_ax_replies(
  IntegerVector ax, IntegerVector na1, IntegerVector na2
) {
  int nx = na1.size();
  int nres = sum(na2);
  IntegerVector replies(nres);

  int offset = 0;
  int ind = 0;
  // Loop through all states
  for (int xrow=0; xrow < nx; xrow++) {
    int a = ax[xrow] - offset -1;
    int a1 = a / na2[xrow];

    // Insert all replies of player 2
    int rep_ind = offset+a1*na2[xrow];
    for (int a2=0; a2<na2[xrow]; a2++) {
      replies[ind] = rep_ind+1;
      ind++;
      rep_ind ++;
    }
    offset += na1[xrow]*na2[xrow];
  }
  return replies;
}

// u_ax is of length nax
// ax is of length nx, it specifies on action profile per state
// changes ax such that player 1 plays a best reply
// Note we assume ax is indexed starting with 1
// [[Rcpp::export]]
IntegerVector c_pl1_best_reply_ax(NumericVector u_ax,
  IntegerVector ax, IntegerVector na1, IntegerVector na2
) {
  int nx = na1.size();
  IntegerVector res(nx);

  double u_br = 0;
  int br = 0;
  int offset = 0;
  // Loop through all states
  for (int xrow=0; xrow < nx; xrow++) {
    int a = ax[xrow] - offset-1; // -1 because ax is 1-based
    int a2 = a % na2[xrow];
    int rep_ind = offset+a2;
    for (int a1=0; a1<na1[xrow]; a1++) {
      double u_cur = u_ax[rep_ind];
      //Rcout << "xrow " << xrow << " a " << a << " a1 "<<a1 << " a2 " << a2 << std::endl <<
      //  " offset " << offset << " rep_ind " << rep_ind << " u_cur " << u_cur << " u_br " << u_br << std::endl;

      if ((a1==0) | (u_cur > u_br)) {
        u_br = u_cur;
        br = rep_ind+1; // Add +1 for R index
      }
      rep_ind += na2[xrow];
    }
    res[xrow] = br;
    offset += na1[xrow]*na2[xrow];
  }
  return res;
}


// [[Rcpp::export]]
IntegerVector c_pl2_best_reply_ax(NumericVector u_ax,
  IntegerVector ax, IntegerVector na1, IntegerVector na2
) {
  int nx = na1.size();
  IntegerVector res(nx);

  int br = 0;
  int offset = 0;
  double u_br = 0;

  // Loop through all states
  for (int xrow=0; xrow < nx; xrow++) {
    int a = ax[xrow] - offset -1;
    int a1 = a / na2[xrow];

    int rep_ind = offset+a1*na2[xrow];
    for (int a2=0; a2<na2[xrow]; a2++) {
      double u_cur = u_ax[rep_ind];
      if ((a2==0) | (u_cur > u_br)) {
        u_br = u_cur;
        br = rep_ind+1; // Add +1 for R index
      }
      rep_ind ++;
    }
    res[xrow] = br;
    offset += na1[xrow]*na2[xrow];
  }
  return res;
}



// [[Rcpp::export]]
IntegerVector c_pl2_ax_best_replies(
  NumericVector u_ax, IntegerVector nai, IntegerVector naj, int nx
) {

  int start_ind = 0;
  IntegerVector br_ax(u_ax.size());
  // Loop through all states
  double u_br = 0;
  int br = 0;
  for (int xrow=0; xrow < nx; xrow++) {
    // Loop through other player's actions
    // in current state
    for (int aj=0; aj<naj[xrow]; aj++) {
      int ind = start_ind + aj*nai[xrow];
      // Loop through own actions
      // To find best reply payoff
      for (int ai=0; ai<nai[xrow]; ai++) {
        double u_cur = u_ax[ind];
        if ((ai==0) | (u_cur > u_br)) {
          u_br = u_cur;
          br = ind+1;
        }
        ind = ind + 1;
      }
      // Loop through own actions to
      // set best reply payoff
      ind = start_ind + aj*nai[xrow];
      for (int ai=0; ai<nai[xrow]; ai++) {
        br_ax[ind] = br;
        ind = ind + 1;
      }
    }
    start_ind = start_ind + naj[xrow]*nai[xrow];
  }
  return br_ax;
}


// [[Rcpp::export]]
NumericVector c_pl1_ax_best_reply_payoffs(
  NumericVector u_ax, IntegerVector nai, IntegerVector naj, int nx
) {

  int start_ind = 0;
  NumericVector br_ax(u_ax.size());
  // Loop through all states
  double u_br = 0;
  for (int xrow=0; xrow < nx; xrow++) {
    // Loop through other player's actions
    // in current state
    for (int aj=0; aj<naj[xrow]; aj++) {
      int ind = start_ind + aj;
      // Loop through own actions
      // To find best reply payoff
      for (int ai=0; ai<nai[xrow]; ai++) {
        double u_cur = u_ax[ind];
        if ((ai==0) | (u_cur > u_br)) {
          u_br = u_cur;
        }
        ind = ind+ naj[xrow];
      }
      // Loop through own actions to
      // set best reply payoff
      ind = start_ind + aj;
      for (int ai=0; ai<nai[xrow]; ai++) {
        br_ax[ind] = u_br;
        ind = ind+naj[xrow];
      }
    }
    start_ind = start_ind + naj[xrow]*nai[xrow];
  }
  return br_ax;
}


// [[Rcpp::export]]
NumericVector c_pl2_ax_best_reply_payoffs(
  NumericVector u_ax, IntegerVector nai, IntegerVector naj, int nx
) {

  int start_ind = 0;
  NumericVector br_ax(u_ax.size());
  // Loop through all states
  double u_br = 0;
  for (int xrow=0; xrow < nx; xrow++) {
    // Loop through other player's actions
    // in current state
    for (int aj=0; aj<naj[xrow]; aj++) {
      int ind = start_ind + aj*nai[xrow];
      // Loop through own actions
      // To find best reply payoff
      for (int ai=0; ai<nai[xrow]; ai++) {
        double u_cur = u_ax[ind];
        if ((ai==0) | (u_cur > u_br)) {
          u_br = u_cur;
        }
        ind = ind + 1;
      }
      // Loop through own actions to
      // set best reply payoff
      ind = start_ind + aj*nai[xrow];
      for (int ai=0; ai<nai[xrow]; ai++) {
        br_ax[ind] = u_br;
        ind = ind + 1;
      }
    }
    start_ind = start_ind + naj[xrow]*nai[xrow];
  }
  return br_ax;
}
