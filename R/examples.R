examples.simple = function() {
  library(RelationalContracts)
  e.seq = seq(0,1,by=0.1);
  g = rel_game("Principal Agent") %>%
    rel_param(delta=0.8, rho=0.5) %>%
    rel_state("x0",A2=list(e=e.seq),
      pi1=~e, pi2=~ -0.5*e*e) %>%
    rel_state("x1",A2=list(e=c(-.1,e.seq)),
      pi1=~e, pi2=~ -0.5*e*e) %>%
    rel_state("x2",A2=list(e=c(-.2,e.seq)),
      pi1=~e, pi2=~ -0.5*e*e) %>%
    rel_state("x3",A2=list(e=e.seq),
      pi1=~e, pi2=~ -0.5*e*e) %>%
    rel_state("x4",A2=list(e=e.seq),
      pi1=~e, pi2=~ -0.5*e*e) %>%
    rel_compile() %>%
    rel_capped_rne(T=2)

  g = rel_capped_rne(g,T=10)

  rne=g$rne
  capped_rne_iterations(g=g,T=1000, rne=filter(g$rne,t==max(rne$t)))

  capped_rne_iterations(g=g,T=2,debug_row=0, rne=filter(g$rne,t==max(rne$t)))
  capped_rne_iterations(g=g,T=2,debug_row=1, rne=filter(g$rne,t==max(rne$t)))



  g = rel_capped_rne(g,T=2, adjusted.delta=0.17, rho=0.7,tie.breaking = "random")


  rne.diagram(g)

}

capped_rne_iterations = function(g,T=1,rne=g$rne, tie.breaking="slack", debug_row=-1, tol=1e-12) {
  restore.point("capped_rne_iterations")

  sdf = g$sdf
  transmats = lapply(1:NROW(sdf), function(row) {
    trans.mat = sdf$trans.mat[[row]]
    if (is.null(trans.mat)) {
      x = sdf$x[row]; na1 = sdf$na1[row]; na2 = sdf$na2[row]
      trans.mat = matrix(1,na1*na2,1)
      colnames(trans.mat) = x
    }
    trans.mat
  })

  new_rne = cpp_capped_rne_iterations(T, sdf=sdf,rne,transmats,delta=g$param$delta, rho=g$param$rho,beta1 = g$param$beta1,
    #find_a_fun = capped.rne.find.actions,
    tie_breaking=tie.breaking, tol=tol, debug_row=debug_row)

  new_rne

  return(new_rne)

  trans_mat = transmats[[1]]
  next_U = rne$U
  dest_rows = 0
  mat_times_vec_rows(trans_mat,next_U, dest_rows)


}


examples.multistage = function() {
  library(RelationalContracts)

  static.A.fun = function(x1,x2, x.seq,e.seq,...) {
    restore.point("static.A.fun")
    A1 = list(b1=c("b",""))
    A2 = quick_df(b2=c("b", rep("", length(e.seq))),e=c(0,e.seq))
    list(A1=A1,A2=A2)
  }

  A.fun = function(x1,x2,stage, x.seq,...) {
    restore.point("static.A.fun")
    A1 = list(a1=x.seq[x.seq>=x1])
    A2 = list(a2=x.seq[x.seq>=x2])
    list(A1=A1,A2=A2)
  }

  vec.static.pi.fun = function(ax.df,...) {
    restore.point("vec.pi.fun")
    ax.df %>%
      transmute(
        x=x,
        pi1= ifelse(b1 == "b" | b2=="b",-x1,e),
        pi2= ifelse(b1 == "b" | b2=="b",-x2,- 1/2 * e^2),
      )
  }

  vec.trans.fun = function(ax.df, final=FALSE,...) {
    restore.point("trans.fun")
    ax.df %>%
      select(x,a1,a2) %>%
      unique() %>%
      transmute(xs=x,xd=paste0(a1, " ",a2),a1=a1,a2=a2, prob=1)
  }


  x.seq = seq(0,1, by=0.1)
  #x.seq = c(0,0.01,0.05,0.1,0.2,0.5,1)
  x.df = as_data_frame(expand.grid(x1=x.seq,x2=x.seq, stringsAsFactors = FALSE)) %>%
    mutate(x= paste0(x1," ", x2))

  g = rel_game("Slowly Intensifying Repeated Principal-Agent") %>%
    rel_param(x.seq=x.seq, e.seq=seq(0,0.1,by=0.01)) %>%
    rel_states(x.df,
      # Static effort stage
      static.A.fun=static.A.fun,
      vec.static.pi.fun = vec.static.pi.fun,
      # Dynamic relationship intensification stage
      A.fun = A.fun,
      pi1 = 0, pi2=0,
      vec.trans.fun=vec.trans.fun
    )

  g = rel_compile(g)
  g = rel_capped_rne(g,T=1, adjusted.delta=0.17, rho=0.7,tie.breaking = "random")


  g$new_rne =rne= capped_rne_multistage_iterations(g,T=100, debug_row=-1, tie.breaking="equal_r")


  rne.diagram(g,eq.field = "new_rne",just.eq.chain = TRUE)
  (rne = g$rne)
}







