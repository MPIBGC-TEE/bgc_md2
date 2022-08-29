library(hier.part)
# additive model
make_hier_a1_func_sum=function(x1,x2){
  a1_from_r1=function(r1){
    s1=r1*x1
    s2=(1-r1)*x2
    y=s1+s2
    df=data.frame(s1=s1,s2=s2)
    res=hier.part(y,df,barplot=F)
    Is=res$IJ$I
    att1=Is[1]
    return(att1)
  }
  return(a1_from_r1)
}

make_cov_a1_func_sum=function(x1,x2){
  a1_from_r1=function(r1){
    s1=r1*x1
    s2=(1-r1)*x2
    y=s1+s2
    vy=cov(y,y)
    c1=cov(s1,y)/vy
    c2=cov(s2,y)/vy
    #print(c1+c2)
    return(c1)
  }
  return(a1_from_r1)
  #return(vecf)
}

vectorize=function(f){
  vecf=function(r1s){
    return(unlist(lapply(r1s,f)))
  }
  return(vecf)
}


n=1000
ts=(0:(n-1))/n*2*pi
x1s=sin(ts)
#x2s=cos(ts)
#x2s=sin(ts+pi/4)
#x2s=sin(ts)
x2s=ts
r1s=seq(0,1,0.01)
make_sub_plot=function(sub_list){
  func1=sub_list$func1
  func2=sub_list$func2
  color=sub_list$color
  a1s_hier=vectorize(make_hier_a1_func_sum(func1(ts),func2(ts)))(r1s)
  a1s_cov=vectorize(make_cov_a1_func_sum(func1(ts),func2(ts)))(r1s)
  t_hier=1
  t_cov=2
  lines(
  	r1s,
	a1s_hier,
  	col=color,
	lty=t_hier,
  	main="attribution of variance to y=r1*x1+(1-r1)*x2",
  )
  lines(
  	r1s,
  	a1s_cov,
  	col=color,
	lty=t_cov
  )
  #lines(
  #	r1s,
  #	r1s,
  #	col='grey'
  #)
  legend(
    0.1,
    0.9,
    legend=c("hier.part","cov.part"),
    lty=c(t_hier,t_cov)
  )
}
plot(
  r1s,
  r1s,
  #ylim=c(
  #  min(a1s_hier,a1s_cov),
  #  max(a1s_hier,a1s_cov)
  #)  
)
func_colors=list(
  a=list(
    func1=sin,
    func2=cos,
    color="green"
  )
  ,
  b=list(
    func1=sin,
    func2=function(x){x},
    color="red"
  )
)
make_sub_plot(func_colors$a)
make_sub_plot(func_colors$b)
