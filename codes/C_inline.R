xx <- faithful$eruptions
fit1 <- density(xx)
fit2 <- replicate(100, {
  x <- sample(xx,replace=TRUE);
  density(x, from=min(fit1$x), to=max(fit1$x))$y
})
fit3 <- apply(fit2, 1, quantile,c(0.025, 0.975))
plot(fit1, ylim=range(fit3))
polygon(c(fit1$x,rev(fit1$x)),
        c(fit3[1,], rev(fit3[2,])),
        col='grey', border=F)
lines(fit1)


## we need a pure C/C++ function as the generated function
## will have a random identifier at the C++ level preventing
## us from direct recursive calls
incltxt <- '
int fibonacci(const int x) {
  if (x == 0) return(0);
  if (x == 1) return(1);
  return fibonacci(x - 1) + fibonacci(x - 2);
  }'

## now use the snipped above as well as one argument conversion
## in as well as out to provide Fibonacci numbers via C++
library(inline)
fibRcpp <- cxxfunction(signature(xs="int"),
                          plugin="Rcpp",
                          incl=incltxt,
                          body='
                          int x = Rcpp::as<int>(xs);
                          return Rcpp::wrap( fibonacci(x) );
                          ')
fibRcpp(20)

myroot <- cppFunction('double myroot(double x) { return ::sqrt(x); }')
myroot(16)

myroot2 <- rcpp(signature(xs="numeric"), 
                body='double x=as<double>(xs); return wrap(::sqrt(x));')
myroot2(16)
