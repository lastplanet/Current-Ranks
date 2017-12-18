	 // All the math for computing the probabilities
	 // pdf of the difference of two poisson random variables, y1 and y2, is a Skellam distribution
	 // the approximate cdf of a Skellam distribution can be computed from the cdf of a Chi-square distrition
	 // From the archives of the skellam R package
	 // p((y1-y2) <=x; x<=0) is pchisq(2*lambda2, -2*x, 2*lambda1)
	 // so 1 - p((y1-y2) <=0) is y1>y2 is 1-pchisq(2*lambda2,0,2*lambda1) (home win)
	 // p((y1-y2) <=-1) is y2>y1 is pchisq(2*lambda2,2,2*lambda1) (away win)

// sources:
// http://www.ciphersbyritter.com/JAVASCRP/BINOMPOI.HTM
// http://www.math.ucla.edu/~tom/distributions/poisson.html
 
	 function LogGamma(Z) {
	with (Math) {
		var S=1+76.18009173/Z-86.50532033/(Z+1)+24.01409822/(Z+2)-1.231739516/(Z+3)+.00120858003/(Z+4)-.00000536382/(Z+5);
		var LG= (Z-.5)*log(Z+4.5)-(Z+4.5)+log(S*2.50662827465);
	}
	return LG
}

function ncchisqdf(x,f,theta) {
    with (Math) {
    	var n=1;
		var lam=theta/2;
		var pois=exp(-lam);
		var v=pois;
		var x2=x/2;
		var f2=f/2;
		var t=pow(x2,f2)*exp(-x2-LogGamma(f2+1));
		var chisq=v*t;
		while (n<=(x-f)/2) {
			pois=pois*lam/n;
			v=v+pois;
			t=t*x/(f+2*n);
			chisq=chisq+v*t;
			n=n+1;
		}
		while (t*x/(f+2*n-x)>.000001) {
			pois=pois*lam/n;
			v=v+pois;
			t=t*x/(f+2*n);
			chisq=chisq+v*t;
			n=n+1;
		}
		return chisq
    }
}

function PoissonPDF( lambda, k ) {
   // by logs
   return  Math.exp( (k * Math.log(lambda)) - lambda - LnFact(k) );
   }

function Fact( x ) {
   // x factorial
   var  t=1;
   while (x > 1)
      t *= x--;
   return t;
   }
	 	 
function LnFact( x ) {
   // ln(x!) by Stirling's formula
   //   see Knuth I: 111
   if (x <= 1)  x = 1;

   if (x < 12)
      return Math.log( Fact(Math.round(x)) );
   else {
      var invx = 1 / x;
      var invx2 = invx * invx;
      var invx3 = invx2 * invx;
      var invx5 = invx3 * invx2;
      var invx7 = invx5 * invx2;

      var sum = ((x + 0.5) * Math.log(x)) - x;
      sum += Math.log(2*Math.PI) / 2;
      sum += (invx / 12) - (invx3 / 360);
      sum += (invx5 / 1260) - (invx7 / 1680);

      return sum;
      }
   }	 

	 function LogGamma(Z) {
	with (Math) {
		var S=1+76.18009173/Z-86.50532033/(Z+1)+24.01409822/(Z+2)-1.231739516/(Z+3)+.00120858003/(Z+4)-.00000536382/(Z+5);
		var LG= (Z-.5)*log(Z+4.5)-(Z+4.5)+log(S*2.50662827465);
	}
	return LG
}

function Gcf(X,A) {        // Good for X>A+1
	with (Math) {
		var A0=0;
		var B0=1;
		var A1=1;
		var B1=X;
		var AOLD=0;
		var N=0;
		while (abs((A1-AOLD)/A1)>.00001) {
			AOLD=A1;
			N=N+1;
			A0=A1+(N-A)*A0;
			B0=B1+(N-A)*B0;
			A1=X*A0+N*A1;
			B1=X*B0+N*B1;
			A0=A0/B1;
			B0=B0/B1;
			A1=A1/B1;
			B1=1;
		}
		var Prob=exp(A*log(X)-X-LogGamma(A))*A1;
	}
	return 1-Prob
}

function Gser(X,A) {        // Good for X<A+1.
    with (Math) {
		var T9=1/A;
		var G=T9;
		var I=1;
		while (T9>G*.00001) {
			T9=T9*X/(A+I);
			G=G+T9;
			I=I+1;
		}
		G=G*exp(A*log(X)-X-LogGamma(A));
    }
    return G
}

function Gammacdf(x,a) {
	var GI;
	if (x<=0) {
		GI=0
	} else if (x<a+1) {
		GI=Gser(x,a)
	} else {
		GI=Gcf(x,a)
	}
	return GI
}

function PoissonCDF( lambda, k ) {
// x < = k
		 return 1-Gammacdf(lambda,k+1);
}
