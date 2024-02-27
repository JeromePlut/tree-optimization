\\ cost of the 2-3 tree
\\ we check that: f(2^e m) = 2^e (m⋅e + f(m)) for m odd
\\ f(2^i+1) = i*(2^i+1)
\\ f(2^i+3) = i*(2^i+3) (i large enough)
\\ f(2^i+5) = i*(2^i+5) (idem)
\\ f(i) is generally close to i*floor(log2(i))
\\ except in intervals: ]7*2^i, 2^(i+3)]

ispow2=n->frac(log2(n))==0
f=n->if(n<=1,0,n+if(ispow2(n/3),3*f(n/3),if(n%2==0,2*f(n/2),f((n-1)/2)+f((n+1)/2))))

