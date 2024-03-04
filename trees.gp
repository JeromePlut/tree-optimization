\\ cost of the 2-3 tree
\\ we check that: f(2^e m) = 2^e (m⋅e + f(m)) for m odd
\\ f(2^i+1) = i*(2^i+1)
\\ f(2^i+3) = i*(2^i+3) (i large enough)
\\ f(2^i+5) = i*(2^i+5) (idem)
\\ f(i) is generally close to i*floor(log2(i))
\\ except in intervals: ]7*2^i, 2^(i+3)]

ispow2=n->frac(log2(n))==0;
T2p=n->if(n<=1,0, n+if(n%2==0, 2*T2p(n/2), T2p((n-1)/2)+T2p((n+1)/2)));
T2s=T2p;
T23p=n->if(n<=1,0,n+if(ispow2(n/3),3*T23p(n/3),if(n%2==0,2*T23p(n/2),T23p((n-1)/2)+T23p((n+1)/2))));
T23s=n->if(n<=1,0,if(ispow2(n/3),2*n+3*T23s(n/3),n+if(n%2==0,2*T23s(n/2),T23s((n-1)/2)+T23s((n+1)/2))));
T3p=n->if(n<=5, [0,2,3,8,10][max(1,n)], n+sum(i=0,2,T3p(floor((n+i)/3))));
T3s=n->if(n<=5, [0,2,6,8,13][max(1,n)], 2*n+sum(i=0,2,T3s(floor((n+i)/3))));

f=T23p;
L=n->ceil(log(n)/log(2));
g2=(b,x)->2^b*(b*x+x-1);

F=(n,M=[;])->{
	for(i=1,#M~,if(n==M[i,1],return(M[i,2])));
	if(n<=1, return(0));
	return(n+F(floor(n/2),M)+F(floor((n+1)/2),M));
};

defun("tabulate", (f=T3p, N=30,l=L)->{
	for(n=1,N,
		b=l(n);
		printf("%3d│%4d‖%5d|%5d\n", b, n, f(n), f(n)-g2(b,n/2^b));
	);
});
defun("tplot", (f=T2p, N=30,l=L)->{
	local(v);
	v=vector(N, n, f(n)-n*l(n));
	print1("[32m");
	forstep(y=vecmax(v), vecmin(v), -1,
		if(y==0,
			print1("[m");
			for(i=1,#v,print1("━"));
			print1("\n[31m");
		);
		for(x=1,#v, if(abs(v[x])>=abs(y), print1("█"), print1(" ")));
		print();
	)
});

tabulate(T2p);
