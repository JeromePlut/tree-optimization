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
T23p=n->if(n<=1,0,if(n==3,3,n+T23p(floor(n/2))+T23p(floor((n+1)/2))));
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
D=3;
L=(n,d=D)->if(n<=0, 0, ceil(log(n)/log(d)));
P=(n,d=D)->if(n <= 1, 0, n+sum(i=0,d-1,P(floor((n+i)/d),d)));
Q=(n,d=D)->{local(e=L(n,d)); if(n<=2*d^(e-1), (e+1)*n-2*d^(e-1), e*n)};
E=(n,d=D)->{if(n<=0,return(1/(1-d)));local(e=L(n,d)); -(d^e)/(d-1)+e*(n)};
Dr=(f,d=D)->(n)->f(n)-sum(i=0,d-1,f(floor((n+i)/d)));
In=(g,d=D)->(n)->if(n<=1,[-g(0)/(d-1),'C][1],g(n)+sum(i=0,d-1,In(g,D)(floor((n+i)/d))));
\\ G=(f=1,d=D)->(n)->{if(n<=d-1,return(0));
\\ 	subst(f,'x,n)+sum(i=1,d-1,G(f,d)(floor((n+i)/d)))};
\\ make_G(f,d=D)=(n)->{ f(n)-sum(i=0,d-1,f(floor((n+i)/d)))};
\\ log_e=(d=D)->(n)->if(n<=0,0,ceil(log(n)/log(d)));
clogD=(n,d=D)->if(n<=1,[l0,0][n+1],ceil(log(n)/log(d)));
flogD=(n,d=D)->if(n<=1,[l0,0][n+1],floor(log(n)/log(d)));
U=(n,d=D)->(n<=1);
V=(n,d=D)->n-(n==1);
J=(n,d=D)->n+(n==0);
\\ G=(d)->(n)->{if(n<=0,return(0));local(e=ceil(log(n)/log(d)));
\\ 	e*d^e-(d^e-1)/(d-1)+e*(n-d^e)};
\\ test(d=3,N=50)={
\\ 	for(n=1,N,local(t);
\\ 	local(g, s, t); g = G(d)(n); s = (n-1)+sum(i=0,d-1,G(d)(floor((n+i)/d)));
\\ 	print([n,g,s,g-s]);
\\ 	);}
\\ (e+2)*(n-d^e)+e*d^e, (e+1)*(n-2*d^e)+2*(e+1)*d^e)};

defun("slopes", (f=P(3), N=300)->{
  my(s = 0, y=0);
  for(n=1, N,
		my(z, t);
		z = f(n+1); t = z-y;
		print([n, y, z]);
\\ 		if(t != s, printf("from %03d: slope = %03d\n", n, t); s=t);
		if(t != s, printf("from %03d: slope = %s (f(%s)=%s)\n", n,
		f(n+1)-f(n), n, f(n)); s=t);
	);
)};
defun("tabulate", (f=T3p, N=25,l=log_e(D))->{
	for(n=1,N,
		b=l(n);
		printf("%3d│%4d‖%5d\n", b, n, f(n));
\\ 		printf("%3d│%4d‖%5d|%5d\n", b, n, f(n), f(n)-g2(b,n/2^b));
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

