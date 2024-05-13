
defun("geom_sum", (d,e,X='X)->{
	local(a,b,c); a=d^e; b=2*a; c=d*a;
	local(T=(u,v,f)->(1-X)*sum(n=u,v-1,X^n*subst(f,'n,n)));
	S1 = T(a,b, (e+2)*n-b ) + T(b, c, (e+1)*n);
	S2 = (e+2)*T(a,b,n) + (e+1)*T(b,c,n) - b*T(a,b,1);
	U1 = (1-X)*S1;
	U2 = (1-X)*S2;
	U3 = (e+2)*((b-1)*X^(b+1)-b*X^b-(a-1)*X^(a+1)+a*X^a);
	U3+= (e+1)*((c-1)*X^(c+1)-c*X^c-(b-1)*X^(b+1)+b*X^b);
	U3+= -b*(X^(b+1)-X^b-X^(a+1)+X^a);
	U4 = (e+1)*(c-1)*X^(c+1) - (e+1)*c*X^c;
	U4+= -X^(b+1);
\\ 	U4+= ((e+2)*(b-1)-(e+1)*(b-1)-b)*X^(b+1);
	U4+= (-e*a+e+2)*X^(a+1);
	U4+= (e*a)*X^a;
	print(U2-U4);
	U2
});

W=b->Mod('w,polcyclo(2*b,'w));
R=(expr,N=1e3)->round(expr*N)/N;
X=(expr,b)->R(subst(lift(expr),'w,exp(I*Pi/b)));

P=b->local(w=W(b));matrix(b+1,b+1,i,j,if(j==b,I^(i-1),if(j==b+1,I^(1-i),if(i>=2&&i<=b,(w^((i-1)*j)-w^-((i-1)*j))/(2*I)))));
Q=b->local(w=Mod('w,polcyclo(2*b,'w)));matrix(b-1,b-1,i,j,(w^(i*j)-w^-(i*j))/(2*I));
M=b->matrix(b+1,b+1,i,j,(i>1)&&(i<=b)&&(abs(i-j)==1));

P1=b->matrix(b+1,b+1,i,j,if(j==b,I^(i-1),if(j==b+1,I^(1-i),if(i>1&&i<=b,sin((i-1)*j*Pi/b)))));
Q1=b->matrix(b+1,b+1,i,j,if(j==1,if(i<b,tan(i*Pi/b),b/2),if(j==b+1,if(i<b,(-1)^(i-1)*tan(i*Pi/b),(-1)^i*b/2*I^b),if(i<b,2*sin(i*(j-1)*Pi/b)))));
D=b->matdiagonal(vector(b+1,j,if(j<b,2*cos(j*Pi/b))));

Z=b->vector(b+1,i,z^(i-1))~;
{ PZ=b->vector(b+1,j,if(j<b,
	1/b*tan(j*Pi/b)/(1-2*cos(j*Pi/b)*z+z^2)*(1-(-1)^j*z^b)*(1+z^2),
	(1+(-1)^j*(I*z)^b)/2))~; }
{ MPZ=b->vector(b+1,j,if(j<b,1/(1-2*x*cos(j*Pi/b))*
	1/b*tan(j*Pi/b)/(1-2*cos(j*Pi/b)*z+z^2)*(1-(-1)^j*z^b)*(1+z^2),
	(1+(-1)^j*(I*z)^b)/2))~; }
F=b->(1-x*M(b))^-1*Z(b);
{F1=b->vector(b+1,j,sum(k=1,b-1,sin((j-1)*k*Pi/b)*
	1/(1-2*x*cos(k*Pi/b))*
	1/b*tan(k*Pi/b)/(1-2*cos(k*Pi/b)*z+z^2)*(1-(-1)^k*z^b)*(1+z^2))
	+ I^(j-1)*(1+(-1)^(j-1)-(I*z)^b+(-1)^(j-1)*(I*z)^b)/2)~;}
\\ 	+ (I^(j-1)+(-I)^(j-1)-I^(j-1)*(I*z)^b+(-I)^(j-1)*(I*z)^b)/2)~;}
R(F(5)-F1(5))
G=b->taylor(substvec(F(b),[x],['X]),z);
H=(b,d)->vector(b+1,j,polcoeff(G(b)[j], d, z))~;
{F0=b->local(f=F(b)); subst(f,z,0);}
{F01=b->R(vector(b+1,j,
	1/b*sum(k=1,b-1,sin((j-1)*k*Pi/b)*tan(k*Pi/b)/(1-2*x*cos(k*Pi/b)))
	+real(I^(j-1))));}
	
\\ {F1=b->vector(b+1,j,
\\ 	1/b*sum(k=1,b-1,sin(j*k*Pi/b)*tan(k*Pi/b)/(1-2*x*cos(k*Pi/b))
\\ 	*(1-(-1)^k*z^b)/(1-2*z*cos(k*Pi/b)+z^2))
\\ 	+(I^j-(-I)^j-I^j*(I*z)^b+(-I)^j*(I*z)^b)/2)~;}
