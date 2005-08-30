double f_powdi(x,i)
double x; int i;
{
int j,k;
double y;
if (i==0) return 1.0;
j=(i>0)?i:-i;
y=1.0;
for (k=0; k<j; k++) y=y*x;
if (i<0) y=1/y;
return y;
}
