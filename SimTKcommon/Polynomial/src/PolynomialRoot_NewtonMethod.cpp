//      max polynomial order is 9         
//      check program with the following polynomials
//	5x3-7x2-40x+100=0  interval [-6.0, 7.0]                          
//	x2+2x-4=0  intervals [0.0, 2.0] and [-5.0, 0.0]
//	3x2-4x-4=0  intervals [-2.0, 0.0] and [0.5, 4.0]                                       

#include <stdio.h>                          
#include <math.h>               
int i,n,count=0,b[10];
double posx, negx, ypos=0.0,yneg=0.0,y=0.0,previous=-999.0, x=0.0, posx0,negx0,temp ;

void PolynomialRootFinder()
	{
	//input order and coefficients of the polynomial...                
	do	{	printf("order of polynomial =");
			scanf("%d",&n);	
		}  
	while (n>10);
	
	printf("coefficients from high to low order\n");
	for (i=n;i>=0;i--) 
		{	printf("b[%2d ]=",i);
			scanf("%d",&b[i]);
		}
	// print the array of coefficients in decreasing order                
	printf("COEFFICIENTS OF THE POLYNOMIAL\n");
	for (i=n;i>=0;i--) printf("b[%2d ]=%4d\n",i,b[i]);

	// input the possible interval of a root [negx,posx]
	do {
			printf("give the upper bound  posx="); 
			scanf("%lf",&posx);
			printf("give the lower bound  negx="); 
			scanf("%lf",&negx);
		}
	while (fabs(posx-negx)<=0.0001);

	// ensure that always negx<posx
	if (negx>posx) 
		{temp=posx; posx=negx; negx=temp;};
	
	// keep the initial values of bounds
		posx0=posx; negx0=negx;

		printf("starting interval [negx = %12.3lf  posx = %12.3lf]\n",negx,posx);
 
	// starting calculations...........................            
	for (i=0;i<=n;i++) 	yneg=yneg+b[i]*pow(negx,i); 
	printf("yneg = %12.4lf\n\n",yneg);

	for (i=0;i<=n;i++) 	ypos=ypos+b[i]*pow(posx,i); 
	printf("ypos = %12.4lf\n\n",ypos);

	if (ypos*yneg>0.0) 
		printf("no root between values  %10.3lf   and %10.3lf\n",negx0,posx0);
	else 	//repeated process
		{	// swap locations posx, negx in case ypos<yneg
			if (ypos<yneg) {temp=posx; posx=negx; negx=temp;
							temp=ypos; ypos=yneg; yneg=temp;};
	
		while ((ypos*yneg<0.0) && (fabs(x-previous)>0.0001) )
		{	count++;	// counts the number of iterations
			previous=x;
	
			//calclulate new value of x and y............................
			x=(posx+negx)/2.0; 
			y=0; 
			for (i=0;i<=n;i++) y=y+b[i]*pow(x,i);
			printf("previous=%10.6lf  x=%10.6lf  y=%12.5lf\n",previous,x,y);
			// check if root is found
			if ((fabs(y)<0.0001) || (fabs(x-previous)<=0.0001))
				{ printf("root is found !!! root=%10.6f  in %4d  iterations\n",x,count); 
					previous=x;  }
			else	{
						if (y>0.0) { posx=x; ypos=y; } 
							else
								if (y<0.0) { negx=x; yneg=y; }
					}
		}; // end while


	}  // end if

}  // end Polynomial Root Finder
