       double precision function zdel(z, d, a, b)
c
c Given M(z), this routine computes 
c      zdel = (b^2 (x-xd)^2) + (a^2 (y-yd)^2) - (a^2 b^2)
c where x is the real part of the complex z
c       y is the imaginary part of the complex z 
c       d is the centre of the ellipse E(d,a,b)
c       xd is the real part of d
c       yd is the imaginary part of d
c       a is the real semi-axis of E(d,a,b) and
c       b is the imaginary semi-axis of E(d,a,b)
c       E(d,a,b) is defined by 
c       (x-xd)^2 + (y-yd)^2 
c       -------   --------- = 1
c         a^2        b^2
c
c if zdel<0 M is inside E(d,a,b)
c if zdel>0 M is outside E(d,a,b)
c if zdel=0 M lies on E(d,a,b)
c
       complex *16 z,d	
       double precision a, b, xd, yd, x, y, d2, one, zero,d3
       parameter (zero = 0.d0, one = 1.d0)
c
       xd=dreal(d)
       yd=dimag(d)
       x=dreal(z)
       y=dimag(z)
       d2 = (x-xd)*(x-xd)
       d3 = (y-yd)*(y-yd)
       zdel  = -one
       if (a.ne.zero) zdel = zdel+(d2/a/a)
       if (b.ne.zero) zdel = zdel+(d3/b/b)
       return
       end
