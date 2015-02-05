       double precision function ddel(x, y, d, a, b)
c
c Given M(x,y), this routine computes 
c      ddel = (b^2 (x-d)^2) + (a^2 y^2) - (a^2 b^2)
c where d is the centre of the ellipse E(d,a,b)
c       a is the real semi-axis of E(d,a,b) and
c       b is the imaginary semi-axis of E(d,a,b)
c       E(d,a,b) is defined by 
c       (x-d)^2 + y^2 
c       -------   --- = 1
c         a^2     b^2
c
c if ddel<0 M is inside E(d,a,b)
c if ddel>0 M is outside E(d,a,b)
c if ddel=0 M lies on E(d,a,b)
c
       double precision a, b, d, x, y, d2, one, zero
       parameter (zero = 0.d0, one = 1.d0)
c
       d2 = (x-d)*(x-d)
       ddel  = -one
       if (a.ne.zero) ddel = ddel+(d2/a/a)
       if (b.ne.zero) ddel = ddel+(y*y/b/b)
       return
       end
