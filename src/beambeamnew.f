c----------------------------------------------------------
c This program calculates the beam-beam-kick with image charges for
c     an elliptical boundary.
c
c      25.08.99 FZ, include offset xof, yof for correct treatment
c               of image charges   
c
       subroutine foliver(x,y,xof,yof,sx,sy,rg,rm1,fr,fi)
       implicit real*8 (a-h,o-z)
c       implicit integer (a-z)
       parameter(nimage=30)
       double precision xt,yt
       real*8 x,y,sx,sy,rg,rm1,fr,fi, efx, efy, r2

       complex*16 zone, zi, z1, z2, z3, z4, efz
       complex*16 z1b, z2b, z3b, z4b

       charge = 1.0d0
       zone = (1.0d0,0.0d0)
       zi = (0.0d0,1.0d0)

c      new variables 
       dx = x-xof
       dy = y-yof
c

       xt = abs(dx)
       yt = abs(dy)     
       if (xt.lt.(1d-7)) xt=1.d-7
       if (yt.lt.(1d-7)) yt=1.d-7
       r2 = xt*xt+yt*yt 
       r = dsqrt(r2)

c       write(*,*) ' in foliver ',dx,dy,sx,sy,rg,rm1
c       write(*,*) xof,yof,x,y 

       pi = acos(-1.0d0)
       res = 2.0d0/sqrt(2.0d0*pi)/r2
       resmod = 2./sqrt(2.*pi)
       sx1 = sx
       sy1 = sy

       if (abs((sx/sy)-1.0d0).lt.1.0d-2) then
          if(r2.le.sx*sy) then
             co = (1.0d0-exp(-r2/(2.0d0*sx*sy)))
             fr = res*co*yt
             fi = res*co*xt
             return
          endif

          if(r2.gt.sx*sy) then
             efx = dx/ r2
             efy = dy/ r2
             efz = efx * zone + efy * zi

             z1 = x/rg * zone + y/rg * zi
             z2 = x/rg * zone - y/rg * zi
             z3 = cdlog(z1 - cdsqrt(z1**2 - zone))
             z4 = cdlog(z2 - cdsqrt(z2**2 - zone))
             rmi = dabs(dreal(z3+z4)/2.0d0)
             rphi = dreal((z4-z3)*zi)/2.0d0

c            introduce the same quantities for the beam
             z1b = xof/rg * zone + yof/rg * zi
             z2b = xof/rg * zone - yof/rg * zi
             z3b = cdlog(z1b - cdsqrt(z1b**2 - zone))
             z4b = cdlog(z2b - cdsqrt(z2b**2 - zone))
             rmib = dabs(dreal(z3b+z4b)/2.0d0)
             rphib = dreal((z4b-z3b)*zi)/2.0d0

             if(rphi*y.lt.0.0d0) rphi = -rphi
             if(rphib*yof.lt.0.0d0) rphib = -rphib

             do 100 i=1,nimage
                charge = 1.0d0
             if ((xof.eq.0.0).and.(yof.eq.0.0)) then
                efz = efz+((-1)**i)*
     *               8.0d0*charge/(dexp(4.0d0*i*rm1)+1.0d0) * 
     *               (cdexp(2.0d0*i*(rmi*zone-rphi*zi)) - 
     *                cdexp(-2.0d0*i*(rmi*zone-rphi*zi)))/
     *               (cdexp((rmi*zone-rphi*zi)) - 
     *                cdexp(-(rmi*zone-rphi*zi)))/ rg
             else 
                efz = efz+2.d0*charge*dexp(-i*rm1)*(
     *             (dexp(i*rmib)+dexp(-i*rmib))*
     *             cos(i*rphib)/(dexp(i*rm1)+dexp(-i*rm1))+
     *             zi*(dexp(i*rmib)-dexp(-i*rmib))*
     *             sin(i*rphib)/(dexp(i*rm1)-dexp(-i*rm1)))*
     *             (cdexp(i*(rmi*zone-rphi*zi))-
     *             cdexp(-i*(rmi*zone-rphi*zi)))/
     *             (cdexp((rmi*zone-rphi*zi))-
     *             cdexp(-(rmi*zone-rphi*zi)))/rg
             endif 

 100         continue
             efx = dreal(efz)
             efy = dreal(-efz*zi)
             fi = abs(efx)*2.0d0/sqrt(2.0d0*pi)
             fr = abs(efy)*2.0d0/sqrt(2.0d0*pi)
             return
          endif

       else  
c      flat beam        

       fr2 = resmod * yt/r * co
       fi2 = resmod * xt/r * co

       call bbkick(pre1,pre2,pim1,pim2,xt,yt,sx1,sy1)
       expon = dexp(-xt**2/(2.d0*sx1**2)-yt**2/(2.d0*sy1**2))
       fr1 = ( pre1 - expon * pre2 )/ sqrt(sx*sx-sy*sy) 
       fi1 = ( pim1 - expon * pim2 ) / sqrt(sx*sx-sy*sy)

       fr = fr2 
       fi = fi2 
c 
       den = sqrt((40.*sx)**2+r*r)
       w1 = 40.*sx/den
       w2 = 1./den

c      image charge: use oliver's routine and subtract round-beam kick
c      and add flat beam kick instead 

          frc = 0.
          fic = 0.

          if(r2.gt.sx*sy) then
             efx = x/ r2
             efy = y/ r2

             efz = efx * zone + efy * zi
             z1 = x/rg * zone + y/rg * zi
             z2 = x/rg * zone - y/rg * zi
             z3 = cdlog(z1 - cdsqrt(z1**2 - zone))
             z4 = cdlog(z2 - cdsqrt(z2**2 - zone))
             rmi = dabs(dreal(z3+z4)/2.0d0)
             rphi = dreal((z4-z3)*zi)/2.0d0

c            introduce the same quantities for the beam
             z1b = xof/rg * zone + yof/rg * zi
             z2b = xof/rg * zone - yof/rg * zi
             z3b = cdlog(z1b - cdsqrt(z1b**2 - zone))
             z4b = cdlog(z2b - cdsqrt(z2b**2 - zone))
             rmib = dabs(dreal(z3b+z4b)/2.0d0)
             rphib = dreal((z4b-z3b)*zi)/2.0d0

             if(rphi*y.lt.0.0d0) rphi = -rphi
             if(rphib*yof.lt.0.0d0) rphib = -rphib

             do 101 i=1,nimage
                charge = 1.0d0 
             if ((xof.eq.0.0).and.(yof.eq.0.0)) then
                efz = efz+((-1)**i)*
     *                8.0d0*charge/(dexp(4.0d0*i*rm1)+1.0d0) * 
     *               (cdexp(2.0d0*i*(rmi*zone-rphi*zi)) - 
     *                cdexp(-2.0d0*i*(rmi*zone-rphi*zi)))/
     *               (cdexp((rmi*zone-rphi*zi)) - 
     *                cdexp(-(rmi*zone-rphi*zi)))/ rg
             else 
                efz = efz+2.d0*charge*dexp(-i*rm1)*(
     *             (dexp(i*rmib)+dexp(-i*rmib))*
     *             cos(i*rphib)/(dexp(i*rm1)+dexp(-i*rm1))+
     *             zi*(dexp(i*rmib)-dexp(-i*rmib))*
     *             sin(i*rphib)/(dexp(i*rm1)-dexp(-i*rm1)))*
     *             (cdexp(i*(rmi*zone-rphi*zi))-
     *             cdexp(-i*(rmi*zone-rphi*zi)))/
     *             (cdexp((rmi*zone-rphi*zi))-
     *             cdexp(-(rmi*zone-rphi*zi)))/rg
             endif 
 101            continue
             efx = dreal(efz)
             efy = dreal(-efz*zi)
             firs = abs(efx)*2.0d0/sqrt(2.0d0*pi)
             frrs = abs(efy)*2.0d0/sqrt(2.0d0*pi)

             cor = (1.0d0-exp(-r2/(2.0d0*sx*sy)))
             frr = res*co*yt
             fir = res*co*xt

             frc = frrs-frr
             fic = firs-fir
          endif


c       arg = sqrt(sx*sy)/5.
c       if (r.lt.arg) then
c       fr = fr1
c       fi = fi1 
c       else
       fr = w1 * fr1 + w2 * fr2 + frc 
       fi = w1 * fi1 + w2 * fi2 + fic 
c       endif       

       return

       endif 

c       write(*,*) "Message from foliver: The beams are not round!"

       end
 
c----------------------------------------------------------
       subroutine ffrank(x,y,xof,yof,sx,sy,fr,fi,ilin)
c      25.08.99 FZ, include offset xof, yof to be 
c               consistent with foliver 

       implicit real*8 (a-h,o-z)
c       implicit integer (a-z)
       double precision sx1,sy1,xt,yt
       real*8 x,y,sx,sy,fr,fi

       pi = 4.*datan(1.0d0)       

       if (ilin.eq.1) then       
       res = 2.0d0/sqrt(2.0d0*pi)/(sx+sy)
       fr = res*abs(y-yof)/sy
       fi = res*abs(x-xof)/sx
       return 
       endif 

       iflag = 0 
       if (sx.lt.sy) iflag = 1 
       if (iflag.eq.0) then
       xt = abs(x-xof)
       yt = abs(y-yof)
       else if (iflag.eq.1) then
       xt = abs(y-yof)
       yt = abs(x-xof)
       endif 

       r = sqrt(xt*xt+yt*yt)
       r2=r*r  
       resmod = 2./sqrt(2.*pi)
c       if (r.gt.(10.0d0*sx)) then
       pi = 4.0d0*atan(1.0d0)
       res = 2.0d0/sqrt(2.0d0*pi)/(r*r)
       if (iflag.eq.0) then
       sx1 = sx
       sy1 = sy
       else if (iflag.eq.1) then       
       sx1 = sy
       sy1 = sx
       endif       

c      round beam 
c       if (abs((sx/sy)-1.0d0).lt.5.0d-2.or.
c     &     r2.gt.(8.*sx1*sy1)) then
       if (abs((sx/sy)-1.0d0).lt.5.0d-2) then
       co = (1.0d0-exp(-r*r/(2.0d0*sx1*sy1)))
       fr = res*co*yt
       fi = res*co*xt
       return

       else  
c      flat beam        

c       write(*,*) ' HERE ', iflag 

       fr2 = resmod * yt/r * co
       fi2 = resmod * xt/r * co

c       write(*,*) ' before bbkick ', xt, yt, sx1, sy1  

       call bbkick(pre1,pre2,pim1,pim2,xt,yt,sx1,sy1)

c       write(*,*) ' after bbkick ', pre1, pre2, pim1, pim2

       expon = dexp(-xt**2/(2.d0*sx1**2)-yt**2/(2.d0*sy1**2))
       fr1 = ( pre1 - expon * pre2 )/ sqrt(sx1*sx1-sy1*sy1) 
       fi1 = ( pim1 - expon * pim2 ) / sqrt(sx1*sx1-sy1*sy1)

       fr = fr2 
       fi = fi2 
c 
       den = sqrt((40.*sx)**2+r*r)
       w1 = 40.*sx/den
       w2 = 1./den

       if (iflag.eq.0) then
       fr = w1 * fr1 + w2 * fr2 
       fi = w1 * fi1 + w2 * fi2 
       else if (iflag.eq.1) then
       fi = w1 * fr1 + w2 * fr2 
       fr = w1 * fi1 + w2 * fi2 
       endif 

       return 
       endif 
       end

c----------------------------------------------------------
c      Oliver 30.06.1997
c      This program calculates the kick from a TEM wave
c      excited in a coaxial transmission line.
       subroutine tem(x,y,wk,time,length,fr,fi)
       implicit real*8 (a-h,o-z)
c       implicit integer (a-z)
       real*8 xt,yt
       real*8 x,y,fr,fi,wk,time,length
       real*8 vl

       parameter(vl=2.9889d8)

       xt = abs(x)
       yt = abs(y)
       r = sqrt(xt*xt+yt*yt) 
       if(r.ne.0.0d0) then
          vcos = x/r
          vsin = y/r 
c       pi = 4.0d0*atan(1.0d0)
          pi = acos(-1.0d0)

          fi = sin(2*pi*mod(wk*vl*time/length,1.0d0))/r 
c          fi = 1.0d0/r
          fr = -vsin * fi
          fi = -vcos * fi
       endif
       if(r.eq.0.0d0) then
          fr = 0.0d0
          fi = 0.0d0
       endif

c       write(*,*) "TEM: ",mod(wk*vl*time/length,1.0d0),fi,1.0d0/r

       return
       end

c----------------------------------------------------------
c      Oliver 30.06.1997
c      This program calculates the kick from a TEM wave
c      excited in a coaxial transmission line.
       subroutine tem6(x,y,wk,time,length,rout,fr,fi,vt0)
       implicit real*8 (a-h,o-z)
c       implicit integer (a-z)
       real*8 xt,yt,vt0
       real*8 x,y,fr,fi,ftmp,wk,time,length
       real*8 vl,ftmp1
       real*8 rr(6), ri(6), x0r(6), y0r(6), vcosr(6), vsinr(6)
       real*8 x0i(6), y0i(6), rbound6, webin(0:500)
       real*8 webinphi(0:500)
       real*8 rin, rout, pi, vcosi(6), vsini(6)

       common / wavguid6 / x0r, y0r, rbound6, webin, webinphi
       
       parameter(vl=2.9889d8)
c       parameter(rin=0.01d0,rout=0.20)

       pi = dacos(-1.0d0)
       fr = 0.0d0
       fi = 0.0d0

c     Electric Field due to the six wires:
       do 10 i=1,6
          rr(i)=sqrt((x-x0r(i))*(x-x0r(i))+(y-y0r(i))*(y-y0r(i)))
          if(rr(i).gt.rbound6) then
             vcosr(i) = (x-x0r(i))/rr(i)
             vsinr(i) = (y-y0r(i))/rr(i) 
c             ftmp = sin(2*pi*mod(wk*vl*time/length,1.0d0))/rr(i) 
c
c     For square wave: (By X .ZHANG)
c     length->wavelength, wk->peroid over pulse duration
c
             ftmp1=mod(time*vl/length, 1.0d0)
             if(ftmp1.le.(1.0d0/wk)) then
                ftmp=(1.0d0+vt0)/rr(i)
               else
                ftmp=vt0/rr(i)
              endif
             fr = fr - vsinr(i) * ftmp
             fi = fi - vcosr(i) * ftmp
          endif
 10    continue

c     Electric Field due to six immage wires:
       do 20 i=1,6
          scale = sqrt(x0r(i)*x0r(i)+y0r(i)*y0r(i))/rout
          if(scale.eq.0.0) goto 20
             x0i(i) = x0r(i)/scale**2
             y0i(i) = y0r(i)/scale**2
             ri(i)=sqrt((x-x0i(i))*(x-x0i(i))+(y-y0i(i))*(y-y0i(i)))
             if(ri(i).gt.rbound6) then
                vcosi(i) = (x-x0i(i))/ri(i)
                vsini(i) = (y-y0i(i))/ri(i) 
c                ftmp = sin(2*pi*mod(wk*vl*time/length,1.0d0))/ri(i) 
c
c     For square wave: (By X. ZHANG)
c     length->wavelength, wk->peroid over pulse duration
c
             ftmp1=mod(time*vl/length, 1.0d0)
             if(ftmp1.le.(1.0d0/wk)) then
                ftmp=(1.0d0+vt0)/ri(i)
               else
                ftmp=vt0/ri(i)
              endif
                fr = fr + vsini(i) * ftmp
                fi = fi + vcosi(i) * ftmp
           endif
 20    continue


       fr = fr/5.4
       fi = fi/5.4
c       write(*,*) "TEM: ",mod(wk*vl*time/length,1.0d0),fi,1.0d0/r

       return
       end
 

c---------------------------------------------------------       
       subroutine ffranklin(x,y,sx,sy,fr,fi)
       implicit real*8 (a-h,o-z)       
       real x,y,sx,sy,fr,fi

       xt = abs(x)
       yt = abs(y)
       r = dsqrt(xt*xt+yt*yt) 
       pi = 4.*datan(1.d0)

       fr = 2./(sy*(sx+sy)) / sqrt(2.*pi) * 
     &      yt
       fi = 2./(sx*(sx+sy)) / sqrt(2.*pi) * 
     &      xt


       return
       end
c---------------------------------------------------------       
c      this routine converts into the first quadrant
c      and then calls Tong's/Talman's programs 
c 
       subroutine bbkick(pre1,pre2,pim1,pim2,x,y,sx,sy)
       implicit double precision(a-h,o-z)
       fac = dsqrt(2.d0*(sx**2-sy**2))
       r = sy/sx 
       u = x 
       v = y 
c       u = x / fac
c       v = y / r / fac
       u1 = dabs(u)
       v1 = dabs(v)

       u1a = u1 / fac
       v1a = v1 / fac 
       call errf(u1a,v1a,wr1,wi1)

       u1b = u1 / fac * r 
       v1b = v1 / fac / r 
       call errf(u1b,v1b,wr2,wi2)

       pre1 = wr1 
       pre2 = wr2
       pim1 = wi1 
       pim2 = wi2
  
       goto 900

       texp1 = 2.d0*dexp(-u**2+(r*v)**2)
       texp2 = 2.d0*dexp(-(r*u)**2+v**2)
       arguvr = 2.d0*u*r*v 

       if(u.gt.0.and.v.gt.0)then 
       call fbbclc(wr1,wr2,wi1,wi2,u,v,r)
       pre1 = wr1 
       pre2 = wr2
       pim1 = wi1 
       pim2 = wi2

       else if(u.gt.0.and.v.lt.0)then 
       call fbbclc(wr1,wr2,wi1,wi2,u,-v,r)
       pre1 = texp1*dcos(arguvr)-wr1 
       pre2 = texp2*dcos(arguvr)-wr2
       pim1 = texp1*dsin(arguvr)+wi1 
       pim2 = texp2*dsin(arguvr)+wi2

       else if(u.lt.0.and.v.gt.0)then 
       call fbbclc(wr1,wr2,wi1,wi2,-u,v,r)
       pre1 = wr1 
       pre2 = wr2
       pim1 = -wi1 
       pim2 = -wi2

       else if(u.lt.0.and.v.lt.0)then 
       call fbbclc(wr1,wr2,wi1,wi2,-u,-v,r)
       pre1 = texp1*dcos(arguvr)-wr1 
       pre2 = texp2*dcos(arguvr)-wr2
       pim1 = -texp1*dsin(arguvr)-wi1 
       pim2 = -texp2*dsin(arguvr)-wi2

       endif 

 900   continue 
       end 
c
c ####################################################################
c
c                                  2   2  2
c                             -(1-r )(v +u )
c      f(u,v,r) = w(u+irv) - e              w(ru+iv)
c
c ####################################################################
c
       subroutine fbbclc(wr1,wr2,wi1,wi2,u,v,r)
       implicit double precision(a-h,o-z)
       call errf(u,r*v,wr1,wi1)
       call errf(r*u,v,wr2,wi2)
       arg=(1.0-r*r)*(v*v+u*u)
       if(arg.gt.100.) then
         arg = 100.
         wr2 = 0.0
         wi2 = 0.0
       endif
       expon=exp(-arg)
c
       fr=wr1-expon*wr2
       fi=wi1-expon*wi2
       return
       end
c
c                 fnctnw.fortran
c
c          this program evaluates the function w(z)
c       (where z = zr +  zi) in the first quadrant of
c       the complex plane (i.e. zr > 0 and zi > 0).
c       three different expressions, pade1, pade2, and
c       asymp, are used, depending on where z lies in
c       the quadrant.
c
       subroutine fcnw(wr,wi,zr,zi)
       implicit double precision(a-h,o-z)
       data x1/4.1e0/,x2/3.6e0/,x3/3.5e0/,x4/2.7e0/,x5/2.2e0/
       data y1/1.275e0/,y2/1.095e0/
       data r2/8.7025e0/
        if(zr-x1)200,30,30
200        eps1=.0625e0*(zr-x3)
           if(zr-x2)300,210,210
210          yc=-1.4e0*(zr-x2)+y2
          if(zi-yc)220,30,30
220       if(zr*zi.lt.eps1)then
          call asymp(wr,wi,zr,zi)
          else
          call pade2(wr,wi,zr,zi)
          endif
300        if(zr-x4)400,310,310
310          yc=-.2e0*(zr-x4)+y1
          if(zi-yc)320,30,30
320       if(zr.ge.x3.and.zr*zi.lt.eps1) then
          call asymp(wr,wi,zr,zi)
          else
          call pade2(wr,wi,zr,zi)
          endif
400        if(zr-x5)500,410,410
410          yc1=-1.4e0*(zr-x4)+y1
          yc2=1.75e0*(zr-x4)+y1
          if(zi-yc1)420,30,30
420          if(zi-yc2)20,10,10
500        if(zr*zr+zi*zi-r2)10,30,30
10        call pade1(wr,wi,zr,zi)
        return
20        call pade2(wr,wi,zr,zi)
        return
30        call asymp(wr,wi,zr,zi)
        return
       end
     
c                 wasymp.fortran
c
c          this program calculates an asymptotic expression of
c       w(z) valid away from the origin.
c
       subroutine asymp(wr,wi,zr,zi)
       implicit double precision(a-h,o-z)
       data a1p/1.94443615e-1/,a2p/7.64384940e-2/,
     1   a3p/1.07825546e-2/,a4p/4.27695730e-4/,a5p/2.43202531e-6/
       data b1/3.42901327e-1/,b2/1.036610830e0/,b3/1.756683649e0/,
     1   b4/2.532731674e0/,b5/3.436159119e0/
       data pi2/1.12837917e0/
       data x1/3.5e0/,x2/4.2e0/
       data eps/.01e0/,check/0.e0/
10        dr1=zr+b1
        d1r=zr-b1
        dr2=zr+b2
        d2r=zr-b2
        dr3=zr+b3
        d3r=zr-b3
        dr4=zr+b4
        d4r=zr-b4
        dr5=zr+b5
        d5r=zr-b5
        de1=dr1*dr1+zi*zi
        d1e=d1r*d1r+zi*zi
        de2=dr2*dr2+zi*zi
        d2e=d2r*d2r+zi*zi
        de3=dr3*dr3+zi*zi
        d3e=d3r*d3r+zi*zi
        de4=dr4*dr4+zi*zi
        d4e=d4r*d4r+zi*zi
        de5=dr5*dr5+zi*zi
        d5e=d5r*d5r+zi*zi
        if(1.e0-check)70,70,20
20          if(zr.lt.x1) goto 60
            eps1=.04e0/(zr-3.29e0)-.034e0
            if(zr*zi.lt.eps1) goto 50
            if(.not.((zr.ge.x2).and.(zr*zi.lt.eps))) goto 60
50          check=1.e0
          wi0=a1p*(dr1/de1+d1r/d1e)+a2p*(dr2/de2+d2r/d2e)+
     1         a3p*(dr3/de3+d3r/d3e)+a4p*(dr4/de4+d4r/d4e)+
     1         a5p*(dr5/de5+d5r/d5e)
          zi0=zi
          zi=0.e0
          go to 10
60        wr=(a1p*(1.e0/de1+1.e0/d1e)+a2p*(1.e0/de2+1.e0/d2e)+
     1      a3p*(1.e0/de3+1.e0/d3e)+a4p*(1.e0/de4+1.e0/d4e)+
     1      a5p*(1.e0/de5+1.e0/d5e))*zi
70        wi=a1p*(dr1/de1+d1r/d1e)+a2p*(dr2/de2+d2r/d2e)+
     1      a3p*(dr3/de3+d3r/d3e)+a4p*(dr4/de4+d4r/d4e)+
     1      a5p*(dr5/de5+d5r/d5e)
        if(1.e0-check)80,80,90
80        wr=exp(-zr*zr)+2.e0*wi*zr*zi0-pi2*zi0
        wi=wi0
        zi=zi0
        check=0.e0
90        return
c
c  ##################################################################
c  ######################################################################
       entry errasymp
c
       write(41,*)'zi,zr,ck=',zi,zr,ck
       return
c
       end
     
c                 wpade1.fortran
c
c          this program calculates a pade approximation of w(z)
c       around the origin.
c
       subroutine pade1(wr,wi,zr,zi)
       implicit double precision(a-h,o-z)
       data c1/-1.25647718e0/,c2/8.25059158e-1/,
     1   c3/-3.19300157e-1/,c4/7.63191605e-2/,
     1   c5/-1.04697938e-2/,c6/6.44878652e-4/
       data d1/-2.38485635e0/,d2/2.51608137e0/,
     1   d3/-1.52579040e0/,d4/5.75922693e-1/,
     1   d5/-1.35740709e-1/,d6/1.85678083e-2/,
     1   d7/-1.14243694e-3/
        u2r=zi*zi-zr*zr
        u2i=-2.e0*zr*zi
        u3r=-u2r*zi-u2i*zr
        u3i=u2r*zr-u2i*zi
        u4r=-u3r*zi-u3i*zr
        u4i=u3r*zr-u3i*zi
        u5r=-u4r*zi-u4i*zr
        u5i=u4r*zr-u4i*zi
        u6r=-u5r*zi-u5i*zr
        u6i=u5r*zr-u5i*zi
        u7r=-u6r*zi-u6i*zr
        u7i=u6r*zr-u6i*zi
        fr=1.e0-c1*zi+c2*u2r+c3*u3r+c4*u4r+c5*u5r+c6*u6r
        fi=c1*zr+c2*u2i+c3*u3i+c4*u4i+c5*u5i+c6*u6i
        dr=1.e0-d1*zi+d2*u2r+d3*u3r+d4*u4r+d5*u5r+d6*u6r+d7*u7r
        di=d1*zr+d2*u2i+d3*u3i+d4*u4i+d5*u5i+d6*u6i+d7*u7i
        de=dr*dr+di*di
         wr=(fr*dr+fi*di)/de
        wi=(fi*dr-fr*di)/de
       return
       end
     
c                 wpade2.fortran
c
c          this program calculates a pade approximation of w(z)
c       around the point z = 3.
c
       subroutine pade2(wr,wi,zr,zi)
       implicit double precision(a-h,o-z)
       data c0r/1.23409804e-4/,c0i/2.01157318e-1/,
     1   c1r/2.33746715e-1/,c1i/1.61133338e-1/,
     1   c2r/1.25689814e-1/,c2i/-4.0422725e-2/,
     1   c3r/8.92089179e-3/,c3i/-1.81293213e-2/
       data d1r/1.19230984e0/,d1i/-1.16495901e0/,
     1   d2r/8.9401545e-2/,d2i/-1.07372867e0/,
     1   d3r/-1.68547429e-1/,d3i/-2.70096451e-1/,
     1   d4r/-3.20997564e-2/,d4i/-1.58578639e-2/
        zr=zr-3.e0
        z2r=zr*zr-zi*zi
        z2i=2.e0*zr*zi
        z3r=z2r*zr-z2i*zi
        z3i=z2r*zi+z2i*zr
        z4r=z3r*zr-z3i*zi
        z4i=z3r*zi+z3i*zr
        fr=c0r+c1r*zr-c1i*zi+c2r*z2r-c2i*z2i+c3r*z3r-c3i*z3i
        fi=c0i+c1r*zi+c1i*zr+c2r*z2i+c2i*z2r+c3r*z3i+c3i*z3r
        dr=1.e0+d1r*zr-d1i*zi+d2r*z2r-d2i*z2i+d3r*z3r-d3i*z3i+
     1      d4r*z4r-d4i*z4i
        di=d1r*zi+d1i*zr+d2r*z2i+d2i*z2r+d3r*z3i+d3i*z3r+d4r*z4i+
     1      d4i*z4r
        de=dr*dr+di*di
        wr=(fr*dr+fi*di)/de
        wi=(fi*dr-fr*di)/de
        zr=zr+3.e0
       return
       end



      SUBROUTINE ERRF(XX, YY, WX, WY)                                      
*----------------------------------------------------------------------*   
* Purpose:                                                             *   
*   Modification of WWERF, double precision complex error function,    *   
*   written at CERN by K. Koelbig.                                     *   
* Input:                                                               *   
*   XX, YY    (real)    Argument to CERF.                              *   
* Output:                                                              *   
*   WX, WY    (real)    Function result.                               *   
*----------------------------------------------------------------------*   
                                                                           
*---- Double precision version.                                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)                   
      PARAMETER         (MWFLT = 2, MREAL = 4)                             
      PARAMETER         (MCWRD = 4)                                        
      PARAMETER         (MCNAM = 16, MWNAM = MCNAM / MCWRD)                
      PARAMETER         (MCFIL = 80, MCRNG = 40, MCSTR = 80)               
                                                                           
      PARAMETER         (CC     = 1.12837 91670 9551D0)                    
      PARAMETER         (ONE    = 1.D0)                                    
      PARAMETER         (TWO    = 2.D0)                                    
      PARAMETER         (XLIM   = 5.33D0)                                  
      PARAMETER         (YLIM   = 4.29D0)                                  
      DIMENSION         RX(33), RY(33)                                     
                                                                           
      X = ABS(XX)                                                          
      Y = ABS(YY)                                                          
                                                                           
      IF (Y .LT. YLIM  .AND.  X .LT. XLIM) THEN                            
        Q  = (ONE - Y / YLIM) * SQRT(ONE - (X/XLIM)**2)                    
        H  = ONE / (3.2D0 * Q)                                             
        NC = 7 + INT(23.0*Q)                                               
        XL = H**(1 - NC)                                                   
        XH = Y + 0.5D0/H                                                   
        YH = X                                                             
        NU = 10 + INT(21.0*Q)                                              
        RX(NU+1) = 0.                                                      
        RY(NU+1) = 0.                                                      
                                                                           
        DO 10 N = NU, 1, -1                                                
          TX = XH + N * RX(N+1)                                            
          TY = YH - N * RY(N+1)                                            
          TN = TX*TX + TY*TY                                               
          RX(N) = 0.5D0 * TX / TN                                          
          RY(N) = 0.5D0 * TY / TN                                          
   10   CONTINUE                                                           
                                                                           
        SX = 0.                                                            
        SY = 0.                                                            
                                                                           
        DO 20 N = NC, 1, -1                                                
          SAUX = SX + XL                                                   
          SX = RX(N) * SAUX - RY(N) * SY                                   
          SY = RX(N) * SY + RY(N) * SAUX                                   
          XL = H * XL                                                      
   20   CONTINUE                                                           
                                                                           
        WX = CC * SX                                                       
        WY = CC * SY                                                       
      ELSE                                                                 
        XH = Y                                                             
        YH = X                                                             
        RX(1) = 0.                                                         
        RY(1) = 0.                                                         
                                                                           
        DO 30 N = 9, 1, -1                                                 
          TX = XH + N * RX(1)                                              
          TY = YH - N * RY(1)                                              
          TN = TX*TX + TY*TY                                               
          RX(1) = 0.5D0 * TX / TN                                          
          RY(1) = 0.5D0 * TY / TN                                          
   30   CONTINUE                                                           
                                                                           
        WX = CC * RX(1)                                                    
        WY = CC * RY(1)                                                    
      ENDIF                                                                
                                                                           
      IF(Y .EQ. 0.) WX = EXP(-X**2)                                        
      IF(YY .LT. 0.) THEN                                                  
        WX =   TWO * EXP(Y*Y-X*X) * COS(TWO*X*Y) - WX                      
        WY = - TWO * EXP(Y*Y-X*X) * SIN(TWO*X*Y) - WY                      
        IF(XX .GT. 0.) WY = -WY                                            
      ELSE                                                                 
        IF(XX .LT. 0.) WY = -WY                                            
      ENDIF                                                                
                                                                           
      END




