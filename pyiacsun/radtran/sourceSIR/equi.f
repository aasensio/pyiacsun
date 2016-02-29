c equi rutina que evalua la presion en equilibrio hidrostatico

	subroutine equi(ntau,tau,t,pe)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'  !solo por kt
	parameter (nex=28)
	real*4 tau(*),t(*),pe(*),x(kt),kap(kt),pg(kt),pgold,kac,d2,d3
	real wgt,abu,ei1,ei2,pp(10),tsi,psi,psg,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
	common/mu/cth
        common/preciso/prec      
        common/pinicial/pcontorno

	g=2.7414e+4*cth		!gravedad cm/s^2 en fotosfera solar   
        negativo=0 
	avog=6.023e23
        do i=1,10
           d1(i)=0
        end do
c calculamos el peso molecular medio pmu
	pmu=0.0
	asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmu=pmu+wgt*abu
	   asum=asum+abu
	end do
	pmu=pmu/asum
	avog=avog/pmu

c calculamos la pe correspondiente al contorno
        i0=1
        taumin=1.e3 
        do i=1,ntau
           if(abs(tau(i)).lt.taumin)then
              taumin=tau(i)
              i0=i
           end if
           x(i)=(10.)**(tau(i)) 
        end do
 
        tcont=t(i0)
        pecont=pe(i0) 
        if(pcontorno.gt.1.e-5)call pefrompg1(tcont,pcontorno,pecont)

c calculamos kap(i0)
        tsi=t(i0)
	psi=pecont
        call gasc(tsi,psi,psg,pp)
        pg(i0)=psg
        pe(i0)=pecont

	call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)

        kap(i0)=kac*avog

c supongo log(pg) lineal con log(tau)  y por tanto 
c pg*kap/(pg'*kap')=x/x'
c y haciendo pg'=pg+2*g*(x'-x)/(kap'+kap) se llega a una ec. de segundo
c grado para kap'

        do i=i0,ntau-1
c          kap(i+1)=kap(i)
           b=(x(i+1)/x(i)-1.d0)*(2.d0*g*x(i)/pg(i) - kap(i))
           kap(i+1)=-b+sqrt(b*b+kap(i)*kap(i)*x(i+1)/x(i)) 
           paso=x(i+1)-x(i)
           dif=1.e10
           n=0
           do while(dif.ge.prec.and.n.le.50)
              pg(i+1)=pg(i)+2.*g*paso/ (kap(i+1)+kap(i))
              nrep=0
              do while(pg(i+1).le.0.and.nrep.le.50)
                 kap(i+1)=2.*kap(i+1)
                 pg(i+1)=pg(i)+2.*g*paso/ (kap(i+1)+kap(i))
                 nrep=nrep+1
              end do 
              if(nrep.eq.50)then
                negativo=1
                return
              end if 
	      n=n+1
              pe(i+1)=pe(i)

              pgold=pg(i+1)
              call pefrompg1(t(i+1),pg(i+1),pe(i+1))
              tsi=t(i+1)
	      psi=pe(i+1)
	      call gasc(tsi,psi,psg,pp)
	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)
              kap(i+1)=kac*avog
              pg(i+1)=pg(i)+2.*g*paso/ (kap(i+1)+kap(i))
              dif=abs(pgold/pg(i+1)-1.)
c              print*,'i+1=',i+1,' n=',n,' dif=',dif,' pg=',pg(i+1)
           end do
           if(n.gt.50)then 
                print*,'WARNING: The maximum number of iterations for computing the electronic pressure'
                print*,'         from hydrostatic equilibrium has been reached.'
                print*,'         The error at tau(',i+1,') >=',dif*100.,' %.'
	   endif
        end do

        do i=i0,2,-1
c          kap(i-1)=kap(i)
           b=(x(i-1)/x(i)-1.d0)*(2.d0*g*x(i)/pg(i) - kap(i))
           kap(i-1)=-b+sqrt(b*b+kap(i)*kap(i)*x(i-1)/x(i)) 
           paso=x(i-1)-x(i)
           n=0
           dif=1.e10
           do while(dif.ge.prec.and.n.le.50)
              pg(i-1)=pg(i)+2.*g*paso/ (kap(i-1)+kap(i))
              pgold=abs(pg(i-1))
              n=n+1
              pe(i-1)=pe(i) 
              call pefrompg1(t(i-1),pg(i-1),pe(i-1))
              tsi=t(i-1)
	      psi=pe(i-1)
	      call gasc(tsi,psi,psg,pp)
	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)
              kap(i-1)=kac*avog
              pg(i-1)=pg(i)+2.*g*paso/ (kap(i-1)+kap(i))
              dif=abs(pgold/pg(i-1)-1.)
c              print*,'i=',i,' n=',n,' dif=',dif,' pg=',pg(i-1)
           end do
           if(n.gt.50)then 
                print*,'WARNING: The maximum number of iterations for computing the electronic pressure'
                print*,'         from hydrostatic equilibrium has been reached.'
                print*,'         The error at tau(',i+1,') >=',dif*100.,' %.'
	   endif
       end do
       return
       end 
c _______________________________________________________________
c equi2 rutina que evalua la presion en equilibrio hidrostatico
c _______________________________________________________________
	subroutine equi2(ntau,tau,t,pe)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'   !solo por kt

	parameter (nex=28)
	real*4 tau(*),t(*),pe(*),x(kt),kap(kt),pg(kt),pgold,kac,d2,d3
	real wgt,abu,ei1,ei2,pp(10),tsi,psi,psg,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
	common/mu/cth
        common/preciso/prec      
        common/pinicial/pcontorno

	g=2.7414e+4*cth		!gravedad cm/s^2 en fotosfera solar   
	avog=6.023e23
        do i=1,10
           d1(i)=0
        end do
c calculamos el peso molecular medio pmu
	pmu=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmu=pmu+wgt*abu
           asum=asum+abu
	end do
        pmu=pmu/asum
	avog=avog/pmu

c calculamos la pe correspondiente al contorno
        do i=1,ntau
           x(i)=10.**(tau(i)) 
        end do
 
        tcont=t(ntau)
        pecont=pe(ntau) 
        if(pcontorno.gt.1.e-5)call pefrompg1(tcont,pcontorno,pecont)

c calculamos kap(ntau)
        tsi=t(ntau)
	psi=pecont
        call gasc(tsi,psi,psg,pp)
        pg(ntau)=psg
        pe(ntau)=pecont

	call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)

        kap(ntau)=kac*avog/asum  !(k /unidad long. )/ro

c supongo log(pg) lineal con log(tau)  y por tanto 
c pg*kap/(pg'*kap')=x/x'
c y haciendo pg'=pg+2*g*(x'-x)/(kap'+kap) se llega a una ec. de segundo
c grado para kap'

        do i=ntau,2,-1
           b=(x(i-1)/x(i)-1.d0)*(2.d0*g*x(i)/pg(i) - kap(i))
           kap(i-1)=-b+sqrt(b*b+kap(i)*kap(i)*x(i-1)/x(i)) 
           if(kap(i-1).le.0)kap(i-1)=kap(i)
           paso=x(i-1)-x(i)
           n=0
           dif=1.e10
           tsi=t(i-1)

           do while(dif.ge.prec.and.n.le.50)
              pg(i-1)=pg(i)+2.*g*paso/ (kap(i-1)+kap(i))
              n=n+1
              pe(i-1)=pe(i)
              pgold=pg(i-1) 

              call pefrompg1(t(i-1),pg(i-1),pe(i-1))
	      psi=pe(i-1)
	      call gasc(tsi,psi,psg,pp)
	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)
              kap(i-1)=kac*avog/asum
              pg(i-1)=pg(i)+2.*g*paso/ (kap(i-1)+kap(i))
              call pefrompg1(t(i-1),pg(i-1),pe(i-1))

	      psi=pe(i-1)
c              pgold=pg(i-1) 
	      call gasc(tsi,psi,psg,pp)
	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d3)
              kap(i-1)=kac*avog/asum
              pg(i-1)=pg(i)+2.*g*paso/ (kap(i-1)+kap(i))
              call pefrompg1(t(i-1),pg(i-1),pe(i-1))

              dif=abs(pgold/pg(i-1)-1.)
c              print*,'i=',i-1,' n=',n,' dif=',dif,' pg=',pg(i-1),'pe=',pe(i-1)
           end do
           if(n.gt.50)then 
                print*,'WARNING: The maximum number of iterations for computing the electronic pressure'
                print*,'         from hydrostatic equilibrium has been reached.'
                print*,'         The error at tau(',i+1,') >=',dif*100.,' %.'
	   endif
       end do
       return
       end 

                        
c____________________________________________________________________
c pefrompg1 evalua la presion electonica p correspondiente a t1 y pg

	subroutine pefrompg1(t,pg,p)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER' !solo por kt
        parameter (nex=28) 
        common/preciso/prec      
 
        epsilon=prec/10.

        dif=1.d0
        n2=0
        do while (dif.gt.epsilon.and.n2.lt.250)
	     n2=n2+1
	     p1=p
             call pe_pg1(t,p,pg)
             dif=abs(p/p1-1.0)
         end do
	return
	end
c____________________________________________________________________

c calcula la presion electronica a partir de la pg y de una estimacion de la pe

      subroutine pe_pg1(t,pe,pg)

      parameter (ncontr=28)
      dimension cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),
     *u0(ncontr),u1(ncontr),u2(ncontr)
	real du0,du1,du2,dcmol(91)

      if(t.lt.500)then
         print*,'pe_pg1: temperature < 500 K '
	 print*,'temperature = 500 K'
	 t=500.
      end if	 
 
      theta=5040./t
      g4=0.
      g5=0.
      if(pe.le.0)then
         pe=1.e-15
         g4=0.
         g5=0.
      else
         call molecb(theta,cmol,dcmol)
         do i=1,2
            call acota(cmol(i),-30.,30.)
         end do
         g4=pe*10.**(cmol(1))
         g5=pe*10.**(cmol(2))
      end if 
c ahora calculo los niveles u0,u1,u2 y sus derivadas
      do 5 i=1,ncontr
      		iii=i
5     		call neldatb(iii,0.,weight,alfai(i),chi1(i),chi2(i))
6     do 4 i=1,ncontr
      	  iii=i
4     	  call nelfctb(iii,t,u0(iii),u1(iii),u2(iii),du0,du1,du2)

      
      g2=saha(theta,chi1(1),u0(1),u1(1),pe)   ! p(h+)/p(h)
      
      g3=saha(theta,0.754,1.,u0(1),pe)        ! p(h)/p(h-) 

      call acota(g3,1.e-30,1.e30)
      g3=1.d0/g3                              ! p(h-)/p(h) 
     
      g1=0.
      do 1 i=2,ncontr
        a=saha(theta,chi1(i),u0(i),u1(i),pe)
        b=saha(theta,chi2(i),u1(i),u2(i),pe)
	c=1.+a*(1.+b)
1	g1=g1+alfai(i)/c*a*(1.+2.*b)
      
        a=1.+g2+g3
        b=2.*(1.+g2/g5*g4)
        c=g5
        d=g2-g3
        e=g2/g5*g4

	call acotasig(a,1.e-15,1.e15)
	call acotasig(d,1.e-15,1.e15)

        c1=c*b**2+a*d*b-e*a**2

        c2=2.*a*e-d*b+a*b*g1
	
        c3=-(e+b*g1)

        f1=0.5*c2/c1

        f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1)
        f5=(1.-a*f1)/b
        f4=e*f5
        f3=g3*f1
        f2=g2*f1
        fe=f2-f3+f4+g1

        call acota(fe,1.e-30,1.e30)
	phtot=pe/fe

        if(f5.gt.1.e-4) goto 2
          const6=g5/pe*f1**2
          const7=f2-f3+g1
	  do 3 i=1,5
      	       f5=phtot*const6
	       f4=e*f5
	       fe=const7+f4
	       phtot=pe/fe
3	  continue
      
2       pe=pg/(1.+(f1+f2+f3+f4+f5+0.1014)/fe)
        if(pe.le.0)pe=1.e-15
  
      return
      end

c ____________________________________________________________________________

	subroutine acota(x,x0,x1)

	if(x.lt.x0)x=x0	
        if(x.gt.x1)x=x1
   
        return
        end 
	subroutine acotasig(x,x0,x1)

        if(x.lt.0)then
           x=-x
           call acota(x,x0,x1)
           x=-x
         else
           call acota(x,x0,x1)
         end if
        return
        end 
