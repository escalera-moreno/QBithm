c
      program qbithm
c
      implicit double precision (a-h,o-z)
c
c     Parameters
c
      integer*8,parameter :: id = 7   ! (2J+1)*(2I+1) (=>2)
      integer*8,parameter :: ig = 3   ! qubit ground state (1<=ig<=id-1)
      integer*8,parameter :: ie = 6   ! qubit excited state (ig+1<=ie<=id)
      integer*8,parameter :: nm = 0   ! number of vibrational modes (=>0)
      integer*8,parameter :: ic = 1   ! number of static magnetic field values (=>1)
      integer*8,parameter :: nd = 1   ! number of static magnetic field directions (=>1)
c
      real*8,parameter :: temp  = 5.00d0  ! temperature (>0) (K)
      real*8,parameter :: sfgw  = 1.00d0  ! scaling factor for gaussian width (>0)
      real*8,parameter :: top   = 5.00d0  ! writing threshold for one-phonon relaxation modes (0<top<=100) (%) 
      real*8,parameter :: ttp   = 5.00d0  ! writing threshold for two-phonon relaxation modes (0<ttp<=100) (%) 
      real*8,parameter :: geme  = 0.000d0 ! additional emission rate (=>0) (1/microsecond)
      real*8,parameter :: gabe  =-1.000d0 ! additional absorption rate (=>0) (1/microsecond) (detailed balance if -1.0d0)  
      real*8,parameter :: tmage = 1.000d-10 ! additional magnetic rate (>0) (1/microsecond)
c
      real*8,parameter    :: gfi  = 2.000d0 ! free-ion Landé factor (>0)
      real*8,parameter    :: bcm  = 1.500d0 ! oscillating magnetic field magnitude (=>0) (mT)
      real*8,parameter    :: firr = 8.99377d0 ! irradiation frequency (=>0) (GHz)
      real*8,parameter    :: alp  = 0.0d0   ! linear polarization angle (0<=alp<360) (º)
      integer*8,parameter :: nang = 1       ! number of rotation directions to integrate (=>1)
c
      real*8,parameter    :: esta = 0.00d0 ! minimum delay/rotation time (=>0) (microsecond)
      real*8,parameter    :: eend = 0.50d0 ! maximum delay/rotation time (>esta) (microsecond)
      integer*8,parameter :: npe  = 201    ! number of magnetization calculated values (=>2)
c
      real*8,parameter    :: ro11  = 0.0d0 ! initial ro(1,1) (0<=ro11<=1)
      real*8,parameter    :: ro22  = 1.0d0 ! initial ro(2,2) with ro(1,1)+ro(2,2)=1 (0<=ro22<=1)
      real*8,parameter    :: ro12r = 0.0d0 ! initial Re(ro(1,2))
      real*8,parameter    :: ro12i = 0.0d0 ! initial Im(ro(1,2))
      integer*8,parameter :: nsa   = 2     ! number of algorithm steps (=>1)
c
c     Constants
c
      real*8,parameter :: pi   = 3.1415926535898d0 ! pi number
      real*8,parameter :: cmub = 9.274009994d-24   ! Bohr magneton (J/T)
      real*8,parameter :: cspe = 2.99792458d10     ! light speed (cm/s)
      real*8,parameter :: hbar = 1.054571817d-34   ! reduced Planck constant (J·s)
      real*8,parameter :: anum = 6.02214076d23     ! Avogadro number (1/mol)
      real*8,parameter :: bol  = 6.950348004d-1    ! Boltzmann constant (cm^-1/K)
c
      integer*8 :: mf,ip,ir,im,is,ia,ij,inpe,inum0,inum3
c
      complex*16 :: opp,tpv,tpda,tpst,tpas
c
      dimension fr(nm),rm(nm),hawi(nm),azwr(ic,nd,4),ene(ic,nd,id)
      dimension rpso(ic,nd,3,2),opp(ic,nd,nm)
      dimension ggg(ic,nd,nang,2),ggr(ic,nd,nang,nsa,2)
      dimension tpv(ic,nd,nm,nm),tpda(ic,nd,nm,nm,ie-ig-1)
      dimension tpst(ic,nd,nm,nm,id-ie),tpas(ic,nd,nm,nm,ig-1)
      dimension wei(nang),ang(nsa,nang),iste(nsa),dura(nsa),rodi(nsa)
      dimension omz(npe,ic+1),omxy(npe,ic+1),cir(4),rmat(4,4),sfr(4)
c
      bcmt=bcm*1.0d-3  ! mT  -> T
      firrh=firr*1.0d9 ! GHz -> Hz
c
      call rdata(ic,id,ig,ie,nm,nd,fr,rm,sfgw,hawi,azwr,
     &ene,rpso,opp,tpv,tpda,tpst,tpas,nsa,iste,dura,rodi)
      azwr(:,:,4)=(azwr(:,:,4)+tmage)*1.0d6 ! 1/mus -> 1/s
      dura(:)=dura(:)*1.0d-9                ! ns    -> s
      call cirrule(pi,nsa,rodi,nang,wei,ang)     
c
      open(7,file='qb.out',status='unknown')
      ii=7;call lab(ii)
      open(8,file='qb.mz.out',status='unknown')
      ii=8;call lab(ii)
      write(8,*) '      time (mus)            magnetization (a.u.)'
      open(9,file='qb.mxy.out',status='unknown')
      ii=9;call lab(ii)
      write(9,*) '      time (mus)            magnetization (a.u.)'
c
      if (nm.gt.0) write(7,*) '   Mode contributions to transition rates
     &'
      write(7,*) 
      do 8 mf=1,ic
        do 1 ip=1,nd
          write(7,103) mf,azwr(mf,ip,1)*(180.0d0/pi),azwr(mf,ip,2)*
     &(180.0d0/pi)
          call fopr(hbar,anum,pi,nm,ic,nd,fr,rm,bol,temp,
     &mf,ip,id,ig,ie,opp,hawi,ene,top,gaop,geop)
          call sopr(hbar,anum,pi,cspe,nm,ic,nd,fr,rm,bol,temp,
     &mf,ip,id,ig,ie,tpv,tpda,tpst,tpas,hawi,ene,ttp,gatp,getp)
          gem=geop+getp+geme*1.0d6
          gabeaux=gabe
          if (gabe.lt.-0.1d0) gabeaux=geme*dexp(-(ene(mf,ip,ie)-
     &ene(mf,ip,ig))/(bol*temp))
          gabs=gaop+gatp+gabeaux*1.0d6
          write(7,104) geop*1.0d-6,getp*1.0d-6,geme
          write(7,105) gaop*1.0d-6,gatp*1.0d-6,gabeaux
          write(7,*) 
          do 2 ir=1,nang
            ggg(mf,ip,ir,1)=gabs;ggg(mf,ip,ir,2)=gem
            do 11 is=1,nsa
              call rabfr(gfi,ic,nd,mf,ip,azwr,nsa,nang,ir,is,ang,rpso,
     &hbar,bcmt,cmub,alp,pi,prafr,qrafr)
              ggr(mf,ip,ir,is,1)=prafr;ggr(mf,ip,ir,is,2)=qrafr
 11         continue
 2        continue
 1      continue
 8    continue
      write(7,*) 
c
      inum0=0;inum3=0
      do 7 ij=1,nsa
        if (iste(ij).eq.0) then
          inum0=inum0+1
        endif
        if (iste(ij).eq.3) then
          inum3=inum3+1
        endif
 7    continue
      if ((inum0.eq.0).and.(inum3.eq.0)) then
        inpe=0
      else
        inpe=npe-1
      endif
c
      write(7,*) '   Time evolution of density matrix'
      write(7,*) 
      do 9 mf=1,ic
        write(7,99) mf
        write(7,*) 
        do 3 im=0,inpe
          if ((inum0.eq.0).and.(inum3.eq.0)) then
            tau=0.0d0
          else 
            tau=esta+dble(im)*((eend-esta)/dble(inpe))
          endif
          taus=tau*1.0d-6 ! mus -> s
          omz(im+1,1)=tau*(inum0+inum3)
          omxy(im+1,1)=tau*(inum0+inum3)
          write(7,100) im+1
          write(7,*) 
          eczeaz=0.0d0;ecxyeaz=0.0d0
          do 4 ip=1,nd
            ecze=0.0d0;ecxye=0.0d0
            do 5 ir=1,nang
              gabs=ggg(mf,ip,ir,1);gem=ggg(mf,ip,ir,2)
              cir(1)=ro11;cir(2)=ro22;cir(3)=ro12r;cir(4)=ro12i
              gap=ene(mf,ip,ie)-ene(mf,ip,ig);delta=(gap*cspe)-firrh
              write(7,101) azwr(mf,ip,1)*(180.0d0/pi),azwr(mf,ip,2)*
     &(180.0d0/pi),(gap*cspe)/1.0d9,delta/1.0d9
              write(7,*) '   step       duration(ns)   rot.dir.(º)   g.
     &R.f(MHz)   ro11      ro22     ro12r     ro12i'
              do 6 ia=1,nsa
                prafr=ggr(mf,ip,ir,ia,1);qrafr=ggr(mf,ip,ir,ia,2)
                if (iste(ia).eq.0) then
                  prafra=0.0d0;qrafra=0.0d0;firrha=0.0d0;time=taus
                elseif (iste(ia).eq.1) then
                  prafra=0.0d0;qrafra=0.0d0;firrha=0.0d0;time=dura(ia)
                elseif (iste(ia).eq.2) then
                  prafra=prafr;qrafra=qrafr;firrha=firrh;time=dura(ia)
                elseif (iste(ia).eq.3) then
                  prafra=prafr;qrafra=qrafr;firrha=firrh;time=taus
                endif
                delta=2.0d0*pi*((gap*cspe)-firrha)
                grf=dsqrt((prafra**2)+(qrafra**2)+(delta**2))/
     &(2.0d0*pi*1.0d6)
                call cmd(rmat,gabs,gem,azwr,prafra,qrafra,
     &pi,cspe,ic,nd,mf,ip,id,ig,ie,ene,firrha)
                call diag(rmat,cir,sfr,time,info,iflag)
                cir(:)=sfr(:)
                if ((iste(ia).eq.0).or.(iste(ia).eq.1)) then
                  write(7,102) ia,time*1.0d9,ang(ia,ir)*(180.0d0/pi),sfr
                elseif ((iste(ia).eq.2).or.(iste(ia).eq.3)) then
                  write(7,106) ia,time*1.0d9,ang(ia,ir)*(180.0d0/pi),
     &grf,sfr
                endif
 6            continue
              ecze=ecze+wei(ir)*(sfr(2)-sfr(1))
              ecxye=ecxye+wei(ir)*2.0d0*dsqrt((sfr(3)**2)+(sfr(4)**2))
              write(7,*) 
 5          continue
            eczeaz=eczeaz+azwr(mf,ip,3)*ecze
            ecxyeaz=ecxyeaz+azwr(mf,ip,3)*ecxye
 4        continue
          omz(im+1,mf+1)=eczeaz
          omxy(im+1,mf+1)=ecxyeaz
 3      continue
 9    continue
c
      do 10 im=0,inpe
        write(8,*) (omz(im+1,mf),mf=1,ic+1)
        write(9,*) (omxy(im+1,mf),mf=1,ic+1)
 10   continue
c
      close(7)
      close(8)
      close(9)
c
  99  format(' ---magnetic field: ',i8)
 100  format(' -magnetization: ',i8)
 101  format('  az(º): ',f10.5,'; ze(º): ',f10.5,'; Gap(GHz): ',f10.5,
     &'; Detuning(GHz): ',f10.5)
 102  format(i8,f19.3,'    ',f10.5,'            -',4f10.5)
 106  format(i8,f19.3,'    ',f10.5,f13.5,4f10.5)
 103  format('  mf: ',i8,'; az(º): ',f10.5,'; ze(º): ',f10.5)
 104  format('  G_em,1ph(1/mus): ',es16.8e3,' G_em,2ph(1/mus): ',
     &es16.8e3,' G_em,add(1/mus): ',es16.8e3)
 105  format('  G_ab,1ph(1/mus): ',es16.8e3,' G_ab,2ph(1/mus): ',
     &es16.8e3,' G_ab,add(1/mus): ',es16.8e3)
c
      end program
c
      subroutine lab(ii)
c
      implicit double precision (a-h,o-z)
c
      write(ii,*)                                             
      write(ii,*)'___________________________________________________'
      write(ii,*) 
      write(ii,*)'.######.....#######................................'
      write(ii,*)'##    ##....##     ##..###...#...#.................'
      write(ii,*)'##    ##....##     ##..# #...#...#.................'
      write(ii,*)'##    ##....########...###..###..#####..###########'
      write(ii,*)'##  # ##....##     ##...#....#...#   #..#    #    #'
      write(ii,*)'##   ####...##     ##...#....#...#   #..#    #    #'
      write(ii,*)'.######..#..#######.....##...##..#   #..#    #    #'
      write(ii,*)'___________________________________________________'
      write(ii,*)                                             
      write(ii,*)                                             
c
      end subroutine
c
      subroutine rdata(ic,id,ig,ie,nm,nd,fr,rm,sfgw,hawi,azwr,
     &ene,rpso,opp,tpv,tpda,tpst,tpas,nsa,iste,dura,rodi)
c
      implicit double precision (a-h,o-z)
c
      integer*8 :: ic,id,ig,ie,nm,nd,nsa
c
      complex*16 :: opp,tpv,tpda,tpst,tpas
c
      dimension fr(nm),rm(nm),hawi0(nm),hawi(nm),azwr(ic,nd,4)
      dimension ene(ic,nd,id),rpso(ic,nd,3,2),opp(ic,nd,nm)
      dimension tpv(ic,nd,nm,nm),tpda(ic,nd,nm,nm,ie-ig-1)
      dimension tpst(ic,nd,nm,nm,id-ie),tpas(ic,nd,nm,nm,ig-1)
      dimension iste(nsa),dura(nsa),rodi(nsa)
c
      if (nm.gt.0) then
        open(1,file='qb.mdata',status='unknown')
        read(1,*)
        do 1 i=1,nm
          read(1,*) fr(i),rm(i),hawi0(i)
 1      continue
        close(1)
        hawi=sfgw*hawi0
      endif
c
      open(2,file='qb.ddata',status='unknown')
      do 9 m=1,ic
        do 2 i=1,nd
          read(2,*)
          read(2,*)
          read(2,*) (azwr(m,i,j),j=1,4)
          read(2,*)
          read(2,*) (ene(m,i,j),j=1,id)
          read(2,*)
          read(2,*) (rpso(m,i,j,1),j=1,3)
          read(2,*)
          read(2,*) (rpso(m,i,j,2),j=1,3)
          if (nm.gt.0) then
            read(2,*)
            read(2,*) (opp(m,i,j),j=1,nm)
            read(2,*)
            do 3 k=1,nm
              read(2,*)
              read(2,*)
              read(2,*) (tpv(m,i,k,j),j=k,nm)
              do 7 j=k+1,nm
                tpv(m,i,j,k)=tpv(m,i,k,j)
 7            continue
              read(2,*)
              do 4 l=1,ie-ig-1
                read(2,*) (tpda(m,i,k,j,l),j=1,nm)
 4            continue
              read(2,*)
              do 5 l=1,id-ie
                read(2,*) (tpst(m,i,k,j,l),j=1,nm)
 5            continue
              read(2,*)
              do 6 l=1,ig-1
                read(2,*) (tpas(m,i,k,j,l),j=1,nm)
 6            continue
 3          continue
          endif
 2      continue
 9    continue
      close(2) 
c
      open(3,file='qb.adata',status='unknown')
      read(3,*)
      read(3,*)
      do 8 i=1,nsa
        read(3,*) iste(i),dura(i),rodi(i)
 8    continue
      close(3)
c
      end subroutine 
c
      subroutine cirrule(pi,nsa,rodi,nang,wei,ang)
c
      implicit double precision (a-h,o-z)
c
      integer*8 :: nsa,nang
c
      dimension rodi(nsa),rodir(nsa),wei(nang),ang(nsa,nang)
c
      rodir(:)=rodi(:)*(pi/180.0d0)
      do 1 i=1,nang
        wei(i)=1.0d0/dble(nang)
        do 2 j=1,nsa
          ang(j,i)=rodir(j)+(2.0d0*pi*dble(i-1))/dble(nang)
 2      continue
 1    continue
c
      end subroutine
c
      subroutine fopr(hbar,anum,pi,nm,ic,nd,fr,rm,bol,temp,
     &mf,ip,id,ig,ie,opp,hawi,ene,top,gaop,geop)
c
      implicit double precision (a-h,o-z)
c
      integer*8 :: nm,ic,nd,mf,ip,id,ig,ie
c
      complex*16 :: opp
c
      dimension fr(nm),rm(nm),hawi(nm),opp(ic,nd,nm),ene(ic,nd,id)
      dimension taop(nm),teop(nm),itaop(nm),iteop(nm)
c
      gap=ene(mf,ip,ie)-ene(mf,ip,ig)
      gaop=0.0d0;geop=0.0d0;naop=0;neop=0
      cte=1.0d23*hbar*anum*dsqrt(pi/2.0d0)
c
      do 1 i=1,nm
        pnu=1.0d0/(dexp(fr(i)/(bol*temp))-1.0d0)
        sqmo=(cdabs(opp(mf,ip,i)))**2
        delta=(dexp(-0.5d0*(((gap-fr(i))/hawi(i))**2)))/hawi(i)
        taop(i)=(pnu/(fr(i)*rm(i)))*sqmo*delta
        teop(i)=((pnu+1.0d0)/(fr(i)*rm(i)))*sqmo*delta
        gaop=gaop+taop(i);geop=geop+teop(i)
 1    continue
c
      if (dabs(gaop).gt.1.0d-14) then
        taop(:)=(taop(:)/gaop)*100.0d0
        do 2 i=1,nm
          if (dabs(taop(i)).gt.top) then
            naop=naop+1;itaop(naop)=i
          endif
 2      continue
      endif 
      if (dabs(geop).gt.1.0d-14) then
        teop(:)=(teop(:)/geop)*100.0d0
        do 3 i=1,nm
          if (dabs(teop(i)).gt.top) then
            neop=neop+1;iteop(neop)=i
          endif
 3      continue
      endif
c
      if (nm.gt.0) then
        write(7,*) ' 1-ph em         : ',(iteop(i),i=1,neop)
        write(7,*) ' 1-ph ab         : ',(itaop(i),i=1,naop)
      endif
c
      gaop=cte*gaop;geop=cte*geop
c
      end subroutine
c
      subroutine sopr(hbar,anum,pi,cspe,nm,ic,nd,fr,rm,bol,temp,
     &mf,ip,id,ig,ie,tpv,tpda,tpst,tpas,hawi,ene,ttp,gatp,getp)
c
      implicit double precision (a-h,o-z)
c
      integer*8 :: nm,ic,nd,mf,ip,id,ig,ie
c
      complex*16 :: scom,tpv,tpda,tpst,tpas
c
      dimension fr(nm),rm(nm),hawi(nm),ene(ic,nd,id)
      dimension tpv(ic,nd,nm,nm),tpda(ic,nd,nm,nm,ie-ig-1)
      dimension tpst(ic,nd,nm,nm,id-ie),tpas(ic,nd,nm,nm,ig-1)
      dimension tabdi(nm,nm),tstok(nm,nm),tabsp(nm,nm)
      dimension temdi(nm,nm),tanst(nm,nm),temsp(nm,nm)
      dimension itabdi(nm*nm,2),itstok(nm*nm,2)
      dimension itabsp(nm*nm,2),itemdi(nm*nm,2)
      dimension itanst(nm*nm,2),itemsp(nm*nm,2)
c
      gap=ene(mf,ip,ie)-ene(mf,ip,ig)
      gabdi=0.0d0;gstok=0.0d0;gabsp=0.0d0
      gemdi=0.0d0;ganst=0.0d0;gemsp=0.0d0
      nabdi=0;nstok=0;nabsp=0;nemdi=0;nanst=0;nemsp=0
      cte=(1.0d46*(hbar**2)*(anum**2))/(4.0d0*cspe*dsqrt(2.0d0*pi))
c
      do 1 i=1,nm
        do 2 j=1,nm
c
          pnui=1.0d0/(dexp(fr(i)/(bol*temp))-1.0d0)
          pnuj=1.0d0/(dexp(fr(j)/(bol*temp))-1.0d0)
c
          sigma=hawi(i)+hawi(j)
          delta1=dexp(-0.5d0*(((gap-fr(i)-fr(j))/sigma)**2))/sigma
          delta2=dexp(-0.5d0*(((gap-fr(i)+fr(j))/sigma)**2))/sigma
          delta3=dexp(-0.5d0*(((gap+fr(i)-fr(j))/sigma)**2))/sigma
c
          scom=dcmplx(0.0d0,0.0d0)
          do 3 k=1,ie-ig-1
            scom=scom+tpda(mf,ip,i,j,k)/(ene(mf,ip,ig)-ene(mf,ip,k+ig)+
     &fr(i))
 3        continue
          scom=scom+0.5d0*tpv(mf,ip,i,j)
          tabdi(i,j)=((pnui)/(fr(i)*rm(i)))*((pnuj)/(fr(j)*rm(j)))*
     &(cdabs(scom)**2)*delta1
          gabdi=gabdi+tabdi(i,j)
c
          scom=dcmplx(0.0d0,0.0d0)
          do 4 k=1,ie-ig-1
            scom=scom+(dconjg(tpda(mf,ip,j,i,k)))/(ene(mf,ip,ie)-
     &ene(mf,ip,k+ig)-fr(i))
 4        continue
          scom=scom+0.5d0*dconjg(tpv(mf,ip,i,j))
          temdi(i,j)=((pnui+1.0d0)/(fr(i)*rm(i)))*((pnuj+1.0d0)/
     &(fr(j)*rm(j)))*(cdabs(scom)**2)*delta1
          gemdi=gemdi+temdi(i,j)
c
          scom=dcmplx(0.0d0,0.0d0)
          do 5 k=1,id-ie
            scom=scom+tpst(mf,ip,i,j,k)/(ene(mf,ip,ig)-ene(mf,ip,k+ie)+
     &fr(i))
 5        continue
          scom=scom+0.5d0*tpv(mf,ip,i,j)
          tstok(i,j)=((pnui)/(fr(i)*rm(i)))*((pnuj+1.0d0)/
     &(fr(j)*rm(j)))*(cdabs(scom)**2)*delta2
          gstok=gstok+tstok(i,j)
c
          scom=dcmplx(0.0d0,0.0d0)
          do 6 k=1,ig-1
            scom=scom+tpas(mf,ip,i,j,k)/(ene(mf,ip,ig)-ene(mf,ip,k)-
     &fr(i))
 6        continue
          scom=scom+0.5d0*tpv(mf,ip,i,j)
          tabsp(i,j)=((pnui+1.0d0)/(fr(i)*rm(i)))*((pnuj)/
     &(fr(j)*rm(j)))*(cdabs(scom)**2)*delta3
          gabsp=gabsp+tabsp(i,j)
c
          scom=dcmplx(0.0d0,0.0d0)
          do 7 k=1,id-ie
            scom=scom+(dconjg(tpst(mf,ip,j,i,k)))/(ene(mf,ip,ie)-
     &ene(mf,ip,k+ie)+fr(i))
 7        continue
          scom=scom+0.5d0*dconjg(tpv(mf,ip,i,j))
          tanst(i,j)=((pnui)/(fr(i)*rm(i)))*((pnuj+1.0d0)/
     &(fr(j)*rm(j)))*(cdabs(scom)**2)*delta3
          ganst=ganst+tanst(i,j)
c
          scom=dcmplx(0.0d0,0.0d0)
          do 8 k=1,ig-1
            scom=scom+(dconjg(tpas(mf,ip,j,i,k)))/(ene(mf,ip,ie)-
     &ene(mf,ip,k)-fr(i))
 8        continue
          scom=scom+0.5d0*dconjg(tpv(mf,ip,i,j))
          temsp(i,j)=((pnui+1.0d0)/(fr(i)*rm(i)))*((pnuj)/
     &(fr(j)*rm(j)))*(cdabs(scom)**2)*delta2
          gemsp=gemsp+temsp(i,j)
c
 2      continue
 1    continue
c
      if (dabs(gabdi).gt.1.0d-14) then
        tabdi(:,:)=(tabdi(:,:)/gabdi)*100.0d0
        do 9 n=1,nm
          do 10 m=1,nm
            if (dabs(tabdi(n,m)).gt.ttp) then
              nabdi=nabdi+1;itabdi(nabdi,1)=n;itabdi(nabdi,2)=m
            endif
 10       continue
 9      continue
      endif
      if (dabs(gstok).gt.1.0d-14) then
        tstok(:,:)=(tstok(:,:)/gstok)*100.0d0
        do 11 n=1,nm
          do 12 m=1,nm
            if (dabs(tstok(n,m)).gt.ttp) then
              nstok=nstok+1;itstok(nstok,1)=n;itstok(nstok,2)=m
            endif
 12       continue
 11     continue
      endif
      if (dabs(gabsp).gt.1.0d-14) then
        tabsp(:,:)=(tabsp(:,:)/gabsp)*100.0d0
        do 13 n=1,nm
          do 14 m=1,nm
            if (dabs(tabsp(n,m)).gt.ttp) then
              nabsp=nabsp+1;itabsp(nabsp,1)=n;itabsp(nabsp,2)=m
            endif
 14       continue
 13     continue
      endif
      if (dabs(gemdi).gt.1.0d-14) then
        temdi(:,:)=(temdi(:,:)/gemdi)*100.0d0
        do 15 n=1,nm
          do 16 m=1,nm
            if (dabs(temdi(n,m)).gt.ttp) then
              nemdi=nemdi+1;itemdi(nemdi,1)=n;itemdi(nemdi,2)=m
            endif
 16       continue
 15     continue
      endif
      if (dabs(ganst).gt.1.0d-14) then
        tanst(:,:)=(tanst(:,:)/ganst)*100.0d0
        do 17 n=1,nm
          do 18 m=1,nm
            if (dabs(tanst(n,m)).gt.ttp) then
              nanst=nanst+1;itanst(nanst,1)=n;itanst(nanst,2)=m
            endif
 18       continue
 17     continue
      endif
      if (dabs(gemsp).gt.1.0d-14) then
        temsp(:,:)=(temsp(:,:)/gemsp)*100.0d0
        do 19 n=1,nm
          do 20 m=1,nm
            if (dabs(temsp(n,m)).gt.ttp) then
              nemsp=nemsp+1;itemsp(nemsp,1)=n;itemsp(nemsp,2)=m
            endif
 20       continue
 19     continue
      endif
c
      if (nm.gt.0) then
        write(7,*) ' 2-ph emdi   (1º): ',(itemdi(k,2),k=1,nemdi)
        write(7,*) '             (2º)  ',(itemdi(k,1),k=1,nemdi)
        write(7,*) ' 2-ph emanst (1º): ',(itanst(k,2),k=1,nanst)
        write(7,*) '             (2º)  ',(itanst(k,1),k=1,nanst)
        write(7,*) ' 2-ph emsp   (1º): ',(itemsp(k,2),k=1,nemsp)
        write(7,*) '             (2º)  ',(itemsp(k,1),k=1,nemsp)
        write(7,*) ' 2-ph abdi   (1º): ',(itabdi(k,1),k=1,nabdi)
        write(7,*) '             (2º)  ',(itabdi(k,2),k=1,nabdi)
        write(7,*) ' 2-ph abstok (1º): ',(itstok(k,1),k=1,nstok)
        write(7,*) '             (2º)  ',(itstok(k,2),k=1,nstok)
        write(7,*) ' 2-ph absp   (1º): ',(itabsp(k,1),k=1,nabsp)
        write(7,*) '             (2º)  ',(itabsp(k,2),k=1,nabsp)
      endif
c
      gatp=cte*(gabdi+gstok+gabsp);getp=cte*(gemdi+ganst+gemsp)
c
      end subroutine
c
      subroutine rabfr(gfi,ic,nd,mf,ip,azwr,nsa,nang,ir,is,ang,rpso,
     &hbar,bcmt,cmub,alp,pi,prafr,qrafr)
c
      implicit double precision (a-h,o-z)
c
      integer*8 :: ic,nd,mf,ip,nsa,nang,ir,is
c
      dimension azwr(ic,nd,4),rpso(ic,nd,3,2),ang(nsa,nang)
c
      fi=azwr(mf,ip,1);th=azwr(mf,ip,2)
      epsi=ang(is,ir);alpr=alp*(pi/180.0d0)
      x1=-dcos(alpr)*dsin(epsi)*dcos(th)*dcos(fi)-
     &dcos(alpr)*dcos(epsi)*dsin(fi)+dsin(alpr)*dsin(th)*dcos(fi)
      x2=-dcos(alpr)*dsin(epsi)*dcos(th)*dsin(fi)+
     &dcos(alpr)*dcos(epsi)*dcos(fi)+dsin(alpr)*dsin(th)*dsin(fi)
      x3=dcos(alpr)*dsin(epsi)*dsin(th)+dsin(alpr)*dcos(th)
      auct=(cmub*gfi*bcmt)/hbar
      prafr=auct*(x1*rpso(mf,ip,1,1)+x2*rpso(mf,ip,2,1)+
     &x3*rpso(mf,ip,3,1))
      qrafr=auct*(x1*rpso(mf,ip,1,2)+x2*rpso(mf,ip,2,2)+
     &x3*rpso(mf,ip,3,2))
c
      end subroutine
c
      subroutine cmd(rmat,gabs,gem,azwr,prafra,qrafra,
     &pi,cspe,ic,nd,mf,ip,id,ig,ie,ene,firrha)
c
      implicit double precision (a-h,o-z)
c
      integer*8 :: ic,nd,mf,ip,id,ig,ie
c
      dimension azwr(ic,nd,4),ene(ic,nd,id),rmat(4,4)
c
      tmag=azwr(mf,ip,4);gap=ene(mf,ip,ie)-ene(mf,ip,ig)
      a11=-gem-tmag/2.0d0;b11=gabs+tmag/2.0d0
      d=(gabs+gem+2.0d0*tmag)/2.0d0;delta=2.0d0*pi*((gap*cspe)-firrha)
c
      rmat(1,1)=a11;rmat(1,2)=b11;rmat(1,3)=-qrafra;rmat(1,4)=-prafra
      rmat(2,1)=-a11;rmat(2,2)=-b11;rmat(2,3)=qrafra;rmat(2,4)=prafra
      rmat(3,1)=qrafra/2.0d0;rmat(3,2)=-qrafra/2.0d0
      rmat(3,3)=-d;rmat(3,4)=delta
      rmat(4,1)=prafra/2.0d0;rmat(4,2)=-prafra/2.0d0
      rmat(4,3)=-delta;rmat(4,4)=-d
c
      end subroutine
c
      subroutine diag(rmat,cir,sfr,time,info,iflag)
c
      implicit double precision (a-h,o-z)
c
      complex*16 :: cvr,cwork,cof
c
      character*1 :: jobvl,jobvr
c
      dimension rmat(4,4),cir(4),sfr(4),wr(4),wi(4)
      dimension vl(1,4),vr(4,4),work(16),evmsq(4),iau(2)
      dimension ipiv(4),cvr(4,4),cwork(16),cof(4)
c
      jobvl='N';jobvr='V';n=4;lda=4;ldvl=1;ldvr=4;lwork=16;info=0
      call dgeev(jobvl,jobvr,n,rmat,lda,wr,wi,vl,ldvl,vr,ldvr,
     &work,lwork,info)
c
      if (info.ne.0) then
        return
      else
        do 1 i=1,4
          evmsq(i)=wr(i)**2+wi(i)**2
 1      continue
        mpos=minloc(evmsq,dim=1);wr(mpos)=0.0d0;wi(mpos)=0.0d0
        iflag=1
        icon=0
        do 2 i=1,4
          if ((wr(i).lt.0.0d0).and.(dabs(wi(i)).lt.1.0d-14)) then
            icon=icon+1;npos=i
          endif
 2      continue
        if (icon.ne.1) then
          iflag=0
          return
        endif
        wi(npos)=0.0d0
        icon=0
        do 3 i=1,4
          if ((wr(i).lt.0.0d0).and.(dabs(wi(i)).gt.1.0d-14)) then
            icon=icon+1;iau(icon)=i
          endif
 3      continue
        if (icon.ne.2) then
          iflag=0
          return
        endif
c
        cvr(:,1)=vr(:,mpos);cvr(:,2)=vr(:,npos)
        cvr(:,3)=dcmplx(vr(:,iau(1)),vr(:,iau(2)))
        cvr(:,4)=dcmplx(vr(:,iau(1)),-vr(:,iau(2)))
c
        m=4;n=4;lda=4
        call zgetrf(m,n,cvr,lda,ipiv,info)
c         
        if (info.ne.0) then
          return
        else
          n=4;lda=4;lwork=16;info=0
          call zgetri(n,cvr,lda,ipiv,cwork,lwork,info)     
          if (info.ne.0) then
            return
          else
            do 4 i=1,4
              cof(i)=cvr(i,1)*cir(1)+cvr(i,2)*cir(2)+
     &cvr(i,3)*cir(3)+cvr(i,4)*cir(4)
 4          continue
            if ((dabs(aimag(cof(1))).gt.1.0d-14).or.
     &(dabs(aimag(cof(2))).gt.1.0d-14).or.
     &(dabs(aimag(cof(3))+aimag(cof(4))).gt.1.0d-14)) then
              iflag=0
              return
            else
              aure=real(cof(3))*dcos(wi(iau(1))*time)-
     &aimag(cof(3))*dsin(wi(iau(1))*time)
              auim=real(cof(3))*dsin(wi(iau(1))*time)+
     &aimag(cof(3))*dcos(wi(iau(1))*time)
              do 5 i=1,4
                sfr(i)=real(cof(1))*vr(i,mpos)+
     &real(cof(2))*dexp(wr(npos)*time)*vr(i,npos)+
     &2.0d0*dexp(wr(iau(1))*time)*(aure*vr(i,iau(1))-auim*vr(i,iau(2)))
 5            continue
            endif
          endif
        endif
      endif
c
      end subroutine 
c
