!Program to calculate the dispersion relation of zigzag nanoribbon graphene
!and other materials 2D in the presence of Coulomb interaction
!Autor: Jorge H.Correa, Marcos Gaussi , Renan Bento
!Institute of physics - Federal Fluminense University
!Institute of physics - Univesidade de Brasilia
!Institute of physics - Federal Fluminense University
!------------------------------------------------------------------------------
 PROGRAM banda_nanofitas
  
  IMPLICIT NONE
   integer, parameter :: npt=500
   integer :: j,l,z,p,q,i
   integer :: fileras,n,m,nmax
   double precision :: a,pi,c,s1,s2,v,rho,ec,passo,delta,rho1,rho2,rho3
   double precision :: kx,ky,hopping,spin,U,spin_orbit,t,k,e1,cp,sup,k1,k2,k3
	complex*16 :: sr,r1,r2,r3,r12,r22,r32,ps
   parameter(fileras=2,pi=dacos(-1.d0),a=1.42d0,n=8*fileras) !modifique o numero de fileiras para aumentar a nanofita
   parameter(t=1.0d0*0.25d0,nmax=500,e1=0.1d0*0.25d0,sr=0.05*0.25d0,cp=0.05d0*0.5d0,sup=0.25d0*0.5d0) !parametros e1 =zeeman sr= spin rashba cp = potencial quimico sup = supercondutor
!t=hopping e1=zeeman sr=spin rashba cp = Chemical potential sup= super conductor
   double precision :: w(n),mf(n),hdel(n) 
   complex*16 :: H(n,n)
   double precision, parameter :: ecmin=-3.15d0,ecmax=3.15d0
	double complex,dimension(n,n) :: h00,h01,h10
	double complex,dimension(n,n) :: hk 
	double complex :: spx,spy,spz
	integer,parameter :: nktot=500
	double precision,dimension(0:nktot,n,3) :: hh
!=========================================================================
  !PARAMETROS PARA A SUBROTINA DE DIAGONALIZAÇÃO
   integer info,lwork,lda
   parameter(lda=n,lwork=2*n-1)
   complex*16 work(lwork)
   double precision rwork(3*n-2)
   external zheev
!=========================================================================
   OPEN(unit=10,file='energy.dat',status='unknown')
   OPEN(unit=11,file='densidade.dat',status='unknown')
	OPEN(unit=12,file='densidadeup005.dat',status='unknown') !saida spin up
	OPEN(unit=13,file='densidadedown005.dat',status='unknown') !saida spin down
	OPEN(unit=14,file='densidadeupdown005.dat',status='unknown') !saida de spin não definido
	OPEN(unit=15,file='densidadetotal005.dat',status='unknown')
	OPEN(unit=16,file='spin.dat',status='unknown')
  passo = (ecmax-ecmin)/float(npt)

 
do p=1,npt+1
   ec = ecmin +passo*float(p-1)
		rho3=0.d0
		rho2=0.d0
		rho1=0.d0
      rho =0.d0
  DO 50 m=-nmax,nmax!0,nmax
   
     kx= DFLOAT(m)*pi/nmax
     C=hopping*2.d0*cos(kx/2.d0)!hopping_strain é com strain,eu só coloco hopping
   !========================================================================
     !HERE IS FORMED THE MATRIX
    DO 100 i=1,n
    DO 200 j=1,n
    H(i,j)=0.0d0
    w(i)=0.0d0
   200 CONTINUE
  100 CONTINUE
		!termos da matriz
				ps=-sup*0.5*Dsin(kx)*(0.0,1.0)
				r1=-2*sr*DCOS(kx*0.5 - 2.d0*PI/3.d0)*(0.0,1.0)
				r2=2*sr*DCOS(kx*0.5 + 2.d0*PI/3.d0)*(0.0,1.0)
				r3=-sr*(0.0,1.0)
				r12=2*sr*DCOS(kx*0.5 + 2.d0*PI/3.d0)*(0.0,1.0)
				r22=-2*sr*DCOS(kx*0.5 - 2.d0*PI/3.d0)*(0.0,1.0)
				r32=SR*(0.0,1.0)
		!Aqui vou por uma sub rotina para montar a matriz 
		Call TERMOS_H(H,t,kx,n,e1,r1,r2,r3,r12,r22,r32,sr,ps,cp) 
!write(*,2) aimag(H)
!2 format(16f5.1)		
!======================================================================== 

!=====================================================================================
!==========CALLING THE ROUTINE TO DIAGONALIZE=========================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%						
		CALL zheev('V','U',n,H,n,W,work,lwork,rwork,info)!eigv
		IF(info/=0)THEN
			WRITE(*,*)"ERRO NO HAMILTONIANO"
			STOP
		ENDIF
!			write(*,3) aimag(H) 
!			3 format(16f5.2)
!stop
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		!aqui precisamos separar por spin 
		DO l=1,n
			WRITE(10,201) kx, W(l)/t!eigv(l)/hopping
		END DO 

		201 FORMAT(f11.8,4x,f11.8)
		j=0
!=============================================================================
		do j=1,n		 
		hh(i,j,2)=0
		hh(i,j,3)=0
			call spinvl2(n,H(:,j),spz)

			hh(i,j,1) = m!*(1.d0/(2.*pi))
			hh(i,j,2) = W(j)
			hh(i,j,3) = spz
			if(hh(i,j,3).GT. 1.d0) hh(i,j,3)= 0.99999d0  
			if(hh(i,j,3).LT. -1.d0) hh(i,j,3)= -0.99999d0  
			write(16,*)  hh(i,j,2) , hh(i,j,3)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
!---SAÍDA DOS DADOS PARA ARQUIVO	
!----------vou por uma condicao para spin up e spin down	
			if	(hh(i,j,3) > 0.5d0 ) then
				k1= ec- hh(i,j,2)
				rho1 = rho1 + (1.d0/(pi*nmax))*delta(k1)

			else if (hh(i,j,3) < -0.5d0 ) then
				k2= ec- hh(i,j,2)
				rho2 = rho2 + (1.d0/(pi*nmax))*delta(k2)
			else 
				k3= ec- hh(i,j,2)
				rho3 = rho3 + (1.d0/(pi*nmax))*delta(k3)
			endif
		 end do

            do l=1,n
            k= ec-w(l)           
            rho = rho + (1.d0/(pi*nmax))*delta(k)
            enddo
 50 continue
			
             write(11,*) ec,rho
				 write(12,*) ec,rho1
				 write(13,*) ec,rho2
				 write(14,*) ec, rho3
				 write(15,*) ec, rho1+rho2+rho3 
				 write(*,1) real(ec)
				 1 format(4f5.1)	
enddo
END PROGRAM banda_nanofitas


double precision function delta(n)
 
 implicit none

 double precision, parameter :: eta=1.d-2
 double precision :: n

 delta = (1.d0/eta)*exp(-n*n/(eta*eta))
end function



SUBROUTINE TERMOS_H(H,t,kx,n,e1,r1,r2,r3,r12,r22,r32,sr,ps,cp) 
integer :: i,j
double precision, INTENT(IN) :: t,kx,e1,cp
complex*16, INTENT(IN) :: sr,r1,r2,r3,r12,r22,r32,ps
COMPLEX(8),INTENT(OUT) :: H(n,n)


!######################POTENCIAL QUIMICO E EFEITO ZEEMAN 

	do i= 1 , N/4 ,2
		do j= 1 , N/4
			if (i==j) then
				H(i,j) =-cp +e1 
				H(i+N/4,j+N/4) = -cp - e1  
				H(i+N/2,j+N/2) = cp -e1 
				H(i+(3*N)/4,j+3*N/4) = cp +e1  
			end if
		end do
	end do	

	do i= 2 , N/4 ,2
		do j= 1 , N/4
			if (i==j) then
				H(i,j) =-cp +e1 
				H(i+N/4,j+N/4) = -cp - e1 
				H(i+N/2,j+N/2) = cp -e1 
				H(i+(3*N)/4,j+3*N/4) = cp +e1 
			end if
		end do
	end do	
!###################### PRIMEIROS VIZINHOS

	do i= 1 , N/4, 2 
		do j= 1 , N/4 
			if (i+1==j) then
				H(i,j) = 2.d0*t*cos(kx/2.d0)
				H(j,i) = 2.d0*t*cos(kx/2.d0)
				H(i+N/4,j+N/4) = 2.d0*t*cos(kx/2.d0)
				H(j+N/4,i+N/4) = 2.d0*t*cos(kx/2.d0)
      		H(i+N/2,j+N/2) =-2.d0*t*cos(kx/2.d0)
				H(j+N/2,i+N/2) = -2.d0*t*cos(kx/2.d0)
				H(i+(3*N)/4,j+3*N/4) =-2.d0*t*cos(kx/2.d0)
				H(j+(3*N)/4,i+3*N/4) = -2.d0*t*cos(kx/2.d0)
			end if
		end do
	end do
	
	do i= 2 , N/4, 2 
		do j= 1 , N/4 
			if (i+1==j) then
				H(i,j) = -t
				H(j,i) = -t
				H(i+N/4,j+N/4) = -t
				H(j+N/4,i+N/4) = -t
  			   H(i+N/2,j+N/2) = t
				H(j+N/2,i+N/2) = t
				H(i+(3*N)/4,j+3*N/4) = t
				H(j+(3*N)/4,i+3*N/4) = t
			end if
		end do
	end do
!##############  RASHBA INTERACTION #####################################

	do i= 1+N/4 , N/2 
		do j= 1 , N/4 , 2
			if (i-1==j+N/4) then
				H(i,j) = r1
				H(j,i) = r1
				H(i+N/2,j+N/2) = -r12
				H(j+N/2,i+N/2) = -r12
			end if
		end do
	end do

	do i= 1+N/4 , N/2 ,2
	do j= 1 , N/4 
		if (i+1==j+N/4) then
		H(i,j) = r2
		H(j,i) = r2
		H(i+N/2,j+N/2) = -r22
		H(j+N/2,i+N/2) = -r22
	end if
	end do
	end do


	do i= 1+N/4 , N/2 
	do j= 2 , N/4 , 2
		if (i-1==j+N/4) then
		H(i,j) = r3
		H(j,i) = r3
		H(i+N/2,j+N/2) = -r32
		H(j+N/2,i+N/2) = -r32
	end if
	end do
	end do


	do i= 2+N/4 , N/2 , 2
	do j= 1 , N/4 
		if (i+1==j+N/4) then
		H(i,j) = -r3
		H(j,i) = -r3
		H(i+N/2,j+N/2) = r32
		H(j+N/2,i+N/2) = r32
	end if
	end do
	end do

!######################POTENCIAL SUPER CONDUTOR###############################    
 	 do i= 1+N/2 , N 
	 do j= 1   , N/4 ,2
		if (i==j+N/2) then
		H(i,j) = ps
		H(j,i) = -ps
		H(i+N/4,j+N/4) = ps
		H(j+N/4,i+N/4) = -ps
	end if
	end do
	end do

 	 do i= 1+N/2 , N 
	 do j= 2   , N/4 ,2
		if (i==j+N/2) then
		H(i,j) = ps
		H(j,i) = -ps
		H(i+N/4,j+N/4) = ps
		H(j+N/4,i+N/4) = -ps
	end if
	end do
	end do

END SUBROUTINE


!###################################################################################
subroutine prodintsq2(vetora,vetorb,n,resultado)

	implicit none

	integer:: n,i
	double complex :: resultado1, auxi
	double complex, dimension(n) :: vetora,vetorb,vetorc
	double complex :: resultado

	

	resultado1=0

	do i=1,n

		auxi=vetora(i)*vetorb(i)

		resultado1=resultado1+auxi

	end do

	resultado=resultado1

end subroutine prodintsq2

subroutine vecconjg2(vetor,n,vetout)

	implicit none

	integer :: i,j,n
	double complex, dimension(n) :: vetor,vetout


	do i=1,n

		vetout(i)=conjg(vetor(i))

	end do



end subroutine vecconjg2

subroutine matvec2(matriz,vetor,N,veout)

	implicit none

	integer :: N !ordem da matriz e tamanho do vetor
	double complex, dimension(N,N) :: matriz
	double complex, dimension(N):: vetor,veout
	double complex :: flag
	integer :: i,j

	veout=0
	
	do i=1,N


		do j=1,N

			flag=matriz(i,j)*vetor(j)

			veout(i)=flag+veout(i)
		

		end do
	

	end do


end subroutine matvec2

subroutine spinvl2(N,vec,spz)!,spx,spy)


	implicit none

        integer :: N

	double complex,dimension(N) :: vec
	double complex,dimension(N) :: vcconj

	double complex,dimension(N) :: xvvaux,zvvaux,yvvaux

	double complex,dimension(N,N) :: spmx,spmy,spmz

	double complex :: spx,spy,spz
     
        double precision :: spz1

	integer :: i,j


	!spmx=0.0
	!spmy=0.0
	spmz=0.0

	do i=1,N/4

		spmz(i,i)= cmplx(1.,0.)
		
		spmz(i+n/4,i+n/4)= cmplx(-1.,0.)
		spmz(i+n/2,i+n/2)= cmplx(1.,0.)
		spmz(i+n/2+n/4,i+n/2+n/4)= cmplx(-1.,0.)
		!spmx(i,n/2+i)= cmplx(1.,0.)

		!spmx(i+n/2,i)= cmplx(1.,0.)


		!spmy(i,n/2+i)= cmplx(0.,-1.)

		!spmy(i+n/2,i)= cmplx(0.,1.)

	end do

	call vecconjg2(vec,N,vcconj)

	!calculando valor medio de sz

	call matvec2(spmz,vec,N,zvvaux)

	call prodintsq2(vcconj,zvvaux,N,spz)


	!calculando valor medio de sy

	!call matvec(spmy,vec,n,yvvaux)

	!call prodintsq(vcconj,yvvaux,n,spy)

	!calculando valor medio de sx

	!call matvec(spmx,vec,n,xvvaux)

	!call prodintsq(vcconj,xvvaux,n,spx)
       	spz1 = dble(spz)



end subroutine spinvl2
