!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!VARIÁVEIS GLOBAIS DO PROGRAMA
module parametros
	INTEGER,PARAMETER :: fileras= 2 !aumenta o tamanho da fita

	INTEGER,PARAMETER :: N=(8*fileras)		!NÚMERO DE ORDEM DA MATRIX

	COMPLEX(8) :: C,S1,S2,R1,R2,R3,V,C1,S12,S22,R12,R22,R32,M!,RI1,RI2,RI3		!TERMO DO COSSENO QUE VAI NA MATRIZ   
!        REAL(8) :: S,F     
!        REAL(8):: V
!	REAL(8) :: hopping
	REAL(8), PARAMETER :: d=1.42d0		!DIST DOIS ÁTOMOS NA REDE
	REAL(8), PARAMETER ::PI=dacos(-1.d0)!PI=3.14159265358979323 !
end module parametros
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


program main
USE parametros
implicit none
external ZHEEV
complex*16, parameter :: Re=(1.0,0.0)
complex*16 :: egvc(N,N) ,H3(N,N), H4(N,N),aux1(N,N)
complex*16 :: egvl(N),  H2(N) !,H3(N)
complex*16 :: Dmat(N,N) , H(N,N) , U(N,N) , ISO(N,N) 
!variaveis Lapack ***********************
character(len=1) :: JOBZ,UPLO
integer :: LDA,LWORK,INFO!,N
real*8, allocatable :: W(:), RWORK(:)
complex*16, allocatable :: A(:,:), WORK(:)
!*********************************************

integer :: i ,j, MT, k , l,MT2
complex*16 :: T1, T2,T3,AT, BT, CT,DT,DT1,SR,CHE, POT,E1,SO,CHER
REAL*8 :: EH, OI, OF, H1, AW, CC, ETA, NT , AH, AW1!, SO
REAL*8 :: AX(4), AY(4), TX(4) , TY(4)   
REAL*8 :: RX , RY

integer, parameter :: npt=1000
double precision :: VN,delta,step,e,ec,passo,rho
double precision, parameter :: emin=-pi,emax=pi
double precision, parameter :: ecmin=-3.d0,ecmax=3.d0
complex :: omega

double complex,dimension(N,N) :: h00,h01,h10
double complex,dimension(N,N) :: hk 
double complex :: spx,spy,spz
integer,parameter :: nktot=1000
double precision,dimension(0:nktot,n,3) :: hh
OPEN(UNIT=50,FILE='banda_teste2.dat',status='unknown')
OPEN(UNIT=60,FILE='banda_spin083.dat',status='unknown') !aquivo para o gnuplot  primeira coluna kx segunda energia terceira spin
!open(unit=10,file='densidade2.dat',status='unknown')
!******************************************************************* 
	OPEN(unit=12,file='bandupNOVO.dat',status='unknown') !saida de spin up
	OPEN(unit=13,file='banddownNOVO.dat',status='unknown') !saida de spin down
	OPEN(unit=14,file='bandspinNOVO.dat',status='unknown') !saida de spin não definido
!*******************************************************************
JOBZ="V"
UPLO="U"
INFO=0
CHER=0.83d0 !potencial quimico 1.33/1.23/0.84/0.73/0.55/0.45/0.05 < onde tem transição
LDA=N
LWORK=N*2-1
allocate(A(N,N),W(N),WORK(LWORK),RWORK(N*3-2))
A=0.0*Re
W=0.0
WORK=0.0*Re
RWORK=0.0
!******************************************************************
	MT2=0
	H4=0.0d0
	ETA=0.005
!	PI=3.14159265358979323
	AH=3.86
	EH=0.D0
        MT=1001
        OI=-PI
        OF=PI
        H1=(OF-OI)/(MT-1)
	step = (emax-emin)/float(npt)
	DO k=1,MT
        AW=(OI+(k-1)*H1)!PI/2
        
!	AW=0
	print*, AW
!**************Elementos da MATRIZ************************


	T1=1.00d0*(1.0,0.0) !primeiros vizinhos
	T3=0.0d0*0.25d0*(1.0,0.0) !terceiros vizinhos 
	SO=0.00*1*(1.0,0.0) !spin orbita intrincico 
	E1=0.1!0.1*(1.0,0.0) !effeito zeeman 
	SR=0.05!0.0*(0.0,1.0) !spin rashba
	POT=0.5*2*Dsin(AW)*(0.0,1.0)!potencial super condutor
	CHE=-CHER*(2,0.0)  !potencial quimico 1.33/1.23/0.84/0.73/0.55/0.45/0.05 < onde tem transição
	
	
	
	C=-2*T1*DCOS(0.5D0*AW)
	DT=-2*T3*DCOS(AW)*(1.0,0.0)
	S2=-2*SO*DSIN(AW)
	S1=2*SO*DSIN(0.5*AW)
	R1=-2*SR*DCOS(AW*0.5 - 2.d0*PI/3.d0)*(0.0,1.0)
	R2=2*SR*DCOS(AW*0.5 + 2.d0*PI/3.d0)*(0.0,1.0)
	R3=-SR*(0.0,1.0)
	

	C1=-2*T1*DCOS(0.5D0*AW)
	DT1=-2*T3*DCOS(AW)*(1.0,0.0)
	S22=-2*SO*DSIN(AW)
	S12=2*SO*DSIN(0.5*AW)
	R12=2*SR*DCOS(AW*0.5 + 2.d0*PI/3.d0)*(0.0,1.0)
	R22=-2*SR*DCOS(AW*0.5 - 2.d0*PI/3.d0)*(0.0,1.0)
	R32=SR*(0.0,1.0)

	
	ISO=0.d0
	H=0.d0
	Dmat=0.d0 
      
    
Call TERMOS_M1(Dmat,AW,POT,CHE,T1, T2,T3,SR,E1,SO,DT)
!	write(*,1) real(Dmat)
!	1 format(16f10.2)
	!	write(*,*)
!Stop
!***********************PREPARANDO************************************************
A=Dmat !Matriz para ser diagonalizada
LWORK=-1
call zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
if (INFO.eq.0) then
 LWORK=int(real(WORK(1)))
 deallocate(WORK)
 allocate(WORK(LWORK))
 WORK=0.0*Re
else if (INFO.ne.0) then
 write(*,*) "From zheev: query for optimal workspace failed: STOP!"
 deallocate(A,W,WORK,RWORK)
 stop
endif

!DIAGONALIZAÇÃO***************************************
A=Dmat

call zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
!Check
if (INFO.ne.0) then
 write(*,*) "Diagonalization failed INFO=",INFO," STOP!"
 deallocate(A,W,WORK,RWORK)
 stop
endif
egvl=W !autovalores da matriz diagonalizada
egvc=A !autovetores da matriz diagonalizada

write(55,*) AW, ( real(egvl(i)), i=1,N) !arquivo de saida fort.55




    		do j=1,n		 

		call spinvl(N,A(:,j),spz)

			!write(50,*) k*(a/(2.*pi)),W(j),spz

			hh(i,j,1) = k!*(1.d0/(2.*pi))
			hh(i,j,2) = W(j)
			hh(i,j,3) = spz
			if(hh(i,j,3).GT. 1.d0) hh(i,j,3)= 0.99999d0  
			if(hh(i,j,3).LT. -1.d0) hh(i,j,3)= -0.99999d0  
			write(50,*) AW,hh(i,j,2),hh(i,j,3) !saida do arquivo com textura de spin
			write(50,*)
		end do

		do j=1,n		 
		hh(i,j,2)=0
		hh(i,j,3)=0
		call spinvl2(N,A(:,j),spz)		 

			!write(50,*) k*(a/(2.*pi)),W(j),spz

			hh(i,j,1) = k!*(1.d0/(2.*pi))
			hh(i,j,2) = W(j)
			hh(i,j,3) = spz
			if(hh(i,j,3).GT. 1.d0) hh(i,j,3)= 0.99999d0  
			if(hh(i,j,3).LT. -1.d0) hh(i,j,3)= -0.99999d0 
			write(60,*) AW ,hh(i,j,2),hh(i,j,3)
			write(60,*)
			if	(hh(i,j,3) > 0.9d0 ) then
				write(12,*) AW,hh(i,j,2)

			else if (hh(i,j,3) < -0.9d0 ) then
				write(13,*) AW,hh(i,j,2)
			else 
				write(14,*) AW, hh(i,j,2)
			endif
			
			
		end do
		!do MT2=1,npt+1
! 		ec = AW
!   		rho =0.d0
!   		do j=1,N
!   		e = egvl(j)
 		!print*,(e)
!    		VN= ec-e
!    		rho = rho + (1.d0/(pi*AW))*delta(VN)
!     		 enddo
! 		write(10,*) ec,rho
		!enddo
do i=1,N
 write(22,*) AW,  real(egvl(i))
 write(66,*) real(egvc(i,:))!saida auto vetores parte real
 write(77,*) aimag(egvc(i,:))!saida auto vetores parte imaginaria


enddo

        ENDDO
deallocate(A,W,WORK,RWORK)
stop
end program main
!******************************************************************
SUBROUTINE TERMOS_M1(Dmat,AW,POT,CHE,T1, T2,T3,SR,E1,SO,DT)


	USE PARAMETROS	
	IMPLICIT NONE
    integer :: i ,j, MT,k
   REAL(8), INTENT(IN)   :: AW
   COMPLEX(8), INTENT(IN) :: POT, CHE, T1, T2,T3,SR,E1,SO,DT
    COMPLEX(8),INTENT(OUT) :: Dmat(N,N)


!######################POTENCIAL QUIMICO E EFEITO ZEEMAN E SPIN ORBITA INTRINSICO

	do i= 1 , N/4 ,2
	do j= 1 , N/4
	!print*, i, j
		if (i==j) then
		Dmat(i,j) =-CHE +E1 + S2
		Dmat(i+N/4,j+N/4) = -CHE - E1  -S2
		Dmat(i+N/2,j+N/2) = CHE -E1 - S22
		Dmat(i+(3*N)/4,j+3*N/4) = CHE +E1  + S22
	end if
	end do
	end do	

	do i= 2 , N/4 ,2
	do j= 1 , N/4
	!print*, i, j
		if (i==j) then
		Dmat(i,j) =-CHE +E1 - S2
		Dmat(i+N/4,j+N/4) = -CHE - E1 +S2
		Dmat(i+N/2,j+N/2) = CHE -E1 + S22
		Dmat(i+(3*N)/4,j+3*N/4) = CHE +E1 - S22
	end if
	end do
	end do	
!###################### PRIMEIROS VIZINHOS

	do i= 1 , N/4, 2 
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j) then
		Dmat(i,j) = -C
		Dmat(j,i) = -C
		Dmat(i+N/4,j+N/4) = -C
		Dmat(j+N/4,i+N/4) = -C
                Dmat(i+N/2,j+N/2) = C1
		Dmat(j+N/2,i+N/2) = C1
		Dmat(i+(3*N)/4,j+3*N/4) = C1
		Dmat(j+(3*N)/4,i+3*N/4) = C1

	end if
	end do
	end do
	

	do i= 2 , N/4, 2 
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j) then
		Dmat(i,j) = -T1 + DT
		Dmat(j,i) = -T1 + DT
		Dmat(i+N/4,j+N/4) = -T1 + DT
		Dmat(j+N/4,i+N/4) = -T1 + DT
                Dmat(i+N/2,j+N/2) = T1 - DT
		Dmat(j+N/2,i+N/2) = T1 - DT
		Dmat(i+(3*N)/4,j+3*N/4) = T1- DT
		Dmat(j+(3*N)/4,i+3*N/4) = T1 - DT

	end if
	end do
	end do

!#######################Terceiros vizinhos##############################

	do i= 1 , N/4, 2 
	do j= 1 , N/4 
	!print*, i, j
		if (i+3==j) then
		Dmat(i,j) = -T3
		Dmat(j,i) = -T3
		Dmat(i+N/4,j+N/4) = -T3
		Dmat(j+N/4,i+N/4) = -T3
                Dmat(i+N/2,j+N/2) = T3
		Dmat(j+N/2,i+N/2) = T3
		Dmat(i+(3*N)/4,j+3*N/4) = T3
		Dmat(j+(3*N)/4,i+3*N/4) = T3

	end if
	end do
	end do



!######################SPIN ORBITA INTRINSICO###############################  

 	do i= 1 , N/4 ,2  
	do j= 1 , N/4 
	!print*, i, j
		if (i+2==j) then
		Dmat(i,j) = S1
		Dmat(j,i) = S1
		Dmat(i+N/4,j+N/4) = -S1
		Dmat(j+N/4,i+N/4) = -S1
                Dmat(i+N/2,j+N/2) = -S12
		Dmat(j+N/2,i+N/2) = -S12
		Dmat(i+(3*N)/4,j+3*N/4) = S12
		Dmat(j+(3*N)/4,i+3*N/4) = S12
	end if
	end do
	end do

	do i= 2 , N/4 ,2  
	do j= 1 , N/4 
	!print*, i, j
		if (i+2==j) then
		Dmat(i,j) = -S1
		Dmat(j,i) = -S1
		Dmat(i+N/4,j+N/4) = S1
		Dmat(j+N/4,i+N/4) = S1
                Dmat(i+N/2,j+N/2) = S12
		Dmat(j+N/2,i+N/2) =S12
		Dmat(i+(3*N)/4,j+3*N/4) = -S12
		Dmat(j+(3*N)/4,i+3*N/4) = -S12
	end if
	end do
	end do

!##############  RASHBA INTERACTION #####################################

	do i= 1+N/4 , N/2 
	do j= 1 , N/4 , 2
	!print*, i, j
		if (i-1==j+N/4) then
		Dmat(i,j) = R1
		Dmat(j,i) = R1
		Dmat(i+N/2,j+N/2) = -R12
		Dmat(j+N/2,i+N/2) = -R12
	end if
	end do
	end do

	do i= 1+N/4 , N/2 ,2
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j+N/4) then
		Dmat(i,j) = R2
		Dmat(j,i) = R2
		Dmat(i+N/2,j+N/2) = -R22
		Dmat(j+N/2,i+N/2) = -R22
	end if
	end do
	end do


	do i= 1+N/4 , N/2 
	do j= 2 , N/4 , 2
	!print*, i, j
		if (i-1==j+N/4) then
		Dmat(i,j) = R3
		Dmat(j,i) = R3
		Dmat(i+N/2,j+N/2) = -R32
		Dmat(j+N/2,i+N/2) = -R32
	end if
	end do
	end do


	do i= 2+N/4 , N/2 , 2
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j+N/4) then
		Dmat(i,j) = -R3
		Dmat(j,i) = -R3
		Dmat(i+N/2,j+N/2) = R32
		Dmat(j+N/2,i+N/2) = R32
	end if
	end do
	end do



	





  

!######################POTENCIAL SUPER CONDUTOR###############################    
 	 do i= 1+N/2 , N 
	 do j= 1   , N/4 ,2
	!print*, i, j
		if (i==j+N/2) then
		Dmat(i,j) = POT
		Dmat(j,i) = -POT
		Dmat(i+N/4,j+N/4) = POT
		Dmat(j+N/4,i+N/4) = -POT
	end if
	end do
	end do

 	 do i= 1+N/2 , N 
	 do j= 2   , N/4 ,2
	!print*, i, j
		if (i==j+N/2) then
		Dmat(i,j) = POT
		Dmat(j,i) = -POT
		Dmat(i+N/4,j+N/4) = POT
		Dmat(j+N/4,i+N/4) = -POT
	end if
	end do
	end do

END SUBROUTINE



subroutine prodintsq(vetora,vetorb,n,resultado)

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

end subroutine prodintsq

subroutine vecconjg(vetor,n,vetout)

	implicit none

	integer :: i,j,n
	double complex, dimension(n) :: vetor,vetout


	do i=1,n

		vetout(i)=conjg(vetor(i))

	end do



end subroutine vecconjg 

subroutine matvec(matriz,vetor,N,veout)

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


end subroutine matvec

subroutine spinvl(N,vec,spz)!,spx,spy)


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

	do i=1,N/2

		spmz(i,i)= cmplx(1.,0.)

		spmz(i+n/2,i+n/2)= cmplx(-1.,0.)


		!spmx(i,n/2+i)= cmplx(1.,0.)

		!spmx(i+n/2,i)= cmplx(1.,0.)


		!spmy(i,n/2+i)= cmplx(0.,-1.)

		!spmy(i+n/2,i)= cmplx(0.,1.)

	end do

	call vecconjg(vec,N,vcconj)

	!calculando valor medio de sz

	call matvec(spmz,vec,N,zvvaux)

	call prodintsq(vcconj,zvvaux,N,spz)


	!calculando valor medio de sy

	!call matvec(spmy,vec,n,yvvaux)

	!call prodintsq(vcconj,yvvaux,n,spy)

	!calculando valor medio de sx

	!call matvec(spmx,vec,n,xvvaux)

	!call prodintsq(vcconj,xvvaux,n,spx)
       	spz1 = dble(spz)



end subroutine spinvl
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

double precision function delta(VN)
 
 implicit none

 double precision, parameter :: eta=1.d-2
 double precision :: VN

 delta = (1.d0/eta)*exp(-VN*VN/(eta*eta))
end function

