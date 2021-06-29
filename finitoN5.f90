!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!VARIÁVEIS GLOBAIS DO PROGRAMA
module parametros
INTEGER,PARAMETER :: fileras=1 !aumenta o tamanho da fita

INTEGER,PARAMETER :: N=(16*fileras) !NÚMERO DE ORDEM DA MATRIX

COMPLEX(8) :: C,S1,S2,R1,R2,R3,V,C1,S12,S22,R12,R22,R32,M!,H1,H2,M,RI1,RI2,RI3		!TERMO DO COSSENO QUE VAI NA MATRIZ   
!        REAL(8) :: S,F     
!        REAL(8):: V
!	REAL(8) :: hopping
!	REAL(8), PARAMETER :: a=3.86d0		!DIST DOIS ÁTOMOS NA REDE
REAL(8), PARAMETER ::PI=dacos(-1.d0)!PI=3.14159265358979323 !
end module parametros
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


program main
USE parametros
implicit none
external ZHEEV
complex*16, parameter :: Re=(1.0,0.0)
complex*16 :: egvc(N,N)
complex*16 :: egvl(N)
complex*16 :: Dmat(N,N) , H(N,N) , U(N,N) , ISO(N,N) , H3(N,N), H4(N,N) 
!variaveis Lapack ***********************
character(len=1) :: JOBZ,UPLO
integer :: LDA,LWORK,INFO!,N
real*8, allocatable :: W(:), RWORK(:)
complex*16, allocatable :: A(:,:), WORK(:)
!*********************************************
integer :: i ,j, MT,k,l
complex*16 :: T1, T2,T3,AT, BT, CT,DT,DT1,SR,CHE, POT,E1,SO
REAL*8 :: EH, OI, OF, H1, AW, CC, ETA, NT , AH, AW1!, SO
INTEGER DATE_TIME (8)

CHARACTER (LEN = 12) REAL_CLOCK (3)
double complex,dimension(N,N) :: h00,h01,h10
double complex,dimension(N,N) :: hk 
double complex :: spx,spy,spz
integer,parameter :: nktot=501!1601
double precision,dimension(0:nktot,N,3) :: hh
!OPEN(UNIT=50,FILE='banda_teste2_FINITOVz0005.dat',status='unknown')
OPEN(UNIT=60,FILE='MINIMO01005004.dat',status='unknown') !saida gnuplot
!******************************************************************* 
JOBZ="V"
UPLO="U"
INFO=0

LDA=N
LWORK=N*2-1
allocate(A(N,N),W(N),WORK(LWORK),RWORK(N*3-2))
A=0.0*Re
W=0.0
WORK=0.0*Re
RWORK=0.0
!******************************************************************


ETA=0.005
!	PI=3.14159265358979323
AH=3.86
EH=0.D0
M=0.05
MT=251
OI=-2
OF=2
H1=(OF-OI)/(MT-1)
DO k=1,MT
AW=(OI+(k-1)*H1)!PI/2
        
!	AW=0
!	print*, AW
!**************Elementos da MATRIZ************************


!	CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &

!	                    REAL_CLOCK (3), DATE_TIME)
!	print*, DATE_TIME(6), DATE_TIME(7)

T1=(1.0,0.0) !primeiros vizinhos
T3=0.0*0.5d0*(1.0,0.0) !terceiros vizinhos 
SO=0.00*(0.0,1.0) !spin orbita intrincico NÃO TIRAR O NEGATIVO
E1=0.05!0.1*(1.0,0.0) !effeito zeeman
SR=0.04!0.0*(0.0,1.0) !spin rashba
POT=0.1*(1.0,0.0)!*DSIN(AW) !potencial super condutor
CHE=2.0*AW!*(1.0,0.0)  !potencial quimico
	
	
	
 C=T1 !primeiros vizinhos
DT=T3!-2*T3*DCOS(AW)*(1.0,0.0) não serve para nada
S2=0!SO!-2*SO*DSIN(AW)*(1.0,0.0) não serve para nada
S1=0!!SO!2*SO*DSIN(0.5*AW)*(1.0,0.0) não serve para nada
R1=SR*0.5D0*((0.0,1.0)*0.5D0 + 0.5D0*SQRT(3.D0)*(1.0,0.0))!2*SR*DCOS(AW*0.5 + 2.d0*PI/3.d0)
R2=SR*0.5D0*((0.0,1.0)*0.5D0 - 0.5D0*SQRT(3.D0)*(1.0,0.0))!-2*SR*DCOS(AW*0.5 - 2.d0*PI/3.d0)
R3=0.5*SR*(0.0,1.0)
	
!	C1=-2*T1*DCOS(0.5D0*AW)!T1 !primeiros vizinhos
	!DT=-2*T3*DCOS(AW)*(1.0,0.0)
S22=0!-2*SO*DSIN(AW)
S12=0!2*SO*DSIN(0.5*AW)
R12=SR*0.5D0*((0.0,-1.0)*0.5D0 + 0.5D0*SQRT(3.D0)*(1.0,0.0))!-2*SR*DCOS(AW*0.5 + 2.d0*PI/3.d0)*(0.0,1.0)
R22=SR*0.5D0*((0.0,-1.0)*0.5D0 - 0.5D0*SQRT(3.D0)*(1.0,0.0))!2*SR*DCOS(AW*0.5 - 2.d0*PI/3.d0)*(0.0,1.0)
R32=-0.5*SR*(0.0,1.0)
	
ISO=0.d0
H=0.d0
Dmat=0.d0
      
    
!Call TERMOS_M1(Dmat,AW,POT,CHE,T1, T2,T3,SR,E1,SO)
!Call TERMOS_M2(Dmat,AW,POT,CHE,T1, T2,T3,SR,E1,SO)
Call  TERMOS_M3(Dmat,AW,POT,CHE,T1, T2,T3,SR,E1,SO) !chamando matriz
!Call TERMOS_M4(Dmat,AW,POT,CHE,T1, T2,T3,SR,E1,SO)
!**********************T2******************************
Print*, AW
!	Print*, "  "
!	write(*,1) int(Dmat)
!	1 format(32i4)
!	Print*, "  "


!	Print*, "  "
!	print*, aimag(Dmat)
!	Print*, "  "
!	print*, int(Dmat)
!    Print*, "  "
!	write(*,1) int(Dmat)
!	1 format(8i4)
!	Print*, "  "
!	Print*, "  "
!	write(*,2) aimag(Dmat)
!	2 format(32f5.2)
!	Print*, "  "
!	Print*, "  "
!	write(*,3) REAL(Dmat)
!	3 format(32f5.2)
!	Print*, "  "

!	print*, real(Dmat)
!	Print*, "  "

!	Print*, "  "
!	do j = N/2+1, N
!	print*, (int(Dmat(i,j)), i=N/2+1,N)
!	end do
!	Print*, "  "


!	Print*, "  "
!	do j = 1, N/4
!	print*, (aimag(Dmat(i,j)), i=1,N/4)
!	end do
!	Print*, "  "
!	do j = 1, N/4
!	print*, (real(Dmat(i,j)), i=1,N/4)
!	end do
 !      Print*, "  "
!***********************PREPARANDO************************************************
A=Dmat
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
egvl=W
egvc=A


!write(88,*) AW

!write(55,*) AW, ( real(egvl(i)), i=1,N)

!write(55,*) AW, ( real(egvl(i)), i=1,N,4), ( real(egvl(i)), i=2,N,4)

!do i=1,N
! write(66,*) real(egvc(i,:))
! write(77,*) aimag(egvc(i,:))
! write(55,*) AW, real(egvl(i))
!enddo
!   do j=1,N
!    call spinvl(N,A(:,j),spz)
		 

!write(50,*) k*(a/(2.*pi)),W(j),spz

!    hh(i,j,1) = k!*(1.d0/(2.*pi))
!    hh(i,j,2) = W(j)
!    hh(i,j,3) = real(spz)
		
!    write(50,*) AW,hh(i,j,2),hh(i,j,3)
!			write(50,*)

!    end do
!write(*,*) A
!stop 

do j=1,n	 
hh(i,j,2)=0
hh(i,j,3)=0


call spinvl2(N,A(:,j),spz)
			!write(50,*) k*(a/(2.*pi)),W(j),spz

hh(i,j,1) = k!*(1.d0/(2.*pi))
hh(i,j,2) = W(j)
hh(i,j,3) = real(spz)

if (hh(i,j,2)<0.1d0 .and. hh(i,j,2)>-0.1d0 ) then
if(hh(i,j,3).GT. 1.d0) hh(i,j,3)= 0.99999d0  
if(hh(i,j,3).LT. -1.d0) hh(i,j,3)= -0.99999d0  
write(60,*) AW ,hh(i,j,2),hh(i,j,3)
write(60,*)
end if	
end do


ENDDO
deallocate(A,W,WORK,RWORK)
stop
end program main
!******************************************************************
!#######################################################################################
SUBROUTINE TERMOS_M3(Dmat,AW,POT,CHE,T1, T2,T3,SR,E1,SO)


	USE PARAMETROS	
	IMPLICIT NONE
    integer :: i ,j, MT,k
   REAL(8), INTENT(IN)   :: AW
   COMPLEX(8), INTENT(IN) :: POT, CHE, T1, T2,T3,SR,E1,SO
    COMPLEX(8),INTENT(OUT) :: Dmat(N,N)

!POTENCIAL QUIMICO E EFEITO ZEEMAN

	do i= 1 , N/4 
	do j= 1 , N/4
	!print*, i, j
		if (i==j) then
		Dmat(i,j) =-CHE +E1
		Dmat(i+N/4,j+N/4) = -CHE - E1
		Dmat(i+N/2,j+N/2) = CHE -E1
		Dmat(i+(3*N)/4,j+3*N/4) = CHE +E1
	end if
	end do
	end do	

!PIMEIROS VIZINHOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i= 1 , N/4, 2 
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j) then
		Dmat(i,j) = -C
		Dmat(j,i) = -C
		Dmat(i+N/4,j+N/4) = -C
		Dmat(j+N/4,i+N/4) = -C
                Dmat(i+N/2,j+N/2) = C
		Dmat(j+N/2,i+N/2) = C
		Dmat(i+(3*N)/4,j+3*N/4) = C
		Dmat(j+(3*N)/4,i+3*N/4) = C

	end if
	end do
	end do

	do i= 2 , N/4, 4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j) then
		Dmat(i,j) = -C
		Dmat(j,i) = -C
		Dmat(i+N/4,j+N/4) = -C
		Dmat(j+N/4,i+N/4) = -C
                Dmat(i+N/2,j+N/2) = C
		Dmat(j+N/2,i+N/2) = C
		Dmat(i+(3*N)/4,j+3*N/4) = C
		Dmat(j+(3*N)/4,i+3*N/4) = C

	end if
	end do
	end do


	do i= 4 , N/4, 4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+3==j) then
		Dmat(i,j) = -C
		Dmat(j,i) = -C
		Dmat(i+N/4,j+N/4) = -C
		Dmat(j+N/4,i+N/4) = -C
                Dmat(i+N/2,j+N/2) = C
		Dmat(j+N/2,i+N/2) = C
		Dmat(i+(3*N)/4,j+3*N/4) = C
		Dmat(j+(3*N)/4,i+3*N/4) = C

	end if
	end do
	end do


	do i= 1 , N/4, 4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+5==j) then
		Dmat(i,j) = -C
		Dmat(j,i) = -C
		Dmat(i+N/4,j+N/4) = -C
		Dmat(j+N/4,i+N/4) = -C
                Dmat(i+N/2,j+N/2) = C
		Dmat(j+N/2,i+N/2) = C
		Dmat(i+(3*N)/4,j+3*N/4) = C
		Dmat(j+(3*N)/4,i+3*N/4) = C

	end if
	end do
	end do
!######################TERCEIROS############################### 
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

	do i= 2 , N/4, 4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+5==j) then
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

	do i= 1 , N/4, 4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+7==j) then
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
!######################POTENCIAL SPIN ORBITA INTRINSICO###############################    


       do i= 1 , N/4 ,4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+2==j) then
		Dmat(i,j) = SO
		Dmat(j,i) = -SO
		Dmat(i+N/4,j+N/4) = -SO
		Dmat(j+N/4,i+N/4) = SO
                Dmat(i+N/2,j+N/2) = -SO
		Dmat(j+N/2,i+N/2) = SO
		Dmat(i+(3*N)/4,j+3*N/4) = SO
		Dmat(j+(3*N)/4,i+3*N/4) = -SO
	end if
	end do
	end do


        do i= 2 , N/4 ,4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+2==j) then
		Dmat(i,j) = SO
		Dmat(j,i) = -SO
		Dmat(i+N/4,j+N/4) = -SO
		Dmat(j+N/4,i+N/4) = SO
                Dmat(i+N/2,j+N/2) = -SO
		Dmat(j+N/2,i+N/2) = SO
		Dmat(i+(3*N)/4,j+3*N/4) = SO
		Dmat(j+(3*N)/4,i+3*N/4) = -SO
	end if
	end do
	end do

        do i= 4 , N/4 ,4  
	do j= 1 , N/4 
	!print*, i, j
		if (i+2==j) then
		Dmat(i,j) = -SO
		Dmat(j,i) = SO
		Dmat(i+N/4,j+N/4) = SO
		Dmat(j+N/4,i+N/4) = -SO
                Dmat(i+N/2,j+N/2) = SO
		Dmat(j+N/2,i+N/2) = -SO
		Dmat(i+(3*N)/4,j+3*N/4) = -SO
		Dmat(j+(3*N)/4,i+3*N/4) = SO
	end if
	end do
	end do


        do i= 1 , N/4 ,2 
	do j= 1 , N/4 
	!print*, i, j
		if (i+4==j) then
		Dmat(i,j) = -SO
		Dmat(j,i) = SO
		Dmat(i+N/4,j+N/4) = SO
		Dmat(j+N/4,i+N/4) = -SO
                Dmat(i+N/2,j+N/2) = SO
		Dmat(j+N/2,i+N/2) = -SO
		Dmat(i+(3*N)/4,j+3*N/4) = -SO
		Dmat(j+(3*N)/4,i+3*N/4) = SO
	end if
	end do
	end do

	do i= 2 , N/4 ,2
	do j= 1 , N/4 
	!print*, i, j
		if (i+4==j) then
		Dmat(i,j) = SO
		Dmat(j,i) = -SO
		Dmat(i+N/4,j+N/4) = -SO
		Dmat(j+N/4,i+N/4) = SO
                Dmat(i+N/2,j+N/2) = -SO
		Dmat(j+N/2,i+N/2) = SO
		Dmat(i+(3*N)/4,j+3*N/4) = SO
		Dmat(j+(3*N)/4,i+3*N/4) = -SO
	end if
	end do
	end do


	do i= 1 , N/4 ,4
	do j= 1 , N/4 
	!print*, i, j
		if (i+6==j) then
		Dmat(i,j) = SO
		Dmat(j,i) = -SO
		Dmat(i+N/4,j+N/4) = -SO
		Dmat(j+N/4,i+N/4) = SO
                Dmat(i+N/2,j+N/2) = -SO
		Dmat(j+N/2,i+N/2) = SO
		Dmat(i+(3*N)/4,j+3*N/4) = SO
		Dmat(j+(3*N)/4,i+3*N/4) = -SO
	end if
	end do
	end do


!######################INTERAÇÃO RASHBA###############################


	do i= 1+N/4 , N/2 
	do j= 1 , N/4 ,4
	!print*, i, j
		if (i-1==j+N/4) then
		Dmat(i,j) = R1
		Dmat(j,i) = R1
		Dmat(i+N/2,j+N/2) = -R12
		Dmat(j+N/2,i+N/2) = -R12
	end if
	end do
	end do


	do i= 1+N/4 , N/2 ,4
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j+N/4) then
		Dmat(i,j) = -R1
		Dmat(j,i) = -R1
		Dmat(i+N/2,j+N/2) = R12
		Dmat(j+N/2,i+N/2) = R12
	end if
	end do
	end do


	do i= 1+N/4 , N/2 
	do j= 2 , N/4 ,4
	!print*, i, j
		if (i-1==j+N/4) then
		Dmat(i,j) = R3
		Dmat(j,i) = R3
		Dmat(i+N/2,j+N/2) = -R32
		Dmat(j+N/2,i+N/2) = -R32
	end if
	end do
	end do


	do i= 2+N/4 , N/2 ,4
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


	do i= 1+N/4 , N/2 
	do j= 3 , N/4 ,4
	!print*, i, j
		if (i-1==j+N/4) then
		Dmat(i,j) = R2
		Dmat(j,i) = R2
		Dmat(i+N/2,j+N/2) = -R22
		Dmat(j+N/2,i+N/2) = -R22
	end if
	end do
	end do


	do i= 3+N/4 , N/2 ,4
	do j= 1 , N/4 
	!print*, i, j
		if (i+1==j+N/4) then
		Dmat(i,j) = -R2
		Dmat(j,i) = -R2
		Dmat(i+N/2,j+N/2) = R22
		Dmat(j+N/2,i+N/2) = R22
	end if
	end do
	end do


	do i= 1+N/4 , N/2 
	do j= 1 , N/4 ,4
	!print*, i, j
		if (i-5==j+N/4) then
		Dmat(i,j) = R2
		Dmat(j,i) = R2
		Dmat(i+N/2,j+N/2) = -R22
		Dmat(j+N/2,i+N/2) = -R22
	end if
	end do
	end do


	do i= 1+N/4 , N/2 ,4
	do j= 1 , N/4 
	!print*, i, j
		if (i+5==j+N/4) then
		Dmat(i,j) = -R2
		Dmat(j,i) = -R2
		Dmat(i+N/2,j+N/2) = R22
		Dmat(j+N/2,i+N/2) = R22
	end if
	end do
	end do


	do i= 1+N/4 , N/2 
	do j= 4 , N/4 ,4
	!print*, i, j
		if (i-3==j+N/4) then
		Dmat(i,j) = R2
		Dmat(j,i) = R2
		Dmat(i+N/2,j+N/2) = -R22
		Dmat(j+N/2,i+N/2) = -R22
	end if
	end do
	end do


	do i= 4+N/4 , N/2 ,4
	do j= 1 , N/4 
	!print*, i, j
		if (i+3==j+N/4) then
		Dmat(i,j) = -R2
		Dmat(j,i) = -R2
		Dmat(i+N/2,j+N/2) = R22
		Dmat(j+N/2,i+N/2) = R22
	end if
	end do
	end do

!######################POTENCIAL SUPER CONDUTOR###############################   
 
	do i= 1+N/2 , N/2+N/4 
	 do j= 2   , N/4 ,3
	!print*, i, j
		if (i+4==j+N/2) then
		Dmat(i,j) = POT
		Dmat(j,i) = -POT
		Dmat(i+N/4,j+N/4) = POT
		Dmat(j+N/4,i+N/4) = -POT
	end if
	end do
	end do

 	 do i= 1+N/2 , N/2+N/4 
	 do j= 1   , N/4 ,3
	!print*, i, j
		if (i-4==j+N/2) then
		Dmat(i,j) = -POT
		Dmat(j,i) = POT
		Dmat(i+N/4,j+N/4) = -POT
		Dmat(j+N/4,i+N/4) = POT
	end if
	end do
	end do

! 	 do i= 1+N/2 , N 
!	 do j= 2   , N/4 ,2
!	!print*, i, j
!		if (i-==j+N/2) then
!		Dmat(i,j) = POT
!		Dmat(j,i) = -POT
!		Dmat(i+N/4,j+N/4) = -POT
!		Dmat(j+N/4,i+N/4) = POT
!	end if
!	end do
!	end do
 
! 	 do i= 1+N/2 , N 
!	 do j= 1   , N/4 ,4
	!print*, i, j
!		if (i-1==j+N/2) then
!		Dmat(i,j) = -POT
!		Dmat(j,i) = POT
!		Dmat(i+N/4,j+N/4) = -POT
!		Dmat(j+N/4,i+N/4) = POT

!	end if
!	end do
!	end do
    
 !   	do i= 1+N/2 , N ,4
!	do j= 1   , N/4
	!print*, i, j
!		if (i+1==j+N/2) then
!		Dmat(i,j) = POT
!		Dmat(j,i) = -POT
!		Dmat(i+N/4,j+N/4) = POT
!		Dmat(j+N/4,i+N/4) = -POT

!	end if
!	end do
!	end do

 !	 do i= 1+N/2 , N 
!	 do j= 3   , N/4 ,4
	!print*, i, j
!		if (i-1==j+N/2) then
!		Dmat(i,j) = POT
!		Dmat(j,i) = -POT
!		Dmat(i+N/4,j+N/4) = POT
!		Dmat(j+N/4,i+N/4) = -POT

!	end if
!	end do
!	end do

 !   	do i= 3+N/2 , N ,4
!	do j= 1   , N/4
!	!print*, i, j
!		if (i+1==j+N/2) then
!		Dmat(i,j) = -POT
!		Dmat(j,i) = POT
!		Dmat(i+N/4,j+N/4) = -POT
!		Dmat(j+N/4,i+N/4) = POT

!	end if
!	end do
!	end do



!ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!#######################################################################################



!Print*, R1
! Print*, "  "
!	write(*,1) int(Dmat)
!	1 format(32i4)
!	Print*, "  "
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
!		spmz(i+n/2,i+n/2)= cmplx(1.,0.)
!		spmz(i+n/2+n/4,i+n/2+n/4)= cmplx(-1.,0.)
		!spmx(i,n/2+i)= cmplx(1.,0.)

		!spmx(i+n/2,i)= cmplx(1.,0.)


		!spmy(i,n/2+i)= cmplx(0.,-1.)

		!spmy(i+n/2,i)= cmplx(0.,1.)
		
!		write(*,1) int(spmz)
!		1 format(16i4) 
!		stop

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

	integer :: i


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
!	write(*,1) int(spz)
!	1 format(16i4)
!	stop

end subroutine spinvl2
