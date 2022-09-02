subroutine fkfun2(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use solver
use system
use const
use results
implicit none
 
real*8 xmsalt
real*8 Penality,testKa,testkeo,testkd
real*8 testneuta,testneutb,testpcka,testpckb
integer*4 ier
real*16 vectfalpha(3),vectfbeta(3),vectfrac(6)
real*8 x(4),f(4)
real*8 potA,potEo,potNa,free_ener
real*8 muAalpha,muAbeta,muEOalpha,muEObeta,fealpha,febeta
real*8 potquimA,elib,potquimEO,potquimNa,muNaalpha,muNabeta
real*8 penalityA,penalityNa,penalityEO
integer i, j
real*8 delta

xmNAalpha=exp(x(1))        !xmNa en alfa
xmNAbeta=exp(x(2))        !xmNa en beta

n_totbeta=exp(x(3))    ! total monomer in beta
ratioEOAbeta=exp(x(4))   !xmpeo/xmA en beta

xmAalpha=n_totalpha/(1.+ratioEOAalpha)/MA ! chain number densities
xmEOalpha=xmAalpha*ratioEOAalpha*MA/MEO

xmAbeta=n_totbeta/(1.+ratioEOAbeta)/MA
xmEObeta=xmAbeta*ratioEOAbeta*MA/MEO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FRACTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vectfalpha(1)=xmAalpha
vectfalpha(2)=xmEOalpha
vectfalpha(3)=xmNaalpha
vectfbeta(1)=xmAbeta
vectfbeta(2)=xmEObeta
vectfbeta(3)=xmNabeta

call fractions(vectfalpha,vectfrac)

!vectfractions-> (1) fEO_aspol,(2)fA_aspol,(3) fEO_asion (4)fA_Asion,(5)fEO_unas,f(6)fA_unas
fEO_aspol_alpha=vectfrac(1)
fA_aspol_alpha=vectfrac(2)
fEO_asion_alpha=vectfrac(3)
fA_asion_alpha=vectfrac(4)
fEO_unas_alpha=vectfrac(5)
fA_unas_alpha=vectfrac(6)

!print*, 'ALFA', 'fEO_aspol:', vectfrac(1), 'fA_aspol:', vectfrac(2)
!print*, 'ALFA', 'fEO_asion:', vectfrac(3), 'fA_asion:', vectfrac(4)
!print*, 'ALFA', 'fEO_unas:', vectfrac(5), 'fA_unas:', vectfrac(6)
!print*,'ALFA, test', fEO_aspol_alpha*xmEOalpha*Meo, fA_aspol_alpha*xmAalpha*Ma

!stop

!testeo fracciones
!testKa=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fA_asion_alpha-fA_aspol_alpha)/fA_asion_alpha )-pKa
!testkeo=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fEO_asion_alpha-fEO_aspol_alpha)/fEO_asion_alpha )-pKEO
!testkd=-log10((Na/1.0d24)*xmnaalpha*vsol*fA_unas_alpha*Ma*xmAalpha*vab*(1.-fEo_asion_alpha-fEO_aspol_alpha)/fEO_aspol_alpha) -pKd
!print*,'testalpha',testKa,testkeo,testkd
!stop

call fractions(vectfbeta,vectfrac)

!vectfractions-> (1) fEO_aspol,(2)fA_aspol,(3) fEO_asion (4)fA_Asion,(5)fEO_unas,f(6)fA_unas

fEO_aspol_beta=vectfrac(1)
fA_aspol_beta=vectfrac(2)
fEO_asion_beta=vectfrac(3)
fA_asion_beta=vectfrac(4)
fEo_unas_beta=vectfrac(5)
fA_unas_beta=vectfrac(6)


!print*, 'BETA', 'fEO_aspol:', vectfrac(1), 'fA_aspol:', vectfrac(2)
!print*, 'BETA', 'fEO_asion:', vectfrac(3), 'fA_asion:', vectfrac(4)
!print*, 'BETA', 'fEO_unas:', vectfrac(5), 'fA_unas:', vectfrac(6)


!testeo fracciones
!testKa=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fA_asion_alpha-fA_aspol_alpha)/fA_asion_alpha )-pKa
!testkeo=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-feO_asion_alpha-fEo_aspol_alpha)/feO_asion_alpha )-pKEO
!testkd=-log10((Na/1.0d24)*xmNabeta*vsol*fA_unas_beta*Ma*xmAbeta*vab*(1-fEo_asion_beta-fEO_aspol_beta)/fEO_aspol_beta) -pKd
!print*,'testbeta',testKa,testkeo,testkd
!stop
!print*,'beta',vectfrac
!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  AUXILIARY CALC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

xSolventalpha=1. -Ma*vpol*vsol*xmAalpha -Meo*vpol*vsol*xmEOalpha&
-vpos*vsol*((fA_asion_alpha+fA_aspol_alpha)*Ma*xmAalpha+fEO_asion_alpha*Meo*xmEOalpha)&
-xmNaalpha*vpos*vsol 

xSolventbeta=1. - Ma*vpol*vsol*xmAbeta -Meo*vpol*vsol*xmEObeta&
-vpos*vsol*((fA_asion_beta+fA_aspol_beta)*Ma*xmAbeta+fEO_asion_beta*Meo*xmEObeta)&
-xmNabeta*vpos*vsol 

!if (xsolventalpha.lt.0.0) xsolventalpha = 1.e-15
!if (xsolventbeta.lt.0.0) xsolventbeta = 1.e-15

xmSolventalpha=xSolventalpha/vsol
xmSolventbeta =xSolventbeta/vsol

packconst=(1./vsol)*(log(xSolventalpha)-log(xSolventbeta) ) ! betapi 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!  CHEMICAL POTS AND FREE ENERGY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


potquimNa=0.
call muNa(potquimNa)
potNa=potquimNa

call muA(potquimA)
potA = potquimA

potquimEO=0.
call muEO(potquimEO)
potEO = potquimEO

elib=0.
call fe(elib)
free_ener=elib

Penality=1
Penality=abs(xmNabeta-xmNaalpha)/(xmNabeta*0.5+xmNaalpha*0.5)
Penality=Penality+abs(n_totalpha-n_totbeta)/(n_totalpha*0.5+n_totbeta*0.5)
Penality=Penality+abs(ratioEOAalpha-ratioEOAbeta)/(ratioeoaalpha*0.5+ratioeoabeta*0.5)

f(1)=free_ener/Penality
f(2)=potNa/Penality
f(3)=potA/Penality
f(4)=potEO/Penality

iter = iter + 1
norma = 0.0

do i = 1, 4
!  print*,'f',f(i)
  norma = norma +(f(i))**2    
enddo

!testneuta=xmNaalpha -xmClalpha-Ma*xmAalpha*fA_unas_alpha+fEO_asion_alpha*Meo*xmEOalpha
!testneutb=xmNabeta -xmClbeta-Ma*xmAbeta*fA_unas_beta+fEO_asion_beta*Meo*xmEObeta

!testpcka=1.-Ma*xmAalpha*vpol*vsol-xmSolventalpha*vsol -Meo*xmEOalpha*vpol*vsol-xmNaalpha*vpos*vsol-xmClalpha*vneg*vsol
!testpcka=testpcka -vneg*vsol*(Ma*xmAalpha*(fA_asion_alpha+fA_aspol_alpha)+Meo*xmEOalpha*fEO_asion_alpha)
 
!testpckb=1.-Ma*xmAbeta*vpol*vsol-xmSolventbeta*vsol -Meo*xmEObeta*vpol*vsol-xmNabeta*vpos*vsol-xmClbeta*vneg*vsol
!testpckb=testpckb -vneg*vsol*(Ma*xmAbeta*(fA_asion_beta+fA_aspol_beta)+Meo*xmEObeta*fEO_asion_beta)
!print*,'testneutra',testneuta,testneutb,testpcka,testpckb


ier = 0.0
return
end subroutine
