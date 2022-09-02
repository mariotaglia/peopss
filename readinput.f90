subroutine readinput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  This routine reads variables from fort.8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use pks
use system
use solver
!use kai

implicit none
integer i
character basura

! read starts here, not that read is performed sequentially! 

read(8,*), basura
read(8,*), Ma    ! Ma  polA

read(8, *), basura
read(8, *), MEO  !   MEO  polEO

read(8, *), basura! vp
read(8, *), vpolcero !

read(8, *), basura! rsal
read(8, *), rsal !

read(8, *), basura
read(8, *), logn_toti, logn_totf, npasosn_tot  !  scan total monomer density in molar

read(8, *), basura
read(8, *), logratioEOAalphai, logratioEOAalphaf, npasosratioEO   ! scan ratio monomer density
!
read(8, *), basura
read(8, *), pKDi, pKDf, pKDp     ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura
read(8, *), pKAi, pKAf, pKAp      ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura
read(8, *), pKEOi, pKEOf, pKEOp     ! polymer-polymer attraction strenght in kBT

read(8, *), basura
read(8, *), chi  ! cutoff for porr sv interaction in lattice sites

read(8,*) basura
read(8,*) justone ! solves only one point

read(8, *), basura
read(8, *), xmNaalphainitial, basura, basura, deltaxmNabetainitial,n_totbetainitial,ratioEOAbetainitial ! read initial guess in the same order as output, ratio and xmtot in alpha are not solved for, so the result is stores in basura

xmNabetainitial = xmNaalphainitial+ deltaxmNabetainitial
end subroutine
