!  ###############################################################################!     
  !     PEO/Pss Molecular Theory Program 
  !    
  !###############################################################################
!use pks
use system
use const
use solver
use results
implicit none
integer i, j, iii
real*8 logratioalpha
real*8 logxmtotalpha
integer k, kk, kkk


print*, 'GIT Version: ', _VERSION

call readinput ! read input from file
call allocation

vab=1.
vpol=vpolcero/vsol

vneg=4./3.*pi*rsal**3/vsol !volume of anion in units of vsol
vpos=4./3.*pi*rsal**3/vsol !volume of cation in units of vsol
yes=0 ! es para  chequear si encuentra o no xalpha, xbeta


do k=1,pKDp
do kk=1,pKAp
do kkk=1,pKEOp

pKD = pKDi + (pKDf-pKDi)*float(k-1)/float(pKDp)
pKA = pKAi + (pKAf-pKAi)*float(kk-1)/float(pKAp)
pKEO = pKEOi + (pKEOf-pKEOi)*float(kkk-1)/float(pKEOp)

KD=10**(-pKD)
KA=10**(-pKA)
KEO=10**(-pKEO)

K0A = (KA)*(vsol*Na/1.0d24) ! thermodynamic constants
K0EO = (KEO)*(vsol*Na/1.0d24)
K0D = (KD)*(vsol*Na/1.0d24)**2 !!!!!!!!!!!!!!!!!!!!!!!

do i = 1, npasosratioEO ! loop over ratio_alpha

  logratioalpha = logratioEOAalphai  + (logratioEOAalphaf-logratioEOAalphai) &
  /float(npasosratioEO)*float(i-1)  !Na

  ratioEOAalpha = 10**(logratioalpha)

   do j=1, npasosn_tot  ! loop over xmtot alpha

      logn_totalpha = (logn_totf-logn_toti)*float(j-1)/float(npasosn_tot) + logn_toti
      n_totalpha= 10**(logn_totalpha)  !Segunda variable que fijamos  xmpoltotalalpha

      iter=0
      call solve

  enddo ! j

enddo ! i

enddo !k
enddo !kk
enddo !kkk


! SAVE RESULTS TO FILE

open (unit=3,file='csal_poltot_mol_alpha.txt',status='replace')

do iii=1,yes
   write (3,*) arraympoltot(1,iii), arraymcsal(1,iii)
end do

open (unit=4,file='csal_poltot_mol_beta.txt',status='replace')

do iii=1,yes
   write (4,*) arraympoltot(2,iii), arraymcsal(2,iii)
end do

open (unit=40,file='cpoltot_ratioEOA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (40,*) arrayratioEOA(1,iii), arraympoltot(1,iii)
end do

open (unit=30,file='cpoltot_ratioEOA_mol_beta.txt',status='replace')

do iii=1,yes
   write (30,*) arrayratioEOA(2,iii), arraympoltot(2,iii)
end do

open (unit=400,file='csal_ratioEOA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (400,*) arrayratioEOA(1,iii), arraymcsal(1,iii)
end do

open (unit=300,file='csal_ratioEOA_mol_beta.txt',status='replace')

do iii=1,yes
   write (300,*) arrayratioEOA(2,iii), arraymcsal(2,iii)
end do

call endall     ! clean up and terminate
end 



subroutine endall
 stop
end subroutine
