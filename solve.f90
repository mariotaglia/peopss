subroutine solve

!use pks 
use system
use results
use solver
use const

implicit none
real*8 tolerancia,criterio,check_Ka_alpha,check_ka_beta,checkb_Kai,KK0check,KKaAcheckplus,kkaBcheckmin
real*8 x1(4)
real*8 x1g(4)
real*8 checkresults
integer ier, i,newt, j, ii, jj
integer flag
integer ngrid,iii

criterio=1E-4!criterio pra la norma
tolerancia=1E-4!criterio pra la diferencia de concentraciones relativa

x1(1)=log(xmNaalphainitial) ! Na in alpha
x1g(1)=x1(1)  

x1(2)=log(xmNabetainitial)     !xNabeta  inicial
x1g(2)=x1(2)

x1(3)=log(n_totbetainitial) ! total monomer density, beta
x1g(3)=x1(3)

x1(4)=log(ratioEOAbetainitial) ! ratioeobeta inicial 
x1g(4)=x1(4)

print*,'Na_alfa', 'n_tot_alfa', 'EO/Na alfa', 'Na_beta-Na_alpha', 'n_tot_beta', 'EO/Na beta'

call call_kinsol(x1, x1g, ier)

print*,'out:',exp(x1(1)), n_totalpha, ratioEOAalpha , exp(x1(2))-exp(x1(1)), exp(x1(3)), exp(x1(4))

!checkresults=0.
!checkresults=abs((x1(1)- x1(2)/(x1(1)+x1(2)))) ! Naalpha - Nabeta
!checkresults=checkresults+abs((n_totalpha- x1(3))/(n_totalpha + x1(3))) ! n_totalpha - n_totbeta
!checkresults=checkresults+abs((log(ratioEOAalpha) - x1(4))/(log(ratioEOAalpha) - x1(4))) ! ratio alfa - ratio beta

!print*, norma
!print*, checkresults
!stop


!if ((norma.lt.criterio).and.( checkresults.gt.tolerancia)) then ! encuentra solucion
if (norma.lt.criterio) then ! encuentra solucion

    print*,'Grid Point OK',yes
    if(justone.eq.1)call endall


! Found a solution, save in arrays

! Solutions for next iteration
      xmNaalphainitial=exp(x1(1))! Ratio inicial alpha 1:1
      xmNabetainitial=exp(x1(2))     !xNabeta  inicial
      n_totbetainitial=exp(x1(3)) ! for next iteration
      ratioEOAbetainitial=exp(x1(4)) !xratioeobeta inicial 

! save arrays

      n_totbeta = exp(x1(3))
      ratioEOAbeta = exp(x1(4))
      xmNaalpha =exp(x1(1)) 
      xmNabeta  =exp(x1(2))    

      yes=yes+1 ! counter of sucessful solutions

      write(9999,*)yes, xmsolventalpha,xmsolventbeta,xmClalpha, xmClbeta
      write(9998,*)yes,fEO_aspol_alpha, fA_aspol_alpha, fEO_asion_alpha,fA_asion_alpha
      write(9997,*)yes,fEO_aspol_beta, fA_aspol_beta, fEO_asion_beta,fA_asion_beta
      write(9996,*)yes, xmNaalpha,xmNabeta, xmNaalpha-xmNabeta
      write(9995,*)yes, packconst, neutralconst
  

      arraymNa(1,yes)=xmNaalpha/Na*1.d24
      arraymNa(2,yes)=xmNabeta/Na*1.d24

      arraymCl(1,yes)=xmClalpha/Na*1.d24
      arraymCl(2,yes)=xmClbeta/Na*1.d24

      arraymcsal(1,yes)=arraymCl(1,yes) ! salt can be calculated from Cl- only
      arraymcsal(2,yes)=arraymCl(2,yes)

      arraymA(1,yes)=MA*xmAalpha/Na*1.e24 ! in units of molar monomers 
      arraymA(2,yes)=MA*xmAbeta/Na*1.e24

      arraymEO(1,yes)=MEO*xmEOalpha/Na*1.e24  ! in units of molar monomers
      arraymEO(2,yes)=MEO*xmEObeta/Na*1.e24

      arraympoltot(1,yes)=arraymA(1,yes)+arraymEO(1,yes)
      arraympoltot(2,yes)=arraymA(2,yes)+arraymEO(2,yes)

      arrayratioEOA(1,yes)=arraymEO(1,yes)/arraymA(1,yes)
      arrayratioEOA(2,yes)=arraymEO(2,yes)/arraymA(2,yes)

endif ! found solution

return
end subroutine

