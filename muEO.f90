
subroutine muEO(potquimEO)
use const
use system
use results
implicit none
real*8 potquimEO

potquimEO=log(xmEOalpha*vsol)-log(xmEObeta*vsol)-chi*Meo*(Ma*(xmAalpha-xmAbeta)&
+Meo*(xmEOalpha-xmEObeta))+Meo*(log(fEO_unas_alpha)-log(fEO_unas_beta))-packconst*Meo*vpol*vsol

end subroutine 
