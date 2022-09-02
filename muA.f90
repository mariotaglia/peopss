
subroutine muA(potquimA)
use const
use system
use results
implicit none
real*8 potquimA

potquimA= log(xmAalpha*vsol)-log(xmAbeta*vsol)- chi*MA*(MA*(xmAalpha-xmAbeta)&
+MEO*(xmEOalpha-xmEObeta))+MA*(log(fA_unas_alpha)-log(fA_unas_beta))&
-packconst*Ma*vpol*vsol+Ma*neutralconst


end  subroutine 
