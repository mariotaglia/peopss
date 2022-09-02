subroutine fe(elib)

use results
use system
use const
implicit none
real*8 Free_energy,elib 

Free_Energy = 0.0
Free_Energy=Free_Energy-xmAalpha-xmEOalpha-xmSolventalpha-xmNaalpha-xmClalpha
Free_Energy=Free_Energy + 0.5*chi*(xmAalpha*Ma+xmEOalpha*MEO)**2+xmAalpha*Ma*fA_aspol_alpha
Free_Energy=Free_Energy+xmAbeta+xmEObeta+xmSolventbeta+xmNabeta+xmClbeta
Free_Energy=Free_Energy+packconst-0.5*chi*(Ma*xmAbeta+Meo*xmEObeta)**2-xmAbeta*Ma*fA_aspol_beta

elib=Free_Energy

return 
end subroutine

