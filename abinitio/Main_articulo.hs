
import HartreeFock
import qualified LinearAlgebra as LA
import Control.Monad(liftM,ap)
import GlobalTypes
import IntegralsEvaluation

main = do
  let sto3g = [(0.444635,0.168856),(0.535328,0.623913),(0.154329,3.42525)]
      sto3gHe = [ (0.182985,0.480844),(0.587136,1.776691),(0.607075,9.753934)]
      r1 = [0.0,0.0,1.4632]
      r2 = [0.0,0.0,0.0]
      base1 = [CGF sto3gHe S]
      base2 = normaCoeff`fmap` [CGF sto3g S]
      listBasis = [base1, base2]
      listCoord = [r1,r2]
      zlist = [2.0,1.0]
  result <- scfHF listCoord listBasis zlist
  print "HF"
  print $ getEnergy result
  print "Final Density"
  print $ getDensity result
  print "Orbital Energy"
  print $ getOrbE result
  print "Final Coeeficients"
  print $ getCoeff result
  print "Final Fock Matrix"
  print $ getFock result


  