
{-|
Module: Science.QuantumChemistry.NumericalTools.EigenValues
Description: Compute all the eigenvalues of the real symmetrical matrix
Copyright: @2016 Felipe Zapata
-}


module Science.QuantumChemistry.NumericalTools.EigenValues (
                   eigenSolve
                   ) where

-- =============================> Standard and third party libraries <===============================
import Control.Arrow ((***))
import Data.Array.Repa as R
import Data.Vector.Storable (convert)
import qualified Data.Vector.Unboxed as U
import qualified Numeric.LinearAlgebra.Data as ND
import Numeric.LinearAlgebra.HMatrix (eigSH')


-- =================> Internal Modules <======================
import Science.QuantumChemistry.GlobalTypes (VecUnbox)
import Science.QuantumChemistry.NumericalTools.JacobiMethod (jacobiP)
import Science.QuantumChemistry.NumericalTools.VectorTools (sortEigenData)


-- | Compute all the eigenvalues and eigenvectors of symmetrical matrix.
eigenSolve :: Array U DIM2 Double -> (VecUnbox, Array U DIM2 Double)
eigenSolve arr = (vs, fromUnboxed (extent arr) mtx )
  where (vs,mtx) = sortEigenData $ hv2VecUnbox *** hm2VecUnbox $ eigSH' $ repa2HM arr 


vecUnbox2HV :: VecUnbox -> ND.Vector Double 
vecUnbox2HV = convert 

hv2VecUnbox :: ND.Vector Double -> VecUnbox
hv2VecUnbox = convert 

repa2HM :: Array U DIM2 Double -> ND.Matrix Double 
repa2HM arr = ND.reshape dim . convert . toUnboxed $ arr 

 where (Z :. dim :. _) = extent arr

hm2Repa :: ND.Matrix Double -> Array U DIM2 Double
hm2Repa mtx = fromUnboxed (ix2 dim dim) . convert . ND.flatten $ mtx
 where dim = ND.rows mtx

hm2VecUnbox :: ND.Matrix Double -> VecUnbox
hm2VecUnbox = convert . ND.flatten
