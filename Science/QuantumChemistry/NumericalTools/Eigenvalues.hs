
module Science.QuantumChemistry.NumericalTools.EigenValues (
                   eigenSolve
                   ) where

import Control.Arrow ((***))
import Data.Array.Repa as R
import Numeric.LinearAlgebra.HMatrix (eigSH')
import qualified Data.Packed.Vector as HV
import qualified Data.Packed.Matrix as HM

import Data.Array.Repa as R
import Data.Vector.Storable (convert)
import qualified Data.Vector.Unboxed as U


-- -----------------> Internal Modules <----------

import Science.QuantumChemistry.GlobalTypes (VecUnbox)
import Science.QuantumChemistry.NumericalTools.VectorTools (sortEigenData)

-- -------------------> <-----------------

eigenSolve :: Array U DIM2 Double -> (VecUnbox, Array U DIM2 Double)
eigenSolve arr = (vs, fromUnboxed (extent arr) mtx )
  where (vs,mtx) = sortEigenData $ hv2VecUnbox *** hm2VecUnbox $ eigSH' $ repa2HM arr 


-- ---------------------> <------------------------------
vecUnbox2HV :: VecUnbox -> HV.Vector Double 
vecUnbox2HV = convert 

hv2VecUnbox :: HV.Vector Double -> VecUnbox
hv2VecUnbox = convert 

repa2HM :: Array U DIM2 Double -> HM.Matrix Double 
repa2HM arr = HM.reshape dim . convert . toUnboxed $ arr 

 where (Z :. dim :. _) = extent arr

hm2Repa :: HM.Matrix Double -> Array U DIM2 Double
hm2Repa mtx = fromUnboxed (ix2 dim dim) . convert . HM.flatten $ mtx
 where dim = HM.rows mtx

hm2VecUnbox :: HM.Matrix Double -> VecUnbox
hm2VecUnbox = convert . HM.flatten
