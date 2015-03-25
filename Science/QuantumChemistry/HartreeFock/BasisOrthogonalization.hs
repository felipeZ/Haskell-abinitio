
module Science.QuantumChemistry.HartreeFock.BasisOrthogonalization  where

import Data.Array.Repa as R
import Data.Array.Repa.Unsafe as R
import Data.Array.Repa.Algorithms.Matrix (mmultP,transpose2P)
import qualified Data.List as DL
import qualified Data.Vector.Unboxed as VU

-- internal modules
import Science.QuantumChemistry.GlobalTypes
import qualified Science.QuantumChemistry.NumericalTools.LinearAlgebra as LA
import Science.QuantumChemistry.NumericalTools.JacobiMethod (jacobiP)



-- ========================> ORTHOGONALIZATION OF THE BASIS <======================================================

-- | Here is diagonalized the Overlap matrix and it is obtained a transformation matrix
--   S^-1/2
symmOrtho :: Monad m => Array U DIM2 Double -> m (Array U DIM2 Double)
symmOrtho arr = do
  (eigVal,eigVecs) <- jacobiP $ arr
  let (Z:.dim:._) = extent arr
      invSqrt = computeUnboxedS . R.map (recip . sqrt) $ fromUnboxed (Z:. dim) eigVal -- For building the S^-1/2 matrix
      diag = LA.vec2Diagonal invSqrt
  eigVecTrans <- transpose2P eigVecs
  mtx1 <- mmultP eigVecs diag
  mmultP mtx1 eigVecTrans
--   LA.mmultP eigenvecs <=< LA.mmultP eigenvecs $ diag 
  
canortho :: Monad m => Array U DIM2 Double -> m (Array U DIM2 Double)
canortho arr = do
  (eigVal,eigVecs) <- jacobiP arr
  let invSqrt = VU.map (recip . sqrt) eigVal -- For building the S^-1/2 matrix
  R.computeUnboxedP $  R.traverse eigVecs id (\f sh@(Z:._:. k) -> (invSqrt VU.! k) * f sh)


