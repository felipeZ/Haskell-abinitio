
module BasisOrthogonalization where

import Data.Array.Repa as R
import Data.Array.Repa.Unsafe as R
import Data.Array.Repa.Algorithms.Matrix 
import qualified Data.List as DL
import qualified Data.Vector.Unboxed as VU

-- internal modules
import GlobalTypes
import qualified LinearAlgebra as LA
import Jacobi (jacobiP)



-- ========================> ORTHOGONALIZATION OF THE BASIS <======================================================

-- | Here is diagonalized the Overlap matrix and it is obtained a transformation matrix
--   S^-1/2
symmOrtho :: Monad m => Matrix -> m (Matrix)
symmOrtho arr = do
  (eigValUnboxed,eigVecs) <- jacobiP $ arr
  let (Z:.dim:._) = extent arr
      invSqrt = fromUnboxed (Z:. dim) $  VU.map (recip . sqrt) $ eigValUnboxed -- For building the S^-1/2 matrix
      diag = LA.vec2Diagonal invSqrt
  eigVecTrans <- transpose2P eigVecs
  mtx1 <- mmultP eigVecs diag
  mmultP mtx1 eigVecTrans
  
canortho :: Monad m => Matrix -> m (Matrix)
canortho arr = do
  (eigVal,eigVecs) <- jacobiP arr
  let  invSqrt = VU.map (recip . sqrt) eigVal -- For building the S^-1/2 matrix
  R.computeUnboxedP $  R.traverse eigVecs id (\f sh@(Z:._:. k) -> (invSqrt VU.! k) * f sh)


