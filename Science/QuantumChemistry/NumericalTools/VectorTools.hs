
module Science.QuantumChemistry.NumericalTools.VectorTools (
                                      diagonalVec
                                     ,sortEigenData
                                     ) where

import Control.Monad (liftM)
import Control.Monad.Primitive (PrimMonad, PrimState)
import Control.Monad.ST (runST)
import Data.Array.Repa.Index
import Data.Array.Repa (toIndex)
import Data.Vector.Algorithms.Merge (sort)
import Data.Vector.Unboxed as U
import Data.Vector.Generic.Mutable (new,unsafeWrite)


-- -----------------------------> Internal Modules <----------------
import Science.QuantumChemistry.GlobalTypes (VecUnbox)

-- --------------------------> <-----------------------------

sortEigenData :: (VecUnbox,VecUnbox) -> (VecUnbox,VecUnbox)
sortEigenData (vals,vss) = runST $   
 do   (sortVals,indexes) <- liftM U.unzip $ sortEigenVals vals
      eigenvecs          <- sortEigenVecs (U.length vals) vss indexes
      return (sortVals,eigenvecs)      

sortEigenVals :: PrimMonad m => VecUnbox -> 
                 m (Vector (Double, Int))
sortEigenVals xs = 
                do vs <- unsafeThaw ts 
                   sort vs 
                   unsafeFreeze vs
   where dim = U.length xs
         is  = U.generate dim id 
         ts  = U.zip xs is

sortEigenVecs :: PrimMonad m =>
                 Int        ->  
                 VecUnbox -> 
                 Vector Int ->
               m VecUnbox
sortEigenVecs dim mtx indexes = 
  do let dimMtx = U.length mtx
     vss <- new dimMtx 
     zipWithM_ (goM vss) (generate dim id) indexes 
     unsafeFreeze vss 

 where goM vss i k = do 
           let xs    = sliceColumn k dim mtx 
           writeColumn i dim vss xs 

sliceColumn ::  Unbox a =>
                  Int ->      -- Column number
                  Int ->      -- square Matrix dimension
                  Vector a -> -- Flatten matrix
                  Vector a 
sliceColumn k dim xs = generate dim $ \i -> xs ! (k + dim*i)
  

writeColumn ::  (Unbox a, PrimMonad m) =>
                     Int ->
                     Int -> 
                     MVector (PrimState m) a ->
                     Vector a -> 
                     m ()                    
writeColumn k dim mvec unvec = zipWithM_ (unsafeWrite mvec) indexes unvec 
  where indexes = generate dim $ \i -> k + dim*i


diagonalVec :: Unbox a => Vector a -> Vector a 
diagonalVec vss = unsafeBackpermute vss diagonalIdxs 
  where dim = floor . sqrt . fromIntegral $ U.length vss
        diagonalIdxs = U.map funIndex $ generate dim id
        funIndex i   = toIndex (ix2 dim dim) (ix2 i i)
        
