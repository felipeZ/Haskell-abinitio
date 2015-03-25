{-# LANGUAGE BangPatterns #-}

module Science.QuantumChemistry.NumericalTools.JacobiMethodST (
               jacobiST
              ) where

import Control.Arrow (second)
import Control.Monad (liftM)
import Control.Monad.Primitive (PrimMonad, PrimState)
import Control.Monad.ST (runST)
import qualified Data.List as DL
import Data.Array.Repa  (extent,fromIndex,fromUnboxed,toIndex,toUnboxed)
import Data.Array.Repa.Index
import Data.Vector.Generic.Mutable (unsafeRead,unsafeWrite)
import Data.Vector.Unboxed as U


-- ----------------------> Internal Modules <---------------

import Science.QuantumChemistry.GlobalTypes (Matrix,VecUnbox)
import Science.QuantumChemistry.NumericalTools.LinearAlgebra (diagonal,identity)
import Science.QuantumChemistry.NumericalTools.VectorTools (diagonalVec,sortEigenData)

-- ============================> Types <=====================

type Step      = Int
type Tolerance = Double

data Parameters = Parameters !Double !Double !Double deriving Show

-- ============================> <======================
   
jacobiST :: Matrix -> (VecUnbox,Matrix)
jacobiST !arr = second (fromUnboxed (ix2 dim dim)) $ sortEigenData $ 
                runST $ do 
                vss <- unsafeThaw $ toUnboxed arr
                ide <- unsafeThaw $ toUnboxed $ identity dim
                loopJacobi vss ide dim 0 (1.0e-9)  
  where (Z:.dim:. _) = extent arr

                  

-- | Loop to carry out the Jacobi Method rotations
loopJacobi :: PrimMonad m => MVector (PrimState m)  Double ->
                             MVector (PrimState m)  Double ->
                             Int       -> 
                             Step      -> 
                             Tolerance ->
                             m (VecUnbox,VecUnbox)
loopJacobi !arr !prr dim step tolerance = 
      if step > 5*dim*dim
         then error "Jacobi method did not converge "
         else do 
           sh      <- maxElemIndex arr dim
           maxElem <- unsafeRead arr $ funIdx sh 
           if maxElem  > tolerance 
              then action sh 
              else  do 
                    eigVals <- liftM diagonalVec $ unsafeFreeze arr
                    eigVecs <- unsafeFreeze prr
                    return (eigVals,eigVecs)

  where funIdx = toIndex (ix2 dim dim)  
        action sh@(Z:.k:.l) = do
           maxElem         <- unsafeRead arr $ funIdx sh 
           [elemLL,elemKK] <- Prelude.mapM (unsafeRead arr . funIdx ) [ix2 l l, ix2 k k]
           let aDiff      = elemLL - elemKK
               parameters = calcParameters maxElem aDiff 
           newA <- rotateA arr sh parameters dim funIdx
           newP <- rotateP prr sh parameters dim funIdx
           loopJacobi newA newP dim (succ step) tolerance           
                                         
rotateA :: PrimMonad m => MVector (PrimState m) Double -> 
                          DIM2          -> 
                          Parameters    -> 
                          Int           ->
                          (DIM2 -> Int) ->
                          m (MVector (PrimState m) Double)
rotateA arr kl@(Z:. k :. l) (Parameters s t tau) dim funIdx = 
 do  maxElem <- unsafeRead arr $ funIdx kl     
     U.foldM'_  (update maxElem) () (generate (dim*dim) id)  
     return arr 

  where update maxElem () idx = do 
            val <- unsafeRead arr idx  
            let sh        = fromIndex (ix2 dim dim) idx
                uwrite    = unsafeWrite arr idx
                gix2 a b  = unsafeRead arr (funIdx $ ix2 a b)
                uwriteijP i j = gix2 i j >>= \x -> uwrite $ val - s*( x + tau*val)
                uwriteijM i j = gix2 i j >>= \x -> uwrite $ val - s*( x - tau*val)
                fun (Z:.n:.m) 
                      | (n,m) == (k,l)           = uwrite $ 0.0
                      | (n,m) == (k,k)           = uwrite $ val - t*maxElem
                      | (n,m) == (l,l)           = uwrite $ val + t*maxElem
                      | n < k && m == k          = uwriteijP n l 
                      | n < k && m == l          = uwriteijM n k 
                      | k < m && m < l && n == k = uwriteijP m l 
                      | k < n && n < l && m == l = uwriteijM k n 
                      | m > l && n == k          = uwriteijP l m 
                      | m > l && n == l          = uwriteijM k m 
                      | otherwise                = return ()                                   
            fun sh
{-# Inline rotateA #-}                     
   
rotateP :: PrimMonad m => MVector (PrimState m) Double ->
                          DIM2       ->
                          Parameters ->
                          Int        -> 
                          (DIM2 -> Int) ->
                          m (MVector (PrimState m) Double)
rotateP prr kl@(Z:. k :. l) (Parameters !s !t !tau) dim funIdx =  
       U.foldM'_  update () (generate (dim*dim) id)  >> return prr

  where update () idx = do 
          val <- unsafeRead prr idx  
          let sh       = fromIndex (ix2 dim dim) idx
              uwrite   = unsafeWrite prr idx 
              gix2 a b = unsafeRead prr (funIdx $ ix2 a b)           
              fun (Z:.x:.y)                      
                 | y == k    = gix2 x l >>= \x -> uwrite $ val - s*(x + tau*val)
                 | y == l    = gix2 x k >>= \x -> uwrite $ val + s*(x - tau*val)
                 | otherwise = return ()
          fun sh

{-# Inline rotateP #-}                     

calcParameters ::  Double -> Double -> Parameters                 
calcParameters !maxElem !aDiff = Parameters s t tau
  where t = if (abs maxElem < abs aDiff *1.0e-36) then maxElem/aDiff
            else let phi = aDiff/(2.0 *maxElem )
                     var = recip (abs phi + (sqrt $ 1 + phi^2))
                 in if phi < 0 then negate var else var
        c = recip $ sqrt(t^2 + 1)
        s = t*c
        tau = s/(1 + c)
        
-- Return the index of the largest off-diagonal element in the array
maxElemIndex  :: PrimMonad m => 
                   MVector (PrimState m) Double ->
                   Int -> 
                   m DIM2
maxElemIndex !arr dim  = liftM (fromIndex sh) $ U.foldM' fun 1 inds

  where inds = U.enumFromN 0 (dim*dim) 
        sh   = ix2 dim dim 
        fun n acc = do let (Z:. i:.j) = fromIndex sh n
                       x <- unsafeRead arr n
                       y <- unsafeRead arr acc
                       return $  if i < j then if (abs x > abs y) then n else acc
                                          else acc                

