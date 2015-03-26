{-# LANGUAGE BangPatterns #-}

module Science.QuantumChemistry.NumericalTools.JacobiMethodST (
               jacobiST
              ,maxElemIndex
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
     uwrite (ix2 k l) 0
     uread (ix2 k k) >>= \akk -> uwrite (ix2 k k) $ akk - t*maxElem
     uread (ix2 l l) >>= \all -> uwrite (ix2 l l) $ all + t*maxElem
     loopArr (\i -> ix2 i k) (\i -> ix2 i l) 0 k
     loopArr (\i -> ix2 k i) (\i -> ix2 i l) (k+1) (l-k-1)
     loopArr (\i -> ix2 k i) (\i -> ix2 l i) (l+1) (dim-l-1)
     return arr 

  where uwrite idx = unsafeWrite arr (funIdx idx)
        uread  = unsafeRead arr . funIdx
        loopArr f g from len = U.forM_ (enumFromN from len ) $ \i ->  
                do  val <- uread $ f i
                    ail <- uread $ g i
                    uwrite (f i) $ val - s*(ail + tau*val)
                    uwrite (g i) $ ail + s*(val - tau*ail)
{-# Inline rotateA #-}                     
   
rotateP :: PrimMonad m => MVector (PrimState m) Double ->
                          DIM2       ->
                          Parameters ->
                          Int        -> 
                          (DIM2 -> Int) ->
                          m (MVector (PrimState m) Double)
rotateP prr kl@(Z:. k :. l) (Parameters !s !t !tau) dim funIdx =  do
    U.forM_ (generate dim id ) $ \i ->  
                do  val <- uread $ (ix2 i k)
                    pil <- uread $ (ix2 i l)
                    uwrite (ix2 i k) $ val - s*(pil + tau*val)
                    uwrite (ix2 i l) $ pil + s*(val - tau*pil)
    return prr 

  where uwrite idx = unsafeWrite prr (funIdx idx)
        uread      = unsafeRead prr . funIdx


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

  where inds = U.generate (dim*dim) id
        sh   = ix2 dim dim 
        fun acc n = do let (Z:. i:.j) = fromIndex sh n
                       x <- unsafeRead arr n
                       y <- unsafeRead arr acc
                       return $  if i < j then if (abs x > abs y) then n else acc
                                          else acc                

