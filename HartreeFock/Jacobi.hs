{-# LANGUAGE BangPatterns #-}

module Jacobi (
               Matrix,
               VecUnbox,
               jacobiP,
               sortEigenV,
               testJacobi
              ) where

import Control.Applicative
import Control.Monad
import qualified Data.List as DL
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import qualified Data.Vector.Unboxed as VU

-- ==============> Internal Imports <==========
import GlobalTypes  (Matrix,Step,Tolerance,VecUnbox)          
import LinearAlgebra (diagonal,identity)

-- ============================> Types <=====================

data Parameters = Parameters !Double !Double !Double deriving Show


-- ============================> <======================
   
jacobiP :: Monad m => Matrix -> m (VecUnbox,Matrix)
jacobiP !arr = sortEigenV =<< loopJacobi arr ide 0 (1.0e-9)                     
  where ide = identity $ extent arr

-- | Loop to carry out the corresponding rotation of the Jacobi Method
loopJacobi :: Monad m => Matrix -> Matrix -> Step -> Tolerance -> m (VecUnbox,Matrix)
loopJacobi !arr !prr step tolerance  
                                   | step > 5*dim^2 = error "Jacobi method did not converge "                             
                                   | otherwise = if abs maxElem > tolerance 
                                                    then action 
                                                    else  liftM2 (,) (diagonal arr) (return prr)
  where (Z:.dim:. _) = extent arr
        mx@(Z:.k:.l) = maxElemIndex arr
        maxElem      = (arr ! mx )
        aDiff        = (arr ! (ix2 l l)) - (arr ! (ix2 k k))
        parameters   = calcParameters (arr ! mx) aDiff
        action       = do            
           newA <- rotateA arr mx parameters
           newP <- rotateP prr mx parameters
           loopJacobi newA newP (succ step) tolerance
           
                                         
rotateA :: Monad m => Matrix -> DIM2 -> Parameters -> m Matrix         
rotateA arr kl@(Z:. k :. l) (Parameters !s !t !tau) = computeUnboxedP $ unsafeTraverse arr id fun

  where fun f sh@(Z:.n:.m)  
                  | (n,m) == (k,l)           = 0.0
                  | (n,m) == (k,k)           = val - t*maxElem
                  | (n,m) == (l,l)           = val + t*maxElem
                  | n < k && m == k          = val - s*(gix2 n l + tau*val)
                  | n < k && m == l          = val + s*(gix2 n k - tau*val)
                  | k < m && m < l && n == k = val - s*(gix2 m l + tau*val)
                  | k < n && n < l && m == l = val + s*(gix2 k n - tau*val)
                  | m > l && n == k          = val - s*(gix2 l m + tau*val)
                  | m > l && n == l          = val + s*(gix2 k m - tau*val)
                  | otherwise = val 
 
               where gix2 a b = f $ ix2 a b
                     val      = f sh
                     maxElem  = f kl                     
{-# Inline rotateA #-}                     
                     
rotateP :: Monad m => Matrix -> DIM2 -> Parameters -> m Matrix         
rotateP prr kl@(Z:. k :. l) (Parameters !s !t !tau) = computeUnboxedP $ unsafeTraverse prr id  fun

  where fun f sh@(Z:.x:.y)                      
                 | y == k    = val - s*(gix2 x l + tau*val)
                 | y == l    = val + s*(gix2 x k - tau*val)
                 | otherwise = val
            where gix2 a b = f $ ix2 a b
                  val      = f sh
{-# Inline rotateP #-}                     

calcParameters ::  Double -> Double -> Parameters                 
calcParameters !maxElem !aDiff = Parameters s t tau
  where t = if (abs maxElem < abs aDiff *1.0e-36) then maxElem/aDiff
            else let phi = aDiff/(2*maxElem )
                     var = recip (abs phi + (sqrt $ 1 + phi^2))
                 in if phi < 0 then negate var else var
        c = recip $ sqrt(t^2 + 1)
        s = t*c
        tau = s/(1 + c)
        
-- Return the index of the largest off-diagonal element in the array
maxElemIndex  :: Matrix -> DIM2
maxElemIndex !arr  = R.fromIndex sh $ VU.foldr fun 1 inds

  where inds = VU.enumFromN 0 (dim^2) 
        sh@(Z:. dim :. _dim) = extent arr
        fun n acc= let sh2@(Z:. i:.j) = R.fromIndex sh n
                       sh3 = R.fromIndex sh acc
                   in if i < j then if (abs (arr ! sh2) > abs (arr ! sh3)) then n else acc
                               else acc                

        
-- ==================> Sort EigenValues and EigenVectors <==============
        
sortEigenV :: Monad m => (VecUnbox,Matrix) -> m (VecUnbox,Matrix)
sortEigenV (vals,mtx) = do 
  let (vs,indexes)       =  unzip .  DL.sort $ zip (VU.toList vals) [0..] 
      sortVals           = VU.fromList vs
      fun f sh@(Z:.i:.j) = f $ ix2 i (indexes !! j)   
  sortMtx  <- computeUnboxedP $ unsafeTraverse mtx id fun
  return (sortVals,sortMtx)


testJacobi :: Monad m => m (VecUnbox,Matrix)
testJacobi = do 
     let mtx = R.fromListUnboxed (ix2 7 7) [1.0,0.2367040913545398,0.0,0.0,0.0,5.362538902922486e-2,5.362532497129294e-2,0.2367040913545398,1.0,0.0,0.0,0.0,0.47299800898159106,0.47299767940412246,0.0,0.0,1.0,0.0,0.0,0.0,0.37004397750859236,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.3924913290390063,-0.1308310371205199,5.362538902922486e-2,0.47299800898159106,0.0,0.0,0.3924913290390063,1.0,0.2328803188393675,5.362532497129294e-2,0.47299767940412246,0.37004397750859236,0.0,-0.1308310371205199,0.2328803188393675,1.0]
     jacobiP mtx


-- respuesta hmatrix
-- [1.9105413182235584,1.3514890369682206,1.0875988669892904,0.9999999999999996,0.8898858674147775,0.41563064419280893,0.3448542662113441]
