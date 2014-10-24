{-# LANGUAGE BangPatterns #-}

-- The HaskellFock SCF Project
-- @2013 Alessio Valentini, Felipe Zapata
-- 

module EigenHmatrix (
                     computeWithHS
                    ,computeWithHP
                    ) where

import Data.Array.Repa.Repr.Vector
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import Data.Packed.Repa
import qualified Data.Vector.Unboxed as VU
import Numeric.LinearAlgebra (eigSH')

 -- =========> Internal Modules <=======


import GlobalTypes (VecUnbox,Matrix)


-- =========================> <===================================

-- | Resolve the EigenValue problem using Hmatrix and return a Data
--   containing the EigenVectors as a Repa Array and the EigenValues
--   as a Vector Unboxed
computeWithHS :: Array U DIM2 Double -> (VecUnbox,Matrix)
computeWithHS x = let 
        solutions = eigSH' $ repaToMatrix $ computeVectorS $ delay x
        vect      = toUnboxed . computeUnboxedS . delay . vectorToRepa $ fst solutions
        matr      = computeUnboxedS . delay . matrixToRepa $ snd solutions
        in (vect,matr)

computeWithHP :: Monad m => Array U DIM2 Double -> m  (VecUnbox,Matrix) 
computeWithHP x = do
        prova <- computeVectorP $ delay x
        let solutions = eigSH' $ repaToMatrix prova
            vect      = delay . vectorToRepa $ fst solutions
            matr      = delay . matrixToRepa $ snd solutions
        vectM <- computeUnboxedP vect
        let vectF = toUnboxed vectM
        matrM <- computeUnboxedP matr
        return (vectF,matrM)
