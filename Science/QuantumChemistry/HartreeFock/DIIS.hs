{-# Language FlexibleContexts,BangPatterns,ViewPatterns #-}


{-|
Module: Science.QuantumChemistry.HartreeFock.DIIS
Description: DIIS acceleration Convergence
Copyright: @2012,2013 Angel Alvarez Adhesive tape
           @2012,2015 Felipe Zapata core SCF machinery 
           @2016 Felipe Zapata
For a complete description, please refer to: 
Pulay, Péter (1980). "Convergence acceleration of iterative sequences. the case of SCF iteration". 
Chemical Physics Letters 73 (2): 393–398.

DIIS ALGORITHM
(1) Construct the error vector according to eq. and transform it to orthogonal basis. Note that
the density matrix in eq. (4) is the density used to construct the Fock matrix. If the largest
element of the error matrix, emax, is less than a threshold (= O.1 Eh),
initiate the DIIS procedure. If the emax is less than another threshold (= 1.0e-9)the SCF
procedure has converged.
-}


module Science.QuantumChemistry.HartreeFock.DIIS (
              DataDIIS(..)
             ,calcErrorMtx
             ,diis
             ,diisDriver
             ,convergeDIIS             
             )   where

import Control.Monad ((<=<))             
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix as R
--import qualified Data.List as DL
import Data.Sequence as S 
import qualified Data.Vector.Unboxed as U

-- Internal imports
import Science.QuantumChemistry.NumericalTools.GaussElimination (gaussElem)
import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.NumericalTools.LinearAlgebra



{- |The Direct Inversion in the iterative space Method is an improved made by Pulay to accelerate
   the direct SCF procedure -}

-- ====================> TYPES <=============

data DataDIIS = DataDIIS (Seq ErrorMatrix) (Seq FlattenFock) Int deriving Show



diisDriver :: Monad m => ErrorMatrix -> FlattenFock -> DataDIIS -> Step -> Switch -> m (FlattenFock,DataDIIS,Switch)
diisDriver e f dat step switch =
  case switch of
       OFF -> if step > 0 && (maxElem e) < 0.1
                       then let newDIIS = updateDIIS dat e f
                            in return (f,newDIIS,ON)
                       else return (f,dat,OFF)
       ON  ->  do
              (newF,currentDIIS) <- diis e f dat
              let newDIIS = updateDIIS currentDIIS e f
              return (newF,newDIIS,ON)


{-| if the GaussElimination becomes ill-conditioned, omit the
first, second, etc., equation until its condition becomes acceptable. -}
diis :: Monad m => ErrorMatrix -> FlattenFock -> DataDIIS -> m (FlattenFock,DataDIIS)
diis !e !f dat@(DataDIIS es fs n) = do
  let dim1 = S.length $ (f <| fs)
      dim2 = succ dim1
      (Z:. h1) = R.extent f
      ctes = genCtesDIIS dim2
  mtx <- genMtxDIIS dim2 (e <| es)
  r   <- gaussElem mtx ctes
  case  r of
       Nothing -> diis e f $ DataDIIS (initS es) (initS fs) n
       Just cs ->  do
           let 
           r <- computeUnboxedP $ R.fromFunction (Z :. h1 :. dim1 ) $
                    \(Z:.i :. k) -> let arr = (f <| fs) `S.index` k
                                        ck  = cs  ! (Z :. k)
                                        val = arr ! (Z :. i)
                                    in ck* val
           newF <- sumP  r
           return (newF,dat)


genMtxDIIS :: Monad m => Int -> Seq ErrorMatrix -> m (Array U DIM2 Double)
genMtxDIIS !dim2 xs@(viewl -> e :< _) =
   computeUnboxedP  $ R.fromFunction (Z:.dim2 :. dim2) $ \(Z:.m :. n ) -> fill m n

  where (Z:. dim1 :. _) = extent e
        finalIx = pred  dim2
        fill i j | i == finalIx && j == finalIx = 0
                 | i == finalIx || j == finalIx  = -1
                 | otherwise = entry i j
        entry i j = let ei = xs `S.index` i
                        ej = xs `S.index` j
                        et = computeUnboxedS $ transpose ei
                    in trace $ mmultS et ej
{-# INLINE genMtxDIIS #-}


genCtesDIIS :: Int -> Array U DIM1 Double
genCtesDIIS !dim2 = R.fromUnboxed (Z:.dim2) $ U.generate dim2 (\n -> if n==(pred dim2) then (-1) else 0)
{-# INLINE genCtesDIIS #-}

updateDIIS :: DataDIIS -> ErrorMatrix ->  FlattenFock -> DataDIIS
updateDIIS dat@(DataDIIS es fs mx) e f =
 if (S.length (e <| es)) < mx then DataDIIS (e <| es) (f <| fs) mx
                              else DataDIIS (initS $ e <| es) ( initS $ f <| fs) mx

{- | The Error Matrix is calculated as the matrix FDS -SDF
where F, D and S are the Fock, Charge Density and Overlap matrices respectively.
the error Vector is the trace of the multiplication of the error matrix for itself -}
calcErrorMtx :: Monad m =>
                   FlattenFock ->
                   FlattenOverlap ->
                   FlattenChargeDensity ->
                   TransformMatrix ->
                   m ErrorMatrix
calcErrorMtx !f !s !d !xmatrix =  do
  ds <- mmultFlattenP d s
  sd <- mmultFlattenP s d
  fds <- mmultFlattenDIM2P f ds
  sdf <- mmultDIM2FlattenP sd f
  e <- computeUnboxedP $ R.zipWith (-) fds sdf
  t <- transpose2P xmatrix
  mmultP t <=< mmultP e $ xmatrix


convergeDIIS :: ErrorMatrix -> Threshold -> Bool
convergeDIIS !e !diisThreshold = if maxE < diisThreshold then True else False
  where maxE = maxElem e

maxElem :: Array U DIM2 Double -> Double
maxElem  = U.maximum . U.map abs . R.toUnboxed

initS :: Seq a -> Seq a 
initS xs = case viewr xs of
                EmptyR    -> empty
                (hs :> _) -> hs      


