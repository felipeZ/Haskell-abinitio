{-# LANGUAGE FlexibleContexts,BangPatterns #-}

module JacobiMethod
    ( jacobiP
    ) where

import qualified Data.List as DL
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import LinearAlgebra as LA



-- authors:       felipe.zapata@edu.uah.es, a.alvarez@uah.es

-- ===============> JACOBI METHOD <==============================
{- Those interested in the theory behind this method can find a good description
here :
"Numerical Methods in Engineering with Python.Jaan Kiusalaas.   . 2005" -}




-- ==============> TYPES <=========================
type Indx = (Int,Int)
type Tolerance = Double
type Step = Int

data Parameters = Parameters {
                   getMaxElem :: !Double
                  ,gett       :: !Double
                  ,gets       :: !Double
                  ,gettau     :: !Double
                  ,getK       :: !Int
                  ,getL       :: !Int
                  } deriving Show
                  
-- ============================> PARALLEL JACOBI <========================================

jacobiP ::  (Monad m,VU.Unbox Double) =>
            Array U DIM2 Double  ->
            m LA.EigenData
jacobiP !arr = let (Z:. dim :. _dim) = extent arr
                   tolerance = 1.0e-9
                 in  jacobi arr (LA.identity dim) 0 tolerance

jacobi :: (Monad m, VU.Unbox Double)
          => Array U DIM2 Double
          -> Array U DIM2 Double
          -> Step
          -> Tolerance
          -> m LA.EigenData
jacobi !arrA !arrP step tol

  | step > 5*dim*dim = error "Jacobi method did not converge "

  | otherwise = case abs maxElem > tol of
       True -> do
               arr1 <- rotateA arrA (matrixA arrA args)
               arr2 <- rotateR arrP (matrixR arrP args)
               jacobi arr1 arr2 (step+1) tol

       False -> return $
                LA.EigenData (diagonalElems arrA) arrP

  where (Z:. dim :. _dim) = extent arrA
        sh@(Z:. k :. l) = maxElemIndex arrA
        maxElem = arrA ! sh
        args = parameters maxElem aDiff k l
        aDiff = toval (l,l) - toval (k,k)
        toval (i,j) = arrA ! (Z:. i :. j)

parameters :: Double -> Double ->Int -> Int -> Parameters
parameters !maxElem !aDiff k l = Parameters maxElem t s tau k l

  where t = if (abs maxElem < abs aDiff *1.0e-36)then maxElem/aDiff
            else let phi = aDiff/(2.0 *maxElem )
                     var = recip (abs phi + (sqrt $ 1.0 + phi**2))
                 in case phi < 0.0 of
                         True -> -var
                         False -> var
        c = 1.0/sqrt(t**2 + 1.0)
        s = t*c
        tau = s/(1.0 + c)


-- rotateA make zero the off-Diagonal component of the matrix that is being diagonalized
-- rotateR is the productory of the matrices that diagonalize the symmetric matrix

rotateA :: (Monad m ,VU.Unbox Double) =>
           Array U DIM2 Double ->
           (Int -> Int -> Double) ->
           m(Array U DIM2 Double)           
rotateA !arr fun =
  computeUnboxedP $ fromFunction (extent arr)
                  $ ( \sh@(Z:. n:. m) ->
                    case n <= m of
                         True -> fun n m
                         False -> arr ! sh)

matrixA :: VU.Unbox Double =>
           Array U DIM2 Double ->
           Parameters ->           
           Int -> Int -> Double
           
matrixA !arr (Parameters !maxElem !t !s !tau !k !l) n m
  | (n,m) == (k,l) = 0.0
  | (n,m) == (k,k) = val - t*maxElem
  | (n,m) == (l,l) = val + t*maxElem
  | n < k && m == k = val - s*(toval (n,l) + tau*val)
  | n < k && m == l = val + s*(toval (n,k) - tau*val)
  | k < m && m < l && n == k = val - s*(toval (m,l) + tau*val)
  | k < n && n < l && m == l = val + s*(toval (k,n) - tau*val)
  | m > l && n == k = val - s*(toval (l,m) + tau*val)
  | m > l && n == l = val + s*(toval (k,m) - tau*val)
  | otherwise = val

  where val = toval (n,m)
        toval (i,j) = arr ! (Z :. i:. j)


rotateR :: (Monad m ,VU.Unbox Double) =>
           Array U DIM2 Double ->
           (Int -> Int -> Double) ->
           m(Array U DIM2 Double)
rotateR !arr fun =
  computeUnboxedP $ fromFunction (extent arr)
                  $ ( \sh@(Z:. n:. m) -> fun n m)
        
matrixR :: VU.Unbox Double =>
           Array U DIM2 Double ->
           Parameters ->
           Int -> Int -> Double
matrixR !arr (Parameters !maxElem !t !s !tau !k !l) n m
  | m == k = val - s*((toval (n,l)) + tau*val)
  | m == l = val + s*((toval (n,k)) - tau*val)
  | otherwise = val

  where val = toval (n,m)
        toval (x,y) = arr ! (Z:. x :. y)

           
-- =========> Functions for creating arrays <===============


  
-- ======================> Auxiliar functions <===================
-- Return the index of the largest off-diagonal element in the array

maxElemIndex !arr  = R.fromIndex sh $ VU.foldr fun 1 inds

  where inds = (VU.enumFromN 0 (dim^2) :: VU.Vector Int)
        sh@(Z:. dim :. _dim) = extent arr
        fun n acc= let sh2@(Z:. i:.j) = R.fromIndex sh n
                       sh3 = R.fromIndex sh acc
                   in case i < j of
                           True ->  if (abs (arr ! sh2) > abs (arr ! sh3)) then n else acc
                           False -> acc


-- maxElemIndex !arr =  R.fromIndex sh $ V.foldr1 fun triang
--   where sh@(Z:. dim :. _dim) = extent arr
--         flatten =  R.toUnboxed arr
--         triang = V.generate (dim-1) (\n -> VU.slice (n*(dim+1)+1) (dim-n-1) inds)
--         inds = (VU.enumFromN 0 (dim^2) :: VU.Vector Int)
--         fun x acc = let n = VU.foldr1 fun2 x
--                     in if (abs (flatten VU.! n) > abs (flatten VU.! acc)) then n else acc
--         fun2 y ac = if (abs (flatten VU.! y) > abs (flatten VU.! ac)) then y else ac

                
diagonalElems :: VU.Unbox a => Array U DIM2 a -> VU.Vector a
diagonalElems ! arr = toUnboxed . computeUnboxedS .
    fromFunction (Z :. dim) $
    (\(Z:. i) -> arr ! (Z:. i :. i))

  where (Z:. dim :. _dim) = extent arr

        
-- -------------------------------------------------------------------------------------        
        
v1:: VU.Unbox Double => Array U DIM2 Double
v1 = LA.list2ArrDIM2 3 [80.0,30.0,0.0,30.0,40.0,0.0,0.0,0.0,60.0]



       