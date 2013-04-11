{-# LANGUAGE FlexibleContexts,BangPatterns #-}

module JacobiMethod
    ( jacobiP
    ) where

import Data.List
import Control.Monad(liftM)
import qualified Data.Map as Map
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import qualified Data.Vector.Unboxed as VU
import LinearAlgebra (EigenValues,EigenVectors,EigenData(..), identity, list2ArrDIM1, list2ArrDIM2,map2indx,maxElemIndex,maxElemArray,progression)

-- ===============> JACOBI METHOD <=========================================
{- Those interested in the theory behind this method can find a good description
here :
"Numerical Methods in Engineering with Python.Jaan Kiusalaas.   . 2005" -}


-- ==============> TYPES <=========================
type Indx = (Int,Int)
type Tolerance = Double
type Step = Int

       
-- ===> Indexes Related Functions < =====

sh2List ::Array U DIM1 Double -> [Int]  
sh2List !arr = [0..(dim-1)]  
  where (Z:.dim) = extent arr  
{-# INLINE sh2List #-}

indx2val :: (Eq a , VU.Unbox Double) => Array U DIM1 Double -> [(a,Int)]-> a -> Double
indx2val !arr pairs indx =  arr ! (Z:. map2indx pairs indx)
{-# INLINE indx2val #-}

-- function for calculating the dimension (as nxn square matrix) store in a DIM1 array
dimTriang :: Int -> Int
dimTriang a = (-1+p) `div` 2
  where p =  floor . sqrt . fromIntegral $ 1 + 8*a
-- ==> Fetching elements related functions <==


toval' :: VU.Unbox Double => Array U DIM1 Double -> (Int,Int) -> Double
toval' !arr indx = indx2val arr pairs2 indx
  where pairs2 = Data.List.zip [(x,y)| x <-[0..dim-1], y <- [0..dim-1], y >= x] [0..]
        dim = dimTriang lenArr
        (Z:. lenArr) = extent arr   
{-# INLINE toval' #-}


-- ==> Filtering Related Functions <===========

diagonalElems :: VU.Unbox Double => Array U DIM1 Double -> [Double]
diagonalElems !arr = Prelude.map (\i -> arr ! (Z:.i)) $ filter (\x -> elem x $ progression n) indexes
  where (Z:.len) = extent arr
        indexes = [0..(len-1)]
        n = dimTriang len
{-# INLINE diagonalElems #-}        
               

-- Filtering off diagonal elements of a one dimensional array
filterDiagonal :: Array U DIM1 Double -> [Int]  
filterDiagonal !arr = filter (\x -> not $ elem x $ progression n) allIndx  
  where allIndx = sh2List arr  
        len = length allIndx  
        n = dimTriang len
{-# INLINE filterDiagonal #-}
 
-- ------------------------------------------------------------------------------------------ --
        
parameters :: Double -> Double -> [Double]
parameters maxElem aDiff = [maxElem,t,s,tau]
  where t = if (abs maxElem < abs aDiff *1.0e-36)then maxElem/aDiff
            else let phi = aDiff/(2.0 *maxElem )
                     var = recip (abs phi + (sqrt $ 1.0 + phi**2))
                 in case phi < 0.0 of
                         True -> -var
                         False -> var
        c = 1.0/sqrt(t**2 + 1.0)
        s = t*c
        tau = s/(1.0 + c)


-- rotate make zero the off-Diagonal component of the matrix that is being diagonalized
-- rotateR is the productory of the matrix wich contains the operator rotating the matrix being diagonalized 
   
rotate :: (Monad m ,VU.Unbox Double) => Array U DIM1 Double -> (Int,Int) -> [Double] -> ((Int,Int) -> Double) -> (Int -> (Int,Int))->  m(Array U DIM1 Double)
rotate !arr (k,l) [maxElem,t,s,tau] toval toindx = computeUnboxedP $ fromFunction (extent arr)
                   $ ( \sh ->  let val = arr ! sh
                                   (Z:. indx) = sh
                                   (n,m) = toindx indx
                               in fun val (n,m) )

  where fun val (n,m)| (n,m) == (k,l) = 0.0
                     | (n,m) == (k,k) = val - t*maxElem
                     | (n,m) == (l,l) = val + t*maxElem
                     | n < k && m == k = val - s*(toval (n,l) + tau*val)
                     | n < k && m == l = val + s*(toval (n,k) - tau*val)
                     | k < m && m < l && n == k = val - s*(toval (m,l) + tau*val)
                     | k < n && n < l && m == l = val + s*(toval (k,n) - tau*val)
                     | m > l && n == k = val - s*(toval (l,m) + tau*val)
                     | m > l && n == l = val + s*(toval (k,m) - tau*val)
                     | otherwise = val 
{-# INLINE rotate #-}
        
rotateR ::  (Monad m,VU.Unbox Double) => Array U DIM2 Double -> (Int,Int)-> [Double] -> m(Array U DIM2 Double)
rotateR !arr (k,l) [maxElem,t,s,tau] = computeUnboxedP $ fromFunction (extent arr)
                    $ ( \sh -> let val = arr ! sh
                                   (Z:. n :. m) = sh
                               in fun val (n,m))
                               
  where fun val (n,m) | m == k = val - s*(toval2 (n,l) + tau*val)
                      | m == l = val + s*(toval2 (n,k) - tau*val)
                      | otherwise = val
                      
        toval2 (x,y) = arr ! (Z:. x :. y)
{-# INLINE rotateR #-}


-- ============================> PARALLEL JACOBI <========================================
 
jacobiP ::  (Monad m,VU.Unbox Double) => Array U DIM1 Double  ->  m EigenData
jacobiP !arr = let (Z:. lenArr) = extent arr
                   dim = dimTriang lenArr
                   tolerance = 1.0e-9
                 in  jacobi' arr (identity dim) 0 tolerance
                 
jacobi' ::  (Monad m,VU.Unbox Double) => Array U DIM1 Double -> Array U DIM2 Double -> Step -> Tolerance -> m EigenData
jacobi' !arrDIM1 !arrDIM2 step tol 
   
    | step > 5*dim*dim = error "Jacobi Method exceeded the default number of steps "

    | step == 0 = do                                       
        arr1 <- rotate arrDIM1 (k,l) pars toval toindx
        arr2 <- rotateR arrDIM2 (k,l) pars
        jacobi' arr1 arr2 1 tol
                                       
                                       
    | otherwise = case abs (maxElemArray arrDIM1) > tol of
                     True -> do
                          arr1 <- rotate arrDIM1 (k,l) pars toval toindx
                          arr2 <- rotateR arrDIM2 (k,l) pars
                          jacobi' arr1 arr2 1 tol

                     False -> return $ EigenData (diagonalElems arrDIM1) arrDIM2
 where (Z:. lenArr) = extent arrDIM1
       dim = dimTriang lenArr
       (k,l) = toindx $ maxElemIndex arrDIM1
       toindx = map2indx $ Prelude.zip [0..] [(x,y)| x <-[0..dim-1], y <- [0..dim-1], y >= x]

       toval = toval' arrDIM1
       aDiff = toval (l,l) - toval (k,k)
       maxElem =  maxElemArray arrDIM1
       pars = parameters maxElem aDiff

       
-- orderEig :: EigenData -> EigenData
-- orderEig eigData = let eigvals = eigenval eigData
--                        eigvec  = eigenvec eigData
--                    in case check of
--                            True -> eigdata
--                            False -> 
--                            
--   where  check = all (\i -> (eigenval !! (succ i)) > (eigenval !! i)) [0..len-2]
--          len = length eigvals 
-- ===========================> SEQUENCIAL JACOBI <================================

-- rotate ::  VU.Unbox Double => Array U DIM1 Double -> (Int,Int) -> [Double] -> ((Int,Int) -> Double) -> (Int -> (Int,Int))->  Array U DIM1 Double
-- rotate arr (k,l) [maxElem,t,s,tau] toval toindx = computeUnboxedS $ fromFunction (extent arr)
--                    $ ( \sh ->  let val = arr ! sh
--                                    (Z:. indx) = sh
--                                    (n,m) = toindx indx
--                                in fun val (n,m) )
-- 
--   where fun val (n,m)| (n,m) == (k,l) = 0.0
--                      | (n,m) == (k,k) = val - t*maxElem
--                      | (n,m) == (l,l) = val + t*maxElem
--                      | n < k && m == k = val - s*(toval (n,l) + tau*val)
--                      | n < k && m == l = val + s*(toval (n,k) - tau*val)
--                      | k < m && m < l && n == k = val - s*(toval (m,l) + tau*val)
--                      | k < n && n < l && m == l = val + s*(toval (k,n) - tau*val)
--                      | m > l && n == k = val - s*(toval (l,m) + tau*val)
--                      | m > l && n == l = val + s*(toval (k,m) - tau*val)
--                      | otherwise = val
-- {-# INLINE rotate #-}

-- rotateR ::  VU.Unbox Double => Array U DIM2 Double -> (Int,Int)-> [Double] -> Array U DIM2 Double
-- rotateR arr (k,l) [maxElem,t,s,tau] = computeUnboxedS $ fromFunction (extent arr)
--                     $ ( \sh -> let val = arr ! sh
--                                    (Z:. n :. m) = sh
--                                in fun val (n,m))
-- 
--   where fun val (n,m) | m == k = val - s*(toval2 (n,l) + tau*val)
--                       | m == l = val + s*(toval2 (n,k) - tau*val)
--                       | otherwise = val
-- 
--         toval2 (x,y) = arr ! (Z:. x :. y)
-- {-# INLINE rotateR #-}

-- jacobi ::  VU.Unbox Double => Array U DIM1 Double -> Tolerance -> Either String (EigenValues, EigenVectors)
-- jacobi arr tol = let (Z:. lenArr) = extent arr
--                      dim = dimTriang lenArr
--                  in  jacobi' arr (identity dim) 0 tol
-- 
-- jacobi' :: VU.Unbox Double => Array U DIM1 Double -> Array U DIM2 Double -> Step -> Tolerance -> Either String (EigenValues, EigenVectors)
-- jacobi' arrDIM1 arrDIM2 step tol | step > 5*dim*dim = Left "Jacobi Method exceeded the default number of steps "
-- 
--                                  | step == 0 = jacobi' (rotate arrDIM1 (k,l) pars toval toindx) (rotateR arrDIM2 (k,l) pars) 1 tol
-- 
--                                  | otherwise = case abs (maxElemArray arrDIM1) > tol of
--                                                     True  -> jacobi' (rotate arrDIM1 (k,l) pars toval toindx )(rotateR arrDIM2 (k,l) pars ) (step+1) tol
--                                                     False -> Right ((diagonalElems arrDIM1),arrDIM2)
-- 
--  where (Z:. lenArr) = extent arrDIM1
--        dim = dimTriang lenArr
--        (k,l) = toindx $ maxElemIndex arrDIM1
--        toindx = map2indx $ Prelude.zip [0..] [(x,y)| x <-[0..dim-1], y <- [0..dim-1], y >= x]
-- 
--        toval = toval' arrDIM1
--        aDiff = toval (l,l) - toval (k,k)
--        maxElem =  maxElemArray arrDIM1
--        pars = parameters maxElem aDiff 

-- =========================================================================================================
       
v1:: VU.Unbox Double => Array U DIM1 Double
v1 = list2ArrDIM1 9 [80.0,30.0,0.0,30.0,40.0,0.0,0.0,0.0,60.0]


        