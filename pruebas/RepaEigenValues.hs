{-# LANGUAGE FlexibleContexts,BangPatterns #-}

module RepaEigenValues
    ( jacobiP
    ) where

import qualified Data.List as DL
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import qualified Data.Vector.Unboxed as VU


-- author:       felipe.zapata@edu.uah.es, a.alvarez@uah.es

-- ===============> JACOBI METHOD <==============================
{- Those interested in the theory behind this method can find a good description
here :
"Numerical Methods in Engineering with Python.Jaan Kiusalaas.   . 2005" -}

{-
Symmetric matrices are represented as a flatten triangular matrix, for instance
the matrix [[1,2,3],[2,4,5],[3,5,6]] is represented as
fromListUnboxed [1,2,3,4,5,6]-}



-- ==============> TYPES <=========================
type Indx = (Int,Int)
type Tolerance = Double
type Step = Int
type EigenValues = [Double]
type EigenVectors = Array U DIM2 Double

data EigenData = EigenData {
             eigenvals :: !EigenValues
           , eigenvec  :: !EigenVectors } deriving (Show)

data Parameters = Parameters {
                   getMaxElem :: !Double
                  ,gett       :: !Double
                  ,gets       :: !Double
                  ,gettau     :: !Double    
                  } deriving Show
                  
-- ============================> PARALLEL JACOBI <========================================

jacobiP ::  (Monad m,VU.Unbox Double) =>
            Array U DIM1 Double  ->
            m EigenData
jacobiP !arr = let (Z:. lenArr) = extent arr
                   dim = dimTriang lenArr
                   tolerance = 1.0e-9
                 in  jacobi arr (identity dim) 0 tolerance

jacobi :: (Monad m, VU.Unbox Double)
          => Array U DIM1 Double
          -> Array U DIM2 Double
          -> Step
          -> Tolerance
          -> m EigenData
jacobi !arrDIM1 !arrDIM2 step tol

  | step > 5*dim*dim = error "Jacobi method did not converge "

  | otherwise = case abs (maxElemArray arrDIM1) > tol of
       True -> do
               arr1 <- rotateA arrDIM1 (k,l) args toval toindx
               arr2 <- rotateR arrDIM2 (k,l) args
               jacobi arr1 arr2 (step+1) tol

       False -> return $
                EigenData (diagonalElems arrDIM1) arrDIM2

 where (Z:. lenArr) = extent arrDIM1
       dim = dimTriang lenArr
       (k,l) = toindx $ maxElemIndex arrDIM1
       args = parameters maxElem aDiff
       aDiff = toval (l,l) - toval (k,k)
       maxElem =  maxElemArray arrDIM1
       toval = toval' arrDIM1
       toindx = map2indx $ DL.zip [0..] [(x,y)| x <-[0..dim-1], y <- [0..dim-1], y >= x]


parameters :: Double -> Double -> Parameters
parameters maxElem aDiff = Parameters maxElem t s tau
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
           Array U DIM1 Double ->
           (Int,Int)  ->
           Parameters   ->
           ((Int,Int) ->
           Double) ->
           (Int -> (Int,Int)) ->
           m(Array U DIM1 Double)
           
rotateA !arr (k,l) (Parameters maxElem t s tau) toval toindx =
  computeUnboxedP $ fromFunction (extent arr)
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
{-# INLINE rotateA #-}

rotateR ::  (Monad m,VU.Unbox Double) =>
            Array U DIM2 Double ->
            (Int,Int) ->
            Parameters  ->
            m(Array U DIM2 Double)
rotateR !arr (k,l) (Parameters maxElem t s tau) =
  computeUnboxedP $ fromFunction (extent arr)
                    $ ( \sh -> let val = arr ! sh
                                   (Z:. n :. m) = sh
                               in fun val (n,m))

  where fun val (n,m) | m == k = val - s*(toval2 (n,l) + tau*val)
                      | m == l = val + s*(toval2 (n,k) - tau*val)
                      | otherwise = val

        toval2 (x,y) = arr ! (Z:. x :. y)
{-# INLINE rotateR #-}

           
-- =========> Functions for creating arrays <===============

-- symmetric matrix represented as a flatten diagonal matrix 
list2ArrDIM1 ::VU.Unbox a => Int -> [a] -> Array U DIM1 a
list2ArrDIM1 dim !list = R.fromListUnboxed (Z:. dim :: DIM1) list

list2ArrDIM2 ::VU.Unbox a => Int -> [a] -> Array U DIM2 a
list2ArrDIM2 dim !list = R.fromListUnboxed (Z:. dim :. dim :: DIM2) list

identity :: VU.Unbox Double  => Int -> Array U DIM2 Double
identity dim =  list2ArrDIM2 dim [dij i j | i <- [1..dim], j <- [1..dim]]
  where dij u v = if u == v then 1.0 else 0.0

       
-- =========> Functions for indexing triangular matrices as flatten arrays <==============

sh2List ::Array U DIM1 Double -> [Int]  
sh2List !arr = [0..(dim-1)]  
  where (Z:.dim) = extent arr  

indx2val :: (Eq a , VU.Unbox Double) => Array U DIM1 Double -> [(a,Int)]-> a -> Double
indx2val !arr pairs indx =  arr ! (Z:. map2indx pairs indx)

map2indx :: Eq a  => [(a,b)]-> a -> b
map2indx pairs indx= case lookup indx pairs of
                          Just p -> p
                          Nothing -> error "Non-existing index"

toval' :: VU.Unbox Double => Array U DIM1 Double -> (Int,Int) -> Double
toval' !arr indx = indx2val arr pairs indx
  where pairs = DL.zip [(x,y)| x <-[0..dim-1], y <- [0..dim-1], y >= x] [0..]
        dim = dimTriang lenArr
        (Z:. lenArr) = extent arr
                          
-- function for calculating the dimension (as nxn square matrix) store in a DIM1 array
dimTriang :: Int -> Int
dimTriang a = (-1+p) `div` 2
  where p =  floor . sqrt . fromIntegral $ 1 + 8*a

-- ==================> Specialised folds <============================

maxElemArray :: Array U DIM1 Double -> Double
maxElemArray !arr  = arr ! (Z:. maxElemIndex arr)

-- Return the index of the largest number in the array
maxElemIndex :: Array U DIM1 Double -> Int
maxElemIndex !arr  = go h t
  where
        t = filterDiagonal arr -- we are interested only in off-Diagonal elements
        h = head t
        go n1 ix@(n2:n3) =case (abs (arr ! (Z:.n1)) > abs (arr ! (Z:.n2))) of
                         True  -> go n1 n3
                         False -> go n2 n3
        go n1 _ = n1

-- ==============> Filter <====================
        
-- Filtering off-diagonal indexes of a one dimensional array
filterDiagonal :: Array U DIM1 Double -> [Int]
filterDiagonal !arr = filter (\x -> not $ elem x $ progression n) allIndx
  where allIndx = sh2List arr
        len = length allIndx
        n = dimTriang len


-- Filtering diagonal elements of a one dimensional array
diagonalElems :: VU.Unbox Double => Array U DIM1 Double -> [Double]
diagonalElems !arr = DL.map (\i -> arr ! (Z:.i)) $ filter (\x -> elem x $ progression n) indexes
  where (Z:.len) = extent arr
        indexes = [0..(len-1)]
        n = dimTriang len
               

progression :: Int -> [Int]
progression n =  tail . DL.map sum . DL.tails. reverse  $ [n,(n-1)..0]
-- ------------------------------------------------------------------------------------------ --

v1:: VU.Unbox Double => Array U DIM1 Double
v1 = list2ArrDIM1 10 [1.0, 2.0, 3.0,4.0, 5.0, 6.0, 7.0, 8.0,9.0,10.0]



       
