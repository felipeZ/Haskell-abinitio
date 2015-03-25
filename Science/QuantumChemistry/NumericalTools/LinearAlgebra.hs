{-# LANGUAGE FlexibleContexts,BangPatterns #-}

-- The HaskellFock SCF Project
-- @2013 Felipe Zapata, Angel Alvarez
-- Linear Algebra General Utilities

module  Science.QuantumChemistry.NumericalTools.LinearAlgebra ( 
       EigenValues
      ,EigenVectors
      ,EigenData(..)
      ,calcCoordCGF
      ,calcIndexShell
      ,diagonal
      ,dimTriang
      ,dotMatrix
      ,flattenZero
      ,fourIndex2Flat
      ,identity
      ,indexDIM2toFlat
      ,indexFlat2DIM2
      ,list2ArrDIM1
      ,list2ArrDIM2
      ,map2indx
      ,map2val
      ,mmultDIM2FlattenP
      ,mmultFlattenDIM2P
      ,mmultFlattenP
      ,mmultFlattenS
      ,scalarxMtx
      ,rowFlattenArray
      ,takeZero
      ,toTriang
      ,trace
      ,triang2DIM2
      ,triang2DIM2S
      ,unitaryTransf
      ,vec2Diagonal
      ,vector2Matrix
      ,zero
      ) where



import Data.List (lookup, tails, transpose, zip5)
import Control.Applicative
import Control.Monad(ap, liftM, mplus)
import qualified Data.Vector.Unboxed as VU
import qualified Data.Map as M
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix as R

-- ============> Internal imports <==============
import Science.QuantumChemistry.GlobalTypes

-- ========================> DATA TYPES <===================
type DIM = Int
type Indx = (Int,Int)
type Tolerance = Double
type Step = Int

-- ----------------------------> <---------------------------


-- |The dot product between two matrices
dotMatrix :: Array U DIM2 Double -> Array U DIM2 Double -> Double
dotMatrix mtx1 mtx2 = sumAllS. computeUnboxedS $ fromFunction (Z:.dim^2) $
                      \(Z:. i) ->  let x = v1 ! (Z:.i)
                                       y = v2 ! (Z:.i)
                                   in  x*y

  where (Z:. dim :._) = extent mtx1
        v1 = computeUnboxedS $ reshape (Z:. dim^2) mtx1
        v2 = computeUnboxedS $ reshape (Z:. dim^2) mtx2
{-# INLINE dotMatrix #-}

        
-- ==============> SIMMETRIC MATRIX MULTIPLICATION <===============

-- |Parallel multiplication of two flatten matrices
mmultFlattenP :: Monad m => Array U DIM1 Double -> Array U DIM1 Double -> m (Array U DIM2 Double)
mmultFlattenP !arr !brr =
   do   computeUnboxedP
         $ fromFunction (Z :. dim :. dim)
         $ \(Z:. i:. j) ->  R.sumAllS
                            $ R.zipWith (*)
                              (rowFlattenArray arr i)
                              (rowFlattenArray brr j)

  where (Z:. h1) = extent arr
        dim =  dimTriang h1
{-# NOINLINE mmultFlattenP #-}

-- |Sequential multiplication of two flatten matrices
mmultFlattenS :: Array U DIM1 Double -> Array U DIM1 Double -> Array U DIM2 Double
mmultFlattenS !arr !brr =
   do   computeUnboxedS
         $ fromFunction (Z :. dim :. dim)
         $ \(Z:. i:.j ) ->   R.sumAllS
                             $ R.zipWith (*)
                               (rowFlattenArray arr i)
                               (rowFlattenArray brr j)

  where (Z:. h1) = extent arr
        dim =  dimTriang h1
{-# NOINLINE mmultFlattenS #-}   

-- | Parallel Square-upper traingular matrices multiplication
mmultDIM2FlattenP :: Monad m => Array U DIM2 Double -> Array U DIM1 Double -> m (Array U DIM2 Double)
mmultDIM2FlattenP !arr !flat =
   do   computeUnboxedP
         $ fromFunction (Z :. dim :. dim)
         $ \ix@(Z:. i:.j) -> R.sumAllS
                             $ R.zipWith (*)
                             (unsafeSlice arr (Any :. (row ix) :. All))
                             (rowFlattenArray flat j)

  where (Z:. dim :. _) = extent arr
{-# NOINLINE mmultDIM2FlattenP #-}

-- | Parallel upper traingular-square matrices multiplication
mmultFlattenDIM2P :: Monad m => Array U DIM1 Double -> Array U DIM2 Double -> m (Array U DIM2 Double)
mmultFlattenDIM2P !flat !brr = do
   trr <- transpose2P brr
   computeP
     $ fromFunction (Z :. dim :. dim)
     $ \ix@(Z:. i:.j) -> R.sumAllS
                         $ R.zipWith (*)
                         (rowFlattenArray flat i)
                         (unsafeSlice trr (Any :. (col ix) :. All))

  where (Z:. dim:. _) = extent brr
{-# NOINLINE mmultFlattenDIM2P #-}

-- |function to take the raw of a symmetric matrix  represented as an upper triangular matrix
rowFlattenArray :: Array U DIM1 Double -> Int -> Array U DIM1 Double
rowFlattenArray !arr !k =
  let (Z:. h1) = extent arr
      dim = dimTriang h1
  in computeUnboxedS $ fromFunction (Z:.dim)
    $  \(Z:.i) -> let x = indexDIM2toFlat dim k i  in arr ! (Z:.x)
{-# INLINE rowFlattenArray #-}

-- | Scalar-matrix multiplication
scalarxMtx :: Array U DIM1 Double -> Double -> Array U DIM1 Double
scalarxMtx !arr !s = R.computeUnboxedS $ R.unsafeTraverse arr id (\f sh -> s * f sh)
{-# INLINE scalarxMtx #-}


diagonal :: Monad m => Matrix -> m VecUnbox        
diagonal mtx = liftM toUnboxed  $ computeUnboxedP $ unsafeBackpermute (ix1 dim) (\(Z:.i) -> ix2 i i) mtx

  where (Z:. dim :._) = extent mtx
        
-- ============> Hartree -Fock Linear Algebra <=========
-- | Transform the given matrix into the orthogonal basis
unitaryTransf :: Monad m => TransformMatrix -> Array U DIM1 Double -> m(Array U DIM2 Double)
unitaryTransf !orthoMtx !mtx = do
                   transp <- transpose2P orthoMtx
                   prod1 <- mmultDIM2FlattenP transp mtx
                   mmultP prod1 orthoMtx


                  
-- ==============> FUNCTIONS FOR GENERATING REPA ARRAYS <=====================

-- Identity matrix
identity :: Int -> Matrix
identity dim = R.computeUnboxedS $ R.fromFunction (ix2 dim dim) $
  \(Z:. x:. y) -> if x==y then 1 else 0

-- zero matrix  
zero :: Int -> Array U DIM2 Double
zero dim = R.computeUnboxedS $ R.fromFunction (Z:. dim:. dim)
           $ (\sh -> 0.0)

-- zero matrix represented as a flatten matrix           
flattenZero :: Int -> Array U DIM1 Double
flattenZero n = R.computeUnboxedS $ R.fromFunction (Z:. dim)
                $ (\sh -> 0.0)
  where dim = (n^2 +n) `div`2              

list2ArrDIM2 ::VU.Unbox a => Int -> [a] -> Array U DIM2 a
list2ArrDIM2 dim !list = R.fromListUnboxed (Z:. dim :. dim :: DIM2) list

list2ArrDIM1 ::VU.Unbox a => Int -> [a] -> Array U DIM1 a
list2ArrDIM1 dim !list = R.fromListUnboxed (Z:. dim :: DIM1) list
        
vec2Diagonal :: Array U DIM1 Double -> Array U DIM2 Double
vec2Diagonal vec = computeUnboxedS $ fromFunction (Z:.dim :. dim)
   (\(Z:. i:. j) -> case i == j of
                         True -> vec ! (Z:. i)
                         False -> 0.0)

  where (Z:.dim) = R.extent vec


-- ===================> Triangular Matrices <====================


takeZero :: Int -> [Int] -> [Int]
takeZero !n !xs | null xs   = [] 
                | n == 0    = [0]
                | otherwise = take n xs 
                 
-- | Flatten Upper triangular  to symmetric square matrix
triang2DIM2 :: (VU.Unbox a, Monad m) => Array U DIM1 a -> m (Array U DIM2 a)
triang2DIM2 ar = R.computeUnboxedP $ R.fromFunction (Z:.dim :. dim)
                 $ \sh@(Z:.i :.j) ->
                      let x = indexDIM2toFlat dim i j 
                      in ar ! (Z:.x)
                    
  where (Z:. len) = R.extent ar
        dim = dimTriang $ len
{-# INLINE triang2DIM2 #-}
        

-- | Flatten Upper triangular  to symmetric square matrix
triang2DIM2S :: VU.Unbox a => Array U DIM1 a -> Array U DIM2 a
triang2DIM2S ar = R.computeUnboxedS $ R.fromFunction (Z:.dim :. dim)
                 $ \sh@(Z:.i :.j) ->
                      let x = indexDIM2toFlat dim i j
                      in ar ! (Z:.x)

  where (Z:. len) = R.extent ar
        dim = dimTriang $ len
        
-- | Symmetric square matrix to flatten triangular matrix        
toTriang :: Monad m => Matrix -> m Flatten
toTriang arr = computeUnboxedP $ unsafeBackpermute (ix1 triang) (indexFlat2DIM2 $ dim) arr              
                
  where (Z:.dim:. _) = extent arr                  
        triang = sum [dim,pred dim .. 1]

-- |Trace of a Matrix        
trace :: Array U DIM2 Double -> Double
trace = trace2S

                               
vector2Matrix :: Int -> [a] -> [[a]]
vector2Matrix !n !list = [fmap (\i -> list !! (xs+i) ) [0..pred n] | xs <- xss]
  where xss = [n*k | k <-[0..pred n] ]
-- ============> Functions to manipulate indexes <================

-- | Function to transform the indexes of a flatten upper triangular matrix to its equivalent
--   symmetric square matrix
indexFlat2DIM2 :: Int -> DIM1 -> DIM2 --
indexFlat2DIM2  !dim  (Z:.i) = ix2 x y
  where x = f i 0
        y = g i x
        f !x !n = if x < (dim-n) then n else f (x-(dim-n)) (succ n)
        g !i !x =  i + x - (calcElem x)
        calcElem x = sum . takeZero x $ iterate pred dim

-- | Function to transform the indexes of a symmetric square matrix to its equivalent 
--   flatten upper triangular matrix  
indexDIM2toFlat :: Int -> Int -> Int -> Int 
indexDIM2toFlat !dim !x !y | x <= y = calcElem x + (y-x)
                           | otherwise = calcElem y + (x-y)
  where calcElem w = sum . takeZero w $ iterate pred dim
        

sumIndexTriang:: Int -> Int
sumIndexTriang !n = (n^2 + n) `div` 2
                   
sumFirstIndex :: Int -> Int                   
sumFirstIndex !n  = n + (sumIndexTriang $ pred n) + (sum [sumIndexTriang n - k | k <- [1..pred n]])                  

totalFirstIndex :: Int -> Int -> Int
totalFirstIndex !n !firstIndex = sum [sumFirstIndex (n-m) | m <- [0..pred firstIndex]]

-- | Function to  retrieve the index of the elem n, in the flatten array (DIM1)
--   representation containing the electronic repulsion integral (ERIs) of the form  <ab|cd>
fourIndex2Flat :: Int -> [Int] -> Int
fourIndex2Flat !n xs@[!i,!j,!k,!l] = a0 + (selectPosition n xs)
  where a0 = totalFirstIndex n i
  
selectPosition :: Int -> [Int] -> Int
selectPosition !n [!i,!j,!k,!l] =  
  case test of
     Just 1 -> l - i
     Just 2 -> (n - i) + nfun (ipred n) (ipred k) (ipred l)
     Nothing   -> (n - i) + (sumIndexTriang $ (pred n) - i) + 
                  (sum [sumIndexTriang (n-i) - m | m <- [1..(pred j - i)]])  +
                  (nfun n k l - nfun n i j)
     
  where test = first3equal `mplus` first2equal
        first3equal = if all (i==) [j,k] then Just 1 else Nothing
        first2equal = if  i == j then Just 2 else Nothing
        ipred x = iterate pred x !! (succ i)
        nfun = indexDIM2toFlat  
  

-- ==========> Lookup Functions <=================================
map2indx :: Eq a  => [(a,b)]-> a -> b
map2indx pairs indx= case lookup indx pairs of
                          Just p -> p
                          Nothing -> error "Boooooooommmmmmmm"

                          
map2val :: Ord k => M.Map k b -> k -> b
map2val mapa key = case M.lookup key mapa of
                        Just val -> val
                        Nothing  -> error "Boooooooommmmmmmm"

calcCoordCGF :: [AtomData] -> Int -> (NucCoord,CGF)
calcCoordCGF atoms i = fun 0 atoms
  where fun !acc (x:xs) =
              let nb = length . getBasis $ x
                  newAcc = nb + acc
                  cgf = (getBasis x) !! (i-acc)
                  r   =  getCoord x
              in  if i < newAcc then (r,cgf) else fun newAcc xs
        fun _ [] = error "there is not such atom"        

calcIndexShell :: [AtomData] -> Int -> Int
calcIndexShell atoms i = fun 0 0 atoms
  where fun !acc !indx (x:xs) =
              let nb = length . getBasis $ x
                  newAcc = nb + acc
              in  if i < newAcc then  indx else fun newAcc (succ indx) xs
        fun _ _ [] = error "there is not such atom"


-- ================> Miscellaneous <==============================
                          
dimTriang :: Int -> Int
dimTriang a = (-1+p) `div` 2
  where p =  floor . sqrt . fromIntegral $ 1 + 8*a


     
