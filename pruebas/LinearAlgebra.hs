{-# LANGUAGE FlexibleContexts,BangPatterns #-}

module LinearAlgebra ( 
       EigenValues
      ,EigenVectors
      ,EigenData(..)
      ,Vec(..)
      ,calcDensity
      ,identity
      ,list2ArrDIM1
      ,list2ArrDIM2
      ,list2Diagonal
      ,list2UpperTriang
      ,map2indx
      ,map2val
      ,maxElemIndex
      ,maxElemArray
      ,mmultP
      ,progression
      ,scalarxMtx
      ,toTriang
      ,tr
      ,transpose2P
      ,triang2DIM2
      ,unitaryTransf
      ,vecCross
      ,vecdot
      ,vecnorm
      ,vecscal      
      ,zero
      ) where

-- ########## MODULE FOR IMPLENTING LinearAlgebra Tools #################


import Data.List (lookup, tails, transpose, zip5)
import Control.Parallel.Strategies
import Control.Applicative
import qualified Data.Vector.Unboxed as VU
import qualified Data.Map as M
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import Control.Monad(ap)

-- ========================> TYPES and CLASSES <===================
type Indx = (Int,Int)
type Tolerance = Double
type Step = Int
type EigenValues = [Double]
type EigenVectors = Array U DIM2 Double
data EigenData = EigenData {
             eigenvals :: !EigenValues
           , eigenvec :: !EigenVectors } deriving (Show)

-- class SupplyEigenData where
--   eigvals :: EigenData -> EigenValues
--   eigvec  :: EigeEigenVectors 


-- ===============> DATA TYPES <=================

newtype Vec a = Vec {runVec :: [a]} deriving Show


-- ===============> INSTANCES <==================

instance Functor Vec where
  fmap f (Vec v) = Vec $ f `fmap` v
{-
class (Functor f) => Applicative f where
pure :: a -> f a
(<*>) :: f (a -> b) -> f a -> f b
-}
instance Applicative Vec where
  pure x = Vec $ (repeat x)
  Vec fs <*> Vec xs = Vec $ Prelude.zipWith (\f x -> f x) fs xs

instance Num a => Num (Vec a) where
  (+)         = liftA2 (+)
  (-)         = liftA2 (-)
  (*)         = liftA2 (*)
  abs         = liftA abs
  signum      = liftA signum
  fromInteger = pure . fromInteger

instance (Fractional a) => Fractional (Vec a) where
   (/)  = liftA2 (/)
   recip v = recip <$> v
   fromRational = pure . fromRational

instance (Floating a) => Floating (Vec a) where
  pi = pure pi
  exp v = liftA exp v
  log v = liftA exp v
  sqrt v = liftA sqrt v
  (**)   = liftA2 (**)
  sin v =  sin <$> v
  cos v =  cos <$> v
  tan v =  tan <$> v
  asin v = asin <$> v
  acos v = acos <$> v
  atan v = atan <$> v
  sinh v = sinh <$> v
  cosh v = cosh <$> v
  tanh v = tanh <$> v
  asinh v = asinh <$> v
  acosh v = acosh <$> v
  atanh v = atanh <$> v

-- ===============> Vec functions <===================

vecdot :: Num a => Vec a -> Vec a -> a
v1 `vecdot` v2 =  sum . runVec $ v1 * v2

vecCross :: Num a => Vec a -> Vec a -> Vec a
v1' `vecCross` v2' = let [[x,y,z],[u,v,t]]= fmap runVec [v1',v2']
                     in Vec $ [(y*t-v*z),(v*z-x*t),(x*v-u*y)]

vecscal :: Num a => a -> Vec a -> Vec a
x `vecscal` vec = (pure x) * vec

vecnorm :: Vec Double -> Double
vecnorm v =  sqrt  $  v `vecdot` v

  
-- ================================== > Functions take from Repa-Examples  <==================
-- | Matrix matrix multiply.
mmultP  :: Monad m
        => Array U DIM2 Double
        -> Array U DIM2 Double
        -> m (Array U DIM2 Double)

mmultP arr brr
 = do   trr      <- transpose2P brr
        let (Z :. h1  :. _)  = extent arr
        let (Z :. _   :. w2) = extent brr
        computeP
         $ fromFunction (Z :. h1 :. w2)
         $ \ix   -> R.sumAllS
                  $ R.zipWith (*)
                        (unsafeSlice arr (Any :. (row ix) :. All))
                        (unsafeSlice trr (Any :. (col ix) :. All))
{-# NOINLINE mmultP #-}


-- | Transpose a 2D matrix.
transpose2P
        :: Monad m
        => Array U DIM2 Double
        -> m (Array U DIM2 Double)

transpose2P arr
 = computeUnboxedP
 $ unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# INLINE transpose2P #-}


-- | Take the row number of a rank-2 index.
row :: DIM2 -> Int
row (Z :. r :. _) = r
{-# INLINE row #-}


-- | Take the column number of a rank-2 index.
col :: DIM2 -> Int
col (Z :. _ :. c) = c
{-# INLINE col #-}



-- ===================> <=============================================
scalarxMtx :: (Monad m,VU.Unbox Double) => Array U DIM2 Double -> Double -> m(Array U DIM2 Double)
scalarxMtx arr s = R.computeUnboxedP $ R.unsafeTraverse arr id (\f sh -> s * f sh)
{-# INLINE scalarxMtx #-}


(//) :: (VU.Unbox a,Shape sh) => Array U sh a -> [(sh,a)] -> Array D sh a
(//) !arr us = fromFunction (extent arr) (\sh -> case lookup sh us of
                                                   Just a  -> a
                                                   Nothing -> index arr sh)


unitaryTransf :: (Monad m,VU.Unbox Double) => Array U DIM2 Double -> Array U DIM2 Double -> m(Array U DIM2 Double)
unitaryTransf orthoMtx arr = do
                   transp <- transpose2P orthoMtx
                   prod1 <- mmultP transp arr
                   mmultP prod1 orthoMtx

                                                  
-- ==============> FUNCTIONS FOR GENERATING REPA ARRAYS <=====================

identity :: VU.Unbox Double  => Int -> Array U DIM2 Double
identity dim =  list2ArrDIM2 dim [dij i j | i <- [1..dim], j <- [1..dim]]
  where dij u v = if u == v then 1.0 else 0.0
{-# INLINE identity #-}

zero :: VU.Unbox Double  => Int -> Array U DIM2 Double
zero dim = R.computeUnboxedS $ R.fromFunction (Z:. dim:. dim)
           $ (\sh -> 0.0)
{-# INLINE zero #-}
           
list2ArrDIM2 ::VU.Unbox a => Int -> [a] -> Array U DIM2 a
list2ArrDIM2 dim !list = R.fromListUnboxed (Z:. dim :. dim :: DIM2) list
{-# INLINE list2ArrDIM2 #-}

list2ArrDIM1 ::VU.Unbox a => Int -> [a] -> Array U DIM1 a
list2ArrDIM1 dim !list = R.fromListUnboxed (Z:. dim :: DIM1) list
{-# INLINE list2ArrDIM1 #-}

list2Diagonal :: VU.Unbox Double  => [Double] -> Array U DIM2 Double
list2Diagonal !list = list2ArrDIM2 dim flatMtx
  where dim = length list
        flatMtx = [dij k l | k <- [0..dim-1], l <- [0..dim-1]]
        dij n m = if n == m then list !! n else 0.0
{-# INLINE list2Diagonal #-}
        
list2UpperTriang :: VU.Unbox Double  => [Double] -> Array U DIM1 Double
list2UpperTriang !list = list2ArrDIM1 dim flatMtx
  where dim = length list
        flatMtx = [dij k l | k <- [0..dim-1], l <- [0..dim-1],l>=k]
        dij n m = if n == m then list !! n else 0.0
{-# INLINE list2UpperTriang #-}



-- ============================> Upper Triangular Matrix Functions <===========================================
mmultUpperP  :: Monad m
        => Array U DIM2 Double
        -> Array U DIM2 Double
        -> m (Array U DIM2 Double)

mmultUpperP arr brr
 = do   trr      <- transpose2P brr
        let (Z :. h1  :. _)  = extent arr
        let (Z :. _   :. w2) = extent brr
        computeP
         $ fromFunction (Z :. h1 :. w2)
         $ \ix   ->
              let (Z:. i :. j) = ix
              in case i <= j of
                     True -> R.sumAllS $ R.zipWith (*)
                        (unsafeSlice arr (Any :. (row ix) :. All))
                        (unsafeSlice trr (Any :. (col ix) :. All))
                     False -> 0.0
{-# INLINE mmultUpperP #-}

triang2DIM2 :: Monad m => Array U DIM1 Double -> m (Array U DIM2 Double)
triang2DIM2 arr = R.computeUnboxedP $ R.fromFunction (Z:.dim :. dim)
                 $ \sh@(Z:.i :.j) ->
                      if i <=j then arr ! (Z:.map2indx pairs (i,j)) else arr ! (Z:.map2indx pairs (j,i))
                    
  where (Z:. len) = R.extent arr
        pairs = zip [(x,y) | x <- [0..dim-1], y<-[0..dim-1],x<=y] [0..len-1]
        dim = dimTriang len                   
{-# INLINE triang2DIM2 #-}

toTriang :: Monad m => Array U DIM2 Double -> m (Array U DIM1 Double)
toTriang arr = R.computeUnboxedP $ R.fromFunction (Z:. dim)
              (\(Z:. i) ->
                 let (k,l) = pairs !! i
                 in arr ! (Z:. k:. l) )

  where (Z:. len :. _) = R.extent arr
        dim = (len^2 + len) `div`2
        pairs = [(x,y) | x <- [0..len-1], y <- [0..len-1], x <=y]

-- ==================> Specialised folds <============================

-- Return the index of the largest number in the array
maxElemIndex :: Array U DIM1 Double -> Int
maxElemIndex !arr  = go h t
  where 
        t = filterDiagonal arr -- we are interested only in offDiagonal elements
        h = head t
        go n1 ix@(n2:n3) =case (abs (arr ! (Z:.n1)) > abs (arr ! (Z:.n2))) of
                         True  -> go n1 n3
                         False -> go n2 n3
        go n1 _ = n1
{-# INLINE maxElemIndex #-}

-- ======================> Indexing Related Functions <======================

sh2List ::Array U DIM1 Double -> [Int]
sh2List !arr = [0..(dim-1)]
  where (Z:.dim) = extent arr
{-# INLINE sh2List #-}

maxElemArray :: Array U DIM1 Double -> Double
maxElemArray !arr  = arr ! (Z:. maxElemIndex arr)
{-# INLINE maxElemArray #-}
                
-- =============> Filtering <============================

filterDiagonal :: Array U DIM1 Double -> [Int]
filterDiagonal !arr = filter (\x -> not $ elem x $ progression n) allIndx
  where allIndx = sh2List arr
        len = length allIndx
        n = dimTriang len
{-# INLINE filterDiagonal #-}

tr ::(Monad m,VU.Unbox Double) => Array U DIM2 Double -> m Double
tr arr = (R.foldAllP (+) 0.0) =<< (R.computeUnboxedP $ R.fromFunction (Z:. dim ) (\(Z:.i) -> arr ! (Z:. i :. i)))
  where (Z:. dim :. _) = R.extent arr
         

-- ==========> Lookup Functions <=================================
map2indx :: Eq a  => [(a,b)]-> a -> b
map2indx pairs indx= case lookup indx pairs of
                          Just p -> p
                          Nothing -> error "Boooooooommmmmmmm"

                          
map2val :: Ord k => M.Map k b -> k -> b
map2val mapa key = case M.lookup key mapa of
                        Just val -> val
                        Nothing  -> error "Boooooooommmmmmmm"


-- =========> HARTREE FOCK AUXILIAR FUNCTIONS <===========

-- this function only sum over half the coefficients corresponding with the occupied
-- molecular orbitals in order to  build up the density matrix
calcDensity :: (Monad m,VU.Unbox Double) => Array U DIM2 Double -> m(Array U DIM2 Double)
calcDensity arr
 = (flip scalarxMtx 2.0) =<<
   do   let (Z :. h1  :. _)  = extent arr
        computeP
         $ fromFunction (Z :. h1 :. h1)
         (\ix@(Z:. k :. l)   ->
                  let v1 = R.slice arr (Any :. k :. All)
                      v2 = R.slice arr (Any :. l :. All)
                      p1 = R.extract (Z:.0) (Z:.(h1 `div`2)) v1
                      p2 = R.extract (Z:.0) (Z:.(h1 `div`2)) v2
                  in R.sumAllS $ R.zipWith (*) p1 p2)
{-# INLINE calcDensity #-}
                        
-- ================> Miscellaneous <===============================
progression :: Int -> [Int]
progression n =  tail . Prelude.map sum . tails. reverse  $ [n,(n-1)..0]

                          
dimTriang :: Int -> Int
dimTriang a = (-1+p) `div` 2
  where p =  floor . sqrt . fromIntegral $ 1 + 8*a                          