{-# Language BangPatterns,FlexibleContexts,TypeFamilies #-}

module Science.QuantumChemistry.NumericalTools.GaussElimination {-(
--                          gaussElem
                        )-}  where

import Control.Arrow((&&&))
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U

--  --------------> Internal Modules <---------
import qualified Science.QuantumChemistry.NumericalTools.LinearAlgebra as LA

    
-- ================> Types <=====================


type Augmented    = Array U DIM2 Double
type Coefficients = Array U DIM2 Double
type Constants    = Array U DIM1 Double
type Triangular   = V.Vector (U.Vector Double)

data GaussElem = GaussElem { getCoeff :: !Coefficients
                            ,getConst :: !Constants
                            } deriving Show
                            
-- ===========================> GAUSS ELIMINATION <====================

gaussElem :: Monad m => Coefficients -> Constants -> m (Maybe (Array U DIM1 Double))
gaussElem mtx vec = eliminationPhase mtx vec (Z:.0:.0) >>= \gauss ->
                    case gauss of
                         Nothing    -> return Nothing
                         Just gauss -> return . Just $ substitution gauss
                    

{-gaussElem :: Monad m => Coefficients -> Constants -> m (Maybe (Array U DIM1 Double))
gaussElem mtx vec = eliminationPhase mtx vec (Z:.0:.0) >>= \gauss ->
                    return $ substitution gauss   -}                 
-- ============> Elimination Phase <==============

eliminationCoeff :: Monad m => Coefficients -> DIM2 -> m Coefficients
eliminationCoeff !mtx pivot@(Z:. k :. _k) = computeUnboxedP . fromFunction (extent mtx) $
  (\sh@(Z:.i:.j) ->
     if i > k && j >= k
        then let alpha = (mtx ! (Z:. i:. k)) / (mtx ! pivot)
             in (mtx ! sh) - alpha*(mtx ! (Z:. k :. j))
     else mtx ! sh)

eliminationConst :: Monad m => Coefficients -> Constants -> DIM2 -> m Constants
eliminationConst !mtx !vec pivot@(Z:. k :. _k) =
  computeUnboxedP . fromFunction (extent vec) $
  (\sh@(Z:. i) ->
          if (i<=k) then vec ! sh
          else
             let alpha = (mtx ! (Z:. i:. k)) / (mtx ! pivot)
             in (vec ! sh) - alpha*(vec ! (Z:. k)) )

eliminationPhase :: Monad m => Coefficients -> Constants -> DIM2 -> m (Maybe GaussElem)
eliminationPhase !mtx !vec sh@(Z:. n:. _n)
   | (extent mtx) > sh  = eliminationCoeff mtx sh >>= \mtx2 ->
                          eliminationConst mtx vec sh >>= \vec2 ->
                          eliminationPhase mtx2 vec2 (Z:. (n+1):. (n+1))
                          
 | otherwise = return $ if (not $ linearIndependence mtx) || anyInfinite
                           then Nothing 
                           else Just $ GaussElem mtx vec
                                                  

  where (Z:. dim:. _dim) = extent mtx
        anyInfinite = U.any ( \x -> isNaN x || isInfinite x) $ R.toUnboxed mtx

linearIndependence :: Coefficients -> Bool
linearIndependence mtx = U.foldl1' (&&) . U.map genVec $ U.generate dim id
  where genVec k = U.any (\x -> abs x > 1e-6) . toUnboxed . computeUnboxedS $ slice mtx (Any :. k :. All)
        (Z:.dim:._) = R.extent mtx
 

-- ==============> Substitution Phase <===============

substitution ::  GaussElem -> Array U DIM1 Double
substitution !ge = R.fromUnboxed (Z :. dim) $
                   V.foldr (backTracking cs) U.empty tvs
  where (mtx,cs) = getCoeff &&& getConst $ ge
        tvs = toTriang (R.toUnboxed mtx) dim
        (Z:.dim) = extent cs

        
backTracking :: Array U DIM1 Double -> U.Vector Double -> U.Vector Double -> U.Vector Double
backTracking cs v acc =
            let b = cs ! (Z:.k)
                k = n-m
                (Z:.n) = extent cs
                m = U.length v
                akk = U.head v
                sumAx = U.sum $ U.zipWith (*) acc (U.tail v)
                val = (b - sumAx) * recip akk
            in val `U.cons` acc

toTriang :: U.Vector Double -> Int -> V.Vector (U.Vector Double)
toTriang !vecU !dim = V.generate dim (\n -> U.slice (n*(dim+1)) (dim-n) vecU)

q1 :: U.Unbox Double => Array U DIM2 Double
q1 =   R.fromListUnboxed (Z:. 3 :. 3 :: DIM2) [4.0, -2.0, 1.0, -2.0, 4.0, -2.0, 1.0, -2.0, 4.0]

v1 :: U.Unbox Double => Array U DIM1 Double
v1 =   R.fromListUnboxed (Z:. 3 :: DIM1) [11.0,-16.0,17.0]

  
