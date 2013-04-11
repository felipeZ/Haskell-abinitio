{-# Language FlexibleContexts,BangPatterns #-}

module HartreeFock (
      CGF(..)
     ,Funtype(..)
     ,HFData(..)
     ,scfHF
     ,funlambda
     ) where

import IntegralsEvaluation
import JacobiMethod (jacobiP)
import qualified LinearAlgebra as LA
import qualified Data.Vector.Unboxed as VU
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import qualified Data.List as DL
import qualified Data.Map as M
import Control.Monad (liftM,(<=<))
import GlobalTypes
import Control.Arrow ((&&&))
import Control.Monad.List(guard)

-- ===============> MODULE HARTREE-FOCK <========================================

-- The Matricial Elements are calculated in the IntegralsEvaluation with these elements are build unboxed vector wihch
-- are subsequently transformed to Repa Arrays

-- This is a prototype for calculating the electronig energy of H2, using a STO-3G minimun basis set
-- the basis simulate a slater determinat with parameter xi = 1.24 (for hydrogen this parameter is 1.0)
-- phi_CGF_1s ( xi = 1,24, STO-3G) = 0.444635 phi_GF_1s (0.168856) + 0.535328  phi_GF_1s (0.623913) + 0.154329 phi_GF_1s (3.42525)
-- references in Modern Quantum Chemistry Attila Szabo and Ostlund N. S.
                  

--  ==================>  TYPES  <======================
type Step = Int 

data HFData = HFData {
           getFock      :: !(Array U DIM1 Double)
         , getCoeff     :: !LA.EigenVectors 
         , getDensity   :: !(Array U DIM2 Double)
         , getOrbE      :: !LA.EigenValues
         , getEnergy    :: !Double} deriving (Show)


-- ========================> ORTHOGONALIZATION OF THE BASIS <======================================================

-- Here is diagonalized the Overlap matrix and it is obtained a transformation matrix

canortho :: (Monad m, VU.Unbox Double) => Array U DIM2 Double -> m (Array U DIM2 Double)
canortho arr = do
  eigData <- jacobiP arr
  let eigVal = LA.eigenvals eigData
      eigVecs = LA.eigenvec  eigData
      invSqrt = VU.map (recip . sqrt) eigVal -- For building the S^-1/2 matrix
  R.computeUnboxedP $  R.traverse eigVecs id (\f sh@(Z:._:. k) -> (invSqrt VU.! k) * f sh)
    
symmOrtho :: (Monad m, VU.Unbox Double) => Array U DIM2 Double -> m (Array U DIM2 Double)
symmOrtho arr = do
  eigData <- jacobiP $ arr
  let eigVal = LA.eigenvals eigData
      eigVecs = LA.eigenvec  eigData
      invSqrt = VU.map (recip . sqrt) eigVal -- For building the S^-1/2 matrix
      diag = LA.vec2Diagonal invSqrt
  eigVecTrans <- LA.transpose2P eigVecs
  mtx1 <- LA.mmultP eigVecs diag
  LA.mmultP mtx1 eigVecTrans

 -- ============================================> TWO ELECTRON INTEGRALS <===========================================

-- the Idea behind the evaluation of the four centers integral it is to create a function which
-- sorts the indexes of the all possible ones, for the two electron integrals. This sorting
-- use the fact that if i,j,k,l are indexes for the four centers then the following
-- symmetry property applies : <ij|kl> = <ji|kl> = <kl|ij>...
-- Then, there are only evaluated those integrals which indexes of the basis set, are the smallest one.
-- Besides, there is provided a function for sorting the keys using the symmetry properties of the Integrals.

sortKeys :: [Int] -> [Int]
sortKeys [i,j,k,l] = let l1 =DL.sort [i,j]
                         l2 =DL.sort [k,l]
                     in case l1 <=l2 of
                             True -> l1 DL.++ l2
                             False-> l2 DL.++ l1
                             
                 
                             
calcIntegrals :: [NucCoord]-> [Basis] -> Nelec -> M.Map [Int] Double
-- calcIntegrals :: VU.Unbox Double => [NucCoord]-> [Basis] -> Nelec -> Array U DIM1 Double
calcIntegrals coords basis nelec = M.fromAscList $ DL.zip cartProd integrals
  where integrals = evalIntbykey coords basis cartProd
        dim = nelec - 1
        cartProd = do
          i <- [0..dim]
          j <- [i..dim]
          k <- [i..dim]
          l <- [k..dim]
          let xs = [i,j,k,l]
          guard (condition xs)
          return $ [i,j,k,l]
        condition = \e -> case compare e $ sortKeys e of
                               EQ -> True
                               otherwise -> False

funlambda =  \e -> case compare e $ sortKeys e of
                        EQ -> True
                        otherwise -> False

                        
-- ==================================> CALCULATED THE FOCK MATRIX <======================================
-- Following the textbook of Szabo and Ostlund the uv element is given by the expression
-- Fuv = Hcore_uv + G_uv where G_uv = Suma_l Suma_s Pls * [<uv|sl> - 0.5*<ul|sv>] and Pls is the charge density matrix
          
-- let the matrix [<uv|sl> - 0.5*<ul|sv>] being called mtxBielec

map2Array :: M.Map [Int] Double
             -> (Int,Int)
             -> Nelec
             -> Array D DIM1 Double
map2Array mapIntegrals (i,l) ne =
  R.fromFunction (Z:.ne)
     (\(Z:.indx) ->
       let coulomb  = LA.map2val mapIntegrals $ sortKeys [a,b,indx,l]
           exchange = LA.map2val mapIntegrals $ sortKeys [a,l,indx,b]
       in coulomb - 0.5* exchange)

  where pairs = [(x,y) | x <- [0..ne-1], y <- [0..ne-1], x<=y ]
        (a,b) = pairs !! i

calcGmatrix :: (Monad m, VU.Unbox Double)
               => Array U DIM2 Double
               -> M.Map [Int] Double
               -> m (Array U DIM1 Double)
calcGmatrix !density !integrals =
  computeUnboxedP $ fromFunction (Z:. dim)
                  (\(Z :. i ) ->  sumAllS $
                  fromFunction (Z :. nelec)
                    (\( Z:. l) ->
                    let vec1 = unsafeSlice density (getRow l)
                        vec2 = map2Array integrals (i,l) nelec
                    in sumAllS . R.zipWith (*) vec1 $ vec2 ))

  where getRow x = (Any :. (x :: Int) :. All)
        (Z:. nelec :. _) = extent density
        dim = (nelec^2 + nelec) `div` 2
{-# INLINE calcGmatrix #-}

fock :: (Monad m, VU.Unbox Double) => Array U DIM1 Double-> Array U DIM2 Double-> M.Map [Int] Double-> m (Array U DIM1 Double)
fock core densityMtx integrals = calcGmatrix densityMtx integrals >>= \mtx ->
                                              computeP $ R.zipWith (+) core mtx
{-# INLINE fock #-}

                                              
-- -- =====================================> SCF PROCEDURE <================================================
-- -- Assuming the Zero matrix the initial Guess for Density Matrix
-- -- as it is show when the step is 0
variationalE ::(Monad m, VU.Unbox Double) => Array U DIM1 Double -> Array U DIM1 Double -> Array U DIM2 Double  -> m Double
variationalE core fockMtx oldDensity = (0.5*)`liftM` do
  sumHF <- (R.computeUnboxedP $ R.zipWith (+) core fockMtx) >>= \arr -> LA.triang2DIM2 arr
  result <- LA.mmultP sumHF oldDensity
  LA.tr result

converge :: VU.Unbox Double => Array U DIM2 Double -> Array U DIM2 Double -> Bool
converge oldDensity newDensity  = if sigma < 1.0e-4 then True else False
  where sigma = sqrt . (0.25*) . R.sumAllS . R.computeUnboxedS . R.map (^2) . R.zipWith (-) oldDensity $ newDensity
{-# INLINE converge #-}
  

diagonalHF :: (Monad m, VU.Unbox Double) => Array U DIM1 Double -> Array U DIM2 Double-> m(HFData)
diagonalHF fock1 xmatrix = do
        fDIM2 <-newFock
        eigData <- jacobiP fDIM2
        f' <- LA.toTriang fDIM2
        let (coeff,orbEs) = LA.eigenvec &&& LA.eigenvals $ eigData
        newCoeff <- LA.mmultP xmatrix coeff
        newDensity <- LA.calcDensity newCoeff
        return $ HFData f' newCoeff newDensity orbEs 0.0

  where newFock = (LA.unitaryTransf xmatrix) <=< LA.triang2DIM2 $ fock1

scfHF :: (Monad m, VU.Unbox Double) => [NucCoord] -> [Basis] -> [ZNumber] -> Nelec -> m (HFData)
scfHF coords basis zlist nelec= do 
        let core = hcore coords basis zlist nelec
            density = LA.zero nelec
            integrals = calcIntegrals coords basis nelec           
        xmatrix <- symmOrtho <=< LA.triang2DIM2 $ mtxOverlap coords basis nelec
        scf core density integrals xmatrix 0 500
       
scf :: (Monad m, VU.Unbox Double)
     => Array U DIM1 Double
     -> Array U DIM2 Double
     -> M.Map [Int] Double
     -> Array U DIM2 Double
     -> Step
     -> Int
     -> m(HFData)
scf !core  !oldDensity !integrals !xmatrix  step maxStep

--        | step == 0   =  do    -- core == Fock => guess Pdensity = Zero
--                 hfData <- diagonalHF core xmatrix
--                 scf core (getDensity hfData) integrals xmatrix 1 maxStep
                                                                 
   | step < maxStep = do
       fockDIM1 <- fock core oldDensity integrals
       hfData <- diagonalHF fockDIM1 xmatrix
       etotal <- variationalE core fockDIM1 oldDensity
       let newHFData = hfData {getFock = fockDIM1, getEnergy = etotal}
           bool =  converge oldDensity . getDensity $ newHFData
       case bool of
            True -> return newHFData
            False -> scf core (getDensity newHFData) integrals xmatrix
                     (step+1) maxStep

   | otherwise =  error "SCF maxium steps exceeded"

-- =============================> GET DATA <==============================

--getInitialData :: IO ([NucCoord],[Basis],[ZNumber])
--getInitialData = undefined


prueba ::VU.Unbox Double => Array U DIM2 Double
prueba = LA.list2ArrDIM2 3 [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]

