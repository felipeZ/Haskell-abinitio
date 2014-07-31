
{-# Language FlexibleContexts,BangPatterns #-}

-- The HaskellFock SCF Project 
-- @2013 Felipe Zapata, Angel Alvarez
-- Hartree-Fock Method
-- The initial starting point for most method is the Hartree-Fock approximation,
-- which established that the ground-state wavefunction   is determined by 
-- products of monoelectronic wavefunctions of the molecular orbitals. 
-- Within this approximation, the electronic field is determined by a self-consistent
-- field method, where every single electron is moving under the mean-potential 
-- created by the other electrons. The wavefunction is built from Slater determinants 
-- as an antisymmetrized product of spin orbitals 

-- For a complete description, please refer to: 
-- Szabo, A. and N. S. Ostlund (1996). Modern Quantum Chemistry. Toronto, Canada, Dover Publications.

-- | Restricted Hartree-Fock Method
module HartreeFock (
      module GlobalTypes
     ,module IntegralsEvaluation
     ,HFData(..)
     ,boysF
     ,calcIntegrals
     ,scfHF
     ,sortKeys
     ) where
     
import Control.Applicative
import Control.Arrow ((&&&))
import Control.Monad.List(guard)
import Control.Monad (liftM,(<=<))
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix 
import qualified Data.List as DL
import qualified Data.Map as M
import qualified Data.Vector.Unboxed as VU
import Text.Printf


-- internal Modules
import BasisOrthogonalization
import Boys(boysF)
import ConcurrencyTools
import DIIS
import GlobalTypes
import IntegralsEvaluation
import Jacobi (jacobiP,testJacobi)
import LinearAlgebra
import Logger


-- ===============> MODULE HARTREE-FOCK <========================================

-- |The Matricial Elements are calculated in the IntegralsEvaluation with these elements are build unboxed vector which
--  are subsequently transformed to Repa Arrays
                  
--  ==================> LOCAL TYPES  <======================



-- =====================================> SCF PROCEDURE <================================================

-- Initialize the required parameters to run the Hartree-Fock procedure
scfHF ::  [AtomData] -> Charge -> (String -> IO ()) -> IO (HFData)
scfHF atoms charge logger = do
        let repulsionN = nuclearRep atoms
            occupied = floor . (/2) . (subtract charge) . sum $  fmap (getZnumber) atoms
            dataDIIS = DataDIIS [] [] 5
        integrals <- calcIntegrals atoms
        core      <- hcore atoms
        s         <- mtxOverlap $ atoms
        xmatrix   <- symmOrtho <=< triang2DIM2 $ s
        density   <- harrisFunctional core xmatrix integrals occupied
        logger $ show v1
        scfDIIS atoms dataDIIS core density s xmatrix integrals repulsionN occupied 0 100 OFF logger

        
-- | Driver to run the DIIS procedure        
scfDIIS :: [AtomData]
     -> DataDIIS
     -> FlattenCore
     -> FlattenChargeDensity
     -> FlattenOverlap
     -> TransformMatrix
     -> Integrals
     -> NuclearRepulsion
     -> OccupiedShells
     -> Step
     -> Int
     -> Switch
     -> (String -> IO())
     -> IO(HFData)
scfDIIS !atoms !dataDIIS !core !oldDensity !overlapMtx !xmatrix !integrals
        !repulsionN !occupied step maxStep switch logger
   | step < maxStep = do
       logger $ printf "Iteration number: %d\n" step
       currentFock <- fock atoms core oldDensity integrals
       errorMtx <- calcErrorMtx currentFock overlapMtx oldDensity xmatrix
       (fockDIIS,newDIIS,newSwitch) <- diisDriver errorMtx currentFock dataDIIS step switch
       (moCoeff, newDensity, orbEnergies) <- diagonalHF fockDIIS xmatrix occupied
       etotal <- fmap (+repulsionN) $ variationalE core fockDIIS newDensity
       let diisThreshold = 1.0e-8
           bool =  if step /= 0 then convergeDIIS errorMtx diisThreshold else False
           newHFData = HFData fockDIIS moCoeff newDensity orbEnergies etotal
       if bool then do
                    logger $ printf "total Energy: %.8f\n" $ getEnergy newHFData
                    return newHFData
               else do
                    logger $ printf "total Energy: %.8f\n" $ getEnergy newHFData
                    scfDIIS atoms newDIIS core newDensity overlapMtx
                            xmatrix integrals repulsionN occupied (step+1) maxStep newSwitch logger 

   | otherwise =  error "SCF maxium steps exceeded"

-- ==========================> <===========================
-- | Initial Density Guess Using the Harris Functional 
harrisFunctional :: Monad m => FlattenCore -> TransformMatrix -> Integrals -> OccupiedShells -> m FlattenChargeDensity
harrisFunctional !fcore !xmatrix !integrals !occupied = do
  let (Z :. dim) = extent fcore
  guess <- computeGuessIntegrals integrals occupied dim
  fockGuess <- computeUnboxedP $ fcore +^ guess
  (_,flattenChargeDensity,_) <- diagonalHF fockGuess xmatrix occupied
  return flattenChargeDensity
{-# INLINE harrisFunctional #-}
  
computeGuessIntegrals :: Monad m =>  Integrals -> OccupiedShells -> Int -> m (Array U DIM1 Double)
computeGuessIntegrals  !integrals !occupied !dim = 
  computeUnboxedP $ fromFunction (Z:.dim) $
     \(Z :. m ) -> let (i,j) = indexFlat2DIM2 dim m
                   in (2*) . sumAllS $ fromFunction (Z :. occupied) $
                        \(Z:.k) -> let flat1 = calcIndex [i,j,k,k]
                                       flat2 = calcIndex [i,k,k,j]
                                       coulomb  = integrals ! (Z:.flat1)
                                       exchange = integrals ! (Z:.flat2)
                                   in coulomb - 0.5* exchange

  where calcIndex = fourIndex2Flat occupied . sortKeys                                   

 -- ============================================> TWO ELECTRON INTEGRALS <===========================================
{-
  | the Idea behind the evaluation of the four centers integral it is to create a function which
   sorts the indexes of the all possible ones, for the two electron integrals. This sorting
   use the fact that if i,j,k,l are indexes for the four centers then the following
   symmetry property applies : <ij|kl> = <ji|kl> = <kl|ij>...
   Then, there are only evaluated those integrals which indexes of the basis set, are the smallest one.
   Besides, it is provided a function for sorting the keys using the symmetry properties of the Integrals.-}

sortKeys :: [Int] -> [Int]
sortKeys [i,j,k,l] = let l1 =DL.sort [i,j]
                         l2 =DL.sort [k,l]
                     in if l1 <=l2 then l1 DL.++ l2 else l2 DL.++ l1
                                              
-- | Compute the electronic integrals                                 
calcIntegrals :: Monad m => [AtomData] -> m (Array U DIM1 Double)
calcIntegrals !atoms = evalIntbykeyStrat atoms cartProd
  where dim = pred . sum . fmap (length . getBasis) $ atoms
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

-- ============> Density Matrix <=================

{- | this function only sum over half the coefficients corresponding with the occupied
     molecular orbitals in order to  build up the density matrix
     Notice that we are assuming a double ocupation in each Shell, therefore
    this matrix differs from the monoelectronic density matrix -}                   
calcDensity :: Monad m => Matrix-> OccupiedShells -> m FlattenChargeDensity
calcDensity !arr !occupied =
   do   let (Z :. dim  :. _)  = extent arr
            flat = (dim^2 + dim)`div`2
        computeP
         $ fromFunction (Z :. flat)
         (\(Z:. x)  ->
                  let (k,l) =  indexFlat2DIM2 dim x
                      v1 = R.slice arr (Any :. k :. All)
                      v2 = R.slice arr (Any :. l :. All)
                      p1 = R.extract (Z:.0) (Z :. occupied) v1
                      p2 = R.extract (Z:.0) (Z :. occupied) v2
                  in (*2) . R.sumAllS $ R.zipWith (*) p1 p2)
{-# INLINE calcDensity #-}
                  
-- ==================================> CALCULATED THE FOCK MATRIX <======================================          
-- | Calculate the current fock matrix as a flatten matrix summing the core Hamiltonian plus the electronic
--   interaction integrals
fock :: Monad m
        => [AtomData]
        -> FlattenCore
        -> FlattenChargeDensity
        -> Integrals
        -> m (Array U DIM1 Double)
fock !atoms !core !densityMtx !integrals =
         calcGmatrix atoms densityMtx integrals >>= \mtx ->
         computeP $ core +^ mtx
{-# INLINE fock #-}

-- |Calculate the average electronic interaction 
calcGmatrix :: Monad m
               => [AtomData]
               -> FlattenChargeDensity
               -> Integrals
               -> m (Array U DIM1 Double)
calcGmatrix !atoms !density !integrals =
  computeUnboxedP $ fromFunction (Z:. flat) $
                  \(Z :. i ) ->  sumAllS $
                  fromFunction (Z :. nshells) $
                    \( Z:. l) ->
                    let vec1 = rowFlattenArray density l
                        vec2 = sliceIntegrals atoms integrals (i,l) nshells
                    in sumAllS $ vec1 *^ vec2 

  where (Z:. flat ) = extent density
        nshells = dimTriang flat
{-# INLINE calcGmatrix #-}
        
-- | Auxiliar function to slice the integrals vector
sliceIntegrals :: [AtomData]
             -> Array U DIM1 Double
             -> (Int,Int)
             -> Nelec
             -> Array D DIM1 Double
sliceIntegrals !atoms !integrals (!i,!l) !nshells = 
  R.fromFunction (Z:.nshells)
     (\(Z:. k) ->
       let flat1 = calcIndex [a,b,k,l]
           flat2 = calcIndex [a,l,k,b]
           coulomb  = integrals ! (Z:.flat1)
           exchange = integrals ! (Z:.flat2)
       in coulomb - 0.5* exchange)

  where (a,b) = indexFlat2DIM2 nshells i
        calcIndex = fourIndex2Flat nshells . sortKeys
{-# INLINE sliceIntegrals #-}        

-- |Diagonalize the Fock matrix
diagonalHF :: Monad m => FlattenFock ->
                         TransformMatrix->
                         OccupiedShells ->
                         m (MOCoefficients,FlattenChargeDensity,EigenValues)
diagonalHF !fock1 !xmatrix !occupied = do
        fDIM2 <- unitaryTransf xmatrix fock1
        (orbEs,coeff) <- jacobiP fDIM2
        newCoeff <- mmultP xmatrix coeff
        newDensity <- calcDensity newCoeff occupied
        return $ (newCoeff, newDensity, orbEs)
                           
                  
-- | Har tree-Fock total energy
variationalE ::Monad m => FlattenCore -> FlattenFock -> FlattenChargeDensity  -> m Double
variationalE !core !fockMtx !oldDensity = do
  sumHF  <- computeUnboxedP $ core +^ fockMtx
  result <- mmultFlattenP oldDensity sumHF
  return $ 0.5 * (trace result)
{-# INLINE variationalE #-}

-- | Function to check convergency
converge :: FlattenChargeDensity -> FlattenChargeDensity -> Bool
converge !oldDensity !newDensity  = if sigma < 1.0e-6 then True else False
  where sigma = sqrt . (0.25*) . sumAllS . computeUnboxedS . R.map (^2) $ oldDensity -^ newDensity
{-# INLINE converge #-}
  

-- =================> Nuclear Repulsion Energy <=====================

-- | Calculate internuclear electronic repulsion
nuclearRep :: [AtomData] -> Double
nuclearRep xs = sum [repulsion atomi atomj | atomi <- xs, atomj <- xs, atomj > atomi]
  where dim = pred . length $ xs
        repulsion at1 at2 = let ([za,zb],[ra,rb]) = fmap getZnumber &&& (fmap getCoord) $ [at1,at2]
                                rab = sqrt . sum . fmap (^2) $ DL.zipWith (-) ra rb
                            in za*zb/rab



