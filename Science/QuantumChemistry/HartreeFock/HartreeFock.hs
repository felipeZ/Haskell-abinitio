
{-# Language FlexibleContexts,BangPatterns #-}

{-|
Module: Science.QuantumChemistry.HartreeFock.HartreeFock
Description: Restricted Hartree-Fock Method
Copyright: @2016 Felipe Zapata

The HaskellFock SCF Project 
The initial starting point for most method is the Hartree-Fock approximation,
which established that the ground-state wavefunction   is determined by 
products of monoelectronic wavefunctions of the molecular orbitals. 
Within this approximation, the electronic field is determined by a self-consistent
field method, where every single electron is moving under the mean-potential 
created by the other electrons. The wavefunction is built from Slater determinants 
as an antisymmetrized product of spin orbitals 

For a complete description, please refer to: 
Szabo, A. and N. S. Ostlund (1996). Modern Quantum Chemistry. Toronto, Canada, Dover Publications.
-}

module Science.QuantumChemistry.HartreeFock.HartreeFock (
      HFData(..)
     ,boysF
     ,calcIntegrals
     ,scfHF
     ,sortKeys
     ) where
     

import Control.Arrow ((&&&))
import Control.Monad.List(guard)
import Control.Monad (liftM,(<=<))
import Control.Parallel.Strategies 
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix as R
import Data.Foldable
import Data.List as L hiding (sum)
import qualified Data.Map.Lazy   as ML
import qualified Data.Map.Strict as M
import qualified Data.Sequence as S
import qualified Data.Vector.Unboxed as VU
import Prelude hiding (sum)
import Text.Printf



-- internal Modules
import Science.QuantumChemistry.ConcurrencyTools.Logger
import Science.QuantumChemistry.HartreeFock.BasisOrthogonalization (symmOrtho)
import Science.QuantumChemistry.HartreeFock.DIIS                   (DataDIIS(..)
                                                                   ,calcErrorMtx
                                                                   ,diis
                                                                   ,diisDriver
                                                                   ,convergeDIIS)
import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.Integrals.IntegralsEvaluation 
import Science.QuantumChemistry.NumericalTools.Boys        (boysF)
import Science.QuantumChemistry.NumericalTools.EigenValues (eigenSolve)
import Science.QuantumChemistry.NumericalTools.LinearAlgebra
import Science.QuantumChemistry.NumericalTools.TableBoys   (Boys,generateGridBoys)



-- ===============> MODULE HARTREE-FOCK <========================================

                  
-- =====================================> SCF PROCEDURE <================================================

-- | Initialize the required parameters to run the Hartree-Fock procedure
scfHF ::  [AtomData] -> Charge -> (String -> IO ()) -> IO HFData
scfHF atoms charge logger = do
            -- internuclear electronic repulsion
        let repulsionN = nuclearRep atoms
             -- Occupied Molecular Orbitals
            occupied   = floor . (/2) . subtract charge . sum 
                                $ fmap getZnumber atoms
            -- DIIS configuration
            dataDIIS   = DataDIIS S.empty S.empty 5
            -- gridBoys dx, where mMax is the maximum order of the boys
            -- function and dx the grid delta
            gridBoys   = generateGridBoys 0.1 
            integrals  = calcIntegrals gridBoys atoms
        core      <- hcore gridBoys atoms
        s         <- mtxOverlap atoms
        xmatrix   <- symmOrtho <=< triang2DIM2 $ s
        density   <- harrisFunctional core xmatrix integrals occupied
        scfDIIS atoms dataDIIS core density s xmatrix integrals repulsionN occupied 0 30 OFF logger

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
       errorMtx    <- calcErrorMtx currentFock overlapMtx oldDensity xmatrix
       (fockDIIS,newDIIS,newSwitch)       <- diisDriver errorMtx currentFock dataDIIS step switch
       (moCoeff, newDensity, orbEnergies) <- diagonalHF fockDIIS xmatrix occupied
       let etotal =  (+repulsionN) $ variationalE core fockDIIS newDensity
           diisThreshold = 1.0e-7
           bool =  (step /= 0) && convergeDIIS errorMtx diisThreshold
           newHFData = HFData fockDIIS moCoeff newDensity orbEnergies etotal
       logger $ printf "New Density:\n%s\n\n" $ show newDensity
       logger $ printf "FockDIIS: \n%s\n\n" $ show fockDIIS
       logger $ printf "total Energy: %.8f\n"  etotal
       if bool then  return newHFData
               else scfDIIS atoms newDIIS core newDensity overlapMtx
                            xmatrix integrals repulsionN occupied (step+1) maxStep newSwitch logger 

  | otherwise =  error "SCF maxium steps exceeded"


-- ==========================> <===========================
-- | Initial Density Guess Using the Harris Functional 
harrisFunctional :: Monad m => FlattenCore -> TransformMatrix -> Integrals -> OccupiedShells -> m FlattenChargeDensity
harrisFunctional !fcore !xmatrix !integrals !occupied = do
  let (Z :. dim) = extent fcore
  guess     <- computeGuessIntegrals integrals occupied dim
  fockGuess <- computeUnboxedP $ fcore +^ guess
  (_,flattenChargeDensity,_) <- diagonalHF fockGuess xmatrix occupied
  return flattenChargeDensity
{-# INLINE harrisFunctional #-}

-- | The force center integral are calculated in this step and subsequently transformed to a Repa Array
computeGuessIntegrals :: Monad m =>  Integrals -> OccupiedShells -> Int -> m (Array U DIM1 Double)
computeGuessIntegrals  !integrals !occupied !dim = 
  computeUnboxedP $ fromFunction (Z:.dim) $
     \idx -> let (Z:.i:.j) = indexFlat2DIM2 dim idx
             in (2*) . sumAllS $ fromFunction (Z :. occupied) $
                  \(Z:.k) -> let flat1 = calcIndex [i,j,k,k]
                                 flat2 = calcIndex [i,k,k,j]
                                 coulomb  = integrals ! (Z:.flat1)
                                 exchange = integrals ! (Z:.flat2)
                             in coulomb - 0.5* exchange

  where calcIndex = fourIndex2Flat occupied . sortKeys                                   

 -- ============================================> TWO ELECTRON INTEGRALS <===========================================

sortKeys :: [Int] -> [Int]
sortKeys [i,j,k,l] = let l1 = L.sort [i,j]
                         l2 = L.sort [k,l]
                     in if l1 <= l2 then l1 L.++ l2 else l2 L.++ l1
                                              
-- | Compute the electronic integrals                                 
{-| the Idea behind the evaluation of the four centers integral it is to create a function which
 sorts the indexes of the all possible ones, for the two electron integrals. This sorting
 use the fact that if i,j,k,l are indexes for the four centers then the following
 symmetry property applies : <ij|kl> = <ji|kl> = <kl|ij>...
 Then, there are only evaluated those integrals which indexes of the basis set, are the smallest one.
 Besides, it is provided a function for sorting the keys using the symmetry properties of the Integrals.
-}
calcIntegrals :: ML.Map Boys Double ->  [AtomData] -> Array U DIM1 Double
calcIntegrals gridBoys atoms = fromListUnboxed (ix1 $ length cartProd) $ parMap rdeepseq funEval cartProd 

  where dim      = pred . sum . fmap (length . getBasis) $ atoms
        funEval  = evalIntbykey gridBoys atoms 
        cartProd = do
          i <- [0..dim]
          j <- [i..dim]
          k <- [i..dim]
          l <- [k..dim]
          let xs = [i,j,k,l]
          guard (condition xs)
          return [i,j,k,l]
        condition e = case compare e $ sortKeys e of
                       EQ        -> True
                       otherwise -> False                               

-- | <ii||kl> == <ij||kk> = 0
-- zeroCondition :: [Int] -> Bool
-- zeroCondition [i,j,k,l] | i == j && k /= i && k /= l && l /= i = True
--                         | k == l && i /= j && i /= k && j /= k = True
--                         | otherwise = False  

-- ============> Density Matrix <=================

{- | this function only sum over half the coefficients corresponding with the occupied
     molecular orbitals in order to  build up the density matrix
     Notice that we are assuming a double ocupation in each Shell, therefore
    this matrix differs from the monoelectronic density matrix -}                   
calcDensity :: Monad m => Array U DIM2 Double -> OccupiedShells -> m FlattenChargeDensity
calcDensity !arr occupied =
   do   let (Z :. dim  :. _)  = extent arr
            flat = (dim^2 + dim)`div`2
        computeP
         $ fromFunction (Z :. flat)
         (\idx  ->
                  let (Z:.k:.l) =  indexFlat2DIM2 dim idx
                      v1 = R.slice arr (Any :. k :. All)
                      v2 = R.slice arr (Any :. l :. All)
                      p1 = R.extract (Z:.0) (Z :. occupied) v1
                      p2 = R.extract (Z:.0) (Z :. occupied) v2
                  in (*2) . R.sumAllS $ p1 *^ p2)
{-# INLINE calcDensity #-}
                  
-- ==================================> CALCULATED THE FOCK MATRIX <======================================          
{-| Calculate the current fock matrix as a flatten matrix summing the core Hamiltonian plus the electronic
 interaction integrals-}
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
calcGmatrix atoms density integrals =
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

  where (Z:.a:.b) = indexFlat2DIM2 nshells $ ix1 i
        calcIndex = fourIndex2Flat nshells . sortKeys
{-# INLINE sliceIntegrals #-}        

-- |Diagonalize the Fock matrix
diagonalHF :: Monad m => FlattenFock ->
                         TransformMatrix->
                         OccupiedShells ->
                         m (MOCoefficients,FlattenChargeDensity,EigenValues)
diagonalHF !fock1 !xmatrix !occupied = do
        fDIM2          <- unitaryTransf xmatrix fock1
        let (orbEs,coeff) = eigenSolve fDIM2
        newCoeff       <- mmultP xmatrix coeff
        newDensity     <- calcDensity newCoeff occupied
        return (newCoeff, newDensity, orbEs)
{-# INLINE diagonalHF #-}                
                   
                  
-- | Hartree-Fock total energy
variationalE :: FlattenCore -> FlattenFock -> FlattenChargeDensity -> Double
variationalE core fockMtx newDensity = (0.5*) $ sumAllS $ fromFunction (ix2 dim dim ) $
 \(Z:.i :. j ) -> let uv = ix1 $ indexDIM2toFlat dim i j 
                      vu = ix1 $ indexDIM2toFlat dim j i
                      pvu = newDensity ! vu 
                      huv = core ! uv
                      fuv = fockMtx ! uv                     
                  in pvu * (huv + fuv) 

 where dim    = dimTriang d 
       (Z:.d) = extent core


-- | Function to check convergency
converge :: FlattenChargeDensity -> FlattenChargeDensity -> Bool
converge !oldDensity !newDensity  = sigma < 1.0e-6 
  where sigma = sqrt . (0.25*) . R.sumAllS . R.computeUnboxedS . R.map (^2) $ oldDensity -^ newDensity
{-# INLINE converge #-}
  

-- =================> Nuclear Repulsion Energy <=====================

-- | Calculate internuclear electronic repulsion
nuclearRep :: [AtomData] -> Double
nuclearRep xs = sum [repulsion atomi atomj | atomi <- xs, atomj <- xs, atomj > atomi]
  where dim = pred . length $ xs
        repulsion at1 at2 = let ([za,zb],[ra,rb]) = fmap getZnumber &&& fmap getCoord $ [at1,at2]
                                rab = sqrt . sum . fmap (^2) $ L.zipWith (-) ra rb
                            in za*zb/rab



