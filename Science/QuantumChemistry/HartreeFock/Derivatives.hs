{-# Language BangPatterns #-}

-- The HaskellFock SCF Project
-- @2013 Felipe Zapata, Angel Alvarez
-- Energy Gradient and Hessian (be patient for the last one)

module Science.QuantumChemistry.HartreeFock.Derivatives  (
                    energyGradient
                    ) where

import Control.Applicative
import Control.Arrow ((&&&),first)
import Control.DeepSeq
import Control.Monad ((<=<),(>=>))
import Control.Monad.List (guard)
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix as R
import qualified Data.List as DL
import qualified Data.Vector.Unboxed as VU

-- =========> Internal Modules <=========
import Science.QuantumChemistry.GlobalTypes as G
import Science.QuantumChemistry.NumericalTools.LinearAlgebra as LA
import Science.QuantumChemistry.HartreeFock.HartreeFock
import Science.QuantumChemistry.HartreeFock.IntegralsEvaluation

-- ===========> Local Types <====================



-- ==================> <====================
{- |First derivate of the Energy with respect to all the Nuclear,
Coordinates. It is composed of four terms, the derivatives of the
Nuclear repulsion, the core Hamiltonian, the overlap Matrix and
the Coulombic and exchange integrals-}

energyGradient :: Monad m => [AtomData] -> HFData -> m Gradient
energyGradient !atoms !hfData = do
        termCore    <- calcCoreTerm      atoms hfData norbital
        termJK      <- calcIntegralsTerm atoms hfData norbital
        termOverlap <- calcOverlapTerm   atoms hfData norbital
        termNuclear <- firstNuclearDerivative atoms
        computeUnboxedP $ termCore +^ termJK -^ termOverlap +^ termNuclear

  where norbital = sum . fmap (length . getBasis) $ atoms         

-- =================> Electron-Electron Integrals Derivatives  <================

calcIntegralsTerm :: Monad m => [AtomData] -> HFData -> Int  -> m Gradient
calcIntegralsTerm !atoms !hfData !norbital = do
  jkDerivatives <- calcJKDerivatives atoms norbital
  computeUnboxedP $ fromFunction (Z:. dim) $
   \(Z:. k) -> let dJKuvrs_dXa = computeUnboxedS $ unsafeSlice jkDerivatives (Any :. k :. All)
               in sumAllS $ fromFunction (Z:. dimTriang :. dimTriang) $
                   \(Z:. m :. n) ->
                      let (Z:.u:.v) = LA.indexFlat2DIM2 norbital $ ix1 m
                          (Z:.r:.s) = LA.indexFlat2DIM2 norbital $ ix1 n
                          uv = calcij u v
                          rs = calcij r s
                          totalkr = product $ DL.zipWith calcKr [u,r,uv] [v,s,rs]
                          [ur,us,vr,vs] = DL.zipWith (indexDIM2toFlat norbital) [u,u,v,v] [r,s,r,s]
                          [duv,dur,dus,dvr,dvs,drs] = fmap (\x -> density ! (Z:.x)) [m,ur,us,vr,vs,n]
                          duvrs = 4*duv*drs - dur*dvs - dus*dvr
                          uvrs = calcIndex [u,v,r,s]
                          dJK_u = dJKuvrs_dXa ! (Z:. 0 :. uvrs)
                          dJK_v = dJKuvrs_dXa ! (Z:. 1 :. uvrs)
                          dJK_r = dJKuvrs_dXa ! (Z:. 2 :. uvrs)
                          dJK_s = dJKuvrs_dXa ! (Z:. 3 :. uvrs)
                      in totalkr * duvrs * (dJK_u + dJK_v + dJK_r + dJK_s)

  where dim = 3 * (length atoms)
        dimTriang = (norbital^2 + norbital) `div`2
        density = getDensity hfData
        calcij i j = let [u,v] = fromIntegral `fmap` [i,j]
                     in floor $ 0.5 * (u * (u-1)) + v
        calcKr i j = 1 - (kronecker i j) / 2
        calcIndex = fourIndex2Flat norbital . sortKeys

calcJKDerivatives :: Monad m => [AtomData] -> Int -> m (Array U DIM3 Double)        
calcJKDerivatives !atoms !norbital = computeUnboxedP . fromFunction (Z :. dimCart :. 4 :. dimJK) $
  \(Z:. k :. l :. m) ->
     let dx = G.toCartLabel $ k `rem` 3
         xs = cartProd !! m
         centers = fmap (LA.calcCoordCGF atoms) xs
         (rl,cgfl) = centers !! l
         derv_cgfl = calcDervCGF dx cgfl
     in dervJKuvrs centers (rl,derv_cgfl) l
                  
  where dimJK = length cartProd
        dimCart   = 3 * (length atoms)
        cartProd = do
          let dim = pred norbital
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

-- rewrite this Ugly function!!!
dervJKuvrs :: [(NucCoord,CGF)] -> (NucCoord,[CGF]) -> Int -> Double
dervJKuvrs [!at1,!at2,!at3,!at4] (!rl,!derv_cgfl) !l =
  case l of
       0 -> sum $ fmap (\cgfDerv -> contracted4Centers [(rl,cgfDerv),at2,at3,at4]) derv_cgfl
       1 -> sum $ fmap (\cgfDerv -> contracted4Centers [at1,(rl,cgfDerv),at3,at4]) derv_cgfl
       2 -> sum $ fmap (\cgfDerv -> contracted4Centers [at1,at2,(rl,cgfDerv),at4]) derv_cgfl
       3 -> sum $ fmap (\cgfDerv -> contracted4Centers [at1,at2,at3,(rl,cgfDerv)]) derv_cgfl

-- =================> Core Hamiltonian Derivatives <=======================

calcCoreTerm :: Monad m => [AtomData] -> HFData -> Int -> m Gradient
calcCoreTerm !atoms !hfData !norbital = do
  coreDerivatives <- calcCoreDerivatives atoms norbital
  computeUnboxedP $ fromFunction (Z:. dim) $
   \(Z:. k) -> let dhuv_dXa = computeUnboxedS $ unsafeSlice coreDerivatives (Any :. k :. All)
               in (2*) . sumAllS $ fromFunction (extent dhuv_dXa) $
                   \sh -> let (Z:.u:.v) = LA.indexFlat2DIM2 norbital sh
                              duv = density ! sh
                              huv = dhuv_dXa ! sh
                              kr = 1- (kronecker u v) / 2
                          in kr * duv * huv

  where dim = 3 * (length atoms)
        dimTriang = (norbital^2 + norbital) `div`2
        density = getDensity hfData
                                  
calcCoreDerivatives :: Monad m => [AtomData] -> Int -> m Matrix
calcCoreDerivatives !atoms !norbital = computeUnboxedP . fromFunction (Z :. dimCart :. dimTriang) $
  \(Z:. k :. m) -> let calcIndex = LA.calcCoordCGF atoms
                       (Z:.i:.j) = LA.indexFlat2DIM2 norbital $ ix1 m
                       [(r1,cgf1),(r2,cgf2)] = fmap calcIndex $ [i,j]
                       dx = G.toCartLabel $ k `rem` 3
                       rA = getCoord $ atoms !!  (k`div`3)
                       [derv_cgf1,derv_cgf2] = fmap (calcDervCGF dx) [cgf1,cgf2]
                       firstTerm  = coreterm1 atoms r1 (r2,cgf2) derv_cgf1
--                        secondTerm = coreterm2 (r1,cgf1) (r2,cgf2) rA dx
                       secondTerm = 0.0
                       thirdTerm  = coreterm1 atoms r2 (r1,cgf1) derv_cgf2
                   in firstTerm + secondTerm + thirdTerm

  where dimTriang = (norbital^2 + norbital) `div`2
        dimCart   = 3 * (length atoms)        
        
coreterm1 :: [AtomData] -> NucCoord -> (NucCoord,CGF) -> [CGF] -> Double
coreterm1 atoms !r1 t2@(!r2,!cgf2) !derv_cgf1 =
  sum $!! do
      cgfDerv <- derv_cgf1
      let t1 = (r1,cgfDerv)
          derv = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)
          sumVij = sum $!! DL.zipWith (\z rc -> ((-z) * t1 <<|Vij rc derv |>> t2)) atomicZ coords
      return $  (t1 <<|Tij|>> t2) + sumVij

  where coords = fmap getCoord atoms
        atomicZ = fmap getZnumber atoms

coreterm2 :: (NucCoord,CGF) -> (NucCoord,CGF) -> NucCoord -> CartesianLabel  -> Double
coreterm2 t1 t2 rA dx | dx == Ax = t1 <<| Vij rA (Dij_Ax 1, Dij_Ay 0, Dij_Az 0) |>> t2
                      | dx == Ay = t1 <<| Vij rA (Dij_Ax 0, Dij_Ay 1, Dij_Az 0) |>> t2
                      | dx == Az = t1 <<| Vij rA (Dij_Ax 0, Dij_Ay 0, Dij_Az 1) |>> t2

-- ===============> Overlap Matrix Derivatives <==================          

calcOverlapTerm :: Monad m =>  [AtomData] -> HFData -> Int -> m Gradient
calcOverlapTerm !atoms !hfdata !norbital = do
  wFlatten <- calcWmatrix hfdata norbital
  overlapDerivatives <- calcOverlapDerivatives atoms norbital
  computeUnboxedP $ fromFunction (Z:. dim) $
   \(Z:. k) -> let dSuv_dXa = computeUnboxedS $ unsafeSlice overlapDerivatives (Any :. k :. All)
               in(4*) . sumAllS $ fromFunction (extent dSuv_dXa) $
                      \sh -> let (Z:.u:.v) = LA.indexFlat2DIM2 norbital sh
                                 wuv = wFlatten ! sh
                                 suv = dSuv_dXa ! sh
                                 kr = 1- (kronecker u v) / 2
                             in kr * wuv * suv
                                    
  where dim = 3 * (length atoms)
{-# INLINE calcOverlapTerm #-}

-- the <<||>> operator is defined as the sijContracted function
-- in the integrals module
-- | derivatives of the overlap matrix with respect to the Nuclear Coordinates
calcOverlapDerivatives :: Monad m => [AtomData] -> Int -> m Matrix
calcOverlapDerivatives !atoms !norbital =  computeUnboxedP . fromFunction (Z :. dimCart :. dimTriang) $
  \(Z:. k :. m) -> let calcIndex = LA.calcCoordCGF atoms
                       (Z:.i:.j) = LA.indexFlat2DIM2 norbital $ ix1 m
                       [(r1,cgf1),(r2,cgf2)] = fmap calcIndex $ [i,j]
                       dx = G.toCartLabel $ k `rem` 3
                       [derv_cgf1,derv_cgf2] = fmap (calcDervCGF dx) [cgf1,cgf2]
                   in sum $ do -- monad List to sum over the result of derivate the CGF
                        derv1 <- derv_cgf1
                        derv2 <- derv_cgf2
                        return $  (r1,derv1) <<||>> (r2,cgf2) + (r1,cgf1) <<||>> (r2,derv2)

  where dimTriang = (norbital^2 + norbital) `div`2
        dimCart   = 3 * (length atoms) 
  
-- | Matrix which elements multiply the derivatives of the overlap matrix
calcWmatrix :: Monad m =>  HFData -> Int -> m FlattenMatrix       
calcWmatrix hfData norbital = do
        let (fockFlatten,cs) = getFock &&& getCoeff $ hfData
            (Z:. h1 :._) = extent cs
            dim = (h1^2 + h1) `div` 2 -- number of elements of the triangular matrix
        fockMtx <- triang2DIM2 fockFlatten
        epsilonMtx <- transpose2P cs >>= (\csT -> mmultP csT <=< mmultP fockMtx $ cs)
        computeUnboxedP $ fromFunction (Z:.dim) $
         \idx -> let (Z:.u:.v) = LA.indexFlat2DIM2 norbital idx
                 in sumAllS $ fromFunction (extent cs) $
                    \sh@(Z:. i :. j) -> let cui = (cs ! (Z:. u:.i))
                                            cvj = (cs ! (Z:. v:.j))
                                            eij = (epsilonMtx ! sh)
                                        in cui * cvj * eij
{-# INLINE calcWmatrix #-}

-- ==========================> Nuclear repulsion derivatives <=========================          

-- | Derivatives of the nuclei-nuclei repulsion energy with respect to the nuclear coordinates                                               
firstNuclearDerivative :: Monad m => [AtomData] -> m Gradient
firstNuclearDerivative atoms =
  computeUnboxedP $ reshape (Z :. totalDIM) $ fromFunction (Z:. numat :. 3) $
    \(Z:.i:.j) -> let atomA = atoms !! i
                      (za,ra) = getZnumber &&& getCoord $ atomA
                  in negate . sum $ do
                     atomB <- atoms
                     guard $ atomB /= atomA                     
                     let (zb,rb) = getZnumber &&& getCoord $ atomB
                         rab = DL.zipWith (-) ra rb
                         rab_mod = sqrt . sum . fmap (^2) $ rab
                         rab3 = rab_mod^3
                         cte = za*zb/rab3
                     return $!!  cte * (rab !! j)

  where  numat = length atoms                                                      
         totalDIM = 3 * numat                                                      
        
-- ============> Auxuliar Functions <================

-- | Calculate the derivatives with respect to the x, y or z component of the 
-- | Gaussian functions which composed the Contracted Gaussian Functions
calcDervCGF :: CartesianLabel -> CGF -> [CGF]
calcDervCGF ax cgf =
  let (angularMomentum,primitives) = G.getfunTyp &&& G.getPrimitives $ cgf
      up   = angularMomentum |+> ax -- increase and decrease operator
      down = angularMomentum |-> ax -- defined in GlobalTypes
      i = getCGFangularExpo angularMomentum ax
      primitivesUp   = fmap (\(coeff,expo) -> (2*expo*coeff,expo)) primitives
      primitivesDown = fmap (first (negate . (*i))) primitives
      giPlus = CGF primitivesUp up
      giDown = CGF primitivesDown down
  in if down == Zero then [giPlus]
                     else [giPlus,giDown]


kronecker :: Int -> Int -> Double
kronecker i j | i == j    = 1
              | otherwise = 0
