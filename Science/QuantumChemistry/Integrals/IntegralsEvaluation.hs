{-# Language BangPatterns #-}

-- The HaskellFock SCF Project
-- @2013 Felipe Zapata, Angel Alvarez
-- Analytical evaluation of the overlap, core, nuclei-electron
-- and electron-electron integrals


module Science.QuantumChemistry.Integrals.IntegralsEvaluation 
  (
   (|>)
  ,(<|)
  ,(|>>)
  ,(<<|)
  ,(<<||>>)
  ,Operator(..)
  ,contracted4Centers
  ,evalIntbykey
  ,hcore
  ,mtxOverlap
  ,normaCoeff
  )
  where

import Control.Applicative
import Control.Arrow ((&&&),first,second)
import Control.DeepSeq
import Control.Monad (liftM,mplus,sequence)
import Control.Monad.List
import Control.Monad.State
import Control.Parallel.Strategies (parMap,rdeepseq)
import Data.Array.Repa         as R hiding (map)
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix (mmultP)
import Data.List as L
import qualified Data.Map as M
import Data.Maybe (fromMaybe)
import Data.Monoid (Monoid(..),mappend,mconcat,mempty)
import qualified Data.Vector.Unboxed as U

-- internal modules 
import Science.QuantumChemistry.NumericalTools.Boys(boysF)
import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.NumericalTools.LinearAlgebra as LA

--  ======================= > MODULE FOR EVALUATING THE ONE AND TWO ELECTRON TYPE INTEGRALS < =============================
{- | This module is based on chapter 9 "Molecular Integral Evaluation" in the book
   "Molecular Electronic-Structure Theory". Authors: Trygve Helgaker, Poul Jorgensen and Jeppe Olsen

   | Remember that the Gaussian primitives are not normalized they must be multiplied for a factor sqr (4*alfa/pi)
   | where alfa is the coefficient of the gaussian function-}

-- parMap strat f xs = L.map f xs `using` parList strat


-- =================================================================================================
-- == the Monad List are used for calculating the elements of the lower triangular matrix         ==
-- == the elements of the triangular matrix is equal to (n^2 + n) / 2 with n the number of atoms  ==
-- =================================================================================================

-- =================================================================================================
-- | The idea for the integral is using to use the monad list for calculating
-- | the cartesian product between basis, CGF and primitive gaussian functions.
-- |  Therefore the integration is of 3 levels.
-- =================================================================================================


-- ==============> TYPES <=============

type CartesianDerivatives = (Derivatives Int,Derivatives Int,Derivatives Int)

type Exponent = Double 

type MapHermite = M.Map HermiteIndex Double

type Pairs = (Double,Double)

type Triplet = (Int,Int,Int)

-- ============> ALGEBRAIC DATA TYPES <=================
    
data Operator    =  Rij          -- ^ Internuclear interaction
                  | Tij          -- ^ Electronic kinetic operator
                  | Vij NucCoord CartesianDerivatives -- ^ Nucleus-electron operator
                 deriving Show
  
                
data HermiteIndex = Rpa {getN :: Int, getijt :: [Int]}
                    | Xpa {getijt :: [Int]}
                    | Ypa {getijt :: [Int]}
                    | Zpa {getijt :: [Int]} 
                    deriving (Show,Eq,Ord)

data EvalHermiteCoeff = EvalHermiteCoeff {
                        pabComp  :: [Double]
                       ,getIndex :: [Int]
                       ,cartComp :: [Int] -> HermiteIndex
                          }            

data HermiteStateCoeff = HermiteStateCoeff  {
               getmapC :: MapHermite
              ,getkeyC :: [[HermiteIndex]]
               } deriving Show

data HermiteStateIntegral = HermiteStateIntegral {
               getmapI :: MapHermite
              ,getkeyI :: [HermiteIndex]
               } deriving Show
               
               
-- ==============> Some Utilities <===================

unfoldWhile :: (b -> Bool) -> (b -> (a, b)) -> b -> [a]   
unfoldWhile p f = L.unfoldr (\x -> guard (not (p x)) >> return (f x))     

safeHead :: [a] -> Maybe a
safeHead (x:xs) = Just x
safeHead []     = Nothing
            
getExponentDerivatives :: CartesianDerivatives -> Int
getExponentDerivatives (Dij_Ax x,Dij_Ay y, Dij_Az z) = x + y + z


-- ===============> OPERATORS <===========
 
 -- | Overlap integral over primitive guassians
(<||>) :: Gauss -> Gauss -> Double
gauss1 <||> gauss2 = sab gauss1 gauss2

-- | Primitive bra in the Dirac notation
(<|) :: Gauss-> Operator -> (Gauss,Operator)
gauss1 <| op = (gauss1, op)

-- | Primitive ket in the dirac notation
(|>) :: (Gauss,Operator) -> Gauss -> Double
(gauss1,op) |> gauss2 = case op of
                             Tij -> tab gauss1 gauss2

-- | Overlap between a set of contracted Gauss functions                            
(<<||>>) :: (NucCoord,CGF) -> (NucCoord,CGF) -> Double
(r1,cgf1) <<||>> (r2,cgf2) = sijContracted r1 r2 cgf1 cgf2

-- | bra for contracted Gaussian functions                             
(<<|) :: (NucCoord,CGF) -> Operator ->  ((NucCoord,CGF),Operator)
b1 <<| op  = (b1,op)

-- | ket for contracted Gaussian functions
(|>>) :: ((NucCoord,CGF),Operator) -> (NucCoord,CGF) -> Double
(b1,op) |>> b2 = case op of
                Tij -> tijContracted b1 b2
                Vij rc derivatives -> vijContracted b1 rc b2 derivatives
                                                    
-- ==============> 2-CENTER OVERLAP INTEGRAL <=================================

-- | overlaping between two primitive gaussian function between arbitrary-l functions  <A|B>

-- S12 = *exp[−a1*b1(Rab)^2/gamma]IxIyIz

-- where Ix = sum f2i (l1,l2,PAx,PBx) ((2i-1)!! / (2gamma)^i sqrt (pi/gamma) for i = 0,l1 + l2 / 2 only the even i
-- f2k is defined above in the module

mtxOverlap :: Monad m => [AtomData] -> m (Array U DIM1 Double) 
mtxOverlap atoms = computeUnboxedP . fromFunction (Z:.dim) $
 (\idx -> let (Z:.i:.j) = LA.indexFlat2DIM2 norbital idx
              [(r1,cgf1),(r2,cgf2)] = map calcIndex $ [i,j]
          in sijContracted r1 r2 cgf1 cgf2)
             
  where norbital = sum . map (length . getBasis) $ atoms
        dim = (norbital^2 + norbital) `div`2
        calcIndex = LA.calcCoordCGF atoms

-- | Overlap matrix entry calculation between two Contracted Gaussian functions
sijContracted :: NucCoord -> NucCoord -> CGF -> CGF -> Double
sijContracted r1 r2 cgf1 cgf2 =
                 if cgf1 == cgf2 && r1 == r2
                                  then  1.0
                                  else  sum $ do
                              g1 <- getPrimitives cgf1
                              g2 <- getPrimitives cgf2
                              let [l1,l2] = map getfunTyp [cgf1,cgf2]
                                  gauss1 = Gauss r1 l1 g1
                                  gauss2 = Gauss r2 l2 g2
                              return (gauss1 <||> gauss2 )


-- | Primitive overlap terms                              
sab :: Gauss -> Gauss -> Double
sab g1@(Gauss r1 shellA (c1,e1)) g2@(Gauss r2 shellB (c2,e2)) =
   (c1*c2) * (product $!! [obaraSaika gamma (s00 !! x) (pa !! x) (pb !! x) (l1 x) (l2 x) |x <- [0..2]])
  
  where [l1,l2] = map funtyp2Index [shellA,shellB]
        [pa,pb] = map (L.zipWith (-) p) [r1,r2]
        p = meanp (e1,e2) r1 r2
        expo = \x2 -> exp $ -mu * x2
        s00 = map ((*cte) . expo . (^2)) $  L.zipWith (-) r1 r2
        cte = sqrt $ pi/ gamma
        gamma = e1+ e2
        mu = e1*e2/gamma

-- | ObaraSaika Scheme to calculate overlap integrals
obaraSaika ::Double -> Double -> Double -> Double -> Int -> Int  -> Double
obaraSaika gamma s00 pax pbx i j = s i j

  where pred2 = pred . pred
        c = recip $ 2.0 * gamma
        s !m !n | m < 0 || n < 0   =  0.0
                | all (== 0) [m,n] =  s00
                | m == 0 =   pbx * (s 0 (pred n)) + c* (n_1 * (s 0 $ pred2 n))
                | n == 0 =   pax * (s (pred m) 0) + c*(m_1 * (s (pred2 i) 0))
                | otherwise  =   pax * (s (pred m) n) + c * (m_1 * (s (pred2 m) n) +
                               (fromIntegral n) * (s (pred m) $ pred n))

          where [m_1,n_1] = map (fromIntegral . pred) [m,n]
         

         
-- ====================> HAMILTONIAN CORE <=======================

-- | Core Hamiltonian Matrix Calculation
hcore :: Monad m => [AtomData] -> m (Array U DIM1 Double)
hcore !atoms  = computeUnboxedP . fromFunction (Z:. dim) $
 (\idx -> let (Z:.i:.j) = LA.indexFlat2DIM2 norbital idx
              [atomi,atomj] = map calcIndex $ [i,j]
              derv = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)
              sumVij = sum $!! L.zipWith (\z rc -> ((-z) * atomi <<|Vij rc derv |>> atomj)) atomicZ coords
             in (atomi <<|Tij|>> atomj) + sumVij)

  where coords    = map getCoord atoms
        atomicZ   = map getZnumber atoms
        norbital  = sum . map (length . getBasis) $ atoms
        dim       = (norbital^2 + norbital) `div`2
        calcIndex = LA.calcCoordCGF atoms
             
-- ==================> Kinetic Energy IntegralsEvaluation  <========================
-- |the kinetic integral for two S-functions is
-- | < A| -0.5*nabla^2 | B > = (3*w - 2*w*w*rab2 ) * <A|B>
-- | where w = e1*e2 / (e1 +e2)

-- | For arbitrary angular momentum functions is
-- | < A| -0.5*nabla^2 | B > = c1*c2 (TijSklSmn + SijTklSmn + SijSklTmn)
-- |where Tij = -2b^2 S(i,j+1) + b(2j+1)S(i,j) - 0.5j(j-1)S(i,j-2)
 
-- | Kinetic Energy operator matrix representation
tijContracted :: (NucCoord,CGF) -> (NucCoord,CGF) -> Double
tijContracted (r1,cgf1) (r2,cgf2) =
          sum $!! do
              g1 <- getPrimitives cgf1
              g2 <- getPrimitives cgf2
              let [l1,l2] = map getfunTyp [cgf1,cgf2]
                  gauss1 = Gauss r1 l1 g1 
                  gauss2 = Gauss r2 l2 g2
              return ( gauss1 <| Tij |> gauss2 )  

-- | Primitive kinetic energy terms          
tab :: Gauss -> Gauss -> Double
tab gA@(Gauss r1 shell1 (c1,e1)) gB@(Gauss !r2 !shell2 (!c2,!e2)) =
   c1*c2 * (sum . map product $! permute [\x -> tx x (j x) (k x),\x -> sx x (j x) (k x),\x -> sx x (j x) (k x)] [0..2])

  where [j,k] = map funtyp2Index [shell1,shell2]
        sx i lang1 lang2 = obaraSaika gamma (s00 !! i) (pa !! i) (pb !! i) lang1 lang2
        tx i lang1 lang2 = let [l1,l2] = map fromIntegral [lang1,lang2]
                           in -2.0 * e2^2 * (sx i lang1 (lang2+2)) +
                              e2*(2*l2 +1)* (sx i lang1 lang2) -
                              0.5*l2*(l2-1)*(sx i lang1 (lang2 -2))
        cte = sqrt (pi/ gamma)
        s00 = map ((*cte) . expo . (^2)) $  L.zipWith (-) r1 r2
        t00 = \x -> e1 - 2*e1^2 *((pa !! x)^2 + (recip $ 2*gamma)) * (s00 !! x )
        expo = \x2 -> exp $ -mu * x2
        [pa,pb] = map (L.zipWith (-) p)  [r1,r2]
        p = meanp (e1,e2) r1 r2
        gamma = e1+ e2
        mu = e1*e2/gamma        
        
permute :: [a -> b] -> [a] -> [[b]]
permute [f,g,h] xs = L.zipWith (L.zipWith ($)) [[f,g,h],[g,f,h],[g,h,f]] $ repeat xs

-- =====================> TWO ELECTRON INTEGRALS <=========================================

-- | Main function to calculate the Coulomb (J) and Interchange Integrals using the
--  indexes of the atomic basis <ab|cd>
evalIntbykey :: [AtomData] -> [Int]->  Double
evalIntbykey atoms keys = contracted4Centers centers 
 where centers = (LA.calcCoordCGF atoms) <$> keys

-- | Calculate electronic interaction among the four center contracted Gaussian functions
contracted4Centers :: [(NucCoord,CGF)] -> Double
contracted4Centers [(ra,cgf1), (rb,cgf2), (rc,cgf3), (rd,cgf4)] = sum cartProd
  where [l1,l2,l3,l4] = map getfunTyp [cgf1,cgf2,cgf3,cgf4]
        cartProd = do
          g1 <- getPrimitives cgf1
          g2 <- getPrimitives cgf2
          g3 <- getPrimitives cgf3
          g4 <- getPrimitives cgf4
          let gauss = L.zipWith3 Gauss [ra,rb,rc,rd] [l1,l2,l3,l4] [g1,g2,g3,g4]
              [f1,f2,f3,f4] = L.map funtype gauss
              [ga,gb,gc,gd] = gauss
--          guard $ schwarz gauss 
          return $ if ( ga == gc) && (gb == gd) 
                      then twoTermsERI ga gb
                      else twoElectronHermite gauss

schwarz :: [Gauss] -> Bool 
schwarz gs@[g1,g2,g3,g4] =  True -- if val < 1e-8 then False else True 
 
 where val = qab * qcd  
       qab = sqrt $ twoTermsERI g1 g2
       qcd = sqrt $ twoTermsERI g3 g4


twoTermsERI :: Gauss -> Gauss -> Double
twoTermsERI ga gb = 
  case L.map funtype gs of 
       [S,S ]    -> sTypeIntegral 0 g1 g2
       -- [S,Px]    -> twoTerms_SP 0 g1 g2  
       -- [S,Py]    -> twoTerms_SP 1 g1 g2  
       -- [S,Pz]    -> twoTerms_SP 2 g1 g2
       -- [Px,Px]   -> twoTerms_PP 0 0 g1 g2
       -- [Px,Py]   -> twoTerms_PP 0 1 g1 g2
       -- [Px,Pz]   -> twoTerms_PP 0 2 g1 g2
       -- [Py,Py]   -> twoTerms_PP 1 1 g1 g2      
       -- [Py,Pz]   -> twoTerms_PP 1 2 g1 g2
       -- [Pz,Pz]   -> twoTerms_PP 2 2 g1 g2
       -- [Dxx,S]   -> twoTerms_DS 0 0 g1 g2 
       -- [Dxy,S]   -> twoTerms_DS 0 0 g1 g2 
       -- [Dxz,S]   -> twoTerms_DS 0 0 g1 g2 
       -- [Dyy,S]   -> twoTerms_DS 0 0 g1 g2 
       -- [Dyz,S]   -> twoTerms_DS 0 0 g1 g2 
       -- [Dzz,S]   -> twoTerms_DS 0 0 g1 g2 
        
-- Dxx | Dxy | Dxz | Dyy | Dyz | Dzz 
       otherwise -> twoElectronHermite [g1,g2,g1,g2]


  where gs@[g1,g2]          = L.sortBy compareAngularM [ga,gb]
        compareAngularM x y = let [f1,f2] = L.map funtype [x,y]
                              in  compare f1 f2 
                               
twoTerms_SP :: Int -> Gauss -> Gauss -> Double
twoTerms_SP i gS gP = (rP !! i) * (sTypeIntegral 1 gS gP)
   where rP = nucCoord gP

twoTerms_PP :: Int -> Int -> Gauss -> Gauss -> Double
twoTerms_PP i j g1 g2 = negate (xa* xb * (sTypeIntegral 2 g1 g2)) + dij 

  where dij     = deltaij i j $ sTypeIntegral 1 g1 g2  
        [ra,rb] = L.map nucCoord [g1,g2]
        xa      = ra !! i
        xb      = rb !! j  


sTypeIntegral ::  Int -> Gauss -> Gauss -> Double 
sTypeIntegral ordIntegral g1@(Gauss r1 angFun1 (c1,e1)) g2@(Gauss r2 angFun2 (c2,e2)) =         
  (g1 <||> g2)^2 * (sqrt $ 2 * sigmaM / pi)
 
  where sigmaM = (e1+e2)^(2*m+1)
        m      = ordIntegral

twoElectronHermite :: [Gauss] -> Double
twoElectronHermite gs = (cte *) . U.sum . U.zipWith (*) coeff1 $!! U.fromList
                        [mcMurchie2 coeff2 rpq alpha abcs tuv | tuv <- coefftuv]

  where coefftuv = genCoeff_Integral [symb1,symb2] derv
        abcs =  genCoeff_Integral [symb3,symb4] derv
        coeff1 = U.unfoldr (calcHermCoeff [rpa,rpb] p) seedC1
        coeff2 = U.unfoldr (calcHermCoeff [rqc,rqd] q) seedC2
        seedC1 = initilized_Seed_Coeff [symb1,symb2] rab mu
        seedC2 = initilized_Seed_Coeff [symb3,symb4] rcd nu
        [ra,rb,rc,rd] = map nucCoord gs
        [symb1,symb2,symb3,symb4] = map funtype  gs
        ps = map gaussP gs
        ([c1,c2,c3,c4],[e1,e2,e3,e4]) = (map fst ) &&& (map snd ) $ ps
        rp = meanp (e1,e2) ra rb
        rq = meanp (e3,e4) rc rd
        [rab,rcd,rpa,rpb,rqc,rqd,rpq] = map restVect  (zip [ra,rc,rp,rp,rq,rq,rp][rb,rd,ra,rb,rc,rd,rq])
        [mu,nu] = (\(a,b) -> a*b* (recip $ a + b)) `fmap` [(e1,e2),(e3,e4)]
        [p,q] = uncurry (+) `fmap` [(e1,e2),(e3,e4)]
        alpha = p*q/(p+q)
        cte = (c1*c2*c3*c4*) . (*(2.0*pi**2.5)) . recip $ (p * q ) * (sqrt $ p + q)
        derv = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)


mcMurchie2 :: VecUnbox -> NucCoord -> Double -> [HermiteIndex] -> HermiteIndex -> Double
mcMurchie2 coeff2 rpq alpha abcs' tuv' = U.sum $!!  U.zipWith3 (\x y z -> x*y*z) sgns coeff2  integrals
  where tuv = getijt tuv'
        abcs = map getijt  abcs'
        sgns = U.fromList $ map (\xs-> (-1.0)^(sum xs)) abcs
        integrals = U.unfoldr (calcHermIntegral rpq alpha) seedI
        seedI = HermiteStateIntegral mapI0 listI
        mapI0 = M.insert k0' f0 M.empty
        k0'= Rpa 0 [0, 0, 0]
        rp2 = sum $ map (^2) rpq
        y = alpha*rp2
        f0 = boysF 0 y
        listI = Rpa 0 `fmap` (map (L.zipWith (+) tuv) abcs )


-- ======================> McMURCHIE -DAVIDSON SCHEME <=========================
-- | The integrals of three and four centers are implemented according to chapter 9
--   of the book "Molecular Electronic-Structure Theory" by Trygve Helgaker,
--   Poul Jorgensen and Jeppe Olsen.

-- | Nuclei-Electron matrix operator presentation        
vijContracted :: (NucCoord,CGF) -> NucCoord -> (NucCoord,CGF) -> CartesianDerivatives -> Double
vijContracted (ra,cgf1) rc (rb,cgf2) derv = sum $!
        do g1 <- getPrimitives cgf1
           g2 <- getPrimitives cgf2
           let [l1,l2] = map getfunTyp [cgf1,cgf2]
               gauss1 =  Gauss ra l1 g1
               gauss2 =  Gauss rb l2 g2
           return $!! vijHermite gauss1 gauss2 rc derv
           
-- | Hermite auxiliar function calculations according to the McMURCHIE -DAVIDSON scheme        
vijHermite :: Gauss -> Gauss -> NucCoord -> CartesianDerivatives -> Double
vijHermite g1 g2 rc derv = ((-1)^sumDervExpo) * cte * (mcMurchie shells [ra,rb,rc] (e1,e2) derv)
  where cte = c1 * c2 * 2.0 * (pi/gamma) 
        gamma = e1+e2     
        [ra,rb] =  map nucCoord [g1,g2]
        shells = map funtype [g1,g2]
        [(c1,e1),(c2,e2)] = map gaussP [g1,g2]
        sumDervExpo = getExponentDerivatives derv

--  McMURCHIE -DAVIDSON scheme of the primitive gaussians       
mcMurchie :: [Funtype] -> [NucCoord]  -> (Exponent,Exponent) -> CartesianDerivatives -> Double
mcMurchie !shells [ra,rb,rc] (e1,e2) derv =
  let coeff = U.unfoldr (calcHermCoeff [rpa,rpb] gamma) seedC 
      rtuv  = U.unfoldr (calcHermIntegral rpc gamma) seedI
      gamma = e1 + e2
      nu = e1*e2/gamma
      rp = meanp (e1,e2) ra rb
      [rab,rpa,rpb,rpc] = map restVect [(ra,rb),(rp,ra),(rp,rb),(rp,rc)]
      seedC = initilized_Seed_Coeff shells rab nu
      seedI = initilized_Seed_Integral shells rpc gamma derv
      
  in  U.sum $ U.zipWith (*) coeff rtuv
        
--  | Hermite Coefficients calculation using the maybe monad  
calcHermCoeff :: [NucCoord] -> Double -> HermiteStateCoeff-> Maybe (Double,HermiteStateCoeff)
calcHermCoeff rpab gamma stCoeff = do
    ls <- safeHead listC
    let  (val,newMap) = updateMap ls
         newSt = stCoeff {getmapC = newMap,getkeyC = tail listC}
    return (val,newSt)

  where listC = getkeyC stCoeff
        oldmap = getmapC stCoeff
        updateMap [xpa,ypa,zpa] = runState (hermEij xpa 0 >>= \e1 ->
                            hermEij ypa 1 >>= \e2 ->
                            hermEij zpa 2 >>= \e3 ->
                            return $ e1*e2*e3 ) oldmap

        hermEij k label = get >>= \st ->
                            let xpab = map (\r -> r !! label) rpab
                                ijt = getijt k
                                Just (x,newMap) = (lookupM k st) `mplus` Just (recursiveHC xpab gamma k st ijt)
                          in put newMap >> return x
        
-- | Hermite analytical integrals       
calcHermIntegral :: NucCoord  -> Double -> HermiteStateIntegral -> Maybe (Double,HermiteStateIntegral)
calcHermIntegral rpc gamma stIntegral =  do
    hi <- safeHead listI
    let Just (val,newMap) = fun hi
        newSt = stIntegral {getmapI = newMap,getkeyI = tail listI}
    return (val,newSt)

  where listI = getkeyI stIntegral
        oldmap = getmapI stIntegral
        fun hi = (lookupM hi oldmap) `mplus` Just (recursiveHI rpc gamma hi oldmap (getijt hi))

-- |Recursive Hermite Coefficients
recursiveHC :: [Double] -> Double -> HermiteIndex ->  MapHermite -> [Int] -> (Double, MapHermite)
recursiveHC pab@[xpa,xpb] gamma hi mapC [i,j,t]

       |(i+j) < t = (0.0,mapC)

       |t == 0 && i>=1 = let ([cij0,cij1],mapO) = updateMap [1,0,0] [1,0,(-1)]
                             val = xpa*cij0 + cij1
                             newMap = M.insert hi val mapO
                         in (val,newMap)

       |t == 0 && j>=1 = let ([cij0,cij1],mapO) = updateMap [0,1,0] [0,1,(-1)]
                             val = xpb*cij0 + cij1
                             newMap = M.insert hi val mapO
                         in (val,newMap)

       |otherwise =      let ([aijt,bijt],mapO) = updateMap [1,0,1] [0,1,1]
                             [i',j',t'] = map fromIntegral [i,j,t]
                             val = recip (2.0*gamma*t') * (i'*aijt + j'*bijt)
                             newMap = M.insert hi val mapO
                         in (val,newMap)

  where updateMap !xs !ys = runState (sequence [calcVal xs, calcVal ys]) mapC
        calcVal !xs = get >>= \st -> let key = sub xs
                                         funAlternative = recursiveHC pab gamma key st
                                         (v,m) = boyMonplus funAlternative key st
                                     in put m >> return v
        sub [a,b,c] = hi {getijt = [i-a,j-b,t-c]}
        
-- |Recursive Hermite Integrals        
recursiveHI :: NucCoord -> Double -> HermiteIndex-> MapHermite -> [Int] -> (Double,MapHermite)
recursiveHI rpc gamma hi mapI [t,u,v]

 | any (<0) [t,u,v] = (0.0,mapI)

 | t >= 1 = let ([boy1,boy2],mapO) = updateMap [2,0,0] [1,0,0]
                val = (fromIntegral t -1)*boy1 + (rpc !! 0)*boy2
                newM = M.insert hi val mapI
            in (val,newM)

 | u >= 1 = let ([boy1,boy2],mapO) = updateMap [0,2,0] [0,1,0]
                val = (fromIntegral u -1)*boy1 + (rpc !! 1)*boy2
                newM = M.insert hi val mapI
            in (val,newM)

 | v >= 1 = let ([boy1,boy2],mapO) = updateMap [0,0,2] [0,0,1]
                val = (fromIntegral v -1)*boy1 + (rpc !! 2)*boy2
                newM = M.insert hi val mapI
            in (val,newM)

 | otherwise =  let arg = (gamma*) . sum . (map (^2))  $ rpc
                    x1 = (-2.0*gamma)^n
                    val = x1 * boysF (fromIntegral n) arg
                    newM = M.insert hi val mapI
                in (val,newM)

  where updateMap !xs !ys = runState (sequence [calcVal xs, calcVal ys]) mapI
        calcVal xs = get >>= \st -> let key = sub xs
                                        funAlternative = recursiveHI rpc gamma key st
                                        (v,m) = boyMonplus funAlternative key st
                                    in put m >> return v
        n = getN hi -- order of the boy function
        sub [a,b,c] = hi {getN = n + 1, getijt = [t-a,u-b,v-c]}
        
-- |Boys Function calculation using the monad plus        
boyMonplus :: ([Int] -> (Double,MapHermite)) -> HermiteIndex-> MapHermite -> (Double,MapHermite)
boyMonplus fun key m = fromMaybe err $ (lookupM key m) `mplus` (Just (res,newM))
  where (res,oldmap) = fun (getijt key)
        newM = M.insert key res oldmap
        err = error "failure in the recursive Hermite Procedure" 
 
lookupM :: Ord k => k -> M.Map k a -> Maybe (a , M.Map k a)
lookupM k m = do 
              val <- M.lookup k m
              return (val,m)                                                             
              
{- we sent to listIndexes two Funtypes and the function return
the exponent for x,y,z according to the l-number of the
orbital that thos Funtypes represent
therefore for ls = [l1,m1,n1,l2,m2,n2]-}

listIndexes :: [Funtype] -> [Int]        
listIndexes shells = funtyp2Index <$> shells <*> [0..2]


initilized_Seed_Coeff :: [Funtype] -> NucCoord -> Double -> HermiteStateCoeff        
initilized_Seed_Coeff shells rab mu = HermiteStateCoeff mapC0 listC
  where mapC0 = L.foldl' (\acc (k, v) -> M.insert k v acc) M.empty $ L.zip k0 e0 
        listC = genCoeff_Hermite shells
        k0 = [Xpa [0, 0, 0], Ypa [0, 0, 0], Zpa [0, 0, 0]] 
        e0 = map (\x -> exp(-mu*(x^2))) rab
        
        
genCoeff_Hermite :: [Funtype] -> [[HermiteIndex]]
genCoeff_Hermite shells = do
  i <-[0..l1+l2]
  j <-[0..m1+m2]
  k <-[0..n1+n2]
  return $ [Xpa [l1,l2,i], Ypa [m1,m2,j], Zpa [n1,n2,k]] 
  where [l1,m1,n1,l2,m2,n2] = listIndexes shells


initilized_Seed_Integral :: [Funtype] -> NucCoord -> Double -> CartesianDerivatives -> HermiteStateIntegral
initilized_Seed_Integral !symbols !rpc !gamma !derv = HermiteStateIntegral mapI0 listI 
  where mapI0 = M.insert k0' f0 M.empty
        k0'= Rpa 0 [0, 0, 0] 
        y= gamma*rpc2
        f0 = boysF 0 y
        rpc2 = sum $ map (^2) rpc
        listI = genCoeff_Integral symbols derv
        
                
genCoeff_Integral :: [Funtype] -> CartesianDerivatives -> [HermiteIndex]
genCoeff_Integral symbols (Dij_Ax e, Dij_Ay f, Dij_Az g) =
  do
  t <-[0..l1+l2+e]
  u <-[0..m1+m2+f]
  v <-[0..n1+n2+g]
  return $ Rpa 0 [t,u,v] 
  where [l1,m1,n1,l2,m2,n2] = listIndexes symbols        

-- ==============> Auxiliar Functions <===============


--(2k -1) !! = (2k)!/(2^k * k!)
-- i = 2k - 1 => k = (i + 1)/ 2
-- | Odd factorial function
facOdd ::Int -> Double
facOdd  i | i `rem`2 == 0  = error "Factorial Odd function required an odd integer as input"
          | otherwise  = case compare i 2 of
                             LT -> 1
                             GT-> let k = (1 + i) `div ` 2
                                  in (fromIntegral $ fac (2*k)) /  (2.0^k * (fromIntegral $ fac k))

-- | Factorial function
fac :: Int -> Int
fac i | i < 0   = error "The factorial function is defined only for natural numbers"
      | i == 0  = 1
      | otherwise = product [1..i]



binomial :: Int -> Int -> Double
binomial l k = fromIntegral $ fac l `div` (fac k * fac (l-k))

-- | Square internuclear distance function
rab2 :: NucCoord -> NucCoord -> Double
rab2 a b = sum . map (^2). L.zipWith (-) a $ b

-- | Mean point between two gaussians
meanp ::(Exponent,Exponent) -> NucCoord -> NucCoord -> [Double]
meanp (e1,e2) ra rb = map (\(a,b) -> (e1*a + e2*b)/(e1+e2)) $ L.zip ra rb

restVect :: (NucCoord,NucCoord) -> NucCoord
restVect (ra,rb) = L.zipWith (-) ra rb

{- |The norm of each gaussian is given by the following equation
    N = sqrt $ ((2l -1)!! (2m-1)!! (2n-1)!!)/(4*expo)^(l+m+n)  * (pi/(2*e))**1.5
    where expo is the exponential factor of the GTO -}
normaCoeff :: CGF -> CGF
normaCoeff b1 = b1 { getPrimitives = newPrimitives}
  where xs = getPrimitives b1
        newPrimitives = map ((uncurry fun ) &&& (snd) ) xs
        fun = \c e -> (c/) . sqrt $ (ang e) * (pi/(2*e))**1.5
        ang x = prod / (4*x)^(sum indexes)
        prod = product $ map (\k -> facOdd (2*k -1)) indexes
        shell = getfunTyp  b1
        indexes = map (LA.map2val mapLAngular) $ L.zip (repeat shell) [0..2]


-- | Transform from the unary angular momentum representation to the corresponding integer value
funtyp2Index :: Funtype -> Int -> Int
funtyp2Index funtyp x =  fromMaybe (error " Unknown Label of the Basis" )
                         $ M.lookup (funtyp,x) mapLAngular

deltaij :: Eq a => a -> a -> Double -> Double
deltaij i j val = if i == j then val else 0  





