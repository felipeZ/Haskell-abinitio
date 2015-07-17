{-# Language BangPatterns, ViewPatterns  #-}

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
import Control.Monad.Memo (Memo, memo, runMemo)
import Control.Monad.State
import Control.Parallel.Strategies (parMap,rdeepseq)
import Data.Array.Repa         as R hiding (map)
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix (mmultP)
import Data.List as L
import qualified Data.Map.Strict as M
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
{- INLINE hcore -}
             
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
          guard $ schwarz gauss 
          return $ twoElectronHermite gauss

schwarz :: [Gauss] -> Bool 
schwarz gs@[g1,g2,g3,g4] =  if val < 1e-8 then False else True 
 
 where val = qab * qcd  
       qab = sqrt $ twoTermsERI g1 g2
       qcd = sqrt $ twoTermsERI g3 g4

-- | integral of the form (ab|ab) only depend on the expansion coefficient of the
-- | hermite expanstion because the derivatives of the boys function are equal to 0
twoTermsERI :: Gauss -> Gauss -> Double
twoTermsERI ga gb =  (*cte) . U.sum $ coeff `deepseq` U.map (*suma) coeff 

  where suma              = U.sum $ U.zipWith (*) sgns coeff
        gs                = [ga,gb]
        coeff             = U.unfoldr (calcHermCoeff [rpa,rpb] p) seedC       
        tuvs              = map getijt $ genCoeff_Integral [symb1,symb2] derv
        seedC             = initilized_Seed_Coeff [symb1,symb2] rab mu
        sgns              = U.fromList . map (\xs-> (-1.0)^(sum xs)) $ tuvs
        [ra,rb]           = map nucCoord gs
        [symb1,symb2]     = map funtype  gs
        ([c1,c2],[e1,e2]) = (map fst ) &&& (map snd ) $ map gaussP gs
        p                 = e1 + e2
        rp                = meanp (e1,e2) ra rb
        [rab,rpa,rpb]     = map restVect  (zip [ra,rp,rp] [rb,ra,rb])
        mu                = e1*e2* (recip $ e1 + e2)
        cte               = (c1^2 * c2^2 ) * (sqrt $ (pi/(2*p))^5)
        derv              = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)


twoElectronHermite :: [Gauss] -> Double
twoElectronHermite gs = (cte *) . U.sum . U.zipWith (*) coeff1 $!! U.fromList
                        [mcMurchie2 coeff2 rpq alpha abcs tuv | tuv <- coefftuv]

  where coefftuv = genCoeff_Integral [symb1,symb2] derv
        abcs     = genCoeff_Integral [symb3,symb4] derv
        coeff1   = U.unfoldr (calcHermCoeff [rpa,rpb] p) seedC1
        coeff2   = U.unfoldr (calcHermCoeff [rqc,rqd] q) seedC2
        seedC1   = initilized_Seed_Coeff [symb1,symb2] rab mu
        seedC2   = initilized_Seed_Coeff [symb3,symb4] rcd nu
        ps      = map gaussP gs
        rp      = meanp (e1,e2) ra rb
        rq      = meanp (e3,e4) rc rd
        [mu,nu] = (\(a,b) -> a*b* (recip $ a + b)) `fmap` [(e1,e2),(e3,e4)]
        [p,q]   = uncurry (+) `fmap` [(e1,e2),(e3,e4)]
        alpha   = p*q/(p+q)
        cte     = (c1*c2*c3*c4*) . (*(2.0*pi**2.5)) . recip $ (p * q ) * (sqrt $ p + q)
        [ra,rb,rc,rd]                 = map nucCoord gs
        [symb1,symb2,symb3,symb4]     = map funtype  gs
        [rab,rcd,rpa,rpb,rqc,rqd,rpq] = map restVect  (zip [ra,rc,rp,rp,rq,rq,rp][rb,rd,ra,rb,rc,rd,rq])
        ([c1,c2,c3,c4],[e1,e2,e3,e4]) = (map fst ) &&& (map snd ) $ ps
        derv = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)


mcMurchie2 :: VecUnbox -> NucCoord -> Double -> [HermiteIndex] -> HermiteIndex -> Double
mcMurchie2 !coeff2 rpq alpha abcs' tuv' = U.sum $ integrals `deepseq` U.zipWith3 (\x y z -> x*y*z) sgns coeff2  integrals
  where tuv       = getijt tuv'
        abcs      = map getijt  abcs'
        sgns      = U.fromList $ map (\xs-> (-1.0)^(sum xs)) abcs
        integrals = U.unfoldr (calcHermIntegral rpq alpha) seedI
        seedI     = HermiteStateIntegral mapI0 listI
        mapI0     = M.insert k0' f0 M.empty
        k0'       = Rpa 0 [0, 0, 0]
        rp2       = sum $ map (^2) rpq
        y         = alpha*rp2
        f0        = boysF 0 y
        listI     = Rpa 0 `fmap` (map (L.zipWith (+) tuv) abcs )



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
mcMurchie !shells [ra,rb,rc] (e1,e2) derv =  U.sum $ coeff `deepseq` rtuv `deepseq` U.zipWith (*) coeff rtuv
 where coeff = U.unfoldr (calcHermCoeff [rpa,rpb] gamma) seedC 
       rtuv  = U.unfoldr (calcHermIntegral rpc gamma) seedI
       gamma = e1 + e2
       nu    = e1*e2/gamma
       rp    = meanp (e1,e2) ra rb
       [rab,rpa,rpb,rpc] = map restVect [(ra,rb),(rp,ra),(rp,rb),(rp,rc)]
       seedC = initilized_Seed_Coeff shells rab nu
       seedI = initilized_Seed_Integral shells rpc gamma derv
      
        
-- | Hermite analytical integrals       
calcHermIntegral ::  NucCoord -> Double -> HermiteStateIntegral -> Maybe (Double,HermiteStateIntegral)
calcHermIntegral rpa alpha stIntegral = do 
    hi <- safeHead listI
    let (val,newMap) = runMemo (calcHermM rpa alpha hi) oldmap
        newSt = stIntegral {getmapI = newMap,getkeyI = tail listI}
    return (val,newSt) 

      where listI = getkeyI stIntegral
            oldmap = getmapI stIntegral
        
-- |Recursive Hermite Integrals        
calcHermM :: NucCoord  -> Double -> HermiteIndex ->  Memo HermiteIndex Double Double
calcHermM rpq@[x,y,z] alpha (Rpa n [t,u,v]) 
     | any (<0) [t,u,v] = return $ 0
   
     | t >= 1 = do  boy1 <- memo (calcHermM rpq alpha) $ Rpa (n+1) [t-2,u,v]
                    boy2 <- memo (calcHermM rpq alpha) $ Rpa (n+1) [t-1,u,v] 
                    return $ (fromIntegral t -1)*boy1 + x*boy2

     | u >= 1 = do  boy1 <- memo (calcHermM rpq alpha) $ Rpa (n+1) [t,u-2,v]
                    boy2 <- memo (calcHermM rpq alpha) $ Rpa (n+1) [t,u-1,v] 
                    return $ (fromIntegral u -1)*boy1 + y*boy2

     | v >= 1 = do  boy1 <- memo (calcHermM rpq alpha) $ Rpa (n+1) [t,u,v-2]
                    boy2 <- memo (calcHermM rpq alpha) $ Rpa (n+1) [t,u,v-1] 
                    return $ (fromIntegral v -1)*boy1 + z*boy2

     | otherwise = do let arg = (alpha*) . sum . (map (^2))  $ rpq 
                      return $ (-2.0*alpha)^n * boysF (fromIntegral n) arg 


calcHermCoeff2 :: [NucCoord] -> Double -> ([Double],[[HermiteIndex]])-> [Double]
calcHermCoeff2 pab@[xspa,xspb] gamma (es0,hss) = map (calcCoeff es0) hss

 where calcCoeff es0 hs = product $ zipWith4 (recursiveHermCoeff2 gamma) xspa xspb es0 hs 

recursiveHermCoeff2 :: Double -> Double -> Double -> Double -> HermiteIndex -> Double
recursiveHermCoeff2 gamma xpa xpb e0 hs@(getijt -> ijt@[i,j,t])

       | all (==0) ijt     = e0
  
       |(i+j) < t && t < 0 = 0 

       |t == 0 && i>=1 = let cij0 = fun $ hs {getijt = [i-1,j,0]}
                             cij1 = fun $ hs {getijt = [i-1,j,1]}
                          in xpa*cij0 + cij1

       |t == 0 && j>=1 = let cij0 = fun $ hs {getijt = [i,j-1,0]}
                             cij1 = fun $ hs {getijt = [i,j-1,1]}
                         in xpb*cij0 + cij1

       |otherwise =      let aijt = fun $ hs {getijt = [i-1,j,t-1]}
                             bijt = fun $ hs {getijt = [i,j-1,t-1]}
                             [i',j',t'] = map fromIntegral [i,j,t]
                         in recip (2.0*gamma*t') * (i'*aijt + j'*bijt)

 where fun = recursiveHermCoeff2 gamma xpa xpb e0
  
-- | Hermite Coefficients calculation using the maybe monad  
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

initilizedSeedCoeff :: [Funtype] -> NucCoord -> Double -> ([Double],[[HermiteIndex]])        
initilizedSeedCoeff shells rab mu =  (e0,listC)
  where listC = genCoeff_Hermite shells
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





