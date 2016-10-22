{-# Language BangPatterns, ViewPatterns  #-}

{-|
Module: Science.QuantumChemistry.Integrals.IntegralsEvaluation 
Description: Analytical evaluation of the overlap, core, nuclei-electron
             and electron-electron integrals.
Command line checkers
Copyright: @2013 Felipe Zapata, Angel Alvarez
           @2016 Felipe Zapata
The HaskellFock SCF Project
-}


module Science.QuantumChemistry.Integrals.IntegralsEvaluation 
  -- (
  --  (|>)
  -- ,(<|)
  -- ,(|>>)
  -- ,(<<|)
  -- ,(<<||>>)
  -- ,Operator(..)
  -- ,contracted4Centers
  -- ,evalIntbykey
  -- ,hcore
  -- ,mtxOverlap
  -- ,normaCoeff
  -- )
  where

-- =============================> Standard and third party libraries <===============================
import Control.Applicative
import Control.Arrow ((&&&),first,second)
import Control.DeepSeq
import Control.Monad (liftM,mplus)
import Control.Monad.List
import Control.Monad.Memo (Memo, memo, runMemo)
import Control.Monad.State
import Control.Parallel.Strategies (parMap,rdeepseq)
import Data.Array.Repa         as R hiding (map)
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix (mmultP)
import Data.Foldable 
import Data.List as L hiding (all,any,product, sum)
import qualified Data.Map.Lazy   as ML
import qualified Data.Map.Strict as M
import Data.Maybe (fromMaybe)
import qualified Data.Vector.Unboxed as U
import Prelude hiding (all,any,product,sum)


-- =================> Internal Modules <======================
import Science.QuantumChemistry.NumericalTools.Boys(boysF)
import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.NumericalTools.LinearAlgebra as LA
import Science.QuantumChemistry.NumericalTools.TableBoys (Boys,boysTaylor)

--  ======================= > MODULE FOR EVALUATING THE ONE AND TWO ELECTRON TYPE INTEGRALS < =============================
{- | This module is based on chapter 9 "Molecular Integral Evaluation" in the book
   "Molecular Electronic-Structure Theory". Authors: Trygve Helgaker, Poul Jorgensen and Jeppe Olsen

   | The Gaussian primitives are not normalized they must be multiplied for a factor sqr (4*alfa/pi)
   | where alfa is the coefficient of the gaussian function-}


-- =================================================================================================
-- == the Monad List are used to calculate the elements of the lower triangular matrix            ==
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
                -- Vij rc derivatives -> vijContracted b1 rc b2 derivatives
                                                    
-- ==============> 2-CENTER OVERLAP INTEGRAL <=================================

-- | overlaping between two primitive gaussian function between arbitrary-l functions  <A|B>

-- S12 = *exp[−a1*b1(Rab)^2/gamma]IxIyIz

-- where Ix = sum f2i (l1,l2,PAx,PBx) ((2i-1)!! / (2gamma)^i sqrt (pi/gamma) for i = 0,l1 + l2 / 2 only the even i
-- f2k is defined above in the module

mtxOverlap :: Monad m => [AtomData] -> m (Array U DIM1 Double) 
mtxOverlap atoms = computeUnboxedP . fromFunction (Z:.dim) $
 (\idx -> let (Z:.i:.j) = LA.indexFlat2DIM2 norbital idx
              [(r1,cgf1),(r2,cgf2)] = calcIndex <$> [i,j]
          in sijContracted r1 r2 cgf1 cgf2)
             
  where norbital = sum . fmap (length . getBasis) $ atoms
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
                              let [l1,l2] = fmap getfunTyp [cgf1,cgf2]
                                  gauss1 = Gauss r1 l1 g1
                                  gauss2 = Gauss r2 l2 g2
                              return (gauss1 <||> gauss2 )


-- | Primitive overlap terms                              
sab :: Gauss -> Gauss -> Double
sab g1@(Gauss r1 shellA (c1,e1)) g2@(Gauss r2 shellB (c2,e2)) =
   (c1*c2) * (product $!! [obaraSaika gamma (s00 !! x) (pa !! x) (pb !! x) (l1 x) (l2 x) |x <- [0..2]])
  
  where [l1,l2] = fmap funtyp2Index [shellA,shellB]
        [pa,pb] = fmap (L.zipWith (-) p) [r1,r2]
        p       = meanp (e1,e2) r1 r2
        expo x2 = exp $ -mu * x2
        s00     = ((*cte) . expo . (^2)) <$>  L.zipWith (-) r1 r2
        cte     = sqrt $ pi/ gamma
        gamma   = e1+ e2
        mu      = e1*e2/gamma

-- | ObaraSaika Scheme to calculate overlap integrals
obaraSaika ::Double -> Double -> Double -> Double -> Int -> Int  -> Double
obaraSaika gamma s00 pax pbx = s 

  where pred2 = pred . pred
        c     = recip $ 2.0 * gamma
        s m n | m < 0 || n < 0   =  0
              | all (== 0) [m,n] =  s00
              | m == 0    =  pbx * s 0 (pred n) + c* n_1 * s 0 (pred2 n)
              | n == 0    =  pax * s (pred m) 0 + c* m_1 * s (pred2 m) 0
              | otherwise =  pax * s (pred m) n + c * (m_1 * (s (pred2 m) n) +
                               (fromIntegral n) * s (pred m) (pred n))

          where [m_1,n_1] = fmap (fromIntegral . pred) [m,n]

                
-- ====================> HAMILTONIAN CORE <=======================

-- | Core Hamiltonian Matrix Calculation
hcore :: Monad m => ML.Map Boys Double -> [AtomData] -> m (Array U DIM1 Double)
hcore gridBoys atoms  = computeUnboxedP . fromFunction (Z:. dim) $
 (\idx -> let (Z:.i:.j)     = LA.indexFlat2DIM2 norbital idx
              [atomi,atomj] =  calcIndex <$> [i,j]
              derv          = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)
              sumVij        = sum $ L.zipWith (\z rc -> ((-z) * vijContracted gridBoys atomi rc atomj derv ) ) atomicZ coords              
             in (atomi <<|Tij|>> atomj) + sumVij)

  where coords    = fmap getCoord atoms
        atomicZ   = fmap getZnumber atoms
        norbital  = sum . fmap (length . getBasis) $ atoms
        dim       = (norbital^2 + norbital) `div`2
        calcIndex = LA.calcCoordCGF atoms
{-# INLINE hcore #-}
             
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
              let [l1,l2] = fmap getfunTyp [cgf1,cgf2]
                  gauss1 = Gauss r1 l1 g1 
                  gauss2 = Gauss r2 l2 g2
              return ( gauss1 <| Tij |> gauss2 )  

-- | Primitive kinetic energy terms          
tab :: Gauss -> Gauss -> Double
tab gA@(Gauss r1 shell1 (c1,e1)) gB@(Gauss !r2 !shell2 (!c2,!e2)) =
   c1*c2 * (sum . fmap product $! permute [\x -> tx x (j x) (k x),\x -> sx x (j x) (k x),\x -> sx x (j x) (k x)] [0..2])

  where [j,k] = fmap funtyp2Index [shell1,shell2]
        sx i  = obaraSaika gamma (s00 !! i) (pa !! i) (pb !! i) 
        tx i lang1 lang2 = let [l1,l2] = fmap fromIntegral [lang1,lang2]
                           in -2.0 * e2^2 *  sx i lang1 (lang2+2) +
                              e2*(2*l2 +1)*  sx i lang1 lang2 -
                              0.5*l2*(l2-1)* sx i lang1 (lang2 -2)
        cte     = sqrt (pi/ gamma)
        s00     = ((*cte) . expo . (^2)) <$>  L.zipWith (-) r1 r2
        t00 x   = e1 - 2*e1^2 *((pa !! x)^2 + recip (2*gamma)) * (s00 !! x )
        expo x2 = exp $ -mu * x2
        [pa,pb] = (L.zipWith (-) p) <$>  [r1,r2]
        p       = meanp (e1,e2) r1 r2
        gamma   = e1+ e2
        mu      = e1*e2/gamma        
        
permute :: [a -> b] -> [a] -> [[b]]
permute [f,g,h] xs = fmap (\ x -> L.zipWith ($) x xs) [[f, g, h], [g, f, h], [g, h, f]]

-- =====================> TWO ELECTRON INTEGRALS <=========================================

-- | Main function to calculate the Coulomb (J) and Interchange Integrals using the
--  indexes of the atomic basis <ab|cd>
evalIntbykey :: ML.Map Boys Double -> [AtomData] -> [Int] ->  Double
evalIntbykey gridBoys atoms keys = contracted4Centers gridBoys centers 
 where centers = LA.calcCoordCGF atoms <$> keys

-- | Calculate electronic interaction among the four center contracted Gaussian functions
contracted4Centers :: ML.Map Boys Double -> [(NucCoord,CGF)] -> Double
contracted4Centers gridBoys [(ra,cgf1), (rb,cgf2), (rc,cgf3), (rd,cgf4)] = sum cartProd
  where [l1,l2,l3,l4] = fmap getfunTyp [cgf1,cgf2,cgf3,cgf4]
        cartProd = do
          g1 <- getPrimitives cgf1
          g2 <- getPrimitives cgf2
          g3 <- getPrimitives cgf3
          g4 <- getPrimitives cgf4
          let gauss = L.zipWith3 Gauss [ra,rb,rc,rd] [l1,l2,l3,l4] [g1,g2,g3,g4]
              [f1,f2,f3,f4] = fmap funtype gauss
              [ga,gb,gc,gd] = gauss
          guard $ schwarz gauss
          return $ twoElectronHermite gridBoys gauss

schwarz :: [Gauss] -> Bool 
schwarz gs@[g1,g2,g3,g4] =  not (val < 1e-8) 
 
 where val = qab * qcd  
       qab = sqrt $ twoTermsERI g1 g2
       qcd = if g1== g3 && g2 == g4
                then qab
                else sqrt $ twoTermsERI g3 g4

-- | integral of the form (ab|ab) only depend on the expansion coefficient of the
-- | hermite expanstion because the derivatives of the boys function are equal to 0
twoTermsERI :: Gauss -> Gauss -> Double
twoTermsERI ga gb =  (*cte) . U.sum $ coeff `deepseq` U.map (*suma) coeff 

  where suma              = U.sum $ U.zipWith (*) sgns coeff
        gs                = [ga,gb]
        coeff             = calcHermCoeff [rpa,rpb] p seedC          
        tuvs              = getijt <$> genCoeff_Integral [symb1,symb2] derv
        seedC             = initilizedSeedCoeff [symb1,symb2] rab mu        
        sgns              = U.fromList . fmap (\xs-> (-1.0)^sum xs) $ tuvs
        [ra,rb]           = fmap nucCoord gs
        [symb1,symb2]     = fmap funtype  gs
        ([c1,c2],[e1,e2]) = fmap fst &&& fmap snd $ fmap gaussP gs
        p                 = e1 + e2
        rp                = meanp (e1,e2) ra rb
        [rab,rpa,rpb]     = fmap restVect  (zip [ra,rp,rp] [rb,ra,rb])
        mu                = e1*e2* recip (e1 + e2)
        cte               = 2 * (c1^2 * c2^2 ) * sqrt ((pi/(2*p))^5)
        -- cte               = 2*(sqrt $ (pi/(2*p))^5)
        derv              = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)


twoElectronHermite :: ML.Map Boys Double -> [Gauss] -> Double
twoElectronHermite gridBoys gs = (cte *) . U.sum . U.zipWith (*) coeff1 $!! U.fromList
                        [mcMurchie2 gridBoys coeff2 rpq alpha abcs tuv | tuv <- coefftuv]

  where coefftuv = genCoeff_Integral [symb1,symb2] derv
        abcs     = genCoeff_Integral [symb3,symb4] derv
        coeff1   = calcHermCoeff [rpa,rpb] p seedC1
        coeff2   = calcHermCoeff [rqc,rqd] q seedC2
        seedC1   = initilizedSeedCoeff [symb1,symb2] rab mu
        seedC2   = initilizedSeedCoeff [symb3,symb4] rcd nu
        ps      = fmap gaussP gs
        rp      = meanp (e1,e2) ra rb
        rq      = meanp (e3,e4) rc rd
        [mu,nu] = (\(a,b) -> a*b* recip (a + b)) <$> [(e1,e2),(e3,e4)]
        [p,q]   = uncurry (+) <$> [(e1,e2),(e3,e4)]
        alpha   = p*q/(p+q)
        cte     = (c1*c2*c3*c4*) . (*(2.0*pi**2.5)) . recip $ (p * q ) * sqrt (p + q)
        [ra,rb,rc,rd]                 = fmap nucCoord gs
        [symb1,symb2,symb3,symb4]     = fmap funtype  gs
        [rab,rcd,rpa,rpb,rqc,rqd,rpq] = fmap restVect  (zip [ra,rc,rp,rp,rq,rq,rp][rb,rd,ra,rb,rc,rd,rq])
        ([c1,c2,c3,c4],[e1,e2,e3,e4]) = fmap fst  &&& fmap snd $ ps
        derv = (Dij_Ax 0, Dij_Ay 0, Dij_Az 0)


mcMurchie2 :: ML.Map Boys Double -> VecUnbox -> NucCoord -> Double -> [HermiteIndex] -> HermiteIndex -> Double
mcMurchie2 gridBoys coeff2 rpq alpha abcs' tuv' = U.sum $ integrals `deepseq` U.zipWith3 (\x y z -> x*y*z) sgns coeff2  integrals
  where tuv       = getijt tuv'
        abcs      = fmap getijt  abcs'
        sgns      = U.fromList $ fmap (\xs-> (-1.0)^sum xs) abcs
        integrals = U.unfoldr (calcHermIntegral gridBoys rpq alpha) seedI
        seedI     = HermiteStateIntegral mapI0 listI
        mapI0     = M.insert k0' f0 M.empty
        k0'       = Rpa 0 [0, 0, 0]
        rp2       = sum $ fmap (^2) rpq
        y         = alpha*rp2
        f0        = boysF 0 y
        listI     = fmap (Rpa 0 . L.zipWith (+) tuv) abcs


-- ======================> McMURCHIE -DAVIDSON SCHEME <=========================
-- | The integrals of three and four centers are implemented according to chapter 9
--   of the book "Molecular Electronic-Structure Theory" by Trygve Helgaker,
--   Poul Jorgensen and Jeppe Olsen.

-- | Nuclei-Electron matrix operator presentation        
vijContracted :: ML.Map Boys Double -> (NucCoord,CGF) -> NucCoord -> (NucCoord,CGF) -> CartesianDerivatives -> Double
vijContracted gridBoys (ra,cgf1) rc (rb,cgf2) derv = sum $
        do g1 <- getPrimitives cgf1
           g2 <- getPrimitives cgf2
           let [l1,l2] = fmap getfunTyp [cgf1,cgf2]
               gauss1 =  Gauss ra l1 g1
               gauss2 =  Gauss rb l2 g2
           return $ vijHermite gridBoys gauss1 gauss2 rc derv
           
-- | Hermite auxiliar function calculations according to the McMURCHIE -DAVIDSON scheme        
vijHermite :: ML.Map Boys Double -> Gauss -> Gauss -> NucCoord -> CartesianDerivatives -> Double
vijHermite gridBoys g1 g2 rc derv = ((-1)^sumDervExpo) * cte * mcMurchie gridBoys shells [ra,rb,rc] (e1,e2) derv
  where cte = c1 * c2 * 2.0 * (pi/gamma) 
        gamma = e1+e2     
        [ra,rb] =  fmap nucCoord [g1,g2]
        shells = fmap funtype [g1,g2]
        [(c1,e1),(c2,e2)] = fmap gaussP [g1,g2]
        sumDervExpo = getExponentDerivatives derv

-- |  McMURCHIE -DAVIDSON scheme of the primitive gaussians       
mcMurchie :: ML.Map Boys Double -> [Funtype] -> [NucCoord]  -> (Exponent,Exponent) -> CartesianDerivatives -> Double
mcMurchie gridBoys shells [ra,rb,rc] (e1,e2) derv =  U.sum $ coeff `deepseq` rtuv `deepseq` U.zipWith (*) coeff rtuv
 where coeff = calcHermCoeff [rpa,rpb] gamma seedC 
       rtuv  = U.unfoldr (calcHermIntegral gridBoys rpc gamma) seedI
       gamma = e1 + e2
       nu    = e1*e2/gamma
       rp    = meanp (e1,e2) ra rb
       [rab,rpa,rpb,rpc] = fmap restVect [(ra,rb),(rp,ra),(rp,rb),(rp,rc)]
       seedC = initilizedSeedCoeff shells rab nu       
       seedI = initilized_Seed_Integral shells rpc gamma derv
      
        
-- | Hermite analytical integrals       
calcHermIntegral :: ML.Map Boys Double -> NucCoord -> Double -> HermiteStateIntegral -> Maybe (Double,HermiteStateIntegral)
calcHermIntegral gridBoys rpa alpha stIntegral = do 
    hi <- safeHead listI
    let (val,newMap) = runMemo (calcHermM gridBoys rpa alpha hi) oldmap
        newSt = stIntegral {getmapI = newMap,getkeyI = tail listI}
    return (val,newSt) 

      where listI    = getkeyI stIntegral
            oldmap   = getmapI stIntegral

        
-- |Recursive Hermite Integrals        
calcHermM :: ML.Map Boys Double -> NucCoord  -> Double -> HermiteIndex ->  Memo HermiteIndex Double Double
calcHermM gridBoys rpq@[x,y,z] alpha (Rpa n [t,u,v]) 
     | any (<0) [t,u,v] = return 0
   
     | t >= 1 = do  boy1 <- memo funHermM $ Rpa (n+1) [t-2,u,v]
                    boy2 <- memo funHermM $ Rpa (n+1) [t-1,u,v] 
                    return $ (fromIntegral t -1)*boy1 + x*boy2

     | u >= 1 = do  boy1 <- memo funHermM $ Rpa (n+1) [t,u-2,v]
                    boy2 <- memo funHermM $ Rpa (n+1) [t,u-1,v] 
                    return $ (fromIntegral u -1)*boy1 + y*boy2

     | v >= 1 = do  boy1 <- memo funHermM $ Rpa (n+1) [t,u,v-2]
                    boy2 <- memo funHermM $ Rpa (n+1) [t,u,v-1] 
                    return $ (fromIntegral v -1)*boy1 + z*boy2

     | otherwise = do let arg = (alpha*) . sum . fmap (^2)  $ rpq 
                      return $ (-2.0*alpha)^n * boysTaylor gridBoys (fromIntegral n) arg 

     where funHermM = calcHermM gridBoys rpq alpha

-- | Hermite Coefficients calculation  
calcHermCoeff :: [NucCoord] -> Double -> ([Double],[[HermiteIndex]])-> VecUnbox
calcHermCoeff pab@[xspa,xspb] gamma (es0,hss) = U.fromList $ fmap (calcCoeff es0) hss

 where calcCoeff es0 hs = product $ zipWith4 (recursiveHermCoeff gamma) xspa xspb es0 hs 

recursiveHermCoeff :: Double -> Double -> Double -> Double -> HermiteIndex -> Double
recursiveHermCoeff gamma xpa xpb e0 hs@(getijt -> ijt@[i,j,t])

       | all (==0) ijt  = e0
  
       | any (<0) ijt  = 0 

       |t == 0 && i>=1 = let cij0 = fun $ hs {getijt = [i-1,j,0]}
                             cij1 = fun $ hs {getijt = [i-1,j,1]}
                          in xpa*cij0 + cij1

       |t == 0 && j>=1 = let cij0 = fun $ hs {getijt = [i,j-1,0]}
                             cij1 = fun $ hs {getijt = [i,j-1,1]}
                         in xpb*cij0 + cij1

       |otherwise =      let aijt = fun $ hs {getijt = [i-1,j,t-1]}
                             bijt = fun $ hs {getijt = [i,j-1,t-1]}
                             [i',j',t'] = fmap fromIntegral [i,j,t]
                         in recip (2*gamma*t') * (i'*aijt + j'*bijt)

 where fun = recursiveHermCoeff gamma xpa xpb e0

              
{- we sent to listIndexes two Funtypes and the function return
the exponent for x,y,z according to the l-number of the
orbital that thos Funtypes represent
therefore for ls = [l1,m1,n1,l2,m2,n2]-}

listIndexes :: [Funtype] -> [Int]        
listIndexes shells = funtyp2Index <$> shells <*> [0..2]


initilizedSeedCoeff :: [Funtype] -> NucCoord -> Double -> ([Double],[[HermiteIndex]])        
initilizedSeedCoeff shells rab mu =  (e0,listC)
  where listC = genCoeff_Hermite shells
        k0 = [Xpa [0, 0, 0], Ypa [0, 0, 0], Zpa [0, 0, 0]] 
        e0 = fmap (\x -> exp(-mu*(x^2))) rab
        
        
genCoeff_Hermite :: [Funtype] -> [[HermiteIndex]]
genCoeff_Hermite shells = do
  i <-[0..l1+l2]
  j <-[0..m1+m2]
  k <-[0..n1+n2]
  return [Xpa [l1,l2,i], Ypa [m1,m2,j], Zpa [n1,n2,k]] 
  where [l1,m1,n1,l2,m2,n2] = listIndexes shells


initilized_Seed_Integral :: [Funtype] -> NucCoord -> Double -> CartesianDerivatives -> HermiteStateIntegral
initilized_Seed_Integral !symbols !rpc !gamma !derv = HermiteStateIntegral mapI0 listI 
  where mapI0 = M.insert k0' f0 M.empty
        k0'= Rpa 0 [0, 0, 0] 
        y= gamma*rpc2
        f0 = boysF 0 y
        rpc2 = sum $ fmap (^2) rpc
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


-- | Odd factorial function
-- | (2k -1) !! = (2k)!/(2^k * k!)
-- | where i = 2k - 1 => k = (i + 1)/ 2
facOdd ::Int -> Double
facOdd  i | even i  = error "Factorial Odd function required an odd integer as input"
          | otherwise  = case compare i 2 of
                             LT -> 1
                             GT-> let k = (1 + i) `div ` 2
                                  in fromIntegral (fac (2*k)) /  (2.0^k * fromIntegral (fac k))

-- | Factorial function
fac :: Int -> Int
fac i | i < 0   = error "The factorial function is defined only for natural numbers"
      | i == 0  = 1
      | otherwise = product [1..i]


-- | Square internuclear distance function
rab2 :: NucCoord -> NucCoord -> Double
rab2 a b = sum . fmap (^2). L.zipWith (-) a $ b

-- | Mean point between two gaussians
meanp ::(Exponent,Exponent) -> NucCoord -> NucCoord -> [Double]
meanp (e1,e2) ra rb = (\(a,b) -> (e1*a + e2*b)/(e1+e2)) <$> L.zip ra rb

restVect :: (NucCoord,NucCoord) -> NucCoord
restVect (ra,rb) = L.zipWith (-) ra rb


-- | Transform from the unary angular momentum representation to the corresponding integer value
funtyp2Index :: Funtype -> Int -> Int
funtyp2Index funtyp x =  fromMaybe (error " Unknown Label of the Basis" )
                         $ M.lookup (funtyp,x) mapLAngular

deltaij :: Eq a => a -> a -> Double -> Double
deltaij i j val = if i == j then val else 0  





