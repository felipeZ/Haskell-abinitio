{-# Language FlexibleContexts,ScopedTypeVariables,TypeOperators,PatternGuards #-}

module IntegralsEvaluation 
  (
  evalIntbykey
  ,hcore
  ,mtxOverlap
  ,normaCoeff
  )
  where

import Data.List as DL
import Control.Parallel.Strategies
import Data.Number.Erf
import Control.Monad.List
import Control.Monad.State
import qualified LinearAlgebra as LA
import qualified Data.Foldable as DF
import qualified Data.Traversable as DT
import qualified Data.Map as M
import qualified Control.Arrow as Arrow
import qualified Data.Vector.Unboxed as VU
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import GlobalTypes
import Control.Monad (liftM,mplus,sequence)
import Data.Maybe (fromMaybe)
import Control.Applicative (Applicative (..),ZipList(..),(<*>),(<$>),getZipList)
import Data.Monoid (Monoid(..),mappend,mconcat,mempty)
import Numeric.GSL.Integration (integrateQNG)
import Control.Arrow ((&&&))

--  ======================= > MODULE FOR EVALUATING THE ONE AND TWO ELECTRON TYPE INTEGRALS < =============================

-- This module is based on chapter 9 "Molecular Integral Evaluation" in the book 
-- "Molecular Electronic-Structure Theory". Authors: Trygve Helgaker, Poul Jorgensen and Jeppe Olsen

--Remember that the Gaussian primitives are not normalized they must be multiplied for a factor sqr (4*alfa/pi)         
--where alfa is the coefficient of the gaussian function

-- parMap strat f xs = DL.map f xs `using` parList strat


-- =================================================================================================
-- == the Monad List are used for calculating the elements of the lower triangular matrix         ==
-- == the elements of the triangular matrix is equal to (n^2 + n) / 2 with n the number of atoms  ==
-- =================================================================================================

-- =================================================================================================
-- == The idea for the integral is using to use the monad list for calculating
-- == the cartesian product between basis, CGF and primitive gaussian functions.
-- == Therefore the integration is of 3 levels.
-- =================================================================================================


-- ==============> TYPES <=============

type Exponent = Double 

type Pairs = (Double,Double)

type MapHermite = M.Map HermiteIndex Double

  
-- ============> ALGEBRAIC DATA TYPES <=================

data Gauss = Gauss {
             nucCoord :: !NucCoord
            ,funtype  :: !Funtype
            ,gaussP   :: !GaussPrimitive
             } deriving Show
  
data Operator =  Rij | Tij | Vij NucCoord
                 deriving Show
  
                
data HermiteIndex = Rpa {getN :: Int, getijt :: ![Int]}
                    | Xpa {getijt :: ![Int]}
                    | Ypa {getijt :: ![Int]}
                    | Zpa {getijt :: ![Int]} 
                    deriving (Show,Eq,Ord)

data EvalHermiteCoeff = EvalHermiteCoeff {
                        pabComp  :: ![Double]
                       ,getIndex :: ![Int]
                       ,cartComp :: [Int] -> HermiteIndex
                          }            

data HermiteStateCoeff = HermiteStateCoeff  {
               getmapC :: !MapHermite
              ,getkeyC :: ![[HermiteIndex]]
               } deriving Show


data HermiteStateIntegral = HermiteStateIntegral {
               getmapI :: !MapHermite
              ,getkeyI :: ![HermiteIndex]
               } deriving Show
                    

data Tree a = EmptyTree | Node a (Tree a) (Tree a) deriving (Show)

-- ================> NEWTYPES <========================

newtype Recursive a = Recursive {runRecursive :: a -> a}

 -- =================> INSTANCES <===================
 
instance Functor Tree where
  fmap f EmptyTree = EmptyTree
  fmap f (Node x l r) = Node (f x) (fmap f l) (fmap f r) 

instance DF.Foldable Tree where
   foldMap f EmptyTree = mempty  
   foldMap f (Node x l r) = f x            `mappend`
                            DF.foldMap f l `mappend`                            
                            DF.foldMap f r

instance DT.Traversable Tree where
   traverse f EmptyTree = pure EmptyTree
   traverse f (Node  k l r) = Node <$> f k <*> DT.traverse f l <*> DT.traverse f r

   
instance Num a => Monoid (Recursive a) where
  mempty = Recursive id
  
  r1 `mappend` r2 = let [f,g] = runRecursive <$> [r1,r2]
                    in Recursive $ f . g 
  
  mconcat xs =  Recursive $ DF.foldr1 (.) $ fmap runRecursive xs
   
-- ==============> Some Utilities <===================
{-
unfoldr      :: (b -> Maybe (a, b)) -> b -> [a]
unfoldr f b  =
  case f b of
   Just (a,new_b) -> a : unfoldr f new_b
   Nothing        -> []-}


unfoldWhile :: (b -> Bool) -> (b -> (a, b)) -> b -> [a]   
unfoldWhile p f = DL.unfoldr (\x -> guard (not (p x)) >> return (f x))     

safeHead :: [a] -> Maybe a
safeHead (x:xs) = Just x
safeHead []     = Nothing
            
            
-- ===============> OPERATORS <===========
 
(<<|>>) :: (NucCoord,CGF) -> (NucCoord,CGF) -> Double
cgf1 <<|>> cgf2 = sijContracted cgf1 cgf2

(<||>) :: Gauss -> Gauss -> Double
gauss1 <||> gauss2 = sab gauss1 gauss2

-- (<||>) :: (NucCoord,GaussPrimitive) -> (NucCoord,GaussPrimitive) -> (NucCoord,GaussPrimitive) -> (NucCoord,GaussPrimitive) -> Double
-- (g1,g2) <||> (g3,g4) = undefined

(<|) :: Gauss-> Operator -> (Gauss,Operator)
gauss1 <| op = (gauss1, op)

(|>) :: (Gauss,Operator) -> Gauss -> Double
(gauss1,op) |> gauss2 = case op of
                             Rij -> 1.0
                             Tij -> tab gauss1 gauss2

                             
(<<|) :: (NucCoord,Basis) -> Operator ->  ((NucCoord,Basis),Operator)
b1 <<| op  = (b1,op)

(|>>) ::  ((NucCoord,Basis),Operator) -> (NucCoord,Basis) -> Double
(b1,op) |>> b2 = case op of
                Rij -> 1.0
                Tij -> tijTotal b1 b2
                Vij rc -> vijTotal b1 rc b2

(|+>) :: Gauss -> Int -> Gauss
a |+> i = up a i

(|->) :: Gauss -> Int -> Gauss
a |-> i = down a i

--  ========> CONSTANTS <===============

overlapZero =  [[Px,Py],[Px,Pz],[Py,Pz],[Dxx,Dyy],[Dxx,Dzz],[Dyy,Dzz]]

            
-- ==============> Auxiliar Functions <===============


--(2k -1) !! = (2k)!/(2^k * k!)
-- i = 2k - 1 => k = (i + 1)/ 2
facOdd ::Int -> Double
facOdd i | i `rem`2 == 0  = error "Factorial Odd function required an odd integer as input"
         | otherwise  = case compare i 2 of
                             LT -> 1
                             GT-> let k = (1 + i) `div ` 2
                                  in (fromIntegral $ fac (2*k)) /  (2.0^k * (fromIntegral $ fac k))
                                                                                        
fac :: Int -> Int
fac i = case compare i 2 of
             LT -> 1
             EQ -> 2
             GT-> DL.foldl1' (*) [i,i-1..1]
                  
                                    
binomial :: Int -> Int -> Double
binomial l k = fromIntegral $ fac l `div` (fac k * fac (l-k))

f2k :: Double -> Double -> Int -> Int -> [Double]
f2k pa pb l1 l2  = [pa' i * pb' j * binomial l1 i * binomial l2 j | i <- [0..l1], j <- [0..l2], (i+j) `rem` 2 == 0]
  where pa' i = pa^(l1-i)
        pb' j = pb^(l2-j)

       
rab2 :: NucCoord -> NucCoord -> Double
rab2 a b = sum . DL.map (^2). DL.zipWith (-) a $ b

meanp ::(Exponent,Exponent) -> NucCoord -> NucCoord -> [Double]
meanp (e1,e2) ra rb = DL.map (\(a,b) -> (e1*a + e2*b)/(e1+e2)) $ zip ra rb

restVect :: (NucCoord,NucCoord) -> NucCoord
restVect (ra,rb) = DL.zipWith (-) ra rb

normaCoeff :: CGF -> CGF
normaCoeff b1 = b1 { getPrimitives = newPrimitives}         
  where list = getPrimitives b1
        newPrimitives = zip newCoeff . snd . unzip $ list 
        newCoeff = DL.map norma list
        norma (c1,e1) = c1 * (2.0*e1/pi)**0.75
  
-- ======================> INTEGRAL EVALUATION OF ARBITRARY ANGULAR MOMENTUM FUNCTIONS <==========

-- Contribution of the Angular part to the integral                               
                               
                                                                                                                       
funtyp2Index :: Funtype -> Int -> Int
funtyp2Index funtyp x =  fromMaybe (error "Label of the Basis function not identify " ) 
                         $ M.lookup (funtyp,x) mapLAngular


ix :: Int -> Int -> Double -> Double -> Double -> Double
ix l1 l2 pax pbx gamma = case l1 + l2 <= 1 of
                              True -> 1.0 
                              False -> sum [fun1 k * coeff !! k | k <- [0..dim]]
                            
  where dim = (l1+l2) `div ` 2
        coeff = f2k pax pbx l1 l2
        fun1 i = facOdd (2*i -1) / (2.0*gamma)^i


angTerm :: Funtype ->  Funtype -> [Double] -> [Double] -> Double -> Double
angTerm symb1 symb2 pa pb gamma = DL.foldr1 (*) [ ix (l1 x) (l2 x) (pa !! x) (pb !! x) gamma |x <- [0..2]]
  where [l1,l2] = DL.map funtyp2Index [symb1,symb2]

        
-- ==============> 2-CENTER OVERLAP INTEGRAL <=================================

-- overlaping between two primitive gaussian function between arbitrary-l functions
-- <A|B>

-- S12 = *exp[−a1*b1(Rab)^2/gamma]IxIyIz

-- where Ix = sum f2i (l1,l2,PAx,PBx) ((2i-1)!! / (2gamma)^i sqrt (pi/gamma) for i = 0,l1 + l2 / 2 only the even i
-- f2k is defined above in the module

mtxOverlap :: VU.Unbox Double => [NucCoord] -> [Basis] -> Nelec -> Array U DIM1 Double
mtxOverlap coords basis nelec = LA.list2ArrDIM1 dim flatMtx
  where flatMtx = parMap rdeepseq sijTotal cartProd
        dim = (nelec^2 + nelec) `div`2
        list = zip coords basis
        cartProd = do
          (i,(r1,b1)) <- zip [1..] list
          (j,(r2,b2)) <- zip [1..] list
          guard (i<=j)
          return (r1,b1,r2,b2)
{-# INLINE mtxOverlap #-}


sijTotal :: (NucCoord,Basis,NucCoord,Basis) -> Double
sijTotal (r1,b1,r2,b2) = DL.foldl1' (+) nonZero
  where nonZero =  do
        g1 <- b1
        g2 <- b2
        let [l1,l2] = DL.map getfunTyp [g1,g2]
        guard ((sort [l1,l2]) `notElem` overlapZero)
        let orb1 = (r1,g1)
            orb2 = (r2,g2)
        return $ ( orb1 <<|>> orb2 )
                           
       
sijContracted :: (NucCoord,CGF) -> (NucCoord,CGF) -> Double
sijContracted (r1,cgf1) (r2,cgf2) =
                 case rab2 r1 r2 > 1.0e-9  of
                      False -> 1.0
                      True -> DL.foldl1' (+) $ do
                              g1 <- getPrimitives cgf1
                              g2 <- getPrimitives cgf2
                              let [l1,l2] = DL.map getfunTyp [cgf1,cgf2]
                                  gauss1 = Gauss r1 l1 g1
                                  gauss2 = Gauss r2 l2 g2
                              return (gauss1 <||> gauss2 )

sab :: Gauss -> Gauss -> Double
sab (Gauss r1 l1 g1@(_c1,e1)) (Gauss r2 l2 g2@(_c2,e2)) = sgg' *  ang * (pi/ gamma)**1.5

  where sgg' = sgg gauss1 gauss2
        gauss1 = (r1,g1)
        gauss2 = (r2,g2)
        ang  = angTerm l1 l2 pa pb gamma
        p = meanp (e1,e2) r1 r2
        [pa,pb] = (DL.zipWith (-) p) `fmap`  [r1,r2]
        gamma = e1 + e2
                              
             
sgg :: (NucCoord,GaussPrimitive) -> (NucCoord,GaussPrimitive) -> Double
sgg (r1,(c1,e1)) (r2,(c2,e2)) = c1*c2 * exp (-e1*e2*r/(e1+e2))
  where r = rab2 r1 r2

             
-- ==================> Kinetic Energy IntegralsEvaluation  <========================
-- the kinetic integral for two S-functions is
-- < A| -0.5*nabla^2 | B > = (3*w - 2*w*w*rab2 ) * <A|B>
-- where w = e1*e2 / (e1 +e2)

--for arbitrary angular momentum functions is
-- < A| -0.5*nabla^2 | B > = Ix Iy Iz
-- where Ix = 0.5* l1*l2 *<-1|-1>x + 2a1b1<+1|+1>x - a1*l2*<+1|-1>x - b1*l1*<-1|+1>x
--       |+1>x = X^(l1+1) Y^m1 Z^n1 exp(-a1*Ra^2)
--       |-1>x = X^(l1-1) Y^m1 Z^n1 exp(-a1*Ra^2)
 

arrayKinetic :: [NucCoord] -> [Basis] ->  Array U DIM1 Double
arrayKinetic coords basis = LA.list2ArrDIM1 dim cartProd
  where n = length basis
        dim = (n^2 + n) `div`2        
        list = zip coords basis
        cartProd = do
          (i,(r1,b1)) <- zip [1..] list
          (j,(r2,b2)) <- zip [1..] list
          guard (i<=j)
          return $ (r1,b1) <<|Tij|>> (r2,b2)        
{-# INLINE arrayKinetic #-}

tijTotal :: (NucCoord,Basis) -> (NucCoord,Basis) -> Double
tijTotal (r1, b1) (r2, b2) = DL.foldl1' (+) [tijContracted (r1,g1) (r2,g2)  | g1 <- b1, g2 <- b2]


tijContracted :: (NucCoord,CGF) -> (NucCoord,CGF) -> Double
tijContracted (r1,cgf1) (r2,cgf2) =
       DL.foldl1' (+) $ do
              g1 <- getPrimitives cgf1
              g2 <- getPrimitives cgf2
              let [l1,l2] = DL.map getfunTyp [cgf1,cgf2]
                  gauss1 = Gauss r1 l1 g1 
                  gauss2 = Gauss r2 l2 g2
              return ( gauss1 <| Tij |> gauss2 )  
 
 
tab :: Gauss -> Gauss -> Double
tab g1@(Gauss r1 symb1 (_c1,e1)) g2@(Gauss r2 symb2 (_c2,e2)) = DL.foldl1' (+) [ tx x (j x) (k x) | x <- [0..2]]
  where [j,k] = DL.map funtyp2Index [symb1,symb2]
        tx i lang1 lang2 =
                     let  [gauss1plus,gauss2plus]   = DL.zipWith (|+>) [g1,g2] $ repeat i 
                          [gauss1minus,gauss2minus] = DL.zipWith (|->) [g1,g2] $ repeat i                      
                          [l1,l2] = DL.map fromIntegral [lang1,lang2] 
                          
                     in  0.5*l1*l2*(gauss1minus  <||> gauss2minus) + 2.0*e1*e2*(gauss1plus  <||> gauss2plus)                                          
                         - e1*l2* (gauss1plus  <||> gauss2minus) - e2*l1* (gauss1minus  <||> gauss2plus)



up :: Gauss -> Int -> Gauss
up gauss i = gauss {funtype=newL} 
  where newL = fromMaybe (error "failure in the integration scheme") $ M.lookup (symb1,i) obaraSaikaUp       
        symb1 = funtype gauss

  
down :: Gauss -> Int -> Gauss
down gauss i = gauss {funtype=newL} 
  where newL = fromMaybe (error "failure in the integration scheme") $ M.lookup (symb1,i) obaraSaikaDown 
        symb1 = funtype gauss



-- ======================> McMURCHIE -DAVIDSON SCHEME <=========================
-- The integrals of three and four centers are implemented according to chapter 9
-- of the book "Molecular Electronic-Structure Theory" by Trygve Helgaker,
-- Poul Jorgensen and Jeppe Olsen.

-- COEFFICIENT RECURSION FORMULA 
 -- nucleus electron attraction

vijTotal :: (NucCoord, Basis) -> NucCoord -> (NucCoord, Basis) -> Double
vijTotal  (ra, b1) rc (rb, b2) =  sum [vijContracted (ra,g1) rc (rb,g2) | g1 <- b1, g2 <- b2]       

        
vijContracted :: (NucCoord,CGF) -> NucCoord -> (NucCoord,CGF) -> Double
vijContracted (ra,cgf1) rc (rb,cgf2) = DL.foldl1' (+) cartProd 
                                       
  where cartProd = do
           g1 <- getPrimitives cgf1
           g2 <- getPrimitives cgf2
           let [l1,l2] = DL.map getfunTyp [cgf1,cgf2]
               gauss1 =  Gauss ra l1 g1
               gauss2 =  Gauss rb l2 g2
           return $ vijHermite gauss1 gauss2 rc
                   
vijHermite :: Gauss -> Gauss -> NucCoord -> Double
vijHermite g1 g2 rc = cte * (mcMurchie symbols [ra,rb,rc] (e1,e2))
  where cte = c1 * c2 * 2.0 * (pi/gamma) 
        gamma = e1+e2     
        [ra,rb] =  nucCoord `fmap` [g1,g2]
        symbols = funtype `fmap` [g1,g2]
        [(c1,e1),(c2,e2)] = gaussP `fmap` [g1,g2]

        
mcMurchie :: [Funtype] -> [NucCoord]  -> (Exponent,Exponent) -> Double    
mcMurchie symbols [ra,rb,rc] (e1,e2) = 
  let coeff = DL.unfoldr (calcHermCoeff [rpa,rpb] gamma) seedC 
      rtuv  = DL.unfoldr (calcHermIntegral rpc gamma) seedI
      gamma = e1 + e2
      nu = e1*e2/gamma
      rp = meanp (e1,e2) ra rb
      [rab,rpa,rpb,rpc] = restVect `fmap` [(ra,rb),(rp,ra),(rp,rb),(rp,rc)]
      seedC = initilized_Seed_Coeff symbols rab nu
      seedI = initilized_Seed_Integral symbols rpc gamma
      
  in DL.foldl1' (+) $ DL.zipWith (*) coeff rtuv 
        
  
calcHermCoeff :: [NucCoord] -> Double -> HermiteStateCoeff-> Maybe (Double,HermiteStateCoeff)
calcHermCoeff rpab gamma stCoeff = do
    ls <- safeHead listC
    let newSt = stCoeff {getmapC = newMap,getkeyC = tail listC}
        (val,newMap) = fun ls 
    return (val,newSt)     
              
  where listC = getkeyC stCoeff
        oldmap = getmapC stCoeff
        fun [xpa,ypa,zpa] = runState (hermEij xpa 0 >>= \e1 ->
                            hermEij ypa 1 >>= \e2 ->  
                            hermEij zpa 2 >>= \e3 ->
                            return $ e1*e2*e3 ) oldmap        
                          
        hermEij k label = get >>= \st -> 
                          let xpab = DL.map (\r -> r !! label) rpab
                              ijt = getijt k
                              Just (x,newMap) = (lookupM k st) `mplus` 
                                   Just (recursiveHC xpab gamma k st ijt)
                          in put newMap >> return x                                              
        
        
calcHermIntegral :: NucCoord  -> Double -> HermiteStateIntegral -> Maybe (Double,HermiteStateIntegral)
calcHermIntegral rpc gamma stIntegral =  do
    hi <- safeHead listI
    let Just (val,newMap) = fun hi
        newSt = stIntegral {getmapI = newMap,getkeyI = tail listI} 
    return (val,newSt)                     
    
  where listI = getkeyI stIntegral
        oldmap = getmapI stIntegral
        fun hi = (lookupM hi oldmap) `mplus` Just (recursiveHI rpc gamma hi oldmap (getijt hi))


recursiveHC :: [Double] -> Double -> HermiteIndex ->  MapHermite -> [Int] -> (Double, MapHermite)
recursiveHC pab@[xpa,xpb] gamma hi mapC [i,j,t]
  
       |(i+j) < t = (0.0,mapC)
 
       |t == 0 && i>=1 = let ([cij0,cij1],mapO) = fun [1,0,0] [1,0,(-1)]                                                                      
                             val = xpa*cij0 + cij1
                             newMap = M.insert hi val mapO                              
                         in (val,newMap)
                                                                    
       |t == 0 && j>=1 = let ([cij0,cij1],mapO) = fun [0,1,0] [0,1,(-1)]
                             val = xpb*cij0 + cij1
                             newMap = M.insert hi val mapO                                                           
                         in (val,newMap)
                                                                                                   
       |otherwise =      let ([aijt,bijt],mapO) = fun [1,0,1] [0,1,1]
                             [i',j',t'] = fromIntegral `fmap`  [i,j,t]
                             val = recip (2.0*gamma*t') * (i'*aijt + j'*bijt)     
                             newMap = M.insert hi val mapO                                                                                        
                         in (val,newMap)
                                                                        
  where fun xs ys = runState (sequence [fun2 xs, fun2 ys]) mapC
        fun2 xs = get >>= \st -> let key = sub xs
                                     funAlternative = recursiveHC pab gamma key st
                                     (v,m) = boyMonplus funAlternative key st                                                                          
                                 in put m >> return v                                     
        sub [a,b,c] = hi {getijt = [i-a,j-b,t-c]}                                                    
        
recursiveHI :: NucCoord -> Double -> HermiteIndex-> MapHermite -> [Int] -> (Double,MapHermite)
recursiveHI rpc gamma hi mapI [t,u,v] 
  
 | any (<0) [t,u,v] = (0.0,mapI)          

 | t >= 1 = let ([boy1,boy2],mapO) = fun [2,0,0] [1,0,0] 
                val = (fromIntegral t -1)*boy1 + (rpc !! 0)*boy2
                newM = M.insert hi val mapI  
            in (val,newM)
                                                 
 | u >= 1 = let ([boy1,boy2],mapO) = fun [0,2,0] [0,1,0] 
                val = (fromIntegral u -1)*boy1 + (rpc !! 1)*boy2    
                newM = M.insert hi val mapI  
            in (val,newM)
                                           
 | v >= 1 = let ([boy1,boy2],mapO) = fun [0,0,2] [0,0,1] 
                val = (fromIntegral v -1)*boy1 + (rpc !! 2)*boy2 
                newM = M.insert hi val mapI  
            in (val,newM)                                           
            
 | otherwise =  let arg = (gamma*) . sqrt . sum . (fmap (^2))  $ rpc
                    x1 = (-2.0*gamma)^n
                    val = x1* boyFunction n arg
                    newM = M.insert hi val mapI  
                in (val,newM)                                                         
      
  where fun xs ys = runState (sequence [fun2 xs, fun2 ys]) mapI
        fun2 xs = get >>= \st -> let key = sub xs
                                     funAlternative = recursiveHI rpc gamma key st
                                     (v,m) = boyMonplus funAlternative key st
                                 in put m >> return v
        n = getN hi -- order of the boy function
        sub [a,b,c] = hi {getN = n + 1, getijt = [t-a,u-b,v-c]}        

        
boyMonplus :: ([Int] -> (Double,MapHermite)) -> HermiteIndex-> MapHermite -> (Double,MapHermite)
boyMonplus fun k m = fromMaybe err $ (lookupM k m) `mplus` (Just (res,newM))
  where (res,oldmap) = fun (getijt k)
        newM = M.insert k res oldmap 
        err = error "failure in the recursive Hermite Procedure"
  

boyFunction :: Int -> Double -> Double 
boyFunction m arg = let boy = \x -> x^(2*m) * exp(-arg*x*x)
                        maxerror = 1E-10 -- maximum allow error
                        (val,numerror) = integrateQNG maxerror boy 0 1
                    in val 

lookupM :: Ord k => k -> M.Map k a -> Maybe (a , M.Map k a)
lookupM k m = do 
              val <- M.lookup k m
              return (val,m)                                                             
              

-- we sent to listIndexes two Funtypes and the function return 
-- the exponent for x,y,z according to the l-number of the
-- orbital that thos Funtypes represent
-- therefore for ls = [l1,m1,n1,l2,m2,n2]

listIndexes :: [Funtype] -> [Int]        
listIndexes symbols = ls
  where ls = funtyp2Index <$> symbols <*> [0..2]

  
initilized_Seed_Coeff :: [Funtype] -> NucCoord -> Double -> HermiteStateCoeff        
initilized_Seed_Coeff symbols rab mu = HermiteStateCoeff mapC0 listC
  where mapC0 = DL.foldl' (\acc (k, v) -> M.insert k v acc) M.empty $ DL.zip k0 e0 
        listC = genCoeff_Hermite symbols 
        k0 = [Xpa [0, 0, 0], Ypa [0, 0, 0], Zpa [0, 0, 0]] 
        e0 = DL.map (\x -> exp(-mu*(x^2))) rab
        
        
genCoeff_Hermite :: [Funtype] -> [[HermiteIndex]]
genCoeff_Hermite symbols = do
  i <-[0..l1+l2]
  j <-[0..m1+m2]
  k <-[0..n1+n2]
  return $ [Xpa [l1,l2,i], Ypa [m1,m2,j], Zpa [n1,n2,k]] 
  where [l1,m1,n1,l2,m2,n2] = listIndexes symbols        


initilized_Seed_Integral :: [Funtype] -> NucCoord -> Double -> HermiteStateIntegral
initilized_Seed_Integral symbols rpc gamma = HermiteStateIntegral mapI0 listI 
  where mapI0 = M.insert k0' f0 M.empty
        k0'= Rpa 0 [0, 0, 0] 
        y= gamma*rpc2
        f0 = if (y > 0.0e-5) then 0.5* sqrt (pi/(y)) * erf(sqrt y) else 1.0
        rpc2 = sum $ DL.map (^2) rpc
        listI = genCoeff_Integral symbols
        
                
genCoeff_Integral :: [Funtype] -> [HermiteIndex]
genCoeff_Integral symbols = do
  t <-[0..l1+l2]
  u <-[0..m1+m2]
  v <-[0..n1+n2]
  return $ Rpa 0 [t,u,v] 
  where [l1,m1,n1,l2,m2,n2] = listIndexes symbols        
  
 
hcore :: [NucCoord] -> [Basis] -> [ZNumber] -> Nelec -> Array U DIM1 Double
hcore coords basis atomicZ nelec =  LA.list2ArrDIM1 dim (cartProd `using` parList rdeepseq)
 
 where dim = (nelec^2 + nelec) `div`2
       list = zip coords basis       
       cartProd = do
          (i,atomi) <- zip [1..] list
          (j,atomj) <- zip [1..] list
          guard (i<=j)
          let sumVij = foldl1' (+) . getZipList $
                       (\z rc -> ((-z) * atomi <<|Vij rc|>> atomj))
                         <$> ZipList atomicZ <*> ZipList coords
          return $ (atomi <<|Tij|>> atomj) + sumVij
                    
           

-- =====================> TWO ELECTRON INTEGRALS <=========================================

-- <XaXb|XcXd> = 2 * <A|B> * <C|D> * sqrt (rho/pi) * Fo (rho(P - Q)^2) where Fo(x)= (pi/x)^½ * erf(x^½) P = e1*Ra + e2*Rb / (e1+e2)  
-- rho = (e1+e2)(e3+e4)/(e1+e2+e3+e4)

evalIntbykey :: [NucCoord]-> [Basis]-> [[Int]]-> [Double]
evalIntbykey coord basis keys = (integrals keys) `using` parList rdeepseq
  where integrals =  DL.map (\xs -> let bs = DL.map (\i -> basis !! i) xs
                                        rs = DL.map (\j -> coord !! j) xs                                       
                                        atoms = getZipList $ AtomData <$> ZipList rs <*> ZipList bs                       
                                    in  fourCenterInt atoms  )
        
        
fourCenterInt ::  [AtomData] -> Double       
fourCenterInt atoms = DL.foldl1' (+) cartProd
  where ([ra,rb,rc,rd],[b1,b2,b3,b4]) = (fmap getCoord) &&& (fmap getBasis) $ atoms
        cartProd = [contracted4Centers [(ra,cgf1),(rb,cgf2),(rc,cgf3),(rd,cgf4)] |cgf1 <- b1, cgf2 <- b2, cgf3 <- b3, cgf4 <- b4]

contracted4Centers :: [(NucCoord,CGF)] -> Double
contracted4Centers [(ra,cgf1), (rb,cgf2), (rc,cgf3), (rd,cgf4)] = DL.foldl1' (+) cartProd

  where [l1,l2,l3,l4] = getfunTyp `fmap` [cgf1,cgf2,cgf3,cgf4]
        cartProd = do               
          g1 <- getPrimitives cgf1
          g2 <- getPrimitives cgf2
          g3 <- getPrimitives cgf3
          g4 <- getPrimitives cgf4
          let gauss = getZipList $ Gauss <$> ZipList[ra,rb,rc,rd] <*> 
                                   ZipList [l1,l2,l3,l4] <*> ZipList [g1,g2,g3,g4]
          return $ twoElectronHermite gauss
  
  
twoElectronHermite :: [Gauss] -> Double 
twoElectronHermite gs = (cte *) . foldl1' (+) . DL.zipWith (*) coeff1 $ 
                        [mcMurchie2 coeff2 rpq alpha abcs tuv| tuv <- coefftuv]

  where coefftuv = genCoeff_Integral [symb1,symb2]
        abcs =  genCoeff_Integral [symb3,symb4]
        coeff1 = DL.unfoldr (calcHermCoeff [rpa,rpb] p) seedC1
        coeff2 = DL.unfoldr (calcHermCoeff [rqc,rqd] q) seedC2
        seedC1 = initilized_Seed_Coeff [symb1,symb2] rab mu
        seedC2 = initilized_Seed_Coeff [symb3,symb4] rcd nu
        [ra,rb,rc,rd] = nucCoord `fmap` gs      
        [symb1,symb2,symb3,symb4] = funtype `fmap` gs      
        ps = gaussP `fmap` gs
        ([c1,c2,c3,c4],[e1,e2,e3,e4]) = (fmap fst ) &&& (fmap snd ) $ ps
        rp = meanp (e1,e2) ra rb
        rq = meanp (e3,e4) rc rd
        [rab,rcd,rpa,rpb,rqc,rqd,rpq] = restVect `fmap` (zip [ra,rc,rp,rp,rq,rq,rp][rb,rd,ra,rb,rc,rd,rq])
        [mu,nu] = (\(a,b) -> a*b* (recip $ a + b)) `fmap` [(e1,e2),(e3,e4)]
        [p,q] = uncurry (+) `fmap` [(e1,e2),(e3,e4)]
        alpha = p*q/(p+q)
        cte = (c1*c2*c3*c4*) . (2.0*pi**2.5 *) . recip $ (p * q ) * (sqrt $ p + q)

        
mcMurchie2 :: [Double] -> NucCoord -> Double -> [HermiteIndex] -> HermiteIndex -> Double
mcMurchie2 coeff2 rpq alpha abcs' tuv' =
  DL.foldl1' (+) . DL.zipWith3 (\x y z -> x*y*z) sgns coeff2  $ integrals 
  
  where tuv = getijt tuv'
        abcs = getijt `fmap` abcs'
        sgns = (\xs-> (-1.0)^(sum xs)) `fmap` abcs
        integrals = DL.unfoldr (calcHermIntegral rpq alpha) seedI 
        seedI = HermiteStateIntegral mapI0 listI 
        mapI0 = M.insert k0' f0 M.empty
        k0'= Rpa 0 [0, 0, 0]
        rp2 = sum $ DL.map (^2) rpq
        y = alpha*rp2
        f0 = if (y > 0.0e-5) then 0.5* sqrt (pi/(y)) * erf(sqrt y) else 1.0
        listI = Rpa 0 `fmap` (DL.map (DL.zipWith (+) tuv) abcs )

