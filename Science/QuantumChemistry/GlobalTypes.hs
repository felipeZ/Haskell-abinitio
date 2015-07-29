{-# Language DeriveFunctor #-}

-- The HaskellFock SCF Project
-- @2013 Felipe Zapata, Angel Alvarez
-- shared types

module Science.QuantumChemistry.GlobalTypes where

import Control.Applicative
import Control.DeepSeq
import Data.Array.Repa          as R
import qualified Data.Foldable as DF
import qualified Data.Map as M
import Data.Maybe (fromMaybe)
import Data.Monoid (Monoid(..),Sum(..),mappend,mconcat,mempty)
import qualified Data.Traversable as DT
import qualified Data.Vector.Unboxed as VU

--  =========================> TYPES <=====================
               
--  | Basis set
type Basis = [CGF]

-- | charge of the system
type Charge = Double

-- | Molecular Orbitals Coefficients
type Coeff = [Double]

-- | EigenValues are represented as an unboxed Vector
type EigenValues = VU.Vector Double

-- | EigenVectors are represented as a DIM2 Repa Array
type EigenVectors = Array U DIM2 Double

type Flatten = Array U DIM1 Double

-- | Coefficient and exponent describing a primitive gaussian function
type GaussPrimitive = (Double,Double)

-- | type Energy derivatives with respect to the cartesian muclear coordinates
type Gradient = Array U DIM1 Double

type Matrix = Array U DIM2 Double

-- | Total number of electrons
type Nelec = Int

-- | Nuclear Coordinates
type NucCoord = [Double]

-- | Nuclei-Nuclei Energy repulsion
type NuclearRepulsion = Double

 --  | Occupation number
type OccupiedShells   = Int


-- | Iteration Step
type Step = Int

type Threshold = Double

type VecUnbox  = VU.Vector Double

-- | Atomic Number
type ZNumber = Double

-- ===================> Array Types <===============
type ErrorMatrix = Array U DIM2 Double
type MOCoefficients       = Array U DIM2 Double
type FlattenChargeDensity = Array U DIM1 Double
type FlattenCore          = Array U DIM1 Double
type FlattenErrorMtx      = Array U DIM1 Double
type FlattenFock          = Array U DIM1 Double
type FlattenMatrix        = Array U DIM1 Double
type FlattenOverlap       = Array U DIM1 Double
type Integrals            = Array U DIM1 Double
type TransformMatrix      = Array U DIM2 Double

-- ========================= ALGEBRAIC DATA TYPES <=================
data AtomData = AtomData {
                 getCoord   :: !NucCoord  -- ^ Cartesian Nuclear Coordinates stored as REPA arrays
                ,getBasis   :: !Basis     -- ^ List of the Contracted Gaussian functions for the current atoms
                ,getZnumber :: !Double    -- ^ Atomic Number
                 } deriving (Show,Eq,Ord)

-- | Symbolic manipulation of X,Y and Z axis improves readability
data CartesianLabel = Ax | Ay | Az deriving (Show,Eq)

-- | Contracted Gaussian function
data CGF = CGF  {
                getPrimitives ::[GaussPrimitive] -- ^ List of primitives
               ,getfunTyp ::Funtype              -- ^ angular momentum
                } deriving (Show,Eq,Ord)

data Derivatives a = Dij_Ax a | Dij_Ay a | Dij_Az a deriving (Show,Functor)

data Gauss = Gauss {
             nucCoord :: !NucCoord        -- ^ Nuclear center where the gaussian function is centered
            ,funtype  :: !Funtype         -- ^ Angular Momentum 
            ,gaussP   :: !GaussPrimitive  -- ^ Coefficient and exponent parameters of the primitive
             } deriving (Show,Eq)                                          

-- | Unary representation of the angular momenta                 
data Funtype = Zero| S | Px | Py | Pz | Dxx | Dxy | Dxz | Dyy | Dyz | Dzz  | Fxxx | Fxxy | Fxxz |
               Fxyz | Fxyy | Fxzz | Fyyy | Fyyz | Fyzz | Fzzz
               deriving (Show,Eq,Enum,Ord)
               

-- |Hartree-Fock Results
data HFData = HFData {
           getFock      :: !FlattenFock           -- ^ Current flatten Fock Matrix
         , getCoeff     :: !EigenVectors          -- ^ Molecular orbital coefficients
         , getDensity   :: !FlattenChargeDensity  -- ^ Charge density considering only occupied orbitals
         , getOrbE      :: !EigenValues           -- ^ Molecular Orbital Energies
         , getEnergy    :: !Double                -- ^ Hartree-Fock Total Energy
                     } deriving (Show)

-- | EigenProblem solution
data EigenData = EigenData {
             eigenvals :: !EigenValues
           , eigenvec :: !EigenVectors } deriving (Show)

           
data Tree a = EmptyTree | Node a (Tree a) (Tree a) deriving (Show)

data Choice a = Ignore | Take a  deriving (Show,Functor)

data Switch = ON | OFF

-- ================> NEWTYPES <========================

newtype Recursive a = Recursive {runRecursive :: a -> a} 


 -- =================> INSTANCES <===================
 
instance Bounded Funtype where
  minBound = S
  maxBound = Fzzz
 
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
  
  
instance Monoid a => Monoid (Choice a) where
  mempty = Ignore
  Ignore `mappend` Ignore = Ignore
  Ignore `mappend` Take x = Take x
  Take x `mappend` Ignore = Take x
  Take x `mappend` Take y = Take $ x `mappend` y

instance NFData a => NFData (Choice a) where
  rnf Ignore   = ()
  rnf (Take x) = rnf x

instance NFData a => NFData (Sum a) where
  rnf (Sum a) = rnf a 
  
-- ===================> <=====================  
choices :: b -> (a -> b) -> Choice a -> b
choices x0 f c =
  case c of
       Ignore -> x0
       Take x -> f x
       
-- ====================> Operators <=======================


-- | Operator to increase the angular momentum
(|+>) :: Funtype -> CartesianLabel ->  Funtype
symb |+> label =  getAngMom newMomentum
                  
  where [x,y,z] = getAxesMom symb
        newMomentum = case label of
                           Ax -> [succ x,y,z]
                           Ay -> [x,succ y,z]
                           Az -> [x,y,succ z]
  

-- | Operator to decrease the angular momentum
(|->) :: Funtype -> CartesianLabel ->  Funtype
symb |-> label =  case label of
                       Ax -> check x [pred x,y,z]
                       Ay -> check y [x,pred y,z]
                       Az -> check z [x,y,pred z]

  where [x,y,z] = getAxesMom symb
        check i xs = if i == 0 then Zero else getAngMom xs

getCGFangularExpo :: Funtype -> CartesianLabel -> Double
getCGFangularExpo symb ax =
  let [x,y,z] =  fmap fromIntegral . getAxesMom $ symb
  in case ax of
          Ax -> x                            
          Ay -> y
          Az -> z       
        
toCartLabel :: Int -> CartesianLabel
toCartLabel k | k == 0 = Ax
              | k == 1 = Ay
              | k == 2 = Az
              | otherwise = error "Unexpected Integer <= toCartLabel"
                  
fromCartLabel :: CartesianLabel -> Int
fromCartLabel k | k == Ax = 0
                | k == Ay = 1
                | k == Az = 2

              
-- ===========> Auxiliar Functions <======

-- | Map from Unary Angular momenta representation to its corresponding integer
mapLAngular :: M.Map (Funtype,Int) Int
mapLAngular = M.fromList $ 
            [((S,0),0),((S,1),0),((S,2),0)
            ,((Px,0),1),((Px,1),0),((Px,2),0)
            ,((Py,0),0),((Py,1),1),((Py,2),0)
            ,((Pz,0),0),((Pz,1),0),((Pz,2),1)
            ,((Dxx,0),2),((Dxx,1),0),((Dxx,2),0)
            ,((Dyy,0),0),((Dyy,1),2),((Dyy,2),0)
            ,((Dzz,0),0),((Dzz,1),0),((Dzz,2),2)
            ,((Dxy,0),1),((Dxy,1),1),((Dxy,2),0)
            ,((Dxz,0),1),((Dxz,1),0),((Dxz,2),1)
            ,((Dyz,0),0),((Dyz,1),1),((Dyz,2),1)
            ,((Fxxx,0),3),((Fxxx,1),0),((Fxxx,2),0)
            ,((Fxxy,0),2),((Fxxy,1),1),((Fxxy,2),0)
            ,((Fxxz,0),2),((Fxxz,1),0),((Fxxz,2),1)
            ,((Fxyz,0),1),((Fxyz,1),1),((Fxyz,2),1)
            ,((Fxyy,0),1),((Fxyy,1),2),((Fxyy,2),0)
            ,((Fxzz,0),1),((Fxzz,1),0),((Fxzz,2),2)
            ,((Fyyy,0),0),((Fyyy,1),3),((Fyyy,2),0)
            ,((Fyyz,0),0),((Fyyz,1),2),((Fyyz,2),1)
            ,((Fyzz,0),0),((Fyzz,2),1),((Fyzz,2),2)
            ,((Fzzz,0),0),((Fzzz,2),0),((Fzzz,2),3)]
            

-- | Map from the Integrals representing the exponents of
--   the gaussian functions (total angular momentum ) to 
--   the unary data representation
mapAngMomentum :: M.Map [Int] Funtype
mapAngMomentum = M.fromList $
              [([0,0,0],S)
              ,([1,0,0],Px),([0,1,0],Py),([0,0,1],Pz)
              ,([2,0,0],Dxx),([1,0,1],Dxz),([1,1,0],Dxy),([0,2,0],Dyy),([0,1,1],Dyz),([0,0,2],Dzz)
              ,([3,0,0],Fxxx),([2,1,0],Fxxy),([2,0,1],Fxxz),([1,1,1],Fxyz),([1,2,0],Fxyy),([1,0,2],Fxzz),([0,3,0],Fyyy),([0,2,1],Fyyz),([0,1,2],Fyzz),([0,0,3],Fzzz) ]

              
-- | Map from the unary data representation to the 
--  Integrals representing the exponents of 
--   the cartesian gaussian functions (total angular momentum )
mapAxesMomentum :: M.Map Funtype [Int]
mapAxesMomentum = M.fromList $
              [(S  ,[0,0,0])
              ,(Px ,[1,0,0])
              ,(Py ,[0,1,0])
              ,(Pz ,[0,0,1])
              ,(Dxx,[2,0,0])
              ,(Dxz,[1,0,1])
              ,(Dxy,[1,1,0])
              ,(Dyy,[0,2,0])
              ,(Dyz,[0,1,1])
              ,(Dzz,[0,0,2])
              ,(Fxxx,[3,0,0])
              ,(Fxxy,[2,1,0])
              ,(Fxxz,[2,0,1])
              ,(Fxyz,[1,1,1])
              ,(Fxyy,[1,2,0])
              ,(Fxzz,[1,0,2])
              ,(Fyyy,[0,3,0])
              ,(Fyyz,[0,2,1])
              ,(Fyzz,[0,1,2])
              ,(Fzzz,[0,0,3])]


getAngMom :: [Int] -> Funtype
getAngMom axis =
  if null axis then Zero
               else fromMaybe (error "there is not such angular momentum label") $ M.lookup axis mapAngMomentum

getAxesMom :: Funtype -> [Int]
getAxesMom symb =
  if symb == Zero then []
                  else fromMaybe (error "there is not such Axis") $ M.lookup symb mapAxesMomentum

-- ===================> Map of the atomic element to nuclear charge <=============
-- | String to Nuclear Charge
atom2charge :: M.Map String Int
atom2charge = M.fromList $
             [("H",1),("He",2),("Li",3),("Be",4),("B",5),
              ("C",6),("N",7),("O",8),("F",9),("Ne",10),
              ("Na",11),("Mg",12)]
