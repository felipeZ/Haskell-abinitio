
module GlobalTypes (
                  AtomData(..)
                 ,Basis
                 ,CGF (..)
                 ,Coeff
                 ,Funtype(..)
                 ,GaussPrimitive
                 ,MO
                 ,Nelec
                 ,NucCoord                 
                 ,Point(..)
                 ,ZNumber
                 ,atom2charge
                 ,mapLAngular
                 ,obaraSaikaUp
                 ,obaraSaikaDown                 
                   ) where

import qualified Data.Map as M                   
                   
-- nuclear Coordinates 
type NucCoord = [Double]

                
-- Coefficient and exponent describing a primitive gaussian function
type GaussPrimitive = (Double,Double)

-- Molecular Orbitals
type MO = Double

-- Molecular Orbitals Coefficients
type Coeff = [Double]

--  Basis set
type Basis = [CGF]

-- Atomic Number
type ZNumber = Double

-- Total number of electrons
type Nelec = Int

-- Contracted Gaussian function (coefficient,exponent)
data CGF = CGF  {
                getPrimitives ::[GaussPrimitive]
               ,getfunTyp ::Funtype
                } deriving (Show)
              
data Point  = Point {
                    getX::Double
                   ,getY::Double
                   ,getZ::Double
                    } deriving (Show,Eq)
                              
data AtomData = AtomData {
                 getCoord :: !NucCoord
                ,getBasis :: !Basis
                 } deriving (Show)

data Funtype = S | Px  | Py | Pz | Dxx | Dxy | Dxz | Dyy | Dyz | Dzz  | Fxxx | Fxxy | Fxxz |
               Fxyz | Fxyy | Fxzz | Fyyy | Fyyz | Fyzz | Fzzz
               deriving (Show,Eq,Enum,Ord)
                         
instance Bounded Funtype where
  minBound = S
  maxBound = Fzzz
  
-- ===========> Auxiliar Functions <======

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
            ,((Fzzz,0),0),((Fzzz,2),0),((Fyzz,2),3)]

obaraSaikaUp :: M.Map (Funtype,Int) Funtype            
obaraSaikaUp = M.fromList $
            [((S,0),Px),((S,1),Py),((S,2),Pz)
            ,((Px,0),Dxx),((Px,1),Dxy),((Px,2),Dxz)
            ,((Py,0),Dxy),((Py,1),Dyy),((Py,2),Dyz)
            ,((Pz,0),Dxz),((Pz,1),Dyz),((Pz,2),Dzz)
            ,((Dxx,0),Fxxx),((Dxx,1),Fxxy),((Dxx,2),Fxxz)
            ,((Dyy,0),Fxyy),((Dyy,1),Fyyy),((Dyy,2),Fyyz)
            ,((Dzz,0),Fxzz),((Dzz,1),Fyzz),((Dzz,2),Fzzz)
            ,((Dxy,0),Fxxy),((Dxy,1),Fxyy),((Dxy,2),Fxyz)
            ,((Dxz,0),Fxxz),((Dxz,1),Fxyz),((Dxz,2),Fxzz)
            ,((Dyz,0),Fxyz),((Dyz,1),Fyyz),((Dyy,2),Fyyz)]  

            
obaraSaikaDown :: M.Map (Funtype,Int) Funtype                        
obaraSaikaDown = M.fromList $        
            [((S,0),S),((S,1),S),((S,2),S)
            ,((Px,0),S),((Px,1),S),((Px,2),S)
            ,((Py,0),S),((Py,1),S),((Py,2),S)
            ,((Pz,0),S),((Pz,1),S),((Pz,2),S)
            ,((Dxx,0),Px),((Dxx,1),Dxx),((Dxx,2),Dxx)
            ,((Dyy,0),Dyy),((Dyy,1),Py),((Dyy,2),Dyy)
            ,((Dzz,0),Dzz),((Dzz,1),Dzz),((Dzz,2),Pz)
            ,((Dxy,0),Py),((Dxy,1),Px),((Dxy,2),Dxy)
            ,((Dxz,0),Pz),((Dxz,1),Dxz),((Dxz,2),Px)
            ,((Dyz,0),Dyz),((Dyz,1),Pz),((Dyy,2),Dyy)] 

-- ===================> Map of the atomic element to nuclear charge <=============
atom2charge :: M.Map String Int
atom2charge = M.fromList $
             [("H",1),("He",2),("Li",3),("Be",4),("B",5),
              ("C",6),("N",7),("O",8),("F",9),("Ne",10),
              ("Na",11),("Mg",12)]