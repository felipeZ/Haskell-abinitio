

-- The thing...
-- A system shell for the HaskellFock SCF Project 
-- @2012,2013 Angel Alvarez Adhesive tape
-- @2012,2013 Felipe Zapata core SCF machinery 

module Science.QuantumChemistry.HsFock.Project where


import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.HartreeFock.HartreeFock

data ProjectData = PD 
    {
          pLabel      :: String       -- "water"
        , pBasisType  :: String       -- "sto-3G"
        , pBasis      :: [Basis]
        , atomList    :: [NucCoord]   -- Atom Cartesian Coordinates   
    } deriving (Show)
