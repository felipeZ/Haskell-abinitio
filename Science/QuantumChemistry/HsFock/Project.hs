
{-|
Module: Science.QuantumChemistry.HsFock.Project where
Description: Data types for the project
Copyright: @2012,2013 Angel Alvarez Adhesive tape
           @2012,2013 Felipe Zapata core SCF machinery 
           @2016 Felipe Zapata
-}


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
