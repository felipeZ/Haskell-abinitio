{-# LANGUAGE DeriveGeneric, OverloadedStrings #-}
-- Basis Parser, a parsec based parser for quantum chemistry basis files
-- @2013-2015 Angel Alvarez, Felipe Zapata
-- 
--  2013/05/05 (Spain, Mother's Day) We kindly appreciate our mothers efforts in bringing us so far...

module Science.QuantumChemistry.ParsecTools.ParserBasis
    (
      Element(..)
     ,GaussShape(..)
     ,parseBasisFile
    ) where

-- ====================> Standard Modules and third party <==============================    
import Data.Attoparsec.ByteString.Char8 
import Data.Serialize 
import GHC.Generics


-- ====================> Internal Modules <=======================
import Science.QuantumChemistry.GlobalTypes 
import Science.QuantumChemistry.ParsecTools.ParseUtils 

-- ================================================ Types =============================================================
type AtomLabel = String

data GaussShape = 
      S0 [(Double,Double)]         -- | A S  type Gaussian primitive (Coefficient,Exponent)
    | SP [(Double,Double,Double)]  -- | A SP type Gaussian primitives contains shared data among S and P types, gaussian functions (Coeficcient S , Coefficcient P, Exponent)
    | P  [(Double,Double)]         -- | A P  type Gaussian primitive (Coefficient,Exponent)
    | D  [(Double,Double)]         -- | A D  type Gaussian primitive (Coefficient,Exponent)
    deriving (Generic, Show)

instance Serialize GaussShape
             
-- Basic elements living along the lines of a basis file
data Element =   Atom  AtomLabel [GaussShape]  deriving (Generic, Show)

instance Serialize Element

-- ====================================== Main processing functions ================================

parseBasisFile :: FilePath -> IO [Element]
parseBasisFile fname = parseFromFile parseBasis fname 
   
parseBasis :: Parser [Element]
parseBasis = undefined

