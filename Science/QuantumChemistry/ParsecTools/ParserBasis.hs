{-# LANGUAGE DeriveGeneric, OverloadedStrings #-}
-- Basis Parser, a parsec based parser for quantum chemistry basis files
-- @2013-2015 Angel Alvarez, Felipe Zapata
-- 
--  2013/05/05 (Spain, Mother's Day) We kindly appreciate our mothers efforts in bringing us so far...

module Science.QuantumChemistry.ParsecTools.ParserBasis
    (
        parseBasisFile
    ) where

-- ====================> Standard Modules and third party <==============================    
import GHC.Generics
import Data.Attoparsec.ByteString.Char8 

-- ====================> Internal Modules <=======================
import Science.QuantumChemistry.GlobalTypes 
import Science.QuantumChemistry.ParsecTools.ParseUtils 

-- ================================================ Types =============================================================
type AtomLabel = String

data GaussShape = 
      S  Double Double               -- | A S  type Gaussian primitive
    | SP Double Double Double Double -- | A SP type Gaussian primitives contains shared data among S and P types, gaussian functions 
    | P  Double Double               -- | A P  type Gaussian primitive
    | D  Double Double               -- | A D  type Gaussian primitive
    deriving (Generic, Show)
    
-- Basic elements living along the lines of a basis file
data Element =   Atom  AtomLabel [GaussShape]  deriving (Generic, Show)

-- ====================================== Main processing functions ================================

parseBasisFile :: FilePath -> IO [Element]
parseBasisFile fname = parseFromFile parseBasis fname 
   
parseBasis :: Parser [Element]
parseBasis = undefined

