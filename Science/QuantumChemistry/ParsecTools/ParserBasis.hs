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

comment :: Parser ()
comment = char '#' *> anyLine' *> endOfLine

basisHeader :: Parser ()
basisHeader = "BASIS \"ao basis\" PRINT" *> endOfLine

parseLabel :: Parser ByteString
parseLabel = takeWhile1 isAlpha_ascii
parsePrimitiveBlock :: Parser GaussShape
parsePrimitiveBlock = 

#Basis SET: (4s,1p) -> [2s,1p]
H    S
     18.7311370              0.03349460       
      2.8253937              0.23472695       
      0.6401217              0.81375733       
H    S
      0.1612778              1.0000000        
H    P
      1.1000000              1.0000000        
#BASIS SET: (4s,1p) -> [2s,1p]
He    S
     38.4216340              0.0237660        
      5.7780300              0.1546790        
      1.2417740              0.4696300        
He    S
      0.2979640              1.0000000        
He    P
      1.1000000              1.0000000
