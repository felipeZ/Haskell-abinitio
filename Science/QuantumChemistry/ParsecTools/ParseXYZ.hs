{-# LANGUAGE OverloadedStrings #-}


module Science.QuantumChemistry.ParsecTools.ParseXYZ where

-- ====================> Standard Modules and third party <==============================
import Control.Applicative 
import Data.Attoparsec.ByteString.Char8
import Data.ByteString.Char8 (unpack)

-- ====================> Internal Modules <=======================
import Science.QuantumChemistry.GlobalTypes (AtomData)
import  Science.QuantumChemistry.ParsecTools.ParseUtils 

-- ==========================> Types <=================================
type Atom = (String,[Double])

-- =========================> <================================       

-- | 
parseFileXYZ :: FilePath -> IO [Atom]
parseFileXYZ file = parseFromFile parserXYZ file 

-- | 
parserXYZ :: Parser [Atom]
parserXYZ = do
   n <- spaces *> decimal  
   count n parseAtomXYZ
     where  parseAtomXYZ = do
              l  <- spaceAscii
              xs <- count 3 spaceDouble
              spaces *> endOfLine
              return (unpack l,xs)

