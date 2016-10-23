{-# LANGUAGE OverloadedStrings #-}

{-|
Module: Science.QuantumChemistry.ParsecTools.ParseXYZ
Description: Read molecular geometries in XYZ format
Copyright: @2016 Felipe Zapata
-}

module Science.QuantumChemistry.ParsecTools.ParseXYZ where

-- ====================> Standard Modules and third party <==============================
import Control.Applicative ((<|>))
import Data.Attoparsec.ByteString.Char8
import Data.ByteString.Char8 (unpack)
import Text.Printf

-- ====================> Internal Modules <=======================
import Science.QuantumChemistry.Error (HSFOCKException(..))
import Science.QuantumChemistry.GlobalTypes (AtomData)
import Science.QuantumChemistry.ParsecTools.ParseUtils 


-- ==========================> Types <=================================
type Atom = (String,[Double])

-- =========================> <================================       

-- | Reads Molecular Geometry in xyz format
parseFileXYZ :: FilePath -> (String -> IO () ) -> IO [Atom]
parseFileXYZ file logger = parseFromFile parserXYZ file

-- | Actual xyz parser
parserXYZ :: Parser [Atom]
parserXYZ = do
   n <- spaces *> decimal
   count 2 anyLine'
   count n parseAtomXYZ
     where  parseAtomXYZ = do
              l  <- spaceAscii
              xs <- count 3 spaceDouble
              endOfLine <|> endOfInput
              return (unpack l,xs)

