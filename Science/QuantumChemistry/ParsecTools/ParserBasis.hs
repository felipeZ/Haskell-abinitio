{-# LANGUAGE DeriveGeneric, LambdaCase , OverloadedStrings, TupleSections #-}
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
import Control.Applicative (liftA)
import Data.Attoparsec.ByteString.Char8 hiding (Number(..))
import qualified Data.ByteString.Char8 as B
import Data.List.Split (chunksOf)
import Data.Serialize 
import GHC.Generics

-- ====================> Internal Modules <=======================
import Science.QuantumChemistry.GlobalTypes 
import Science.QuantumChemistry.ParsecTools.ParseUtils 

-- ================================================ Types =============================================================
type AtomLabel = B.ByteString

data GaussShape = 
      S0 [(Double,Double)]         -- ^ A S  type Gaussian primitive (Coefficient,Exponent)
    | SP [(Double,Double,Double)]  -- ^ A SP type Gaussian primitives contains shared data among S and P types, gaussian functions (Coeficcient S , Coefficcient P, Exponent)
    | P  [(Double,Double)]         -- ^ A P  type Gaussian primitive (Coefficient,Exponent)
    | D  [(Double,Double)]         -- ^ A D  type Gaussian primitive (Coefficient,Exponent)
    deriving (Generic, Show)

instance Serialize GaussShape
             
data Element =   Atom  AtomLabel [GaussShape]  deriving (Generic, Show)

instance Serialize Element

-- ====================================== Main processing functions ================================

-- | Parse the primitives (Exponent and Coefficients) for several atomic elements.
-- The Basis are download from <https://bse.pnl.gov/ bse> as plain text.
parseBasisFile :: FilePath -> IO [Element]
parseBasisFile fname = parseFromFile parseBasis fname 
   
parseBasis :: Parser [Element]
parseBasis = parseIntro *> basisHeader *> many1 parseOneBasis

-- | Basis set start
basisHeader :: Parser ()
basisHeader = "BASIS \"ao basis\" PRINT" *> endOfLine

-- | Skips the comments explaining the basis set
parseIntro :: Parser ()
parseIntro = many1 comment *> pure ()

parseOneBasis :: Parser Element
parseOneBasis = comment *> liftA funElement (many1 parsePrimitiveBlock)
 where funElement :: [(AtomLabel, GaussShape)] -> Element
       funElement xs =
                 let l = fst . head $ xs
                 in  Atom l (map snd xs)
     
parseLabel :: Parser B.ByteString
parseLabel = takeWhile1 isAlpha_ascii

-- | Primitive Exponents and coefficients for an atomic element
parsePrimitiveBlock :: Parser (AtomLabel, GaussShape)
parsePrimitiveBlock = do
    atomName   <- parseLabel <* spaces
    shape      <- parseLabel <* spaces
    primitives <- many1 spaceDouble <* endOfLine
    funShape atomName primitives shape 

-- | Pack The GaussShape in triples if the S and P Share the exponents
--   Otherwise makes pairs for coefficients and exponents
funShape :: AtomLabel -> [Double] -> B.ByteString -> Parser (AtomLabel, GaussShape)
funShape atomName primitives = return . (atomName, ) . fun 
  where ps  = pairs primitives
        ts  = triple primitives 
        fun = \case
           "SP" ->  SP ts
           "S"  ->  S0 ps
           "P"  ->  P  ps
           "D"  ->  D  ps
           
comment :: Parser ()
comment = char '#' *> anyLine' *> endOfLine

triple ::  [Double] -> [(Double,Double,Double)]
triple = map fun . chunksOf 3
  where fun [x,y,z] = (x,y,z)
  
pairs :: [Double] -> [(Double,Double)]
pairs = map fun . chunksOf 2
  where fun [x,y] = (x,y)
  
  
