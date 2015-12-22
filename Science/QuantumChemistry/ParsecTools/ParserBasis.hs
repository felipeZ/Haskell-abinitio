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
parseBasis = parseIntro *> many1 parseOneBasis

-- | Skips the comments explaining the basis set
--   Till the basis set starts
parseIntro :: Parser ()
parseIntro = skipTill "BASIS \"ao basis\" PRINT" *> endOfLine

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
    [atomName,shape] <-  count 2 (parseLabel <* spaces)
    primitives <- many1 spaceDouble
    anyLine'
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
comment = char '#' *> anyLine'

-- | For convention primitives are stored as (Coefficient, Exponent). 
-- #BASIS SET: (10s,4p,1d) -> [3s,2p,1d]
-- C    S
--    3047.5249000              0.0018347        
--     457.3695100              0.0140373        
--     103.9486900              0.0688426        
--      29.2101550              0.2321844        
--       9.2866630              0.4679413        
--       3.1639270              0.3623120        
-- C    SP
--       7.8682724             -0.1193324              0.0689991        
--       1.8812885             -0.1608542              0.3164240        
--       0.5442493              1.1434564              0.7443083        
-- C    SP
--       0.1687144              1.0000000              1.0000000        
-- C    D
--       0.8000000              1.0000000        
--
-- We operate with primitives using the (Coefficient, Exponent) representation
triple ::  [Double] -> [(Double,Double,Double)]
triple = map fun . chunksOf 3
  where fun [e,c1,c2] = (c1,c2,1)

-- |         
pairs :: [Double] -> [(Double,Double)]
pairs = map fun . chunksOf 2
  where fun [e,c] = (c,e)
  
  
