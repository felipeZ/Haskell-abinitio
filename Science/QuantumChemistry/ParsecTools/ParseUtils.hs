{-|
Module: Science.QuantumChemistry.ParsecTools.ParseUtils 
Description: Attoparsec utilities.
Copyright: @2016 Felipe Zapata
-}


module Science.QuantumChemistry.ParsecTools.ParseUtils where

-- ====================> Standard Modules and third party <==============================
import Control.Applicative
import Control.Exception (throwIO)
-- import Data.Attoparsec.ByteString
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString.Char8 as B


-- ====================> Internal Modules <=======================
import Science.QuantumChemistry.Error (HSFOCKException(..))

-- ======================================================================

-- | Similar to the Parsec utility
parseFromFile :: Parser a -> FilePath -> IO a
parseFromFile p file = do
 xs <- B.readFile file
 case parseOnly p xs of
      Left  msg  -> throwIO ParseError
      Right rs   -> return rs 

-- | Skip the content till the pattern is found
skipTill :: B.ByteString -> Parser ()
skipTill pattern = skipWhile (/= head (B.unpack pattern)) *>
  ( (string pattern *> pure () )  <|> (anyChar *> skipTill pattern))

-- | Return line content
anyLine :: Parser B.ByteString
anyLine = takeTill (== '\n')

-- | Discard content until end of line
anyLine' :: Parser ()
anyLine' = skipWhile (/= '\n') *> endOfLine 

-- | Parse an integer number preceded by one or more spaces
spaceDecimal :: Parser Int
spaceDecimal = spaces *> decimal

-- | Parse a real number preceded by one or more spaces
spaceDouble :: Parser Double
spaceDouble = spaces *> double

-- | Parse as many ascii characters as possible, preceding by 1 or more spaces
spaceAscii :: Parser B.ByteString
spaceAscii =  spaces *> takeWhile1 isAlpha_ascii

spaces :: Parser ()
spaces = skipWhile isSpace


