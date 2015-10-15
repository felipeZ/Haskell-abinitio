

module  Science.QuantumChemistry.ParsecTools.ParseUtils where

-- ====================> Standard Modules and third party <==============================
import Control.Applicative
import Control.Exception (throwIO)
-- import Data.Attoparsec.ByteString
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString as B


-- ====================> Internal Modules <=======================
import Science.QuantumChemistry.Error (HSFOCKException(..))

-- ======================================================================
-- | 
parseFromFile :: Parser a -> FilePath -> IO a
parseFromFile p file = do
 xs <- B.readFile file
 case parseOnly p xs of
      Left  msg  -> throwIO ParseError
      Right rs   -> return rs 



anyLine :: Parser B.ByteString
anyLine = takeTill  (== '\n')
 
anyLine' :: Parser ()
anyLine' = skipWhile (/= '\n') *> endOfLine 

spaceDecimal :: Parser Int
spaceDecimal = spaces *> decimal

spaceDouble :: Parser Double
spaceDouble = spaces *> double

-- | Parse as many ascii characters as possible, preceding by 1 or more spaces
spaceAscii :: Parser B.ByteString
spaceAscii =  spaces *> takeWhile1 isAlpha_ascii

spaces :: Parser ()
spaces = skipWhile isSpace
