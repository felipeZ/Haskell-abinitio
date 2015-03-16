

module  Science.QuantumChemistry.ParsecTools.ParsecText where

import Data.Char
import Text.Parsec 
import Text.Parsec.ByteString (Parser)

oneOfStrings :: [String] -> Parser String
oneOfStrings listOfStrings = choice $ map (try . string) listOfStrings

-- Parse a text line as a whole into a string
anyLine :: Parser String
anyLine = manyTill anyChar newline     -- whatever chars we find till we hit a newline

-- Parser a line containing some string
stringLine :: String -> Parser String
stringLine str = do
    spaces
    result <- string str
    manyTill anyChar newline
    return result

blankLines :: Parser  [String]
blankLines = many1 blankLine

-- Parse a blank line as a whole into a string
blankLine :: Parser   String
blankLine = manyTill space newline     -- whatever spaces we find till we hit a newline
