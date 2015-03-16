
module  Science.QuantumChemistry.ParsecTools.ParsecNumbers where

import Data.Char
import Control.Monad
import Text.Parsec 
import Text.Parsec.Bytestring (Parser)


-- 
-- Parse a real as <Integer>'.'<Decimal>['E'<Integer>] 
-- Usually : "['+','-']"[0-9]"."{8,[0-9]}"E"["-","+"][0-9][0-9]
realNumber :: Parser  Double
realNumber = try $ do
    integerPart    <- integerNumber
    fractionalPart <- fractionNumber
    exponentPart   <- option 1.0 exponentNumber
    let baseNumber = case integerPart > 0 of
         True -> (fromInteger integerPart) + fractionalPart
         False ->(fromInteger integerPart) - fractionalPart
    return (baseNumber * exponentPart)

-- Parse a fractional part , including leading '.'
fractionNumber :: Parser  Double
fractionNumber = do
    char '.'                                   -- leading dot
    digits     <- many1 digit                  -- digits after the dot
    let result = foldr op 0.0 digits
    seq result (return result)                 -- force computation before returning 
    where
        op d f = (f + fromIntegral (digitToInt d))/10.0

-- Parse an exponent part 
exponentNumber :: Parser  Double
exponentNumber = do
    oneOf "eE"
    f          <- sign
    e          <- decimalNumber
    let result = power (f e)
    seq result (return result)
    where
        power e  | e < 0      = 1.0/power(-e)
                 | otherwise  = fromInteger (10^e)


-- Parse an int as "-/+"{*,[0-9]} 
intNumber :: Parser  Int
intNumber = do
    number <- integerNumber
    return $ fromIntegral number

-- Parse an integer as "-/+"{*,[0-9]} 
integerNumber :: Parser  Integer
integerNumber = do
    f <- sign
    n <- decimalNumber
    return (f n)
 
-- Parse sign '+-'
sign :: Parser  (Integer -> Integer)
sign = 
    (char '-' >> return negate) 
    <|> (char '+' >> return id)
    <|> return id

-- Parse a number in base 16 (computer geeks like this)
hexadecimalNumber :: Parser  Integer
hexadecimalNumber = oneOf "xX" >> naturalNumber 16 hexDigit

-- Parse a number in base 10
decimalNumber :: Parser  Integer
decimalNumber = naturalNumber 10 digit

-- Parse one or more digits into a suitable natural number (as an integer in base BaseDigit)
naturalNumber :: Integer -> Parser  Char -> Parser Integer
naturalNumber base baseDigit = do
    digits <- many1 baseDigit
    let n = foldl (\x d -> base*x + toInteger (digitToInt d)) 0 digits
    seq n (return n)

