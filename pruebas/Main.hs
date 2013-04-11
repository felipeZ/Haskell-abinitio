{-# LANGUAGE BangPatterns #-}

module Main where

import Prueba
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.ByteString.Lex.Lazy.Double as L
import qualified Data.Vector.Unboxed             as VU
import  Data.Array.Repa.Unsafe          as R
import  Data.Array.Repa                 as R
import Control.Exception (evaluate)
--import System.Environment


-- Fill a new vector from a file containing a list of numbers.
parse = VU.unfoldr step
  where
     step !s = case L.readDouble s of
        Nothing       -> Nothing
        Just (!k, !t) -> Just (k, L.tail t)

readDouble  x = read x :: Double
main = do
    s   <- L.readFile "symm.txt"
    let r = parse $ s
    x <- jacobiP $ R.fromUnboxed ((Z :. 100 :. 100) :: DIM2) r
    evaluate x
 