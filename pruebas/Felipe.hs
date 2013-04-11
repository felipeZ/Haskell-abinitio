{-# LANGUAGE BangPatterns #-}

module Main where

import System.IO
import Control.DeepSeq
-- import Prueba
import Prueba
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.ByteString.Lex.Lazy.Double as L
import qualified Data.Vector.Unboxed             as VU
import  Data.Array.Repa.Unsafe          as R
import  Data.Array.Repa                 as R
--import System.Environment


-- Fill a new vector from a file containing a list of numbers.
parse = VU.unfoldr step
  where
     step !s = case L.readDouble s of
        Nothing       -> Nothing
        Just (!k, !t) -> Just (k, L.tail t)

readDouble  x = read x :: Double
main = do
--     h <- openFile "symm.txt" ReadMode
--     s <- hGetContents h
--     s `deepseq` hClose h
-- --    v <- ((map readDouble). words) `fmap` readFile "symm.txt"
--     let v = ((map readDouble). words) $ s
--     let r = (v `seq` list2ArrDIM2 100 v)
    s   <- L.readFile "symm.txt"
    let r = parse $ s
    x <- jacobiP $ R.fromUnboxed ((Z :. 100 :. 100) :: DIM2) r        
    print . show $ x
