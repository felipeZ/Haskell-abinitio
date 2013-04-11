
module Main where

import System.IO
import Control.DeepSeq
-- import Prueba
import Prueba

readDouble  x = read x :: Double
main = do
    h <- openFile "symm.txt" ReadMode
    s <- hGetContents h
    s `deepseq` hClose h
--    v <- ((map readDouble). words) `fmap` readFile "symm.txt"
    let v = ((map readDouble). words) $ s
    let r = (v `deepseq` list2ArrDIM2 100 v)
    x <- jacobiP r
    let result = show x
    result `deepseq` return ()
