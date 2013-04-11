{-# LANGUAGE BangPatterns #-}

module Main where
import qualified Data.ByteString.Lazy.Char8      as L
import qualified Data.ByteString.Lex.Lazy.Double as L
import qualified Data.Vector.Unboxed             as VU
import Data.Array.Repa.Unsafe          as R
import Data.Array.Repa                 as R
import Data.Time.Clock
import Control.Exception (evaluate)
import Criterion.Main
import System.Environment (getArgs)
-- import System.Environment

import Prueba

   -- Fill a new vector from a file containing a list of numbers.
parse = VU.unfoldr step
  where
     step !s = case L.readDouble s of
        Nothing       -> Nothing
        Just (!k, !t) -> Just (k, L.tail t)

readDouble  x = read x :: Double

jacobiBench inputData = do
    x <- jacobiP $ R.fromUnboxed ((Z :. 100 :. 100) :: DIM2) inputData
    evaluate x

main = do
    args <- getArgs
    putStrLn $ "Loading Jacobi algorithm test data symm.txt file..."
    source   <- L.readFile "symm.txt"
    putStrLn $ "Parsing Jacobi data..."
    let sourceData = parse source
    case "bench" `elem` args of
         True -> do
            putStrLn $ "Criterion benchmark:"
            defaultMain [
                bench "Jacobi" (jacobiBench sourceData) 
                ]
         False -> do
             putStrLn $ "Normal benchmark:"
             result <-jacobiBench sourceData
             return ()

   -- mai-     source   <- L.readFile "symm.txt"
--     let sourceData = parse $ source
--     start <- getCurrentTime
--     jacobiBench sourceData
--     end <- getCurrentTime
--     putStrLn $ "Jacobi algorithm on symm.txt file took " Prelude.++ show (diffUTCTime end start)
