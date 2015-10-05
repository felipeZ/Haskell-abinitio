
module Science.QuantumChemistry.ConcurrencyTools.ConcurrencyUtilities where

import Control.Concurrent
import Control.Concurrent.Async


concurrentlyAll :: [IO a] -> IO [a]
concurrentlyAll =  foldr conc (return []) 
  where conc a acc = do
                  (a,as) <- concurrently a acc
                  return (a:as)
