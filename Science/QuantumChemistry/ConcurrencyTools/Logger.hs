
 {-|
Module: Science.QuantumChemistry.ConcurrencyTools.Logger where
Description: Restricted Hartree-Fock Method
Copyright: @2012,2013 Angel Alvarez Adhesive tape
           @2012,2015 Felipe Zapata core SCF machinery 
           @2016 Felipe Zapata

Function to log events Taking as model the logger Implementation of 
Simon Marlow in "Parallel and Concurrent Programming in Haskell.
-}
module Science.QuantumChemistry.ConcurrencyTools.Logger where

import Control.Concurrent
import Control.Concurrent.Chan


-- ==========> <==========

data Logger = Logger (Chan LogCommad)

data LogCommad = Message String | Stop (MVar ())


-- ===========> <=========

initLogger :: FilePath -> IO Logger
initLogger name = do 
  c <- newChan 
  let l = Logger c
  forkIO (logger l name)
  return l 

logMessage :: Logger -> String -> IO ()
logMessage (Logger c) s  = writeChan c $ Message s

logStop :: Logger -> IO ()
logStop (Logger c) = do
  s <- newEmptyMVar
  writeChan  c (Stop s)
  takeMVar s

logger :: Logger -> FilePath -> IO ()
logger (Logger c) name = loop
  where 
    loop = do
      cmd <- readChan c
      case cmd of
           Message msg -> do 
                        appendFile name msg
                        loop
           Stop s -> do
                    print "logger: stop"
                    putMVar s ()
                    
action = do
  l <- initLogger "log"
  logMessage l "prueba1\n"
  logMessage l "launch!\n"
  logStop l
      

