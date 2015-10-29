{-# LANGUAGE  RecordWildCards #-}
-- Command line checkers..
-- @2015 Felipe Zapata core SCF machinery 

module Science.QuantumChemistry.HsFock.OptsCheck (checkOpts) where

-- =============================> Standard and third party libraries <===============================
import Control.Monad (unless, when)
import System.Directory ( doesDirectoryExist, doesFileExist )
import System.Exit (exitFailure)
import Text.Printf

-- =================> Internal Modules <=====================
import Science.QuantumChemistry.GlobalTypes (HSFOCK(..))

-- ==============================<>==============================

-- | Sanity check of arguments
checkOpts :: HSFOCK -> IO ()
checkOpts opts@HSFOCK{..}= do
  unless  scf       $ printExit "We only Serve HartreeFock SCF (for the moment)"
  when (null basis) $ printExit "You must provide a basis set name"
  when (charge < 0) $ printExit "Charge must be a natural number"
  when (multi <= 0) $ printExit "Multiplicity must be 1 (singlet), 2(doblet), etc."
  --todo check that basis set exist
  [xyzBool,outBool]   <- mapM doesFileExist [xyz,outFile]
  unless  xyzBool . printExit $ printf "Coordinates File %s does not exist" xyz 
  when (outBool)  . printExit $ printf "file %s is going to be overwritten\n" outFile 

   where printExit s = putStrLn s *> exitFailure 

checkOpts opts@BasisConfig{..}= do
   fileBool <- doesFileExist basisPath
   when fileBool ( printf "The file:%s\n already exists, it is going to overwritten\n" basisPath)
